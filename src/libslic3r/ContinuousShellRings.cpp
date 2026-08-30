#include "ContinuousShellRings.hpp"

#include "AABBTreeLines.hpp"
#include "EdgeGrid.hpp"
#include "MarchingSquares.hpp"

#include <tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <limits>

// The field raster, registered with the in-tree marching-squares implementation. Only the
// topology of the level set comes from the raster; every vertex it produces is then Newton-
// refined onto the exact level set, so the grid may be coarse.
namespace Slic3r { struct CsURaster { std::vector<float> v; size_t rows = 0, cols = 0; }; }

namespace marchsq {
template<> struct _RasterTraits<Slic3r::CsURaster>
{
    using Rst       = Slic3r::CsURaster;
    using ValueType = float;
    static float  get(const Rst &r, size_t row, size_t col) { return r.v[row * r.cols + col]; }
    static size_t rows(const Rst &r) { return r.rows; }
    static size_t cols(const Rst &r) { return r.cols; }
};
} // namespace marchsq

namespace Slic3r {

using LinesDist = AABBTreeLines::LinesDistancer<Line>;

struct CsRingField::Impl
{
    LinesDist all;    // contour + holes: region membership, by even-odd ray cast
    LinesDist out;    // contour only:    d_out and its nearest point
    LinesDist bores;  // every hole:      d_in and its nearest point
    LinesDist medial;  // bore-less regions: the medial axis, which is what d_in closes onto
    CsURaster raster;

    // Gridded distance for the RASTER only. An exact tree query costs about 2 us, and the
    // raster wants two of them in every cell - which measured 43 ms on a 20 mm section and
    // 900 ms on a 100 mm one, most of the field's cost. EdgeGrid's Danielsson-propagated SDF
    // answers the same question by bilinear lookup, and the raster is a topology carrier, so
    // a fraction of a cell of error in it costs nothing: every vertex it yields is corrected
    // against the exact trees afterwards.
    EdgeGrid::Grid sdf_out, sdf_bores;
    // EdgeGrid::Contour keeps raw pointers into the caller's Points, so the polygons it was
    // built from have to outlive the grid. Handing `create` a temporary reads freed memory.
    Polygons sdf_out_src, sdf_bores_src;
};

// The raster is a topology carrier, not a geometry carrier: every vertex it yields is Newton-
// corrected onto the exact level set afterwards, so a coarse grid costs topology resolution
// and not accuracy. It must still be fine enough that the thinnest neck of the region does not
// close up between samples. These two caps are what stop a plate-sized region from spending
// the whole slicing budget here - a 50 mm annulus at a 0.1 mm pitch is a megapixel per layer.
static constexpr size_t CS_FIELD_MAX_SIDE  = 1024;
static constexpr size_t CS_FIELD_MAX_CELLS = 250000;

bool CsRingField::build(const ExPolygon &interior, double grid_pitch_mm, std::string *fail,
                        const Points &sites)
{
    auto bail = [&](const char *why) {
        if (fail) *fail = std::string("continuous shell ring field: ") + why;
        m_valid = false;
        return false;
    };

    m_valid = false;
    m_region = interior;
    m_sites  = sites;
    m_impl = std::make_shared<Impl>();

    if (interior.contour.points.size() < 3)
        return bail("the interior region has no contour");
    if (grid_pitch_mm <= 0.)
        return bail("grid pitch must be positive");

    const BoundingBox bb = get_extents(interior);
    if (! bb.defined || bb.size().x() <= 0 || bb.size().y() <= 0)
        return bail("the interior region has no extent");

    // One pixel of margin on every side so a level set that runs close to the contour still
    // has a cell outside it to close against.
    double pitch_s = scale_(grid_pitch_mm);
    const double side = double(std::max(bb.size().x(), bb.size().y()));
    if (side / pitch_s > double(CS_FIELD_MAX_SIDE))
        pitch_s = side / double(CS_FIELD_MAX_SIDE);
    const double cells = (double(bb.size().x()) / pitch_s) * (double(bb.size().y()) / pitch_s);
    if (cells > double(CS_FIELD_MAX_CELLS))
        pitch_s *= std::sqrt(cells / double(CS_FIELD_MAX_CELLS));
    m_pitch = pitch_s;

    m_cols = size_t(std::ceil(double(bb.size().x()) / pitch_s)) + 3;
    m_rows = size_t(std::ceil(double(bb.size().y()) / pitch_s)) + 3;
    m_origin = Vec2d(double(bb.min.x()) - pitch_s, double(bb.min.y()) - pitch_s);

    // Distance trees. `all` answers membership over the whole ExPolygon (even-odd, so a
    // point inside a hole reads as outside the region); `out` and `bores` carry the two
    // distances the field is built from.
    Lines contour_lines = interior.contour.lines();
    Lines all_lines     = contour_lines;
    Lines bore_lines;
    m_bore_of_line.clear();
    m_n_bores = interior.holes.size();
    for (size_t h = 0; h < interior.holes.size(); ++h) {
        Lines hl = interior.holes[h].lines();
        for (const Line &l : hl) {
            bore_lines.emplace_back(l);
            m_bore_of_line.push_back(h);
        }
        all_lines.insert(all_lines.end(), hl.begin(), hl.end());
    }
    m_impl->all   = LinesDist(std::move(all_lines));
    m_impl->out   = LinesDist(std::move(contour_lines));
    m_impl->bores = LinesDist(std::move(bore_lines));

    // With no bore the family has to close onto something, and the hypothesis said the centroid.
    // The centroid is right only when the region is round about it: for a disc, d_out + d_centroid
    // is exactly the radius everywhere, so the extent is constant and the rings come out evenly.
    // For anything elongated or curved it is not the thickness at all but "distance from the wall
    // plus distance to a far-away point", and that inflates the extent enormously - on a crescent
    // it measured 2.70..14.71 mm, a 5.4-fold spread on a wall whose half-thickness never exceeds
    // 5.5 mm, and the ring count that came out of it was 12 where 4 is called for.
    //
    // The MEDIAL AXIS is what the centroid was standing in for. Closing onto it makes
    // d_out + d_medial the local half-thickness everywhere - the quantity that actually decides
    // how many rings fit - and it degenerates to the centroid for a disc, where the medial axis
    // IS the centre. Sites cannot substitute: no point placement bounds the extent's maximum
    // without collapsing its minimum, because the extent at a site is just d_out there.
    if (m_n_bores == 0 && m_sites.empty()) {
        // max_width filters the axis by the LOCAL WIDTH it stands for, so it has to admit the
        // whole region: bounded at a couple of millimetres - the value the sub-bead fallback
        // uses, for its own quite different purpose - it rejects every part of a wall thicker
        // than that and hands back nothing, which silently reverts to the centroid.
        const BoundingBox mb = get_extents(interior);
        const double max_w = 2.0 * double(std::max(mb.size().x(), mb.size().y()));
        Polylines med = interior.medial_axis(0.05 * scale_(1.), max_w);
        Lines ml;
        for (const Polyline &pl : med)
            for (const Line &l : pl.lines())
                if (l.length() > 0)
                    ml.push_back(l);
        if (ml.size() >= 2)
            m_impl->medial = LinesDist(std::move(ml));
    }
    if (m_n_bores == 0 && m_sites.empty() && m_impl->medial.get_lines().empty()) {
        const Point c = interior.contour.centroid();
        m_centroid = c.cast<double>();
        // A centroid outside its own contour (a horseshoe) would put the whole family on the
        // wrong side of the wall; that region needs re-sectioning, not a ring family.
        if (m_impl->all.outside(c) != -1)
            return bail("the region's centroid lies outside it - re-section before filling");
    }

    // Gridded distance for the raster. Half the raster pitch keeps the bilinear error well
    // under a cell, which is all the topology needs.
    {
        const coord_t res = std::max<coord_t>(coord_t(std::llround(pitch_s * 0.5)), 1);
        m_impl->sdf_out_src = Polygons{ interior.contour };
        m_impl->sdf_out.create(m_impl->sdf_out_src, res);
        m_impl->sdf_out.calculate_sdf();
        if (m_n_bores > 0) {
            m_impl->sdf_bores_src.clear();
            for (const Polygon &h : interior.holes)
                m_impl->sdf_bores_src.push_back(h);
            m_impl->sdf_bores.create(m_impl->sdf_bores_src, res);
            m_impl->sdf_bores.calculate_sdf();
        }
    }

    m_impl->raster.rows = m_rows;
    m_impl->raster.cols = m_cols;
    m_impl->raster.v.assign(m_rows * m_cols, 0.f);
    m_wgrid.assign(m_rows * m_cols, 0.f);

    // The extent statistics decide the ring count, so they are taken EXACTLY - the gridded
    // SDF is good enough to place a level set that gets corrected afterwards, but it reads a
    // 6.000 mm wall as 6.206 and that error would land in n. A percentile does not need every
    // cell, so the exact queries run on a subsample and cost a fraction of a full pass.
    const size_t stride = std::max<size_t>(1, size_t(std::llround(
        std::sqrt(double(m_rows * m_cols) / 4000.0))));
    std::vector<std::vector<float>> ext_rows(m_rows);
    std::vector<std::pair<float, Point>> deepest_rows(m_rows, { 0.f, Point(0, 0) });

    // Membership by SCANLINE, not by a ray cast per cell. The even-odd test that
    // `LinesDistancer::outside` performs walks the tree twice per point, and at a raster cell
    // per 0.1 mm that dominates everything else the field does - it measured 30 us a cell,
    // which is a second per layer on a 20 mm part.
    //
    // The crossings come from plain scan conversion over the edge list rather than from the
    // AABB tree: a full-width horizontal query descends most of the tree, so per row it costs
    // more than the two distance queries it was meant to save. A linear pass over the edges is
    // rows x edges in total, which is microseconds at these sizes. The half-open rule
    // (y0 <= y < y1) counts each vertex once, so a scanline through a vertex needs no
    // special case.
    struct Edge { double x0, y0, x1, y1; };
    auto edges_of = [](const Polygons &ps) {
        std::vector<Edge> e;
        for (const Polygon &poly : ps)
            for (size_t i = 0; i < poly.points.size(); ++i) {
                const Point &a = poly.points[i];
                const Point &b = poly.points[(i + 1) % poly.points.size()];
                if (a.y() != b.y())
                    e.push_back({ double(a.x()), double(a.y()), double(b.x()), double(b.y()) });
            }
        return e;
    };
    Polygons all_polys{ interior.contour };
    for (const Polygon &h : interior.holes)
        all_polys.push_back(h);
    const std::vector<Edge> edges_all = edges_of(all_polys);
    const std::vector<Edge> edges_out = edges_of(Polygons{ interior.contour });

    auto crossings = [](const std::vector<Edge> &edges, double y, std::vector<double> &xs) {
        xs.clear();
        for (const Edge &e : edges) {
            const double lo = std::min(e.y0, e.y1), hi = std::max(e.y0, e.y1);
            if (y < lo || y >= hi)
                continue;
            xs.push_back(e.x0 + (y - e.y0) * (e.x1 - e.x0) / (e.y1 - e.y0));
        }
        std::sort(xs.begin(), xs.end());
    };
    auto inside_spans = [](const std::vector<double> &xs, double x) {
        // even-odd: inside iff the number of crossings strictly left of x is odd
        return (std::upper_bound(xs.begin(), xs.end(), x) - xs.begin()) % 2 == 1;
    };

    tbb::parallel_for(tbb::blocked_range<size_t>(0, m_rows), [&](const tbb::blocked_range<size_t> &range) {
        std::vector<double> xs_all, xs_out;
        for (size_t r = range.begin(); r < range.end(); ++r) {
            const double y = m_origin.y() + double(r) * pitch_s;
            crossings(edges_all, y, xs_all);
            if (m_n_bores > 0)
                crossings(edges_out, y, xs_out);
            for (size_t c = 0; c < m_cols; ++c) {
                const Vec2d p(m_origin.x() + double(c) * pitch_s, y);
                double u = 0.;
                float  ext = 0.f;
                const bool inside = inside_spans(xs_all, p.x());
                if (inside) {
                    const Point ps(coord_t(std::llround(p.x())), coord_t(std::llround(p.y())));
                    const double d_out = std::abs(double(m_impl->sdf_out.signed_distance_bilinear(ps)));
                    if (float(d_out) > deepest_rows[r].first)
                        deepest_rows[r] = { float(d_out), ps };
                    double d_in = m_n_bores > 0
                        ? std::abs(double(m_impl->sdf_bores.signed_distance_bilinear(ps)))
                        : (!m_impl->medial.get_lines().empty()
                               ? m_impl->medial.distance_from_lines<false>(ps)
                               : (m_sites.empty() ? (p - m_centroid).norm()
                                                  : std::numeric_limits<double>::max()));
                    for (const Point &st : m_sites)
                        d_in = std::min(d_in, (p - st.cast<double>()).norm());
                    const double W = d_out + d_in;
                    ext = float(W * SCALING_FACTOR);
                    u = (W > 1e-9) ? std::clamp(d_out / W, 0., 1.) : 1.;
                } else {
                    // Outside the region. Inside the contour means inside a bore, where the
                    // field must read 1 so the level set wraps it; beyond the contour, 0.
                    u = (m_n_bores > 0 && inside_spans(xs_out, p.x())) ? 1. : 0.;
                }
                m_impl->raster.v[r * m_cols + c] = float(u);
                m_wgrid[r * m_cols + c]          = ext;
                if (ext > 0.f && (r % stride) == 0 && (c % stride) == 0) {
                    double eu = 0., eext = 0.;
                    int    ebore = -1;
                    this->sample_inside(p, eu, eext, ebore);
                    ext_rows[r].push_back(float(eext));
                }
            }
        }
    });

    m_extent.clear();
    for (const std::vector<float> &row : ext_rows)
        m_extent.insert(m_extent.end(), row.begin(), row.end());
    float deepest = 0.f;
    for (const auto &w : deepest_rows)
        if (w.first > deepest) { deepest = w.first; m_deepest = w.second; }

    m_valid = true;
    return true;
}

void CsRingField::sample(const Vec2d &p, double &u, double &extent, int &bore) const
{
    const Point ps(coord_t(std::llround(p.x())), coord_t(std::llround(p.y())));
    bore   = -1;
    extent = 0.;

    if (m_impl->all.outside(ps) != -1) {
        // Outside the region. Inside the contour means inside a bore, where the field must
        // read 1 so the level set wraps the bore; beyond the contour it must read 0. Keeping
        // both saturated makes u continuous across each boundary.
        u = (m_impl->out.outside(ps) == -1) ? 1. : 0.;
        return;
    }
    this->sample_inside(p, u, extent, bore);
}

// The field at a point already known to lie inside the region.
void CsRingField::sample_inside(const Vec2d &p, double &u, double &extent, int &bore) const
{
    const Point ps(coord_t(std::llround(p.x())), coord_t(std::llround(p.y())));
    bore = -1;

    const auto [d_out, io, po] = m_impl->out.distance_from_lines_extra<false>(ps);
    (void) io; (void) po;

    double d_in = 0.;
    if (m_n_bores > 0) {
        const auto [d, li, pi] = m_impl->bores.distance_from_lines_extra<false>(ps);
        (void) pi;
        d_in = d;
        if (li < m_bore_of_line.size())
            bore = int(m_bore_of_line[li]);
    } else if (!m_impl->medial.get_lines().empty()) {
        d_in = m_impl->medial.distance_from_lines<false>(ps);
    } else if (m_sites.empty()) {
        d_in = (p - m_centroid).norm();
    } else {
        d_in = std::numeric_limits<double>::max();
    }
    for (size_t k = 0; k < m_sites.size(); ++k) {
        const double d = (p - m_sites[k].cast<double>()).norm();
        if (d < d_in) { d_in = d; bore = int(m_n_bores + k); }
    }

    const double W = d_out + d_in;
    extent = W * SCALING_FACTOR;
    u = (W > 1e-9) ? std::clamp(d_out / W, 0., 1.) : 1.;
}

bool CsRingField::refine(Vec2d &p, double t, double &extent, int &bore) const
{
    extent = 0.;
    bore = -1;
    for (int it = 0; it < 4; ++it) {
        const Point ps(coord_t(std::llround(p.x())), coord_t(std::llround(p.y())));
        const auto [d_out, io, po] = m_impl->out.distance_from_lines_extra<false>(ps);
        (void) io;

        double d_in = 0.;
        Vec2d  pin  = m_centroid;
        if (m_n_bores > 0) {
            const auto [d, li, pi] = m_impl->bores.distance_from_lines_extra<false>(ps);
            d_in = d;
            pin  = pi;
            // The cell a ring point belongs to is which bore is nearest, and that is defined
            // wherever the point is - including a hair inside the bore wall, where a
            // membership test would refuse to answer and hand back "no bore".
            if (li < m_bore_of_line.size())
                bore = int(m_bore_of_line[li]);
        } else if (!m_impl->medial.get_lines().empty()) {
            const auto [d, li, pi] = m_impl->medial.distance_from_lines_extra<false>(ps);
            (void) li;
            d_in = d; pin = pi;
        } else if (m_sites.empty()) {
            d_in = (p - m_centroid).norm();
        } else {
            d_in = std::numeric_limits<double>::max();
        }
        for (size_t k = 0; k < m_sites.size(); ++k) {
            const Vec2d sp = m_sites[k].cast<double>();
            const double d = (p - sp).norm();
            if (d < d_in) { d_in = d; pin = sp; bore = int(m_n_bores + k); }
        }

        const double W = d_out + d_in;
        extent = W * SCALING_FACTOR;
        if (W < 1e-9 || d_out < 1e-9 || d_in < 1e-9)
            return this->bisect(p, t, extent, bore);   // on a boundary: Newton has no footing

        const double u   = d_out / W;
        const double err = u - t;
        if (std::abs(err) < 1e-4)
            return true;

        // grad u = (d_in * grad d_out - d_out * grad d_in) / W^2, with each grad d the unit
        // vector pointing away from that boundary's nearest point.
        const Vec2d g_out = (p - po) / d_out;
        const Vec2d g_in  = (p - pin) / d_in;
        const Vec2d grad  = (d_in * g_out - d_out * g_in) / (W * W);
        const double g2   = grad.squaredNorm();
        if (g2 < 1e-24)
            return this->bisect(p, t, extent, bore);   // the two boundaries pull the same way

        Vec2d step = -(err / g2) * grad;
        // Where the nearest-point directions are not opposed the gradient is small and the
        // Newton step overshoots; half the local extent is as far as a correction can ever
        // legitimately need to move.
        const double cap = 0.5 * W;
        if (step.norm() > cap)
            step *= cap / step.norm();
        p += step;
    }
    return this->bisect(p, t, extent, bore);
}

// The corrector that does not need a gradient. u runs from 0 on the outer boundary to 1 on
// the bore, so on the segment joining a point's two nearest boundary points it is bracketed
// and bisection always lands on the level set.
//
// This is not a rare fallback. u is continuous across a Voronoi cut - that is the property the
// whole subsectioning rests on - but it is NOT differentiable there, because the nearest bore
// swaps from one to the other. Newton has nothing to work with along a cut, at a crescent's
// tangency, or on a medial axis; bisection has no such trouble.
bool CsRingField::bisect(Vec2d &p, double t, double &extent, int &bore) const
{
    const Point ps(coord_t(std::llround(p.x())), coord_t(std::llround(p.y())));
    const auto [d_out, io, po] = m_impl->out.distance_from_lines_extra<false>(ps);
    (void) d_out; (void) io;
    Vec2d pin = m_centroid;
    double best = m_sites.empty() ? std::numeric_limits<double>::max() : (p - m_centroid).norm();
    if (!m_impl->medial.get_lines().empty()) {
        const auto [d, li, pi] = m_impl->medial.distance_from_lines_extra<false>(ps);
        (void) li;
        pin = pi; best = d;
    }
    if (m_n_bores > 0) {
        const auto [d, li, pi] = m_impl->bores.distance_from_lines_extra<false>(ps);
        (void) li;
        pin = pi; best = d;
    }
    for (const Point &s : m_sites) {
        const double d = (p - s.cast<double>()).norm();
        if (d < best) { best = d; pin = s.cast<double>(); }
    }

    Vec2d a = po, b = pin;                       // u(a) = 0, u(b) = 1
    if ((b - a).squaredNorm() < 1e-12)
        return false;
    double u = 0., ext = 0.;
    int    lab = -1;
    for (int it = 0; it < 30; ++it) {
        const Vec2d m = 0.5 * (a + b);
        this->sample_inside(m, u, ext, lab);
        if (std::abs(u - t) < 1e-4) {
            p = m; extent = ext; bore = lab;
            return true;
        }
        (u < t ? a : b) = m;
    }
    p = 0.5 * (a + b);
    this->sample_inside(p, u, ext, lab);
    extent = ext;
    bore   = lab;
    return std::abs(u - t) < 1e-3;
}

double CsRingField::extent_percentile(double pct) const
{
    if (! m_valid)
        return 0.;
    std::vector<float> w;
    w.reserve(m_extent.size() / 4);
    for (const float e : m_extent)
        if (e > 0.f)
            w.push_back(e);
    if (w.empty())
        return 0.;
    const size_t k = std::min(w.size() - 1,
                              size_t(std::llround(std::clamp(pct, 0., 1.) * double(w.size() - 1))));
    std::nth_element(w.begin(), w.begin() + k, w.end());
    return double(w[k]);
}

double CsRingField::extent_at(const Vec2d &p) const
{
    if (! m_valid || m_wgrid.empty())
        return 0.;
    const double fc = (p.x() - m_origin.x()) / m_pitch;
    const double fr = (p.y() - m_origin.y()) / m_pitch;
    const long   c0 = long(std::floor(fc)), r0 = long(std::floor(fr));
    if (c0 < 0 || r0 < 0 || size_t(c0) + 1 >= m_cols || size_t(r0) + 1 >= m_rows)
        return 0.;
    // Interpolate over the cells that HAVE a value, not over the zeros outside the region. A
    // ring that runs along the interior's own boundary - which the outermost and innermost
    // lattice rings nearly do - has outside cells among its four corners, and averaging their
    // zeros in halves the extent it reports. That put a uniform 6 mm wall at the bead-width
    // floor along both its faces.
    const double tx = fc - double(c0), ty = fr - double(r0);
    const double wt[4] = { (1. - tx) * (1. - ty), tx * (1. - ty), (1. - tx) * ty, tx * ty };
    const long   rr[4] = { r0, r0, r0 + 1, r0 + 1 };
    const long   cc[4] = { c0, c0 + 1, c0, c0 + 1 };
    double acc = 0., wsum = 0.;
    for (int k = 0; k < 4; ++k) {
        const double v = double(m_wgrid[size_t(rr[k]) * m_cols + size_t(cc[k])]);
        if (v > 0.) { acc += wt[k] * v; wsum += wt[k]; }
    }
    return wsum > 1e-9 ? acc / wsum : 0.;
}

double CsRingField::fraction_narrower_than(double w_mm) const
{
    if (! m_valid || m_extent.empty())
        return 0.;
    size_t below = 0, seen = 0;
    for (const float e : m_extent)
        if (e > 0.f) { ++seen; below += (double(e) < w_mm) ? 1 : 0; }
    return seen ? double(below) / double(seen) : 0.;
}

size_t CsRingField::ring_count_for_max_gap(double max_gap_mm, double pct) const
{
    if (! m_valid || max_gap_mm <= 0.)
        return 0;
    // The WIDEST place in the region, not a percentile of it. A percentile leaves the widest
    // quartile free to open past the bound, and the widest place is exactly where a lattice
    // stops welding - the failure this bound exists to prevent.
    const double W = this->extent_percentile(pct);
    const double need = W / max_gap_mm;          // gaps required to stay within the bound
    if (need <= 1.0 + 1e-9)
        return 0;                                // the skin already spans it
    return size_t(std::ceil(need - 1e-9)) - 1;
}

bool CsRingField::level_set(double t, std::vector<CsRing> &out, std::string *fail) const
{
    out.clear();
    if (! m_valid) {
        if (fail) *fail = "continuous shell ring field: level set asked of an unbuilt field";
        return false;
    }
    if (! (t > 0. && t < 1.)) {
        if (fail) *fail = "continuous shell ring field: level must lie strictly inside (0, 1)";
        return false;
    }

    const std::vector<marchsq::Ring> rings = marchsq::execute(m_impl->raster, float(t), {2, 2});

    for (const marchsq::Ring &ring : rings) {
        if (ring.size() < 8)
            continue;                            // a few pixels of noise, not a ring

        // Correct every vertex onto the exact level set. This is the whole reason the raster
        // is allowed to be coarse, and it is the bulk of the extraction cost, so it runs wide.
        std::vector<Vec2d>  pos(ring.size());
        std::vector<double> ext(ring.size(), 0.);
        std::vector<int>    lab(ring.size(), -1);
        for (size_t i = 0; i < ring.size(); ++i)
            pos[i] = Vec2d(m_origin.x() + double(ring[i].c) * m_pitch,
                           m_origin.y() + double(ring[i].r) * m_pitch);
        std::vector<char> ok(ring.size(), 0);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, pos.size()),
                          [&](const tbb::blocked_range<size_t> &range) {
            for (size_t i = range.begin(); i < range.end(); ++i)
                ok[i] = this->refine(pos[i], t, ext[i], lab[i]) ? 1 : 0;
        });

        // A vertex the correction could not place is a vertex that is not on the level set.
        // Where the two boundaries pull the same way - past the tangency of a crescent tip,
        // say - the field has no usable gradient and a few vertices legitimately fail; those
        // are dropped. A ring that loses a real share of itself that way is not a ring, and
        // saying so is better than emitting a loop with a chord cut across it.
        const size_t bad = std::count(ok.begin(), ok.end(), char(0));
        if (bad * 50 > ring.size()) {
            if (fail)
                *fail = "continuous shell ring field: the level set at t=" + std::to_string(t) +
                        " could not be placed - the field has no gradient over " +
                        std::to_string(100 * bad / ring.size()) + "% of it";
            out.clear();
            return false;
        }

        CsRing cr;
        cr.t = t;
        cr.poly.points.reserve(pos.size());
        cr.extent.reserve(pos.size());
        cr.bore.reserve(pos.size());
        for (size_t i = 0; i < pos.size(); ++i) {
            if (! ok[i])
                continue;
            const Point q(coord_t(std::llround(pos[i].x())), coord_t(std::llround(pos[i].y())));
            // Consecutive raster cells can correct onto the same point; a repeated vertex is a
            // zero-length segment downstream.
            if (! cr.poly.points.empty() && cr.poly.points.back() == q)
                continue;
            cr.poly.points.emplace_back(q);
            cr.extent.push_back(float(ext[i]));
            cr.bore.push_back(lab[i]);
        }
        if (cr.poly.points.size() > 2 && cr.poly.points.front() == cr.poly.points.back()) {
            cr.poly.points.pop_back();
            cr.extent.pop_back();
            cr.bore.pop_back();
        }
        if (cr.poly.points.size() < 8)
            continue;
        // Every ring of the family is an outward-facing loop, whatever it encircles, so the
        // downstream stack sees one consistent winding.
        if (! cr.poly.is_counter_clockwise()) {
            std::reverse(cr.poly.points.begin(), cr.poly.points.end());
            std::reverse(cr.extent.begin(), cr.extent.end());
            std::reverse(cr.bore.begin(), cr.bore.end());
        }
        out.emplace_back(std::move(cr));
    }

    if (out.empty()) {
        if (fail) *fail = "continuous shell ring field: a level set came back empty";
        return false;
    }
    return true;
}

bool CsRingField::family(size_t n, CsRingFamily &out, std::string *fail) const
{
    out.rings.clear();
    if (! m_valid) {
        if (fail) *fail = "continuous shell ring field: family asked of an unbuilt field";
        return false;
    }
    if (n == 0)
        return true;                             // no interior rings is a legitimate answer

    for (size_t j = 0; j < n; ++j) {
        const double t = double(j + 1) / double(n + 1);
        std::vector<CsRing> level;
        if (! this->level_set(t, level, fail)) {
            // With no bore the field closes onto a POINT, so the innermost level sets shrink
            // until there is nothing left to trace. That is the family ending, not a failure -
            // u rises monotonically toward the centroid, so once a level is empty every level
            // inside it is too. With a bore the level sets wrap it and never vanish, so there
            // an empty level really is one.
            if (m_n_bores == 0 && ! out.rings.empty()) {
                if (fail) fail->clear();
                return true;
            }
            return false;
        }
        for (CsRing &r : level)
            out.rings.emplace_back(std::move(r));
    }
    return true;
}

} // namespace Slic3r
