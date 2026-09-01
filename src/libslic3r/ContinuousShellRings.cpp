#include "ContinuousShellRings.hpp"

#include "AABBTreeLines.hpp"
#include "ClipperUtils.hpp"
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
    if ((m_n_bores == 0 && m_sites.empty()) || m_close_medial) {
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
                        ? std::min(std::abs(double(m_impl->sdf_bores.signed_distance_bilinear(ps))),
                                   m_close_medial && !m_impl->medial.get_lines().empty()
                                       ? m_impl->medial.distance_from_lines<false>(ps)
                                       : std::numeric_limits<double>::max())
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

    this->build_phi();

    m_valid = true;
    return true;
}

// The azimuth field, carried inward from the outer boundary along the gradient flow of u.
//
// Held as cos/sin rather than as an angle so the seam - where the arc length wraps from 1 back
// to 0 - averages correctly instead of tearing. Cells are visited in order of increasing u, so
// every cell is filled from the cells it came from and one pass suffices.
void CsRingField::build_phi()
{
    m_phi_c.assign(m_rows * m_cols, 0.f);
    m_phi_s.assign(m_rows * m_cols, 0.f);
    std::vector<char> known(m_rows * m_cols, 0);

    // Cumulative arc length along the outer boundary, so a nearest point on it can be turned
    // into a position around it.
    const Lines &cl = m_impl->out.get_lines();
    if (cl.empty())
        return;
    std::vector<double> arc(cl.size() + 1, 0.0);
    for (size_t i = 0; i < cl.size(); ++i)
        arc[i + 1] = arc[i] + cl[i].length();
    const double total = arc.back();
    if (total < 1e-9)
        return;

    const double seed_band = 2.0 * m_pitch;
    std::vector<std::pair<float, size_t>> order;
    order.reserve(m_rows * m_cols / 2);
    for (size_t r = 0; r < m_rows; ++r)
        for (size_t c = 0; c < m_cols; ++c) {
            const size_t k = r * m_cols + c;
            if (m_wgrid[k] <= 0.f)
                continue;                          // outside the region
            const Vec2d p(m_origin.x() + double(c) * m_pitch, m_origin.y() + double(r) * m_pitch);
            const Point ps(coord_t(std::llround(p.x())), coord_t(std::llround(p.y())));
            const auto [d, li, np] = m_impl->out.distance_from_lines_extra<false>(ps);
            if (d < seed_band && li < cl.size()) {
                const double a = (arc[li] + (np - cl[li].a.cast<double>()).norm()) / total;
                m_phi_c[k] = float(std::cos(2.0 * M_PI * a));
                m_phi_s[k] = float(std::sin(2.0 * M_PI * a));
                known[k] = 1;
            }
            order.emplace_back(m_impl->raster.v[k], k);
        }
    std::sort(order.begin(), order.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });

    static const int DR[8] = { -1, 1, 0, 0, -1, -1, 1, 1 };
    static const int DC[8] = { 0, 0, -1, 1, -1, 1, -1, 1 };
    for (const auto &it : order) {
        const size_t k = it.second;
        if (known[k])
            continue;
        const long r = long(k / m_cols), c = long(k % m_cols);
        double sc = 0.0, ss = 0.0;
        int n = 0;
        for (int d = 0; d < 8; ++d) {
            const long rr = r + DR[d], cc = c + DC[d];
            if (rr < 0 || cc < 0 || size_t(rr) >= m_rows || size_t(cc) >= m_cols)
                continue;
            const size_t g = size_t(rr) * m_cols + size_t(cc);
            if (known[g] && m_impl->raster.v[g] < it.first) {
                sc += double(m_phi_c[g]); ss += double(m_phi_s[g]); ++n;
            }
        }
        if (n) {
            const double h = std::hypot(sc, ss);
            if (h > 1e-9) { m_phi_c[k] = float(sc / h); m_phi_s[k] = float(ss / h); known[k] = 1; }
        }
    }
}

double CsRingField::phi_at(const Vec2d &p) const
{
    if (! m_valid || m_phi_c.empty())
        return 0.;
    const double fc = (p.x() - m_origin.x()) / m_pitch;
    const double fr = (p.y() - m_origin.y()) / m_pitch;
    const long   c0 = long(std::floor(fc)), r0 = long(std::floor(fr));
    if (c0 < 0 || r0 < 0 || size_t(c0) + 1 >= m_cols || size_t(r0) + 1 >= m_rows)
        return 0.;
    const double tx = fc - double(c0), ty = fr - double(r0);
    const double wt[4] = { (1. - tx) * (1. - ty), tx * (1. - ty), (1. - tx) * ty, tx * ty };
    const long   rr[4] = { r0, r0, r0 + 1, r0 + 1 };
    const long   cc[4] = { c0, c0 + 1, c0, c0 + 1 };
    double sc = 0., ss = 0.;
    for (int k = 0; k < 4; ++k) {
        const size_t g = size_t(rr[k]) * m_cols + size_t(cc[k]);
        if (m_wgrid[g] > 0.f) { sc += wt[k] * double(m_phi_c[g]); ss += wt[k] * double(m_phi_s[g]); }
    }
    if (std::hypot(sc, ss) < 1e-9)
        return 0.;
    double a = std::atan2(ss, sc) / (2.0 * M_PI);
    return a < 0. ? a + 1.0 : a;
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
        if (m_close_medial && !m_impl->medial.get_lines().empty()) {
            const double dm = m_impl->medial.distance_from_lines<false>(ps);
            if (dm < d_in) { d_in = dm; bore = -1; }   // the skeleton's, not any bore's
        }
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
            if (m_close_medial && !m_impl->medial.get_lines().empty()) {
                const auto [dm, lm, pm] = m_impl->medial.distance_from_lines_extra<false>(ps);
                (void) lm;
                if (dm < d_in) { d_in = dm; pin = pm; bore = -1; }
            }
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
        // Under medial closure the axis stands earlier in this chain and holds when it is
        // nearer; without it the medial lines are empty for a bored region and the bore
        // wins unconditionally, exactly as before.
        if (!(m_close_medial && !m_impl->medial.get_lines().empty()) || d < best) {
            pin = pi; best = d;
        }
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

        // SMOOTH, THEN RE-PROJECT.
        //
        // Every vertex is corrected onto the level set independently, so each lands with its own
        // residual and the polyline zig-zags at the sub-vertex scale. Measured on an eccentric
        // annulus with the weave switched off entirely, 1% of joints turned more than 30 degrees
        // where a smooth curve sampled this finely turns about 2 - so the rings were already
        // tangled before anything displaced them, and every attempt to smooth the WEAVE was
        // treating a symptom. A Laplacian pass pulls each vertex toward its neighbours and a
        // correction puts it back on the level set: the noise goes and the curve stays exact.
        for (int pass = 0; pass < 3; ++pass) {
            std::vector<Vec2d> sm(pos.size());
            for (size_t i = 0; i < pos.size(); ++i) {
                const Vec2d &a = pos[(i + pos.size() - 1) % pos.size()];
                const Vec2d &b = pos[(i + 1) % pos.size()];
                sm[i] = 0.25 * a + 0.5 * pos[i] + 0.25 * b;
            }
            tbb::parallel_for(tbb::blocked_range<size_t>(0, sm.size()),
                              [&](const tbb::blocked_range<size_t> &range) {
                for (size_t i = range.begin(); i < range.end(); ++i)
                    if (ok[i])
                        ok[i] = this->refine(sm[i], t, ext[i], lab[i]) ? 1 : 0;
            });
            pos.swap(sm);
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

// Even arc-length resample of a closed loop, in a consistent rotational direction.
static std::vector<Vec2d> cs_resample_ccw(const Polygon &poly, size_t N)
{
    Polygon p = poly;
    if (! p.is_counter_clockwise())
        p.reverse();
    const size_t m = p.points.size();
    std::vector<double> cum(m + 1, 0.0);
    for (size_t i = 0; i < m; ++i)
        cum[i + 1] = cum[i] + (p.points[(i + 1) % m] - p.points[i]).cast<double>().norm();
    const double tot = cum[m];
    std::vector<Vec2d> out(N);
    if (tot < 1e-9) {
        std::fill(out.begin(), out.end(), p.points.empty() ? Vec2d(0., 0.)
                                                           : p.points.front().cast<double>());
        return out;
    }
    size_t k = 0;
    for (size_t i = 0; i < N; ++i) {
        const double s = tot * double(i) / double(N);
        while (k + 1 < m && cum[k + 1] < s) ++k;
        const double seg = std::max(1e-12, cum[k + 1] - cum[k]);
        const double f   = std::clamp((s - cum[k]) / seg, 0.0, 1.0);
        const Vec2d a = p.points[k].cast<double>();
        const Vec2d b = p.points[(k + 1) % m].cast<double>();
        out[i] = a + f * (b - a);
    }
    return out;
}

// The two boundaries resampled to a common parametrisation, and the start offset that pairs
// them the way they face each other. Equal arc length is the correspondence; where it STARTS is
// free, and minimising the total squared distance is what makes two circles pair equal angles -
// which is what makes the blend come out circular.
bool CsRingField::morph_pair(std::vector<Vec2d> &A, std::vector<Vec2d> &B, size_t &best,
                             std::string *fail) const
{
    if (! m_valid || m_region.holes.size() != 1) {
        if (fail) *fail = "continuous shell ring field: a morph needs exactly one bore";
        return false;
    }
    const size_t N = std::max<size_t>(256, std::max(m_region.contour.points.size(),
                                                    m_region.holes.front().points.size()) * 2);
    A = cs_resample_ccw(m_region.contour, N);
    B = cs_resample_ccw(m_region.holes.front(), N);
    best = 0;
    double bestd = std::numeric_limits<double>::max();
    const size_t step = std::max<size_t>(1, N / 720);
    for (size_t d = 0; d < N; d += step) {
        double acc = 0.0;
        for (size_t i = 0; i < N; i += std::max<size_t>(1, N / 128))
            acc += (A[i] - B[(i + d) % N]).squaredNorm();
        if (acc < bestd) { bestd = acc; best = d; }
    }
    return true;
}

bool CsRingField::morph_family(size_t n, CsRingFamily &out, std::string *fail) const
{
    out.rings.clear();
    std::vector<Vec2d> A, B;
    size_t best = 0;
    if (! this->morph_pair(A, B, best, fail))
        return false;
    if (n == 0)
        return true;
    const size_t N = A.size();

    for (size_t j = 0; j < n; ++j) {
        const double t = double(j + 1) / double(n + 1);
        CsRing cr;
        cr.t = t;
        cr.poly.points.reserve(N);
        cr.extent.reserve(N);
        cr.bore.assign(N, 0);
        for (size_t i = 0; i < N; ++i) {
            const Vec2d q = (1.0 - t) * A[i] + t * B[(i + best) % N];
            cr.poly.points.emplace_back(coord_t(std::llround(q.x())), coord_t(std::llround(q.y())));
            cr.extent.push_back(float(this->extent_at(q)));
        }
        // A morph folds where the two boundaries are dissimilar. Refusing the whole family is
        // the point: the caller then uses the level sets, which never fold.
        if (union_ex(Polygons{ cr.poly }).size() != 1) {
            if (fail)
                *fail = "continuous shell ring field: the morph folds at t=" + std::to_string(t);
            out.rings.clear();
            return false;
        }
        out.rings.emplace_back(std::move(cr));
    }
    return true;
}

bool CsRingField::morph_spiral(size_t n, double lobes, double weave, double closure_mm,
                               double layer_phase, Polygon &out, std::string *fail) const
{
    out.points.clear();
    std::vector<Vec2d> A, B;
    size_t best = 0;
    if (! this->morph_pair(A, B, best, fail))
        return false;
    if (n < 3) {
        if (fail) *fail = "continuous shell ring field: too few rings for a spiral";
        return false;
    }
    const size_t N = A.size();
    const double denom = double(n + 1);

    // The point at RING INDEX v and azimuth index i, with the weave folded in as a modulation of
    // the index itself. v runs 0 at the outer boundary to n+1 at the bore, so one unit of v is
    // one ring gap - W(phi)/(n+1) millimetres, which VARIES around the wall.
    //
    // Deliberately not called a slot. In the marched path a slot is one bead spacing, a fixed
    // distance, and the entire reason for this file is that consecutive marched offsets are not
    // that distance apart. Reusing the word for a dimensionless index would put both meanings in
    // one codebase.
    //
    // Amplitude is half a ring gap less what a bead needs, so two arms one gap apart just touch;
    // the pi per unit of v in the phase is what puts them in anti-phase to do it.
    // The return type is spelled out on purpose. With `auto`, Eigen deduces the type of
    // `(1-t)*a + t*b` as an expression TEMPLATE holding references to `a` and `b` - which are
    // locals of this lambda and are destroyed on the way out. The caller then evaluates a
    // dangling expression, and every point comes back at the origin.
    auto at = [&](double v, size_t i) -> Vec2d {
        const Vec2d a = A[i % N], b = B[(i + best) % N];
        const double W = (a - b).norm() * SCALING_FACTOR;          // mm across the wall here
        const double gap = W / denom;
        // Half a ring gap less what a bead needs, and never more than 0.45 of a gap: at exactly
        // half, two anti-phase arms would meet at a point, and the clamp is what keeps a positive
        // distance between them everywhere so the spiral cannot close on itself.
        const double amp = (gap > 1e-9)
            ? std::min(0.45, weave * std::max(0.0, 0.5 - closure_mm / gap)) : 0.0;
        const double phi = double(i % N) / double(N);
        const double vv = v + amp * std::sin(2.0 * M_PI * lobes * phi + M_PI * v + layer_phase);
        const double t  = std::clamp(vv / denom, 0.0, 1.0);
        return (1.0 - t) * a + t * b;
    };

    // Arm A takes the odd ring indices, arm B the even ones, each gaining two per revolution, so
    // together they cover every ring once. Walk A inward, U-turn one ring at the bore, walk B
    // back out; its far end lands one ring from A's start, which closes the loop.
    // Each arm runs a WHOLE number of revolutions, so both finish at the same azimuth and the
    // U-turn between them is the one ring that separates them. Letting the arms simply run until
    // they reach the bore leaves them half a revolution apart, and the U-turn becomes a chord
    // straight across the middle - which the self-intersection guard duly caught.
    const size_t rev_a = (n >= 1) ? (n - 1) / 2 : 0;
    const size_t rev_b = (n >= 2) ? (n - 2) / 2 : 0;
    if (rev_a == 0 || rev_b == 0) {
        if (fail) *fail = "continuous shell ring field: too few rings for two spiral arms";
        return false;
    }
    std::vector<Vec2d> arm_a, arm_b;
    for (size_t k = 0; k <= rev_a * N; ++k)
        arm_a.push_back(at(1.0 + 2.0 * double(k) / double(N), k));
    for (size_t k = 0; k <= rev_b * N; ++k)
        arm_b.push_back(at(2.0 + 2.0 * double(k) / double(N), k));
    if (arm_a.size() < N || arm_b.size() < N) {
        if (fail) *fail = "continuous shell ring field: the spiral did not complete a turn";
        return false;
    }
    out.points.reserve(arm_a.size() + arm_b.size());
    for (const Vec2d &q : arm_a)
        out.points.emplace_back(coord_t(std::llround(q.x())), coord_t(std::llround(q.y())));
    for (size_t k = arm_b.size(); k-- > 0; )
        out.points.emplace_back(coord_t(std::llround(arm_b[k].x())),
                                coord_t(std::llround(arm_b[k].y())));
    // NO self-intersection test here, and that is deliberate. The arms are MEANT to touch -
    // that is what a weld is - and a polygon touching itself comes back from union_ex as several
    // pieces, so the test that guards the ring family rejects every spiral that is working. What
    // keeps the arms apart is arithmetic instead: they sit one ring apart in anti-phase, so
    // their separation is 1 - 2*amp ring gaps and amp is clamped below 0.45. The correspondence
    // underneath has already been checked by the ring family, which does fold detectably and is
    // computed first.
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
