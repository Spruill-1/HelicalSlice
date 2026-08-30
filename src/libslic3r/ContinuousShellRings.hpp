#pragma once

// Interpolated ring generation for the continuous shell.
//
// The rings of a continuous-shell interior used to come from repeated `offset_ex(band,
// -spacing)`. Offsets do not stay parallel: measured on real slices, consecutive
// undecimated offsets sit 1.00 bead slots apart at the median and 11.58 at the maximum, so
// "one slot" is not a distance. Every pathology the mode has fought - the medial pocket,
// the tip voids, the solid core, the gap cap that counts slots - is that one fact.
//
// Here the rings are level sets of the normalised distance field
//
//     u(p) = d_out(p) / (d_out(p) + d_in(p))
//
// where d_out is the distance to the region's outer boundary and d_in the distance to the
// NEAREST bore (or to the centroid, for a region with no bore). u is 0 on the outer
// boundary and 1 on a bore, so the ring at level t interpolates between them. Level sets of
// a continuous function never cross, always close, and are spaced W(p)/(n+1) apart where
// W = d_out + d_in is the local extent - evenly, by construction, in the fraction.
//
// Taking d_in to the NEAREST bore is also what subsections a multi-bore region: the locus
// where two bores are equidistant is a Voronoi cut, and because u is built from the bore
// distances themselves, u is CONTINUOUS across that cut. Neighbouring cells' level sets at
// the same t meet there and continue, so a cut is not an artificial face and the rings need
// not terminate at one. Each ring point carries the index of the bore nearest to it, which
// is the cell it belongs to.

#include "ExPolygon.hpp"
#include "Point.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace Slic3r {

// One level set of the field. `extent` and `bore` are parallel to `poly.points`.
struct CsRing
{
    Polygon            poly;
    double             t = 0.;      // the level this ring was extracted at, in (0, 1)
    std::vector<float> extent;      // W(p) = d_out + d_in at each point, mm
    std::vector<int>   bore;        // index of the nearest bore at each point, -1 if none
};

struct CsRingFamily
{
    std::vector<CsRing> rings;      // t ascending: outermost first, innermost last
};

// The normalised distance field over one peeled interior region, and the level sets of it.
//
// Two-phase by design: the field is what tells you how wide the region actually is, and the
// ring count depends on that, so `build` first, ask `extent_percentile` what n to use, then
// extract. Rebuilding the field per level would pay for the distance trees n times over.
class CsRingField
{
public:
    // `interior` is the region left after the surface skin has been peeled: its contour is
    // the outer boundary of the infill zone and its holes are the bores. Fails loudly
    // rather than returning something degenerate.
    // `sites` are VIRTUAL bores: extra sources for d_in, used to re-section a region whose ring
    // count will not fit everywhere. Adding one where the region is widest pulls d_in in there,
    // so the widest extent drops and the count with it - and because the field partitions by
    // NEAREST source, the new cut is an equidistant set exactly like a real bore's, so u stays
    // continuous across it and the level sets still meet and continue. Cutting the geometry
    // instead would put an artificial face there for the rings to terminate against.
    bool build(const ExPolygon &interior, double grid_pitch_mm, std::string *fail,
               const Points &sites = {});

    bool   valid() const { return m_valid; }
    size_t n_bores() const { return m_n_bores; }

    // Percentile of the local extent W over the interior, mm. `extent_percentile(0.75)` is
    // what picks the ring count: biasing high closes voids where the region is wide and
    // lets the thin places coalesce into the degradation ladder.
    double extent_percentile(double pct) const;

    // Where the region is DEEPEST - furthest from its outer wall. That is where a re-sectioning
    // site belongs. Placing one at the widest EXTENT instead is actively harmful and was measured
    // so: the extent at a site is just d_out there, so a site dropped at the argmax of
    // d_out + d_in lands somewhere d_out is small, the extent collapses to nearly nothing around
    // it, and the spread the site was meant to tighten explodes - 3.3x to 231x on an eccentric
    // annulus, taking the collapsed fraction from 9.3% to 13.2% with ONE site. At the deepest
    // point the extent is locally at its largest, so no such singularity is created.
    Point deepest_point() const { return m_deepest; }

    // Ring count such that NO gap anywhere exceeds max_gap_mm. n rings leave n+1 gaps across
    // the extent and the gap is W(p)/(n+1), so the widest place in the region is what sets the
    // count: anywhere thinner is crowded, which is the degradation ladder's business, not a
    // reason to let the widest place open past what a weave can bridge. 0 - no room for a
    // buried ring at all - is a legitimate answer, not a failure.
    // `pct` picks which extent the count answers to: 1.0, the default, is the widest place, so
    // no gap anywhere exceeds max_gap. Lowering it answers to a typical extent instead - the
    // thin parts then stop crowding, and the wide parts open past the weld distance in exchange.
    // That trade is the one-count alternative to subsectioning, and it is a trade, not a fix.
    size_t ring_count_for_max_gap(double max_gap_mm, double pct = 1.0) const;

    // The local extent W(p) in mm, by bilinear lookup of the gridded field - O(1), and good to
    // a few percent, which is far finer than the 0.05 mm buckets a bead width is quantised to.
    // This is what an adaptive bead width is a function of, and it must be a FUNCTION of the
    // point rather than an array parallel to the rings: the ring assembly resamples and reorders
    // points, and an index-parallel width array cannot survive that.
    double extent_at(const Vec2d &p) const;

    // Fraction of the region whose local extent is below `w_mm`. What decides whether a ring
    // count needs the degradation ladder is not whether the narrowest corner collapses - some
    // corner nearly always does - but how much of the wall it takes with it.
    double fraction_narrower_than(double w_mm) const;

    // One level can come back as several loops: past the saddle where two bores separate,
    // the level set of the field genuinely is one loop per bore.
    bool level_set(double t, std::vector<CsRing> &out, std::string *fail) const;
    bool family(size_t n, CsRingFamily &out, std::string *fail) const;

    // Exposed for tests and for the callers that need the field itself: the value of u, the
    // local extent W, and the nearest-bore index at an arbitrary point (scaled coords).
    // `u` is 0 outside the region and 1 inside a bore, so the field is continuous across
    // both boundaries.
    void sample(const Vec2d &p, double &u, double &extent, int &bore) const;

    // Same, for a point already known to be inside the region - the field build knows
    // that by scanline and must not pay for the ray cast again.
    void sample_inside(const Vec2d &p, double &u, double &extent, int &bore) const;

private:
    // Newton-correct p onto the level set u = t with the closed-form gradient. False means
    // the point sits where the field has no usable gradient (on a boundary, or where both
    // boundaries pull the same way) and was left as close as it could be put.
    // Also reports the local extent (mm) and the nearest bore at the corrected point.
    bool refine(Vec2d &p, double t, double &extent, int &bore) const;

    // Gradient-free corrector, used wherever the field has a kink - along a Voronoi cut
    // above all, where the nearest bore swaps and u is continuous but not smooth.
    bool bisect(Vec2d &p, double t, double &extent, int &bore) const;

    bool                m_valid = false;
    double              m_pitch = 0.;      // scaled
    Vec2d               m_origin{0., 0.};  // scaled, position of raster cell (0, 0)
    size_t              m_rows = 0, m_cols = 0;
    std::vector<float>  m_extent;          // exact W on a subsample, mm - what the percentiles use
    std::vector<float>  m_wgrid;           // gridded W over the whole raster, mm - what extent_at reads
    ExPolygon           m_region;
    size_t              m_n_bores = 0;
    std::vector<size_t> m_bore_of_line;    // line index in the bore distancer -> bore index
    Vec2d               m_centroid{0., 0.};// inner target when there is no bore or site, scaled
    Points              m_sites;           // virtual bores, scaled
    Point               m_deepest{0, 0};  // argmax of d_out: furthest from the outer wall

    struct Impl;                            // holds the distance trees
    std::shared_ptr<Impl> m_impl;
};

} // namespace Slic3r
