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

    // Close the family onto the medial axis AS WELL as the bores: d_in becomes the minimum
    // of the two sources, so the extent is the local half-thickness even where a region's
    // thin parts sit far from its bores. u stays continuous - min over sources is the same
    // property the virtual sites rest on. Set BEFORE build, and the flag survives rebuilds.
    // The caller keeps this field only when the measured collapse share improves: on a
    // region whose bores already carry the thickness, the closure halves every wall-to-bore
    // corridor instead of relieving anything, so it has to earn its keep per region.
    void set_close_medial(bool v) { m_close_medial = v; }
    bool close_medial() const { return m_close_medial; }

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

    // The AZIMUTH field: the coordinate that runs around the region while u runs across it.
    // Seeded on the outer boundary with its own normalised arc length and carried inward along
    // the gradient flow of u - which is the direction the nearest-point correspondence runs, so
    // two rings at the same phi are the two points that actually touch.
    //
    // This exists because the alternatives were measured and both fail. Inheriting phase from
    // the ring outside by nearest point gives a staircase with half its steps running backwards
    // by the second ring. Giving each ring its own arc-length ramp is perfectly smooth and
    // drifts up to 0.79 of a lobe by the ninth, which turns anti-phase into in-phase. The
    // transported field holds correspondence to about 0.01 of a lobe at the median.
    //
    // Returns a fraction in [0, 1). Zero-length regions and points off the field return 0.
    double phi_at(const Vec2d &p) const;

    // Fraction of the region whose local extent is below `w_mm`. What decides whether a ring
    // count needs the degradation ladder is not whether the narrowest corner collapses - some
    // corner nearly always does - but how much of the wall it takes with it.
    double fraction_narrower_than(double w_mm) const;

    // One level can come back as several loops: past the saddle where two bores separate,
    // the level set of the field genuinely is one loop per bore.
    bool level_set(double t, std::vector<CsRing> &out, std::string *fail) const;
    bool family(size_t n, CsRingFamily &out, std::string *fail) const;

    // SPATIALLY-VARYING family: the uniform count answers to the WIDEST place, which puts the
    // inter-bore saddle - where the level sets spread as |grad u|->0 - right at the weld limit,
    // the strip the plate peeled along. This starts from the n uniform levels and, wherever two
    // consecutive level sets sit physically farther apart than `saddle_gap_mm` (measured as the
    // max distance from one to the next, which is largest exactly at the saddle), splices in
    // enough intermediate levels to bring that corridor under the target. The extra levels are
    // full level sets of the same field - they thread the wide corridor and only crowd elsewhere,
    // where the weave and adaptive width already absorb crowding - so the saddle densifies without
    // cutting the section or adding a stitch. `cap` bounds the total ring count against a
    // pathological field. `max_extra` budgets how many intermediate levels may be added, spent
    // WIDEST-corridor-first (0 = unlimited: bring every corridor under the target); a small budget
    // closes only the worst saddle(s), a light densification that keeps a multi-saddle section from
    // printing solid. saddle_gap_mm <= 0 falls back to the uniform family.
    bool family_saddle(size_t n, double saddle_gap_mm, size_t cap, size_t max_extra,
                       CsRingFamily &out, std::string *fail) const;

    // DISTANCE-anchored rings (arc-capable work, phase 1): rings at a fixed DISTANCE from each
    // boundary - offsets at (k+1/2)*spacing - rather than a fixed FRACTION of the local extent.
    // Where the wall is thick the offsets nest as usual; where it thins, the deep offsets simply
    // do not exist, so the ring count follows the local thickness and the thin flank stops
    // crowding. Wall offsets (each ExPolygon's contour) and bore offsets (its holes) are both
    // returned, meeting at the medial. `t` is a synthetic monotone level (offset depth normalised)
    // only so the downstream group-by-level still orders them outer->inner.
    bool distance_family(coord_t spacing, CsRingFamily &out, std::string *fail) const;

    // Rings by MORPH: the point at fraction t is the straight blend of the two boundary points
    // at equal normalised arc length. Between two circles this returns circles, which the level
    // sets do not - they return Cartesian ovals, off by up to 0.357 mm on a 7.5 mm wall, and
    // they put a genuine 17-degree corner in whichever ring crosses the outer boundary's medial
    // axis, because d_out has a local maximum there and level sets are singular at critical
    // points. A morph has no critical points and no such corner.
    //
    // Needs a correspondence, which is why it is not the general answer: where the two
    // boundaries are dissimilar the blend can fold. Every ring is checked and the whole family
    // is refused if any of them does, so the caller falls back to the level sets rather than
    // printing a folded ring.
    //
    // Single-bore regions only - a bore-less one has no second curve to morph to.
    bool morph_family(size_t n, CsRingFamily &out, std::string *fail) const;

    // The same family as ONE double-arm Fermat spiral, already woven.
    //
    // A layer of nested closed rings needs the nozzle to step inward once per ring, and every
    // one of those steps is a near-right-angle radial hop: measured, 19 to 23 of them per layer
    // against about 21 loops, and they are what is left in the picture once the weave is gentle.
    // A spiral has two - the U-turn at the bore and the join back to its own start - because it
    // never stops going round.
    //
    // In morph coordinates this is almost free. A ring is (1-t)A + tB at fixed t; the spiral is
    // the same expression with t advancing continuously, gaining one ring per half revolution on
    // each arm. The weave is then a modulation of t rather than a displacement along a normal:
    // amplitude is half a RING GAP less the bead allowance, and the phase carries a pi per ring,
    // so an arm and its neighbour one ring away are anti-phase and their excursions kiss.
    //
    // "Ring gap" and not "slot" on purpose: in the marched path a slot is one bead spacing, a
    // fixed distance, and this whole file exists because consecutive marched offsets are not
    // that distance apart. A ring gap is W(phi)/(n+1) and varies around the wall.
    //
    // `weave` scales the amplitude, `lobes` sets its frequency, `closure_mm` is how close two
    // beads' centres may come. Returns one CLOSED polygon: a double-arm spiral joined at its
    // ends is a single loop, which is what lets it travel the rest of the pipeline as a ring.
    bool morph_spiral(size_t n, double lobes, double weave, double closure_mm,
                      double layer_phase, Polygon &out, std::string *fail) const;

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

    void build_phi();

    bool morph_pair(std::vector<Vec2d> &A, std::vector<Vec2d> &B, size_t &best,
                    std::string *fail) const;

    bool                m_valid = false;
    double              m_pitch = 0.;      // scaled
    Vec2d               m_origin{0., 0.};  // scaled, position of raster cell (0, 0)
    size_t              m_rows = 0, m_cols = 0;
    std::vector<float>  m_extent;          // exact W on a subsample, mm - what the percentiles use
    std::vector<float>  m_wgrid;           // gridded W over the whole raster, mm - what extent_at reads
    std::vector<float>  m_phi_c, m_phi_s;  // the azimuth field as cos/sin, so its seam averages
    ExPolygon           m_region;
    size_t              m_n_bores = 0;
    bool                m_close_medial = false; // d_in = min(bores, medial axis); see setter
    std::vector<size_t> m_bore_of_line;    // line index in the bore distancer -> bore index
    Vec2d               m_centroid{0., 0.};// inner target when there is no bore or site, scaled
    Points              m_sites;           // virtual bores, scaled
    Point               m_deepest{0, 0};  // argmax of d_out: furthest from the outer wall

    struct Impl;                            // holds the distance trees
    std::shared_ptr<Impl> m_impl;
};

} // namespace Slic3r
