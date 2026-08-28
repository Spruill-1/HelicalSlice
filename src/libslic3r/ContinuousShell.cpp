#include "ContinuousShell.hpp"
#include <cstdio>
#include <cstdlib>

#include "Arachne/BeadingStrategy/BeadingStrategyFactory.hpp"
#include "ClipperUtils.hpp"
#include "KDTreeIndirect.hpp"
#include "AABBTreeLines.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace Slic3r {

namespace {

// Sparse WOVEN infill for a SOLID continuous-shell cross-section. Past the solid outer shell
// (wall_loops rings) the fill loops are thinned to a sparse spacing and (optionally) perturbed
// into a 3D interlocking weave - all carried by the SAME continuous bead. Driven by config:
//   continuous_shell_infill_density  0 => solid (disabled); else target interior density
//   continuous_shell_infill_weave    0 => plain concentric; 100 => full radial+Z weave
//   wall_loops                       => solid shell ring count
struct WovenInfill {
    bool   enabled     = false;
    size_t shell_loops = 1;
    size_t stride      = 1;         // keep every Nth interior loop (1 = solid)
    double lobes       = 6.0;       // angular lobe count of the weave
    double weave       = 0.0;       // 0 = plain concentric, 1 = full weave (scales amp/zamp)
    double dphase      = 0.9;       // per-layer phase advance (rad)
    double zamp        = 0.0;       // Z amplitude, fraction of layer_height
    double zphase      = 1.5708;    // quadrature offset (rad) between the radial and Z waves
    bool   serpentine  = false;     // radial-traverse interior instead of concentric rings
    // Where this layer must stay SOLID, when that is only part of the cross-section.
    //
    // force_solid is a whole-layer boolean: if any part of the section has air within the shell
    // count, the ENTIRE section prints solid. For a genuine top or bottom shell that is right -
    // the whole face is exposed. For an INTERNAL surface it is not: the floor of a blind bore
    // needs solid over the bore's footprint, and the wall around it is buried and could stay
    // sparse. That layer currently prints solid edge to edge.
    //
    // Empty means "no mask" - either the whole section is solid (force_solid) or none of it is.
    const ExPolygons *solid_mask = nullptr;
    double solid_stagger = 0.0;     // radial inward shift (scaled) of a SOLID cap's interior rings this layer
};
WovenInfill woven_infill_params(const PrintRegionConfig *config)
{
    WovenInfill w;
    if (config == nullptr)
        return w;
    const double density = config->continuous_shell_infill_density.value / 100.0;
    if (density <= 0.0 || density >= 1.0)
        return w;                  // 0% or >=100% => keep the interior fully solid
    w.enabled     = true;
    w.stride      = std::max<size_t>(2, size_t(std::llround(1.0 / std::clamp(density, 0.05, 0.99))));
    w.shell_loops = std::max<size_t>(1, size_t(config->wall_loops.value));
    // Weave FREQUENCY knob: it sets the lobe count (spatial frequency) of the weave, NOT its amplitude.
    // A partial-amplitude weave barely interconnects (pointless), so the weave always uses full inter-
    // lock (anti-phase lobes just touch). 100% = 6 lobes (baseline); scales linearly; 0% = no weave
    // (plain concentric sparse rings). Higher = more, smaller lobes (finer lattice).
    const double freq = config->continuous_shell_infill_weave.value;   // percent; 0 turns the weave off
    if (freq > 0.0) {
        w.weave = 1.0;                                                  // full-amplitude interlock (lobes touch)
        w.lobes = std::clamp(std::round(6.0 * freq / 100.0), 2.0, 32.0);
        w.zamp  = 0.45;                                                 // < 0.5 so the Z excursion stays within ~half a layer
        // Per-layer phase shift = how far the weave rotates between layers (the inter-layer weave). 100%
        // is a half-lobe shift (pi) - a sharp two-layer criss-cross; lower is a gentler helical weave.
        w.dphase = std::clamp(config->continuous_shell_infill_weave_shift.value / 100.0, 0.0, 2.0) * M_PI;
    } else {
        w.weave = 0.0;                                                  // plain concentric sparse rings
        w.zamp  = 0.0;
    }
    w.serpentine = config->continuous_shell_interior_serpentine.value;
    return w;
}

ExtrusionPath make_path(const Polyline3 &polyline, const Flow &flow, double layer_height)
{
    ExtrusionPath path(erExternalPerimeter);
    path.polyline = polyline;
    path.mm3_per_mm = flow.mm3_per_mm();
    path.width = flow.width();
    path.height = float(layer_height);
    return path;
}

// Resample a closed polygon to exactly N points spaced evenly by arc length. Even
// spacing keeps the seam anchoring, the weave phase, and the arc-length Z ramp well
// conditioned (no point-density bias). N is chosen PER RING by the caller: these
// points are the final toolpath (z_contoured beads are exempt from downstream
// simplification), so surface rings are sampled finely and buried rings coarsely.
std::vector<Point> resample_closed(const Polygon &poly, size_t N)
{
    std::vector<Point> out;
    const size_t n = poly.points.size();
    if (n < 3 || N == 0)
        return out;

    std::vector<double> cum(n + 1, 0.);
    for (size_t i = 0; i < n; ++i)
        cum[i + 1] = cum[i] + (poly.points[(i + 1) % n] - poly.points[i]).cast<double>().norm();
    const double total = cum[n];
    if (total <= 0.)
        return out;

    out.reserve(N);
    size_t seg = 0;
    for (size_t k = 0; k < N; ++k) {
        const double target = (double(k) / double(N)) * total;
        while (seg + 1 < n && cum[seg + 1] < target)
            ++seg;
        const double seg_len = cum[seg + 1] - cum[seg];
        const double u = seg_len > 0. ? (target - cum[seg]) / seg_len : 0.;
        const Point &a = poly.points[seg];
        const Point &b = poly.points[(seg + 1) % n];
        out.emplace_back(a.x() + coord_t(std::llround((b.x() - a.x()) * u)),
                         a.y() + coord_t(std::llround((b.y() - a.y()) * u)));
    }
    return out;
}

// The woven infill, the solid stagger, and the Z-wave all parametrize each ring by the
// ANGLE of its points about one centre - the perturbation phase is sin(lobes*theta), and
// the shell-clamp profiles bin "the shell radius at angle theta". That is only meaningful
// when every ring is star-shaped about that centre: each ray from the centre crosses the
// ring exactly once. On a folded section (a crescent, an L whose centroid sits outside
// the material) a bin averages multiple crossings, and clamping to that fiction can drag
// a bead half the section width off its ring (colliding with its neighbours). Detect the
// fold: traversing the ring, theta must advance monotonically through one full turn -
// whichever orientation - with no backtracking beyond discretization noise.
bool star_shaped_about(const Polygon &poly, const Vec2d &centre)
{
    const size_t n = poly.points.size();
    if (n < 3)
        return false;
    Vec2d d0 = poly.points.front().cast<double>() - centre;
    if (d0.norm() < 1.0)
        return false;                       // ring passes (numerically) through the centre
    double pos = 0., neg = 0.;
    double prev = std::atan2(d0.y(), d0.x());
    for (size_t i = 1; i <= n; ++i) {
        const Vec2d d = poly.points[i % n].cast<double>() - centre;
        if (d.norm() < 1.0)
            return false;
        const double th = std::atan2(d.y(), d.x());
        double inc = th - prev;
        inc = std::atan2(std::sin(inc), std::cos(inc));   // wrap to (-pi, pi]
        (inc >= 0. ? pos : neg) += inc;
        prev = th;
    }
    constexpr double tol = 0.1;             // rad of slack for discretization noise
    const double winding = pos + neg;
    if (std::abs(winding - 2. * M_PI) < tol)
        return -neg < tol;                  // CCW ring: negligible backtracking
    if (std::abs(winding + 2. * M_PI) < tol)
        return pos < tol;                   // CW ring (a hole contour)
    return false;                           // centre outside the ring, or a multi-fold
}

// Per-angle mean-radius profile of a ring about a centre: NB angular bins, empty bins
// forward-filled over a doubled sweep. One shared implementation for every angular
// parametrization. The woven-infill shell clamps and the solid-cap stagger use
// sample_nearest (the historical nearest-bin lookup - their output must not change);
// the per-angle bead plan (C-full) uses sample_linear, which interpolates between bin
// centres so ring geometry stays C0 across bin boundaries.
struct AngularProfile {
    std::vector<double> r;

    static AngularProfile build(const Polygon &ring, const Vec2d &centre, int NB = 360)
    {
        AngularProfile ap;
        ap.r.assign(NB, 0.0);
        std::vector<int> cnt(NB, 0);
        for (const Point &q : ring.points) {
            const Vec2d d = q.cast<double>() - centre;
            int b = int((std::atan2(d.y(), d.x()) + M_PI) / (2.0 * M_PI) * NB);
            b = ((b % NB) + NB) % NB;
            ap.r[b] += d.norm();
            ++cnt[b];
        }
        for (int b = 0; b < NB; ++b)
            if (cnt[b]) ap.r[b] /= cnt[b];
        double last = 0.0;                        // fill empty bins from the last seen
        for (int b = 0; b < 2 * NB; ++b) {
            const int i = b % NB;
            if (ap.r[i] > 0.0) last = ap.r[i]; else if (last > 0.0) ap.r[i] = last;
        }
        return ap;
    }

    double sample_nearest(double th) const
    {
        const int NB = int(r.size());
        int b = int((th + M_PI) / (2.0 * M_PI) * NB);
        b = ((b % NB) + NB) % NB;
        return r[b];
    }

    double sample_linear(double th) const
    {
        const int NB = int(r.size());
        const double u  = (th + M_PI) / (2.0 * M_PI) * NB - 0.5;   // bin-centre coordinates
        const double fl = std::floor(u);
        const double frac = u - fl;
        int b0 = int(fl);
        b0 = ((b0 % NB) + NB) % NB;
        const int b1 = (b0 + 1) % NB;
        return r[b0] * (1.0 - frac) + r[b1] * frac;
    }
};

// ---- Per-angle bead plan (adaptive beading, C-full) ---------------------------------------------
// For an annular star-shaped chain the ring count and per-ring widths are PLANNED from the
// measured wall thickness per angle, via Arachne's pure/stateless BeadingStrategy stack,
// instead of emerging from fixed-spacing band marching. Rings become iso-index curves of
// the plan; per-point widths are evaluated as a function of angle at emission time.

// Count law (squeeze bias). The printed T4 taper showed the failure regimes are asymmetric:
// SQUEEZED beads (narrower than nominal) consolidate with only pinprick corner voids, while
// STRETCHED beads develop a periodic void line along the medial. So the ring count is not
// chosen at the equidistant thickness: an inner bead may never stretch past (1+eps) of the
// nominal spacing - the count flips up early and the wall spends the transition band
// squeezed. eps = 0.10 caps stretch at 10% (the T4 void regime began at ~29%) while a wall
// just past an integer tiling (t = 4.05 spacings) still takes a 1% stretch instead of a
// gratuitous 19% squeeze. Stateless (pure in t): no hysteresis, deterministic flip layers.
// The two outermost beads are pinned at nominal by RedistributeBeadingStrategy, so the law
// applies to the INNER region; counts <= 2 delegate to the stack (thin-wall semantics).
// The added bead is TRANSITIONAL while narrow (it is being born; the structural
// min_bead_width does not bind it, same philosophy as the scarf floor), so it is only
// rejected below the termination floor w_floor_sp - otherwise the count law would be
// unable to flip before the stretch reached the void regime it exists to prevent.
coord_t squeeze_biased_count(coord_t t_sp, const Arachne::BeadingStrategy &strategy,
                             coord_t w_opt_sp, coord_t w_floor_sp)
{
    const coord_t opt = strategy.getOptimalBeadCount(t_sp);
    if (opt <= 2)
        return opt;
    constexpr double eps = 0.10;                 // max inner-bead stretch before adding a ring
    const coord_t inner = t_sp - 2 * w_opt_sp;
    coord_t n = 2 + coord_t(std::ceil(double(inner) / ((1.0 + eps) * double(w_opt_sp))));
    n = std::clamp(n, opt, opt + 1);
    while (n > opt) {
        const Arachne::BeadingStrategy::Beading b = strategy.compute(t_sp, n);
        coord_t narrowest = std::numeric_limits<coord_t>::max();
        for (const coord_t w : b.bead_widths)
            if (w > 0)
                narrowest = std::min(narrowest, w);
        if (b.bead_widths.size() == size_t(n) && narrowest >= w_floor_sp)
            break;
        --n;                                     // the extra ring would be unprintable even as a tip
    }
    return n;
}

// Blend two beadings of the SAME thickness across a ring-count step k-1 -> k, with
// front/back index pairing: outer slots pair from index 0, bore slots pair from the
// back, and the unpaired middle slot of the k-bead branch is the appearing bead, its
// width scaled by the blend factor rho (0 = k-1 side, 1 = k side). This is what
// SkeletalTrapezoidation::interpolate CANNOT do for a ring topology: its index-wise
// pairing would ramp the BORE bead instead of the medial one (un-pinning the bore
// surface) and it drops the appearing bead entirely at equal thickness. Both inputs
// sum to the same thickness, so any convex blend conserves the covered width exactly
// (minus the sub-floor sliver at the appearing bead's tips).
Arachne::BeadingStrategy::Beading blend_ring_beadings(const Arachne::BeadingStrategy::Beading &lo,
                                                      const Arachne::BeadingStrategy::Beading &hi,
                                                      double rho, size_t mid)
{
    Arachne::BeadingStrategy::Beading out = hi;
    const size_t k = hi.bead_widths.size();
    for (size_t i = 0; i < k; ++i) {
        if (i == mid) {
            out.bead_widths[i] = coord_t(std::llround(rho * double(hi.bead_widths[i])));
            continue;                              // location stays at the k-branch value
        }
        const size_t j = i < mid ? i : i - 1;      // paired slot in the k-1 branch
        out.bead_widths[i] = coord_t(std::llround((1.0 - rho) * double(lo.bead_widths[j]) + rho * double(hi.bead_widths[i])));
        out.toolpath_locations[i] = coord_t(std::llround((1.0 - rho) * double(lo.toolpath_locations[j]) + rho * double(hi.toolpath_locations[i])));
    }
    return out;
}

// The per-chain plan: angular radius profiles of the region's two boundaries plus, per
// angular bin and radial slot, the planned extrusion width (mm, 0.05-bucketed; 0 = the
// slot is absent at that angle) and the radial placement fraction of the local extent.
// Slot 0 is the outer surface, slot n_slots-1 the bore surface. Queries interpolate
// between bin centres so ring geometry and widths are C0 around the circumference.
struct ChainBeadPlan {
    Vec2d centre{0., 0.};
    int   NB = 720;
    AngularProfile r_out, r_bore;                 // radial profiles, scaled units
    double edge_margin = 0.;                      // (scaled_width - scaled_spacing)/2
    size_t n_slots = 0;
    std::vector<std::vector<float>>  slot_width;  // [NB][n_slots], extrusion mm
    std::vector<std::vector<double>> slot_frac;   // [NB][n_slots], fraction of radial extent
    // Interior spiral (phase B): when set, the buried rings print as one continuous
    // migrating arm whose comb has the SAME bead count and pitch as ring mode (no
    // stretch; the squeeze doctrine holds). The arm geometry - tangent merges onto the
    // surface rings, merge-zone clearance widths - is evaluated analytically at
    // generation time in spiralize_rings from the ordinary slot tables via the
    // fractional-slot queries below. The arm phase anchors at spiral_theta0 (this
    // layer's seam angle). ramp_frac is the fraction of each revolution spent migrating
    // between columns: 1 = a pure helical thread, lower values hold each column at a
    // constant radius and concentrate the migration into one short ramp (see
    // spiralize_rings for why this is the ring<->spiral continuum).
    bool   spiral = false;
    double spiral_theta0 = 0.;
    double ramp_frac = 1.0;
    float  w_term = 0.30f;                            // transitional/merge width floor, mm

    // Bin-centre linear interpolation with wraparound; shared by the queries below.
    template<typename Table, typename Get>
    double interp(const Table &tab, size_t slot, double th, Get get) const
    {
        const double u  = (th + M_PI) / (2.0 * M_PI) * NB - 0.5;
        const double fl = std::floor(u);
        const double f  = u - fl;
        int b0 = int(fl);
        b0 = ((b0 % NB) + NB) % NB;
        const int b1 = (b0 + 1) % NB;
        return get(tab[b0][slot]) * (1.0 - f) + get(tab[b1][slot]) * f;
    }
    float width_at(size_t slot, double th) const
    {
        const double w = interp(slot_width, slot, th, [](float v) { return double(v); });
        return std::round(float(w) / 0.05f) * 0.05f;          // Arachne's width buckets
    }
    double frac_at(size_t slot, double th) const
    {
        return interp(slot_frac, slot, th, [](double v) { return v; });
    }
    // Slot-table queries at a FRACTIONAL slot coordinate (the pairwise ring morph):
    // lerp across adjacent slot columns as well as across angular bins.
    double frac_at_f(double slot_coord, double th) const
    {
        const double fl = std::floor(slot_coord);
        const double f  = slot_coord - fl;
        const size_t j0 = size_t(std::clamp<double>(fl, 0., double(n_slots - 1)));
        const size_t j1 = std::min<size_t>(j0 + 1, n_slots - 1);
        const double a  = interp(slot_frac, j0, th, [](double v) { return v; });
        const double b  = interp(slot_frac, j1, th, [](double v) { return v; });
        return a * (1.0 - f) + b * f;
    }
    double width_at_f(double slot_coord, double th) const
    {
        const double fl = std::floor(slot_coord);
        const double f  = slot_coord - fl;
        const size_t j0 = size_t(std::clamp<double>(fl, 0., double(n_slots - 1)));
        const size_t j1 = std::min<size_t>(j0 + 1, n_slots - 1);
        const double a  = interp(slot_width, j0, th, [](float v) { return double(v); });
        const double b  = interp(slot_width, j1, th, [](float v) { return double(v); });
        return a * (1.0 - f) + b * f;
    }
};

// Build the plan for one annular chain, or decline (nullptr + decline_note) when the
// section is out of the plan's scope so the caller keeps the fixed-width C-lite rings.
// In-scope inconsistencies set *hard_fail instead: the layer must fail loudly rather
// than silently print a wall the plan itself says is wrong.
std::unique_ptr<ChainBeadPlan> build_chain_plan(const ExPolygon &region, const Flow &flow,
                                                const ContinuousShellWidthLimits &limits,
                                                double layer_height,
                                                bool spiral_mode, double junction_ramp, double seam_theta,
                                                std::string *decline_note, std::string *hard_fail)
{
    const coord_t w_sp  = flow.scaled_spacing();
    const coord_t w_ext = flow.scaled_width();
    const double  min_bead_mm    = limits.min_bead_width_mm   > 0. ? limits.min_bead_width_mm   : 0.85 * flow.nozzle_diameter();
    const double  min_feature_mm = limits.min_feature_size_mm > 0. ? limits.min_feature_size_mm : 0.25 * flow.nozzle_diameter();
    const double  trans_len_mm   = limits.wall_transition_length_mm > 0. ? limits.wall_transition_length_mm : 1.0 * flow.nozzle_diameter();
    (void) layer_height;

    auto plan = std::make_unique<ChainBeadPlan>();
    const int NB = plan->NB;

    // Plan centre = centroid of the region's outer contour; every angular parametrization
    // below requires both boundaries star-shaped about it (same gate as the weave).
    Vec2d centre(0., 0.);
    if (region.contour.points.empty() || region.holes.size() != 1)
        { if (decline_note) *decline_note = "not a simple annulus"; return nullptr; }
    for (const Point &p : region.contour.points)
        centre += p.cast<double>();
    centre /= double(region.contour.points.size());
    if (!star_shaped_about(region.contour, centre) || !star_shaped_about(region.holes.front(), centre))
        { if (decline_note) *decline_note = "section is not star-shaped"; return nullptr; }
    plan->centre = centre;

    // Radial profiles from densified copies of the boundaries (4 points per bin, so a
    // valid star section can never leave a bin empty).
    auto densify = [&](const Polygon &poly) -> Polygon {
        std::vector<Point> pts = resample_closed(poly, size_t(4 * NB));
        return pts.empty() ? poly : Polygon(Points(pts.begin(), pts.end()));
    };
    const Polygon outer_d = densify(region.contour);
    const Polygon bore_d  = densify(region.holes.front());
    plan->r_out  = AngularProfile::build(outer_d, centre, NB);
    plan->r_bore = AngularProfile::build(bore_d, centre, NB);
    for (int b = 0; b < NB; ++b)
        if (!(plan->r_out.r[b] > 0.) || !(plan->r_bore.r[b] > 0.) || plan->r_out.r[b] <= plan->r_bore.r[b])
            { if (decline_note) *decline_note = "degenerate radial profile"; return nullptr; }

    // Wall thickness per bin, in the strategy's spacing domain. The thickness is the
    // Euclidean nearest distance from the outer-boundary point to the bore boundary - a
    // pure radial difference overstates it by 1/cos(alpha) where a ray meets the wall
    // obliquely (+41% on a square annulus diagonal) and would plan phantom beads. Minus
    // both half edge margins: beads are placed with the outer bead's centre a half
    // EXTRUSION width inside each surface, but widths add up in the spacing domain.
    plan->edge_margin = 0.5 * double(w_ext - w_sp);
    std::vector<double> t_sp(NB, 0.);
    double t_sp_max = 0.;
    for (int b = 0; b < NB; ++b) {
        const double th = -M_PI + (double(b) + 0.5) * (2.0 * M_PI / NB);
        const Vec2d  p  = centre + plan->r_out.r[b] * Vec2d(std::cos(th), std::sin(th));
        double d2 = std::numeric_limits<double>::max();
        for (const Point &q : bore_d.points)
            d2 = std::min(d2, (q.cast<double>() - p).squaredNorm());
        t_sp[b] = std::sqrt(d2) - 2.0 * plan->edge_margin;
        t_sp_max = std::max(t_sp_max, t_sp[b]);
        if (t_sp[b] < scale_(min_feature_mm))
            { if (decline_note) *decline_note = "wall below minimum feature size"; return nullptr; }
    }

    // The standalone strategy stack, mirroring WallToolPaths' construction: spacings as
    // the preferred widths, split/add thresholds derived from the EXTRUSION width, a cap
    // far above any real count (so LimitedBeadingStrategy's zero-width marker beads never
    // appear). The interior distribution radius is effectively unbounded: the two surface
    // beads stay pinned at nominal (RedistributeBeadingStrategy) and ALL interior beads
    // share the squeeze uniformly - the consolidation pattern the printed T4 sections
    // showed - instead of dumping it onto the medial bead alone.
    const double split_thr = std::clamp(2. * min_bead_mm / double(flow.width()) - 1., 0.01, 0.99);
    const double add_thr   = std::clamp(min_bead_mm / double(flow.width()), 0.01, 0.99);
    // Cap comfortably above any thickness the plan actually measured, so
    // LimitedBeadingStrategy's marker-bead machinery never engages even on
    // count-law probes at optimal+1.
    coord_t max_beads = coord_t(std::llround(t_sp_max / double(w_sp))) + 8;
    if (max_beads & 1)
        ++max_beads;
    const Arachne::BeadingStrategyPtr strategy = Arachne::BeadingStrategyFactory::makeStrategy(
        w_sp, w_sp,
        coord_t(std::llround(scale_(trans_len_mm))),
        float(M_PI / 4.0) /* inert inside the strategies */,
        /*print_thin_walls=*/true,
        coord_t(std::llround(scale_(min_bead_mm))),
        coord_t(std::llround(scale_(min_feature_mm))),
        split_thr, add_thr, max_beads,
        /*outer_wall_offset=*/0,
        /*inward_distributed_center_wall_count=*/9999);

    // Termination floor for transitional (appearing/retiring) beads, and its spacing-
    // domain equivalent for the count law. Below ~1.2x layer height a TPU bead
    // die-swells and gaps (the printed T5 scarf verdict), so no tip may dip under it.
    const float   w_term      = std::max(0.30f, 0.6f * float(min_bead_mm));
    const coord_t w_floor_sp  = std::max<coord_t>(1, coord_t(std::llround(scale_(double(w_term)))) - (w_ext - w_sp));

    // Per-bin bead count.
    std::vector<coord_t> n_bin(NB, 0);
    coord_t n_min = std::numeric_limits<coord_t>::max(), n_max = 0;
    for (int b = 0; b < NB; ++b) {
        n_bin[b] = squeeze_biased_count(coord_t(std::llround(t_sp[b])), *strategy, w_sp, w_floor_sp);
        n_min = std::min(n_min, n_bin[b]);
        n_max = std::max(n_max, n_bin[b]);
    }
    if (n_min < 1)
        { if (decline_note) *decline_note = "wall too thin for a planned bead"; return nullptr; }
    if (n_max - n_min > 1) {
        // The count law picks a count per bin, so a wall whose thickness merely DRIFTS can be
        // assigned counts two apart even where a single count would serve the whole section
        // comfortably. Declining those sends the wall to marching, which cannot fill a wall of
        // varying thickness: it lays rings at fixed spacing from each face, so the leftover at
        // the medial is a function of azimuth and no scalar correction can close it. Measured
        // on a 5.5/6.5 mm wall (a 0.5 mm bore offset - an utterly ordinary geometry): 12 rings
        // and a 1.54 mm void up the thick side, 9% of the section's material missing.
        //
        // So test what actually matters - whether ONE count covers every angle inside the same
        // +-eps stretch budget the count law uses - and if it does, use it uniformly. That runs
        // the section down the well-tested constant-count path (rho == 1, no transition window)
        // rather than asking the k/k-1 blend to span two steps, which it cannot represent.
        // The width extremes sit at the thinnest and thickest bins, so two probes suffice.
        const coord_t w_ceil_sp = coord_t(std::llround(1.10 * double(w_sp)));   // the law's eps = 0.10
        const auto [it_lo, it_hi] = std::minmax_element(t_sp.begin(), t_sp.end());
        const coord_t t_lo = coord_t(std::llround(*it_lo)), t_hi = coord_t(std::llround(*it_hi));
        auto count_fits = [&](coord_t k, coord_t t) {
            const Arachne::BeadingStrategy::Beading b = strategy->compute(t, k);
            if (b.bead_widths.size() != size_t(k))
                return false;
            for (const coord_t w : b.bead_widths)
                if (w < w_floor_sp || w > w_ceil_sp)
                    return false;
            return true;
        };
        coord_t k_uniform = 0;
        for (coord_t k = n_min; k <= n_max && k_uniform == 0; ++k)
            if (count_fits(k, t_lo) && count_fits(k, t_hi))
                k_uniform = k;
        if (k_uniform == 0)
            // Genuinely out of budget: no single count is printable across this wall, and the
            // blend cannot cross two steps. DECLINE (fall through to marching) rather than
            // hard-fail - refusing to slice a part that used to slice is a worse failure than
            // the underfill, and per-angle bead COUNTS are the real fix (see W1b in the handoff).
            { if (decline_note) *decline_note = "wall thickness crosses more than one bead-count threshold"; return nullptr; }
        std::fill(n_bin.begin(), n_bin.end(), k_uniform);
        n_min = n_max = k_uniform;
    }

    // Blend factor per bin. Where the count is uniform this is constant 1 and the medial
    // slot is a full ring. Where t(theta) crosses a count threshold, the bins at the
    // higher count form angular windows; each window edge gets a transition ramp of the
    // configured arc length, anchored mostly on the squeezed (higher-count) side per the
    // strategy's own anchor formula, and the blend inside a ramp is a smoothstep of arc
    // position. The appearing (medial) bead is emitted only where its blended width
    // clears the termination floor - below ~1.2x layer height a TPU bead die-swells and
    // gaps (the printed T5 verdict), so a transitional tip must not dip that low.
    // Two mechanisms lower the blend below 1 (a mid-blend medial bead, floored below):
    //  * ANGULAR: t(theta) crosses a count threshold (an eccentric wall) - the bins at
    //    the higher count form windows with ramped edges where the geometry demands the
    //    extra bead; through the thin sector the medial thins toward the floor and its
    //    neighbours squeeze to compensate.
    //  * VERTICAL: t is uniform but sits just past a flip (a tapering wall) - the whole
    //    circle blends at the coverage fraction c of the flip band (0.25 spacings of
    //    thickness), so the newborn ring is BORN as a full thin ring at the termination
    //    floor and WIDENS to plan width across the band's many layers, instead of two
    //    band-derived rings appearing inside one layer. No bead ever terminates mid-air
    //    (the printed T5 verdict: tips gap) and no chord ever crosses the bore - the
    //    count changes topology only at the single flip layer, as a one-floor-width step.
    //    Stateless: c is a pure function of THIS layer's thickness.
    coord_t k = n_max;                              // slot count of the plan
    std::vector<double> rho(NB, 1.0);
    std::vector<char> hi;
    if (n_min != n_max) {
        hi.assign(NB, 0);
        for (int b = 0; b < NB; ++b)
            hi[b] = char(n_bin[b] == k);
    } else if (k > 2) {
        // Find the flip threshold just below the (uniform) thickness by bisecting the
        // count law itself - robust to every clamp and floor inside it. Sitting more
        // than one flip band above it means the medial is a settled full ring.
        double t_mean = 0.;
        for (int b = 0; b < NB; ++b)
            t_mean += t_sp[b];
        t_mean /= double(NB);
        const double dt_band = 0.25 * double(w_sp);
        const double below   = t_mean - dt_band;
        if (below > 0 && squeeze_biased_count(coord_t(std::llround(below)), *strategy, w_sp, w_floor_sp) < k) {
            double a = below, b2 = t_mean;          // the flip threshold sits in (a, b2]
            for (int it = 0; it < 24; ++it) {
                const double m = 0.5 * (a + b2);
                if (squeeze_biased_count(coord_t(std::llround(m)), *strategy, w_sp, w_floor_sp) < k)
                    a = m;
                else
                    b2 = m;
            }
            const double c = std::clamp((t_mean - b2) / dt_band, 0.0, 1.0);
            if (c < 1.0)
                std::fill(rho.begin(), rho.end(), c);
        }
    }
    if (!hi.empty()) {
        // Mean radius of the wall midline, for arc <-> angle conversion.
        double r_mid = 0.;
        for (int b = 0; b < NB; ++b)
            r_mid += 0.5 * (plan->r_out.r[b] + plan->r_bore.r[b]);
        r_mid /= double(NB);
        const double bin_rad   = 2.0 * M_PI / NB;
        const double L_rad     = std::min(M_PI / 2, scale_(trans_len_mm) / r_mid);
        const double F_rad     = std::min(M_PI / 4, scale_(std::max(0., limits.wall_transition_filter_deviation_mm)) / r_mid);
        const int    L_bins    = std::max(1, int(std::llround(L_rad / bin_rad)));
        const int    F_bins    = int(std::llround(F_rad / bin_rad));
        // Window hygiene on the cyclic high-count set: close gaps shorter than a ramp,
        // then drop islands shorter than two ramps plus the filter deviation. Runs are
        // found by linear scan with wraparound handling.
        auto morph = [&](char target, int min_len) {
            std::vector<char> next = hi;
            int b = 0;
            while (b < NB) {
                if (hi[b] != target) { ++b; continue; }
                int len = 1;
                while (len < NB && hi[(b + len) % NB] == target)
                    ++len;
                if (len >= NB)
                    break;                                   // uniform: nothing to do
                if (b == 0 && hi[NB - 1] == target) { ++b; continue; }  // handled from its true start
                if (len < min_len)
                    for (int i = 0; i < len; ++i)
                        next[(b + i) % NB] = char(1 - target);
                b += len;
            }
            hi = std::move(next);
        };
        morph(0, L_bins);                                    // merge windows across short gaps
        morph(1, 2 * L_bins + F_bins);                       // drop windows too short for a bead
        int hi_count = 0;
        for (int b = 0; b < NB; ++b)
            hi_count += hi[b];
        if (hi_count == 0) {
            k = n_max - 1;                                   // the higher count filtered away: uniform k-1
        } else if (hi_count < NB) {
            // Anchor: fraction of the ramp placed on the lower-count side. With the
            // squeeze-biased law the flip sits near the lower optimal thickness, so the
            // anchor lands close to 1 and the ramp lives on the squeezed side.
            const double lower_opt = double(strategy->getOptimalThickness(k - 1));
            const double upper_opt = double(strategy->getOptimalThickness(k));
            double t_flip = 0.;                              // mean t at the window edges
            int    n_edges = 0;
            for (int b = 0; b < NB; ++b)
                if (hi[b] != hi[(b + 1) % NB])
                    { t_flip += t_sp[b]; ++n_edges; }
            if (n_edges > 0) t_flip /= n_edges;
            const double anchor = std::clamp(1.0 - (t_flip - lower_opt) / std::max(1.0, upper_opt - lower_opt), 0.0, 1.0);
            // rho from the signed bin distance to the nearest window edge: ramps rise
            // across [edge - anchor*L, edge + (1-anchor)*L] into the window.
            std::vector<int> dist(NB, NB);                   // +into-window bins, -outside
            for (int b = 0; b < NB; ++b) {
                if (!hi[b]) continue;
                int din = NB;
                for (int d = 0; d <= NB / 2; ++d)
                    if (!hi[((b + d) % NB)] || !hi[((b - d % NB) + NB) % NB]) { din = d; break; }
                dist[b] = din;
            }
            for (int b = 0; b < NB; ++b) {
                if (hi[b]) continue;
                int dout = NB;
                for (int d = 0; d <= NB / 2; ++d)
                    if (hi[((b + d) % NB)] || hi[((b - d % NB) + NB) % NB]) { dout = d; break; }
                dist[b] = -dout;
            }
            for (int b = 0; b < NB; ++b) {
                const double x = (double(dist[b]) + anchor * double(L_bins)) / double(L_bins);
                const double u = std::clamp(x, 0.0, 1.0);
                rho[b] = u * u * (3.0 - 2.0 * u);            // smoothstep
            }
        }
        // else: hygiene filled the whole circle -> uniform k, rho stays 1.
    }

    // Fill the per-bin slot tables and check the plan's own invariants. These are hard
    // failures: the stack guarantees sum(widths) + left_over == t, so a violation here
    // means the plan and the geometry disagree and the layer must not print.
    // The transitional (mid) bead is FLOORED at the termination width - it prints as a
    // full thin ring rather than terminating (tips gap; the printed T5 verdict) - and
    // the floor deficit is taken back from the interior beads so the covered thickness
    // is conserved exactly and everything stays in the squeeze regime.
    plan->n_slots = size_t(k);
    const size_t mid = size_t(k / 2);                        // the appearing/retiring slot
    plan->slot_width.assign(NB, std::vector<float>(plan->n_slots, 0.f));
    plan->slot_frac.assign(NB, std::vector<double>(plan->n_slots, 0.));
    const double overlap_delta_mm = double(w_ext - w_sp) * SCALING_FACTOR;
    // Floor the transitional (mid) bead of a beading and take the deficit back from the
    // interior beads (the pinned surfaces only when there is no interior). Shared by the
    // ring tables and the spiral turn tables below.
    auto floor_transitional = [&](Arachne::BeadingStrategy::Beading &bd, size_t count) {
        const size_t midc = count / 2;
        if (count <= 2 || bd.bead_widths[midc] >= w_floor_sp)
            return;
        const coord_t orig_mid = bd.bead_widths[midc];
        const coord_t deficit  = w_floor_sp - orig_mid;
        // Captured BEFORE any width changes: this is the outer flank of bead 0, the surface
        // the beading pins, and the fixed point the re-tiling below works from.
        const coord_t flank = bd.toolpath_locations.front() - bd.bead_widths.front() / 2;
        std::vector<size_t> recipients;
        for (size_t i = 1; i + 1 < count; ++i)
            if (i != midc)
                recipients.push_back(i);
        if (recipients.empty())
            for (size_t i = 0; i < count; ++i)
                if (i != midc)
                    recipients.push_back(i);
        // Take the deficit only from beads that can AFFORD it. A recipient rebated below
        // w_term is exactly the bead the floor exists to prevent - below ~1.2x layer height a
        // bead die-swells and gaps (the printed T5 verdict) - so paying for the floor by
        // pushing a neighbour under it trades one unprintable bead for another. PlateF's lobed
        // tube showed the cost: a 0.05 mm width ladder running down to 0.25 mm, 164 mm of
        // extrusion per layer in one azimuth sector, and a column of voids in the cut.
        //
        // Where the budget is short the transitional bead is lifted only as far as the
        // interior can pay for. It is the bead that is being born or retired, and the one the
        // design already treats as special; the alternative puts several structural beads
        // under the floor to lift one transient bead over it.
        coord_t headroom = 0;
        for (const size_t i : recipients)
            headroom += std::max<coord_t>(0, bd.bead_widths[i] - w_floor_sp);
        const coord_t pay = std::min(deficit, headroom);
        if (pay <= 0)
            return;                       // nothing affordable: leave the beading untouched
        coord_t remaining = pay;
        for (int pass = 0; pass < 2 && remaining > 0; ++pass) {
            coord_t room_total = 0;
            for (const size_t i : recipients)
                room_total += std::max<coord_t>(0, bd.bead_widths[i] - w_floor_sp);
            if (room_total <= 0)
                break;
            const coord_t want = remaining;
            for (const size_t i : recipients) {
                if (remaining <= 0)
                    break;
                const coord_t room = std::max<coord_t>(0, bd.bead_widths[i] - w_floor_sp);
                if (room <= 0)
                    continue;
                // Proportional to each bead's headroom, clamped to it; the second pass mops up
                // integer rounding against whatever headroom is left.
                coord_t take = coord_t(double(want) * double(room) / double(room_total));
                take = std::clamp<coord_t>(take, 1, std::min(room, remaining));
                bd.bead_widths[i] -= take;
                remaining -= take;
            }
        }
        bd.bead_widths[midc] = orig_mid + (pay - remaining);
        // Arachne places beads EDGE TO EDGE, so rebating widths without moving the
        // centrelines does not conserve coverage the way the totals suggest: the widened
        // transitional bead overlaps its neighbours by deficit/2 on each side, and every
        // rebated recipient opens share/2 of void on each of ITS sides. At count 3 the
        // recipient list falls back to {0, count-1} - the PINNED SURFACE beads - so both
        // free skins of an ordinary tube lose ~share/2 of their flank, every layer of a
        // transition band. Re-tile from the outer flank of bead 0, which is the surface
        // the beading pins: a narrowed surface bead then keeps its flank ON the surface
        // and moves its centreline inward, instead of pulling the flank in and
        // under-extruding the skin.
        coord_t acc = flank;
        for (size_t i = 0; i < count; ++i) {
            bd.toolpath_locations[i] = acc + bd.bead_widths[i] / 2;
            acc += bd.bead_widths[i];
        }
    };
    for (int b = 0; b < NB; ++b) {
        const coord_t t_b = coord_t(std::llround(t_sp[b]));
        Arachne::BeadingStrategy::Beading bd = strategy->compute(t_b, k);
        if (bd.bead_widths.size() != plan->n_slots || bd.toolpath_locations.size() != plan->n_slots)
            { if (hard_fail) *hard_fail = "Continuous shell adaptive beading: bead count mismatch in the plan."; return nullptr; }
        if (rho[b] < 1.0)
            bd = blend_ring_beadings(strategy->compute(t_b, k - 1), bd, rho[b], mid);
        if (std::abs(double(bd.left_over)) > 0.01 * t_sp[b])
            { if (hard_fail) *hard_fail = "Continuous shell adaptive beading: plan left uncovered wall thickness."; return nullptr; }
        floor_transitional(bd, plan->n_slots);
        const double t_geom = t_sp[b] + 2.0 * plan->edge_margin;
        coord_t prev_loc = std::numeric_limits<coord_t>::min();
        for (size_t i = 0; i < plan->n_slots; ++i) {
            const coord_t w = bd.bead_widths[i];
            const coord_t loc = bd.toolpath_locations[i];
            if (w < w_floor_sp / 2 || loc <= prev_loc
                || (prev_loc != std::numeric_limits<coord_t>::min() && loc - prev_loc < w_sp / 4))
                { if (hard_fail) *hard_fail = "Continuous shell adaptive beading: implausible bead placement in the plan."; return nullptr; }
            prev_loc = loc;
            const float w_mm = float(double(w) * SCALING_FACTOR + overlap_delta_mm);
            plan->slot_width[b][i] = std::round(w_mm / 0.05f) * 0.05f;
            plan->slot_frac[b][i]  = (plan->edge_margin + double(loc)) / t_geom;
        }
    }

    // ---- Interior spiral (phase B) ---------------------------------------------------
    // The buried rings become one double-arm spiral whose comb (numerically validated)
    // has the SAME bead count and pitch as ring mode: both arms span the full interior
    // surface-to-surface at pitch 2 per revolution with a constant one-slot offset, so
    // internal comb gaps equal one pitch at every azimuth and the arm ends tangent-
    // merge onto the surface rings (no terminations in air). All geometry is analytic
    // and evaluated in spiralize_rings from the ordinary slot tables; nothing more to
    // precompute here beyond the flag, the phase anchor, and the merge width floor.
    if (spiral_mode) {
        if (k < 6) {
            if (decline_note)
                *decline_note = "interior too thin for a spiral (needs at least 4 buried rings)";
        } else {
            plan->spiral        = true;
            plan->spiral_theta0 = seam_theta;
            plan->ramp_frac     = std::clamp(junction_ramp, 0.05, 1.0);
            plan->w_term        = w_term;
        }
    }
    return plan;
}

// Generate the plan's rings: iso-index curves sampled at the bin centres. The radial
// position interpolates between the boundary profiles by the slot's planned fraction,
// which projects to the exact planned perpendicular depth on straight wall segments.
// spiralize_rings resamples every ring to its final toolpath density afterwards.
Polygons plan_rings(const ChainBeadPlan &plan)
{
    Polygons rings;
    rings.reserve(plan.n_slots);
    for (size_t slot = 0; slot < plan.n_slots; ++slot) {
        Polygon ring;
        ring.points.reserve(plan.NB);
        for (int b = 0; b < plan.NB; ++b) {
            const double th = -M_PI + (double(b) + 0.5) * (2.0 * M_PI / plan.NB);
            const double ro = plan.r_out.r[b];
            const double rb = plan.r_bore.r[b];
            const double r  = ro - plan.slot_frac[b][slot] * (ro - rb);
            ring.points.emplace_back(coord_t(std::llround(plan.centre.x() + r * std::cos(th))),
                                     coord_t(std::llround(plan.centre.y() + r * std::sin(th))));
        }
        rings.push_back(std::move(ring));
    }
    return rings;
}

// Deterministic per-layer pseudo-random value in [-1, +1]. Deterministic so
// re-slicing the same model yields identical G-code.
double seam_jitter_unit(size_t layer_index)
{
    uint64_t x = uint64_t(layer_index) + 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x =  x ^ (x >> 31);
    const double unit = double(x >> 11) / double(uint64_t(1) << 53); // [0, 1)
    return unit * 2.0 - 1.0;
}

// Accumulated random walk of the per-layer dither: the running sum of seam_jitter_unit
// over all prior layers. Each STEP is small (so the per-layer seam move stays a short
// inter-layer join), but the SUM wanders around the part over the whole print. Scaled up
// it scatters the seam instead of leaving one coherent helical line; scaled to zero the
// seam is a clean helix. Deterministic so re-slicing is stable. The prefix sums are
// memoized per thread (callers run inside TBB'd make_perimeters), keeping the whole
// print O(L) instead of O(L^2).
double seam_jitter_walk(size_t layer_index)
{
    thread_local std::vector<double> prefix{0.0};      // prefix[k] = sum of units 1..k
    while (prefix.size() <= layer_index)
        prefix.push_back(prefix.back() + seam_jitter_unit(prefix.size()));
    return prefix[layer_index];
}

// Trace a radially-ordered set of nested contours (outer at index 0 -> innermost) as a
// single continuous polyline using a CONNECTED FERMAT SPIRAL: every contour is emitted
// as a complete loop, ordered inward along even indices and back outward along odd
// indices with a U-turn at the innermost. Both endpoints therefore land on the two
// outermost contours, so the outer surface (and, for an annulus, the bore) is a full
// clean loop with no radial seam - unlike an inward morph, which lays the unavoidable
// once-per-revolution one-bead radial step onto the surface as a visible seam groove.
//
// Every contour is anchored to a single ABSOLUTE seam angle about the cross-section
// centre. `reverse` flips the contour order per layer so consecutive layers start/end
// on the same contour (adjacent points -> tiny joins). `seam_step_rad` (a fixed
// per-layer ANGLE) and `seam_jitter_rad` (a bounded random angle) walk that anchor a
// little each layer so the residual junction spreads into a shallow helix / scatter.
//
// The seam MUST accumulate an angle, not an arc length. An earlier version anchored
// the seam by absolute arc (layer_index * spacing) divided by the layer's reference
// radius; where a part tapers toward the top (a doming cap) that radius shrinks while
// the accumulated arc is already tens of revolutions, so the seam angle swung by
// radians per layer. Consecutive layers then started/ended far apart, and because
// continuous-shell mode drops the inter-layer travel move, each gap was laid down as
// a long extruding chord straight across the part (visible stringing under a cap).
Polyline spiralize_rings(const Polygons &rings_in, bool reverse, double seam_step_rad, double seam_jitter_rad,
                         double ring_spacing, size_t layer_index, bool has_bore, std::vector<char> *internal,
                         coord_t interlock_amp = 0, std::vector<coord_t> *z_off = nullptr,
                         double connector_spread = 0.0, size_t infill_shell_loops = 0,
                         const ExPolygon *region = nullptr,
                         const std::vector<float> *ring_widths = nullptr, std::vector<float> *pt_width = nullptr,
                         float nominal_width = 0.f,
                         coord_t scarf_len = 0, float scarf_floor = 0.f,
                         const ChainBeadPlan *plan = nullptr, coord_t layer_h_scaled = 0,
                         // Which plan SLOT each ring came from. Identity while nothing has been
                         // decimated; after sparse infill drops rings it is the only thing that
                         // still ties a surviving ring to its own per-angle widths.
                         const std::vector<size_t> *ring_slot = nullptr)
{
    Polyline spiral;
    if (rings_in.empty())
        return spiral;
    const Polygons &rings = rings_in;    // original radial order: outer (0) -> inner

    // Per-ring resample density. These points ARE the printed toolpath: z_contoured beads
    // are exempt from downstream simplification and arc fitting (ExtrusionPath::simplify),
    // so every sample emitted here becomes one G-code segment. Target a half-bead chord so
    // detail down to the bead scale survives - the old fixed 256-point cap put the corners
    // of a large ring up to ~0.8 mm off their true position. The two FREE-SURFACE rings get
    // a generous cap (their fidelity is the printed finish; ~2k points matches what normal
    // slicing emits for a large perimeter). Buried rings only shape density and the weave,
    // so they keep a tight cap: extra points there are pure G-code bloat.
    const double chord = std::max(0.5 * ring_spacing, 1.0);
    auto ring_samples = [&](const Polygon &r, bool surface) -> size_t {
        const size_t cap = surface ? 2048 : 512;
        return std::clamp<size_t>(size_t(std::llround(r.length() / chord)), size_t(64), cap);
    };

    // Cross-section centre (centroid of the largest-area ring = the outer contour) and
    // its mean radius. The seam is anchored by absolute angle about this centre.
    Vec2d centre(0., 0.);
    double r_ref = 1.;
    {
        size_t best = 0;
        double best_area = -1.;
        for (size_t i = 0; i < rings.size(); ++i) {
            const double a = std::abs(rings[i].area());
            if (a > best_area) { best_area = a; best = i; }
        }
        const Polygon &outer = rings[best];
        if (!outer.points.empty()) {
            for (const Point &p : outer.points)
                centre += p.cast<double>();
            centre /= double(outer.points.size());
            double rr = 0.;
            for (const Point &p : outer.points)
                rr += (p.cast<double>() - centre).norm();
            r_ref = std::max(1.0, rr / double(outer.points.size()));
        }
    }

    // Seam angle for this layer: a fixed per-layer angular precession (the helix) plus an
    // ACCUMULATED random walk (the scatter). Both are ANGLES (radians), independent of the
    // reference radius, so a tapering cross-section cannot make the seam jump. The per-layer
    // CHANGE stays small either way, so consecutive layers still join with a short hop; the
    // random part accumulates, so a high jitter scatters the seam around the whole part
    // rather than leaving one coherent helical line. jitter=0 -> a clean helix. Layer 0
    // anchors at angle 0.
    const double seam_theta = double(layer_index) * seam_step_rad
                            + seam_jitter_rad * seam_jitter_walk(layer_index);

    // Per-ring cut rotation: rotate each successive contour's seam a little so the
    // Fermat connectors and the innermost U-turn spiral inward instead of stacking into
    // one radial cut line (both a residual junction and an aligned internal weak seam).
    // Bounded to ~one bead per 2-ring Fermat skip and capped so the innermost rings do
    // not over-spin. This is WITHIN a layer, so unlike the seam precession it never
    // accumulates across layers (no r_ref instability).
    //
    // `connector_spread` (A2, user knob) ADDS to the base twist: it fans each ring's
    // handoff further around the circumference so the once-per-layer surface->interior
    // crossings no longer STACK at one angle into a deep radial groove (the visible spiral
    // seam) but distribute into shallow, spread-out steps. At 100% the handoffs spread
    // across a full revolution; at 0% only the tight base twist remains (current behaviour).
    //
    // Restricted to ANNULAR chains (has_bore). The knob exists to break up the deep spiral
    // seam groove on a TUBE WALL, which is inherently an annulus. On a SOLID cross-section
    // there is no surface->bore seam to spread; fanning the inner connectors instead swings
    // them across the shape's OWN concavities (the notches between the arms of a '+'/star, a
    // dogbone waist), laying stringing chords through empty space outside the part. Solid
    // chains therefore keep only the tight base twist. This also leaves a tube's solid
    // closing cap (has_bore=false there) with tight connectors, which it wants anyway.
    const double base_twist   = std::clamp(0.5 * ring_spacing / r_ref, 0.0, 0.2);
    const double spread_twist = (connector_spread > 0.0 && rings.size() > 1 && has_bore)
                              ? connector_spread * (2.0 * M_PI / double(rings.size())) : 0.0;
    const double ring_twist   = base_twist + spread_twist;

    // Radius-aware seam anchors. A straight hop between two rings laid at a spread seam angle
    // chords across the hollow bore where the rings are small (the inner rings); the more
    // connector_spread, the wider that chord, stringing through the void. So build each ring's
    // anchor as a CUMULATIVE angle whose per-ring step is clamped by the local radius: a hop at
    // radius r spanning angle phi has closest approach r*cos(phi/2), which must stay outside the
    // innermost (bore) radius. Outer rings (large r) still spread the full twist - preserving
    // A2's surface-seam distribution - while inner rings near the bore barely rotate, so no hop
    // crosses the void.
    const size_t n_ring = rings.size();
    std::vector<double> ring_r(n_ring, r_ref);
    for (size_t i = 0; i < n_ring; ++i) {
        double rr = 0.; size_t cnt = 0;
        for (const Point &p : rings[i].points) { rr += (p.cast<double>() - centre).norm(); ++cnt; }
        if (cnt) ring_r[i] = std::max(1.0, rr / double(cnt));
    }
    const double r_bore = *std::min_element(ring_r.begin(), ring_r.end());
    std::vector<double> ring_angle(n_ring, seam_theta);
    {
        double theta = seam_theta;
        for (size_t i = 0; i < n_ring; ++i) {
            if (i > 0) {
                // The 2-ring Fermat hop spans ~2 steps, so cap each step at acos((bore+bead)/r):
                // the hop's chord then stays at least one bead-spacing clear of the bore.
                const double lim  = std::clamp((r_bore + ring_spacing) / std::max(ring_r[i], r_bore + ring_spacing), 0.0, 1.0);
                const double safe = std::acos(lim);
                theta += std::min(ring_twist, safe);
            }
            ring_angle[i] = theta;
        }
    }

    auto anchor_to_angle = [&](std::vector<Point> &pts, double target) {
        size_t best = 0;
        double best_d = std::numeric_limits<double>::max();
        for (size_t k = 0; k < pts.size(); ++k) {
            const Vec2d d = pts[k].cast<double>() - centre;
            double a = std::atan2(d.y(), d.x()) - target;
            a = std::atan2(std::sin(a), std::cos(a)); // shortest wrap to [-pi, pi]
            if (std::abs(a) < best_d) { best_d = std::abs(a); best = k; }
        }
        if (best > 0)
            std::rotate(pts.begin(), pts.begin() + best, pts.end());
    };

    std::vector<std::vector<Point>> grid;
    grid.reserve(rings.size());
    for (size_t ri = 0; ri < rings.size(); ++ri) {
        // Every ring must sweep the same rotational direction, otherwise morphing a
        // CCW ring (outer contour) to a CW ring (a hole contour, per Slic3r's
        // ExPolygon winding convention) maps each point to its antipode and the
        // interpolated chords cut straight across the bore (the "leaf" artifact).
        Polygon oriented = rings[ri];
        oriented.make_counter_clockwise();
        const bool surface_ring = ri == 0 || (has_bore && ri + 1 == rings.size());
        std::vector<Point> pts = resample_closed(oriented, ring_samples(oriented, surface_ring));
        if (pts.empty())
            return {};
        anchor_to_angle(pts, ring_angle[ri]);
        grid.push_back(std::move(pts));
    }

    // A ring is a free SURFACE if it traces a real region boundary: the outermost
    // ring is always the outer surface, and (for an annulus) the innermost ring is
    // the bore surface. Interior rings are buried. `internal` (if requested) tags
    // each emitted point so the caller can over-extrude only provably-buried beads.
    const size_t n_rings = grid.size();
    auto is_surface = [&](size_t j) -> bool {
        return j == 0 || (has_bore && j + 1 == n_rings);
    };
    // When sparse woven infill is active the loop list is already [outer shell | sparse infill |
    // (bore shell, annular)]. A ring is INFILL (thinned + woven, so tagged internal and given the
    // Z-wave) when it is past the outer shell and, for an annulus, before the bore shell. This
    // supersedes is_surface for the internal tag so BOTH solid skins stay flat and un-woven.
    auto is_infill = [&](size_t j) -> bool {
        if (infill_shell_loops == 0 || j < infill_shell_loops) return false;
        return has_bore ? (j + infill_shell_loops < n_rings) : (j + 1 < n_rings);
    };
    auto flag = [&](bool is_internal) { if (internal) internal->push_back(is_internal ? char(1) : char(0)); };
    auto zoff = [&](coord_t z) { if (z_off) z_off->push_back(z); };
    // Per-point extrusion width (adaptive bead width): each point carries its ring's
    // planned width, nominal when the chain has no explicit plan.
    auto width_of = [&](size_t idx) -> float {
        return (ring_widths != nullptr && idx < ring_widths->size()) ? (*ring_widths)[idx] : nominal_width;
    };
    auto pw = [&](float w) { if (pt_width) pt_width->push_back(w); };
    // Per-POINT width: an adaptive bead plan supplies width as a function of the point's
    // angle (so ring resampling can never desynchronize a width array); chains without a
    // plan keep the per-ring width.
    auto point_w = [&](size_t idx, const Point &p) -> float {
        // Position in `rings_in` is NOT the plan slot once decimation has dropped rings, so go
        // through the slot map whenever one is supplied.
        const size_t slot = (ring_slot != nullptr && idx < ring_slot->size()) ? (*ring_slot)[idx] : idx;
        if (plan != nullptr && slot < plan->n_slots) {
            const Vec2d d = p.cast<double>() - plan->centre;
            return plan->width_at(slot, std::atan2(d.y(), d.x()));
        }
        return width_of(idx);
    };

    // Brick interlock: nudge buried wall rings up/down in Z so radially-adjacent beads
    // stop sharing a Z level (a continuous horizontal interface is the wall's weak plane).
    // Even rings rise +amp, odd rings drop -amp, so any two adjacent rings differ by 2*amp
    // (= half a bead height at 100%). Only ANNULAR walls are staggered, and the two rings
    // that carry the layer-to-layer join are pinned to the smooth ramp so the bead stays
    // continuous: the outer surface (ring 0), the bore surface (innermost), and ring 1 -
    // the Fermat spiral's other free end. The interior rings then form a static radial
    // corrugation that interlocks vertically with the neighbouring layer's beads.
    auto ring_z_off = [&](size_t idx) -> coord_t {
        if (interlock_amp == 0 || !has_bore || infill_shell_loops > 0) return 0;
        if (is_surface(idx) || idx == 1)     return 0;
        return (idx % 2 == 0) ? interlock_amp : -interlock_amp;
    };

    if (grid.size() == 1) {
        for (size_t k = 0; k < grid[0].size(); ++k) { spiral.points.push_back(grid[0][k]); flag(false); zoff(0); pw(point_w(0, grid[0][k])); }
        spiral.points.push_back(grid[0][0]); flag(false); zoff(0); pw(point_w(0, grid[0][0]));
        return spiral;
    }

    // ---- Interior single-arm spiral assembly (phase B, v5) --------------------------
    // The two grid rings are the closed surface loops; the buried interior is ONE
    // spiral arm at a pitch of one slot per revolution, generated analytically from
    // the plan. Arm tips tangent-merge onto the surface rings at clearance widths
    // floored at the transitional width - no bead terminates in air. The traversal
    // direction is FIXED (outer -> bore, every layer): a per-layer reversal
    // reintroduced the documented "reverse-join doubling" (paired layers with
    // two-layer voids). The inter-layer join is one radial diagonal across the wall,
    // through material, at the precessing seam azimuth.
    //
    // The intrinsic Z-stagger (v4 sine -> v5 staircase). Three printed revisions
    // (double-arm, single-arm reversing, single-arm fixed) all delaminated along
    // helical interfaces; the peel survived fixing both cold arm-to-arm joints
    // (v2) and reverse-join doubling (v3), so neither was the root cause. What
    // ring mode has and every spiral traversal lacks is RADIAL INTERDIGITATION:
    // the Fermat ring ORDER scrambles the arc-length ramp fraction across radial
    // neighbours (the natural bricklaying the interlock A/B print showed
    // tear-proofs the wall), while a monotone spiral offsets radial neighbours by
    // a near-constant h/k - smooth laminate planes the printed tubes cleaved
    // along. v4 restored interdigitation with z = (h/4)*sin(2*pi*g*v) and the
    // printed tube did NOT peel (mechanism confirmed), but its section showed
    // organized void rows: a sinusoid makes the interdigitation depth OSCILLATE
    // along each interface, so every interface sweeps through the void-producing
    // intermediate and coplanar configurations somewhere on the circumference.
    // Ring mode's interfaces hold one CONSTANT depth all the way around (both
    // neighbours' ramp fractions advance together), each interface a different
    // depth. v5 held one golden level per revolution and printed with pinprick
    // packing in part - but split along remaining weak interfaces, because the
    // levels were designed on the staircase ALONE. The physical offset between
    // adjacent revolutions is staircase PLUS the arc-length ramp's own tilt
    // (~0.13 h per revolution), which shifted every designed 0.38/0.62 h gap to
    // an effective ~0.25 h: uniformly T-junction grade, the explicit interlock's
    // documented tear geometry. v6 designs the TOTAL offset instead:
    //  * Level targets follow the golden low-discrepancy sequence in total
    //    offset; the emitted staircase subtracts the accumulated ramp tilt
    //    (estimated per revolution from the plan radii) and wraps modulo one
    //    layer height - offsets are nest-equivalent mod h - so every physical
    //    gap is a true golden gap, effective 0.38 h from coplanar.
    //  * Steps stay SIMULTANEOUS at the seam azimuth. (Scattering them per
    //    revolution was considered and is wrong: on a single continuous arm,
    //    revolution j+1 before its own step still holds revolution j's level,
    //    so scattered steps make every interface exactly coplanar along the
    //    whole arc between the two step azimuths.) The cost of simultaneity is
    //    that an opposite-sign gap pair interpolates through zero inside the
    //    0.08-revolution step window - a short coplanar ribbon at the seam
    //    column, the already-sacrificial line, far milder than ring mode's own
    //    k-junction column which the wall demonstrably tolerates.
    //  * The stagger tapers out only within a quarter revolution of each
    //    tangent merge (surfaces stay pinned to the ramp); layer 0 is flat.
    if (plan != nullptr && plan->spiral && has_bore && grid.size() == 2 && plan->n_slots >= 6) {
        const double S   = double(plan->n_slots - 1);        // slot span, surface to surface
        const double th0 = seam_theta;                       // both junctions ride the seam column
        const double g_conj = 0.6180339887498949;            // golden-ratio conjugate
        const int    n_revs = int(S);                        // n_slots is integral, so S is too
        const bool   ramped = plan->ramp_frac < 0.999;
        const double f      = std::clamp(plan->ramp_frac, 0.05, 1.0);
        const int    n_col  = std::max(1, n_revs - 1);       // buried columns, slots 1 .. S-1
        // Stagger amplitude. A PURE SPIRAL needs the staircase unconditionally: its
        // interior is one helical thread and interdigitation is the only thing resisting
        // radial cleavage. A RAMPED arm does not - its ring-like columns are topologically
        // blocked from radial extraction, AND the arc-length ramp already offsets radially
        // adjacent columns by its own rise per revolution (measured 0.17 of a layer on the
        // T5 tube, against designed offsets of 0.19/0.19/0.31 - so the staircase was
        // contributing +0.02/+0.02/+0.14, i.e. nothing on two interfaces and a
        // least-overlapping outlier on the third). So for a ramped arm the staircase rides
        // the existing brick-interlock knob, whose default of 0 means "keep the natural
        // nesting" - the same verdict that knob already carries for ring mode.
        // The wrap modulus below is one full LAYER, not twice the amplitude: offsets are
        // nest-equivalent mod the layer height, whatever amplitude is in use.
        const coord_t wamp = (layer_index == 0) ? 0
                           : (ramped ? interlock_amp : coord_t(layer_h_scaled / 4));
        const int    n_arm_rev = ramped ? n_col : n_revs;    // revolutions the arm actually makes
        std::vector<double> lvl(size_t(n_revs) + 1, 0.0);    // staircase level per revolution (scaled)
        if (wamp > 0) {
            const double h_sp = double(layer_h_scaled);      // one layer height, scaled
            const double ro   = plan->r_out.sample_linear(th0);
            const double rb   = plan->r_bore.sample_linear(th0);
            double L_est = 2.0 * M_PI * (ro + rb);           // the two surface rings...
            std::vector<double> rev_len(size_t(n_arm_rev), 0.0);
            for (int j = 0; j < n_arm_rev; ++j) {            // ...plus each arm revolution
                // A ramped revolution sits ON its column's radius; a spiral revolution
                // sweeps between two, so its mean radius is the midpoint. Getting this
                // wrong mis-estimates the arc-length ramp's rise per revolution, and the
                // tilt compensation then leaves the interdigitation depths off target -
                // which is exactly the quantity that sets column-to-column overlap.
                const double vq = ramped ? (1.0 + double(j)) : (double(j) + 0.5);
                const double r  = ro - plan->frac_at_f(std::min(vq, S), th0) * (ro - rb);
                rev_len[size_t(j)] = 2.0 * M_PI * r;
                L_est += rev_len[size_t(j)];
            }
            double tilt = 0.0;                               // accumulated ramp rise up to level j
            for (int j = 0; j <= n_revs; ++j) {
                double fr = (double(j) + 1.0) * g_conj;
                fr -= std::floor(fr);
                const double target = (2.0 * fr - 1.0) * double(wamp);     // TOTAL offset target
                const double x      = target - tilt;
                lvl[size_t(j)]  = x - h_sp * std::floor(x / h_sp + 0.5);   // wrap to [-h/2, h/2)
                // Level j is worn by the revolution that follows it: revolution j for the
                // spiral, revolution j-1 for a ramped arm (whose levels index columns 1..n).
                const int ri = ramped ? (j - 1) : j;
                if (ri >= 0 && ri < n_arm_rev)
                    tilt += h_sp * rev_len[size_t(ri)] / L_est;
            }
        }
        auto arm_wave = [&](double v) -> coord_t {
            const double env = std::clamp((std::min(v, S - v) - 0.25) * 2.0, 0.0, 1.0);
            const int    j   = std::clamp(int(std::floor(v)), 0, n_revs);
            const double u   = v - std::floor(v);            // revolution fraction past the seam azimuth
            constexpr double W = 0.08;                       // step ramp width, fraction of a revolution
            const double prev = (j > 0) ? lvl[size_t(j) - 1] : lvl[0];
            const double t    = std::clamp(u / W, 0.0, 1.0);
            const double s    = t * t * (3.0 - 2.0 * t);     // smoothstep into this revolution's level
            return coord_t(std::llround(env * (prev + (lvl[size_t(j)] - prev) * s)));
        };
        auto sp_point = [&](double v, double th) -> Point {
            const double ro = plan->r_out.sample_linear(th);
            const double rb = plan->r_bore.sample_linear(th);
            const double r  = ro - plan->frac_at_f(v, th) * (ro - rb);
            return Point(coord_t(std::llround(plan->centre.x() + r * std::cos(th))),
                         coord_t(std::llround(plan->centre.y() + r * std::sin(th))));
        };
        auto pitch_mm = [&](double th) -> double {
            const double ro = plan->r_out.sample_linear(th);
            const double rb = plan->r_bore.sample_linear(th);
            return (plan->frac_at(1, th) - plan->frac_at(0, th)) * (ro - rb) * SCALING_FACTOR;
        };
        auto arm_width = [&](double v, double th) -> float {
            double w = plan->width_at_f(std::clamp(v, 1.0, S - 1.0), th);
            const double d = std::min(v, S - v);             // slots to the nearest surface centre
            if (d < 1.0)
                w = std::min(w, std::max(double(plan->w_term), d * pitch_mm(th)));   // tangent-merge clearance
            return std::round(float(w) / 0.05f) * 0.05f;
        };
        // ---- Column schedule: the ring <-> spiral continuum -------------------------
        // tau is the revolution counter (angle / 2pi); v is the fractional slot the arm
        // occupies. A PURE SPIRAL sets v = tau: the arm migrates continuously, so every
        // element of the interior sits at a different radius from its neighbour along
        // the path and the whole thread can UNWIND - printed tubes peel radially apart
        // (v1-v6). A RAMPED schedule instead holds v at a whole slot for most of a
        // revolution and concentrates the migration into one short ramp of `ramp_frac`:
        // each column is then a constant-radius arc of (1 - ramp_frac) of a turn, well
        // past a half turn, so it cannot be pulled out radially without breaking - ring
        // mode's topological blocking - while the join ridge is spread over the ramp arc
        // instead of stacking into a seam column.
        //
        // Ramp azimuths MUST advance monotonically with the slot index. They are spread
        // evenly across one full turn from the first buried column to the last (so the
        // ramps disperse instead of stacking). If they were scattered non-monotonically,
        // two revolutions could occupy the SAME slot at one azimuth while the next slot
        // went uncovered - a void. Monotone, exactly one slot is doubled at any azimuth
        // (the junction overfill, as in ring mode) and no slot is ever empty.
        // Where the transitions sit is not a free parameter: DISPERSING them costs material,
        // one revolution per full turn of dispersion, because a column held from one
        // transition to the next spans 1 + (delta azimuth) turns. Measured, dispersed over
        // the whole circle: +15% extrusion against ring mode at ramp 6, 12 AND 20 alike -
        // flat in the ramp fraction, which is the signature of a structural extra
        // revolution rather than the transit arcs. Paying for it by thinning the arm was
        // printed and made the peel WORSE (thinner beads = less contact between columns,
        // worst at the bore where the merge clearance thins the arm too).
        //
        // So the transitions are ALIGNED at the seam azimuth, which makes the schedule
        // volumetrically exact: one revolution per buried column, two surface rings, the
        // same bead count ring mode lays. The arm therefore covers the BURIED columns only
        // (slots 1 .. S-1) and reaches the surfaces by a one-pitch hop at the seam, and
        // each column holds a constant radius for (1 - ramp) of a turn before handing off.
        // All the columns hand off inside the same narrow azimuth band - the seam column,
        // which is already the sacrificial line - so what the ramp buys over ring mode is a
        // junction spread into a shallow ramp instead of an abrupt radial jump.
        auto v_of_tau = [&](double tau) -> double {
            if (!ramped)
                return tau;
            const int    j = std::clamp(int(std::floor(tau)), 0, n_col - 1);
            const double u = tau - double(j);                  // position within revolution j
            if (j >= n_col - 1)
                return 1.0 + double(j);                        // innermost column: hold, then hop to the bore
            const double t = std::clamp((u - (1.0 - f)) / f, 0.0, 1.0);
            return 1.0 + double(j) + t * t * (3.0 - 2.0 * t);  // hold, then hand off at the seam
        };
        // The Z stagger follows the SAME schedule: interpolating the per-column levels at
        // fractional v means radius and height step together, in one localized event per
        // column, instead of the pure spiral's separate step window.
        auto ramp_wave = [&](double v, double tau) -> coord_t {
            // The end taper must run along the PATH (tau), not along the slot coordinate.
            // During a hold v is constant at the column's own slot, so a v-based taper does
            // not ease in at all - it zeroes the stagger for the WHOLE of the first and last
            // buried column, while their neighbours keep theirs. That leaves those two
            // interfaces at a mismatched offset (measured ~0.48 of a layer, barely half
            // overlapping) at one fixed radius, which prints as a continuous void line.
            // Tapering over a quarter revolution of path still kills the vertical step at
            // the entry and exit hops, without flattening a whole column.
            const double env_s = std::clamp((std::min(v, S - v) - 0.25) * 2.0, 0.0, 1.0);
            const double env_e = std::clamp(std::min(tau, double(n_col) - tau) / 0.25, 0.0, 1.0);
            const double env   = std::min(env_s, env_e);
            const int    j     = std::clamp(int(std::floor(v)), 0, n_revs - 1);
            const double q     = std::clamp(v - double(j), 0.0, 1.0);
            return coord_t(std::llround(env * (lvl[size_t(j)] + (lvl[size_t(j) + 1] - lvl[size_t(j)]) * q)));
        };
        // Sampling density per revolution follows the buried-ring chord policy; a short
        // ramp needs enough samples of its own to stay smooth.
        const double r_mid = 0.5 * (plan->r_out.sample_linear(th0) + plan->r_bore.sample_linear(th0));
        size_t n_rev = std::clamp<size_t>(size_t(std::llround(2.0 * M_PI * r_mid / chord)), size_t(64), size_t(512));
        if (ramped)
            n_rev = std::max(n_rev, size_t(std::ceil(24.0 / f)));
        const double tau_end = ramped ? double(n_col) : S;     // ramped: one turn per buried column
        const int    n_arm   = std::max(2, int(std::ceil(tau_end * double(n_rev))));
        std::vector<Point>   arm_pts;
        std::vector<float>   arm_ws;
        std::vector<coord_t> arm_z;
        arm_pts.reserve(size_t(n_arm) + 1);
        arm_ws.reserve(size_t(n_arm) + 1);
        arm_z.reserve(size_t(n_arm) + 1);
        for (int i = 0; i <= n_arm; ++i) {
            const double tau = tau_end * double(i) / double(n_arm);
            const double v   = v_of_tau(tau);
            const double th  = th0 + tau * 2.0 * M_PI;
            arm_pts.push_back(sp_point(v, th));
            arm_ws.push_back(arm_width(v, th));
            arm_z.push_back(ramped ? ramp_wave(v, tau) : arm_wave(v));
        }
        // Surface widths come from the RING tables' surface columns (the truncated
        // 2-ring grid does not index the plan's slot space directly).
        auto surf_w = [&](size_t slot, const Point &p) -> float {
            const Vec2d d = p.cast<double>() - plan->centre;
            return plan->width_at(slot, std::atan2(d.y(), d.x()));
        };
        anchor_to_angle(grid[0], th0);                       // both tips land at th0 (S integer revs)
        anchor_to_angle(grid[1], th0);

        // No width compensation for the ramp's redundant material. A ramp lays extra
        // material in the interstice between columns, costing roughly the ramp fraction
        // (+14% extrusion at ramp 20 against ring mode). Scaling the arm's widths down to
        // pay for it was TRIED AND PRINTED: thinning every buried bead ~15% shrinks the
        // contact between neighbouring columns, and the tube peeled WORSE than the pure
        // spiral - worst at the bore-adjacent interface, where the merge clearance already
        // thins the arm, so the two reductions compound. Interface contact is the thing
        // this whole workstream is trying to buy; trading it for volume is backwards.
        // Keep the planned widths and keep the ramp SHORT instead: redundancy scales with
        // the ramp arc, so a small ramp holds each column for most of a turn (the
        // topological property) at a few percent of extra material.
        auto emit_ring = [&](size_t gi, size_t slot) {
            for (const Point &p : grid[gi]) { spiral.points.push_back(p); flag(false); zoff(0); pw(surf_w(slot, p)); }
        };
        emit_ring(0, 0);                                     // outer -> in -> bore, every layer
        for (size_t i = 0; i < arm_pts.size(); ++i) {
            spiral.points.push_back(arm_pts[i]);
            flag(true);
            zoff(arm_z[i]);
            pw(arm_ws[i]);
        }
        emit_ring(1, plan->n_slots - 1);
        return spiral;
    }

    // Ordering of the nested contours into one path. The right primitive depends on
    // whether the innermost terminus is a free SURFACE or a buried core:
    //
    //  * Annulus (has_bore): the loop list runs outer-surface -> medial -> bore-surface,
    //    so both path ends are free surfaces. We use the connected-Fermat interleave
    //    (even contours inward, odd contours back outward; Zhao et al.) which lands both
    //    ends at the outer vicinity and keeps the bore a clean mid-path loop - no radial
    //    seam groove on either surface.
    //
    //  * Solid core (!has_bore): the contours nest all the way to a buried centre. Here
    //    the Fermat interleave produces a "tent" height profile - its two arms put the
    //    same radius at heights up to a full layer apart (weak diagonal bonds, uneven
    //    density). A SIMPLE sequential spiral (outer -> core) is monotone in BOTH path
    //    and depth, giving a clean gentle ramp with strong full-face lateral bonds; the
    //    only seam is the buried core point. Per-layer reversal alternates the buried
    //    core join with a (precessed) outer join.
    std::vector<size_t> order;
    order.reserve(n_rings);
    if (has_bore) {
        for (size_t i = 0; i < n_rings; i += 2)     // inward arm: 0, 2, 4, ...
            order.push_back(i);
        const size_t start_odd = (n_rings % 2 == 0) ? (n_rings - 1) : (n_rings - 2);
        for (size_t i = start_odd; ; i -= 2) {      // outward arm: ..., 5, 3, 1
            order.push_back(i);
            if (i < 2) break;
        }
    } else {
        for (size_t i = 0; i < n_rings; ++i)        // simple spiral: outer -> core
            order.push_back(i);
    }
    // Per-layer reversal keeps consecutive SOLID layers' start/end on the same contour
    // (adjacent points -> tiny joins) because a solid simple-spiral has its two ends at
    // DIFFERENT radii (outer surface vs buried core). An ANNULUS Fermat instead has both
    // ends near the OUTER surface (ring 0 and ring 1 are radially adjacent), so it already
    // joins layer-to-layer with a tiny hop WITHOUT reversal - and reversing it is exactly
    // what makes each ring's arc-position (hence its Z under the arc-length ramp) alternate
    // every layer, pairing beads up into a once-per-two-layers vertical void
    // ("reverse-join doubling"). So reverse only the solid case; keep the annulus stable.
    if (reverse && !has_bore)
        std::reverse(order.begin(), order.end());

    {
        size_t total = 0;
        for (const std::vector<Point> &c : grid) total += c.size();
        spiral.points.reserve(total + 2);
    }
    std::vector<size_t> ring_starts;                 // emitted start index of each ring, path order
    std::vector<char>   ring_surface;                // is the emitted ring a free surface? (path order)
    std::vector<float>  ring_base_w;                 // the ring's own (pre-scarf) width, path order
    ring_starts.reserve(order.size());
    ring_surface.reserve(order.size());
    ring_base_w.reserve(order.size());
    for (const size_t idx : order) {
        ring_starts.push_back(spiral.points.size());
        ring_surface.push_back(char(is_surface(idx) ? 1 : 0));
        ring_base_w.push_back(width_of(idx));
        const bool internal_ring = (infill_shell_loops > 0) ? is_infill(idx) : !is_surface(idx);
        const coord_t zo = ring_z_off(idx);
        std::vector<Point> &c = grid[idx];
        // Anchoring by absolute angle assumes rings are star-shaped about the centre.
        // On a concave section (a notched band, a lobe that does not enclose the
        // centre) consecutive rings' anchors can land on OPPOSITE sides of the notch,
        // turning the ring-to-ring hop into a chord across AIR. When the hop leaves
        // the region, re-anchor this ring POSITIONALLY to the point nearest the
        // previous ring's end: rings are one bead apart everywhere, so the nearest
        // point is radially adjacent and the hop stays local. Containment (not hop
        // length) is the trigger, so legitimate in-material hops - square-annulus
        // corners, through-waist crossings - keep their angle anchors and existing
        // parts' output is untouched.
        if (region != nullptr && !spiral.points.empty()) {
            const Point prev = spiral.points.back();
            if ((c.front() - prev).cast<double>().norm() > 2.0 * ring_spacing &&
                !region->contains(Line(prev, c.front()))) {
                size_t bestk = 0; double bd = std::numeric_limits<double>::max();
                for (size_t k = 0; k < c.size(); ++k) {
                    const double d = (c[k] - prev).cast<double>().squaredNorm();
                    if (d < bd) { bd = d; bestk = k; }
                }
                if (bestk > 0)
                    std::rotate(c.begin(), c.begin() + bestk, c.end());
            }
        }
        for (size_t k = 0; k < c.size(); ++k) { spiral.points.push_back(c[k]); flag(internal_ring); zoff(zo); pw(point_w(idx, c[k])); }
    }

    // In-plane scarf (option A). Every ring-to-ring junction is a one-spacing radial step
    // whose diagonal crosses track the neighbouring rings already cover - an abrupt overfill
    // ridge, the once-per-layer visible seam when it involves ring 0. Soften it by tapering
    // the bead width toward a floor over the last `scarf_len` of the outgoing ring and back
    // up over the first `scarf_len` of the incoming ring; the junction diagonal itself (the
    // segment keyed on the incoming ring's first point) is laid at the floor. WIDTH/FLOW
    // ONLY - no centreline moves - because SpiralVase's forced-on XY smoothing pulls points
    // toward the PREVIOUS layer (whose junction sits at a precessed angle) and would smear
    // any deliberate centreline shift; flow modulation passes through it untouched.
    // Tapered widths are quantized to Arachne's 0.05mm bucket so the taper becomes a handful
    // of sub-paths, not one per point. Windows never cross their own ring's far end (a short
    // ring gets a shorter taper), and overlapping windows keep the smaller width. Woven
    // chains are excluded: their inter-ring diagonals are structural struts across open
    // sparse space, not redundant track. Floor: half the ring's width, but never below the
    // caller's printable floor.
    if (scarf_len > 0 && pt_width != nullptr && infill_shell_loops == 0 && ring_starts.size() > 1) {
        auto taper = [&](size_t begin, size_t end, size_t join_end, bool towards_end, float base_w) {
            // Points [begin, end) belong to one ring; the junction sits at index join_end
            // (one past the ring's last point when towards_end, at the ring's first point
            // otherwise). Walk away from the junction accumulating XY arc length. The taper
            // is computed from the ring's ORIGINAL width, never the current point width -
            // on a short ring the head and tail windows overlap, and re-tapering an
            // already-tapered point would compound below the floor. Overlaps keep min().
            // The window is capped at a quarter of the ring's own arc so a small inner ring
            // is never thinned wall-to-wall by an aggressive scarf length: at most half of
            // any ring (a quarter from each end) is ever below full width.
            double ring_arc = 0.;
            for (size_t j = begin + 1; j < end; ++j)
                ring_arc += (spiral.points[j] - spiral.points[j - 1]).cast<double>().norm();
            double acc = 0.;
            const double L = std::min(double(scarf_len), 0.25 * ring_arc);
            if (L <= 0.)
                return;
            const float floor_w = std::max(scarf_floor, 0.5f * base_w);
            size_t i = towards_end ? (join_end - 1) : begin;
            size_t prev = i;
            while (i >= begin && i < end) {
                acc += (spiral.points[i] - spiral.points[prev]).cast<double>().norm();
                if (acc >= L)
                    break;
                float tapered = floor_w + (base_w - floor_w) * float(acc / L);
                tapered = std::round(tapered / 0.05f) * 0.05f;              // 0.05mm buckets
                tapered = std::max(tapered, floor_w);
                float &w = (*pt_width)[i];
                w = std::min(w, tapered);
                prev = i;
                if (towards_end) { if (i == begin) break; --i; } else { ++i; }
            }
        };
        // Free-surface rings are never tapered: narrowing a surface bead recedes its
        // OUTER flank and would trade the overfill ridge for a visible relief groove.
        // Their junctions still soften via the diagonal (floor width, keyed on the
        // incoming ring's first point) and the buried partner's taper.
        for (size_t k = 1; k < ring_starts.size(); ++k) {
            const size_t b = ring_starts[k];                                 // incoming ring's first point
            if (!ring_surface[k - 1])
                taper(ring_starts[k - 1], b, b, /*towards_end=*/true, ring_base_w[k - 1]);   // outgoing ring's tail
            if (!ring_surface[k]) {
                const size_t ring_end = (k + 1 < ring_starts.size()) ? ring_starts[k + 1] : spiral.points.size();
                taper(b, ring_end, b, /*towards_end=*/false, ring_base_w[k]);                // incoming ring's head
            } else if (b < pt_width->size()) {
                // Even into a surface ring, the junction DIAGONAL (the segment keyed on the
                // ring's first point) is redundant track and stays at the floor; only the
                // ring's own subsequent points keep full surface width.
                float &w = (*pt_width)[b];
                w = std::min(w, std::max(scarf_floor, 0.5f * ring_base_w[k]));
            }
        }
    }

    return spiral;
}

// Collect concentric offset bands by marching inward from `band`, appending the nested
// sequence as ONE chain to `chains`. When the region splits during inward offset - a
// SOLID waist (e.g. a dogbone) or a pinched/branching annular wall (a non-uniform tube
// wall that necks off, common at fine spacing) - it ends the current chain and recurses
// so each piece becomes its OWN chain (filled contiguously, joined by a single bridge).
// `vanished` records, parallel to `chains`, whether each chain ended because the next
// inward offset came up EMPTY (the band was consumed at the medial axis) - only then is
// the leftover a measurable medial deficit the caller can close by widening beads; a
// chain ended by a pinch or a split has material continuing elsewhere instead.
// Returns false only for topologies we still do not model (a solid band that grows a
// hole) or on a band/point budget overflow, so the caller can hard-stop.
// Record a way this layer took a weaker path than the options asked for.
//
// APPENDS rather than only filling an empty slot. Every site used to guard on the note being
// empty, so a layer that declined twice reported once - and which one it reported depended on
// the order the code happened to run in. A test written against that behaviour passed while
// asserting nothing: it saw a pre-existing note from an unrelated decline and was satisfied.
inline void add_note(std::string *note, const std::string &msg)
{
    if (note == nullptr) return;
    if (! note->empty()) *note += "; ";
    *note += msg;
}

bool march_bands(ExPolygon band, coord_t spacing, bool has_hole,
                 std::vector<std::vector<ExPolygon>> &chains, std::vector<char> &vanished,
                 std::vector<char> &outer_pinch, size_t &point_count, bool &multi_hole,
                 bool outer_real = true)
{
    constexpr size_t max_chains = 4096;
    constexpr size_t max_points = 400000;
    std::vector<ExPolygon> chain;
    bool ended_by_vanish = false;
    for (;;) {
        // A wall band must keep exactly one hole; when it pinches we hit the medial.
        if (has_hole && band.holes.size() != 1) {
            if (band.holes.empty()) {
                // ...but a wall of VARYING thickness pinches on its thinnest azimuth first, while
                // thicker azimuths still have bands to give. Stopping here abandons that material:
                // an eccentric bore (4 mm wall one side, 8 mm the other) kept only the five bands
                // each face could give before the thin side closed and left a 0.9 mm void ring up
                // the middle of the thick side - 31% of the part's material missing, silently.
                // The remainder is no longer annular, so hand it on as a SOLID region and let it
                // march to its own medial.
                if (! chain.empty()) {
                    chains.push_back(std::move(chain));
                    vanished.push_back(0);
                    outer_pinch.push_back(outer_real ? 0 : 1);
                }
                // The remainder starts where this chain pinched out - an artificial end.
                return march_bands(band, spacing, false, chains, vanished, outer_pinch,
                                   point_count, multi_hole, /*outer_real=*/false);
            }
            // Bore split into several holes: not a topology we model. This ABANDONS the material
            // the chain would have covered from here inward, so it must be said out loud - a
            // silent version of exactly this cost an eccentric bore 31% of its material with
            // nothing in the log.
            multi_hole = true;
            break;
        }
        // A solid region that grows a hole is a topology we do not model.
        if (!has_hole && !band.holes.empty())
            return false;

        chain.push_back(band);
        point_count += band.contour.points.size();
        if (!band.holes.empty())
            point_count += band.holes.front().points.size();
        if (chains.size() > max_chains || point_count > max_points)
            return false;

        ExPolygons next = offset_ex(ExPolygons{ band }, -float(spacing));
        if (next.empty()) {
            ended_by_vanish = true;
            break;                           // vanished at the medial axis
        }
        if (next.size() == 1) {              // still one piece -> keep marching
            band = std::move(next.front());
            continue;
        }
        // The chain split into several pieces (a solid waist, or a pinched annular
        // wall). End this chain and recurse into each piece as its own chain, using
        // each piece's OWN hole status (a piece may keep the bore or close up solid).
        if (!chain.empty()) {
            chains.push_back(std::move(chain));
            vanished.push_back(0);
            outer_pinch.push_back(outer_real ? 0 : 1);
        }
        for (const ExPolygon &comp : next)
            if (!march_bands(comp, spacing, !comp.holes.empty(), chains, vanished, outer_pinch,
                             point_count, multi_hole, /*outer_real=*/false))
                return false;
        return true;
    }
    if (!chain.empty()) {
        chains.push_back(std::move(chain));
        vanished.push_back(char(ended_by_vanish ? 1 : 0));
        outer_pinch.push_back(outer_real ? 0 : 1);
    }
    return true;
}

// One connected piece of a (possibly split) cross-section: its concentric loops in
// radial order (outer contours outer->medial, then, if annular, hole contours
// medial->bore) plus whether that piece is annular (so spiralize_rings knows its
// innermost loop is a bore surface, not a buried centre).
struct ContourChain {
    Polygons loops;
    bool     has_bore = false;
    // True when this chain's OUTER end is where a previous chain pinched out, not a real face
    // of the object. A wall of varying thickness closes on its thin azimuth first, and the
    // marching hands the remainder on as a new chain - so that chain's outermost rings sit in
    // the MIDDLE of the wall. Giving them the solid shell allowance printed a solid ring buried
    // mid-wall: on the eccentric test object it took wall coverage to 63% against 44% for a
    // plain wall, and was visible in the section as a solid band between two foams.
    bool     outer_pinch = false;
    bool     woven_infill = false;   // true once this chain's interior was thinned to sparse infill
    // False when the chain's rings are not star-shaped about the perturbation centre, so
    // every angular-parametrized perturbation (radial weave, solid stagger, and this
    // chain's share of the Z-wave) was/must be declined. Sparse DECIMATION is unaffected -
    // dropping rings needs no parametrization - so a declined chain still prints plain
    // concentric sparse rings.
    bool     star_ok = true;
    // Adaptive bead width (C-lite): per-ring extrusion width in mm, parallel to `loops`.
    // Empty means every ring prints at the nominal flow width. Only chains whose medial
    // deficit was closed (widened innermost pair) or whose sub-bead centreline was
    // width-fitted carry explicit widths.
    std::vector<float> ring_widths;
    // Adaptive beading (C-full): when set, `loops` are the plan's iso-index rings and
    // per-point widths come from the plan (ring_widths stays empty); the C-lite medial
    // deficit correction is skipped - the plan covers the full wall by construction.
    std::shared_ptr<ChainBeadPlan> plan;
    // The plan SLOT each loop came from, parallel to `loops`. Non-empty exactly when `plan`
    // is set. Position in `loops` used to BE the slot index, which is why sparse infill and
    // adaptive beading could not coexist: decimation drops rings, and every surviving ring
    // after the first gap would have read a neighbour's per-angle widths. Carrying the slot
    // explicitly lets the two compose - decimation compacts this in lockstep with `loops`,
    // and a surviving ring keeps its OWN slot's widths.
    std::vector<size_t> ring_slot;
};

// RADIAL-TRAVERSE INTERIOR.
//
// The concentric lattice has no radial members. Its only radial connections are the places where
// two undulating rings happen to touch, so between them lies a continuous cylindrical void - and
// a peel needs exactly one such surface. Printed de-lid parts came apart into concentric sections
// along those voids at every density tried, and three successive fixes only moved the void to
// whichever interface tied least (measured: 0.8% of one interface's circumference tied, against
// 74% two rings away). The void is a property of the generator, not of its tuning.
//
// A serpentine sweeping across the wall has the opposite property: every traverse crosses every
// ring, so a crack anywhere in the depth must sever traverses rather than merely separate. On the
// de-lid geometry at density 25 that is 24 members at every radius against 0 to 1 at the ring
// lattice's worst, for the same path length - both lay the same area of bead.
//
// It is built in RING-INDEX space, not in angle about a centroid. Angle would have confined this
// to star-shaped single-bore sections, which is a small share of what the slicer meets: the
// marching already produces the rings for concave outlines, multi-bore walls, and solid sections
// that close on a medial terminus, and a serpentine is just a path that crosses those rings
// instead of following one. So the traverse is a function of NORMALISED ARC LENGTH along the
// rings, and every topology the marching can produce is carried for free.
//
// It closes on itself after one circuit, so it enters the chain as a single ring and every piece
// of downstream machinery - chaining, junctions, the arc-length Z ramp, the seam - is unchanged.

// Equal-arc-length resample of `ring` to `S` points, starting at the point nearest `anchor` and
// wound the same way as `reference_ccw`. Correspondence between successive rings is what the
// traverse interpolates along, so both the start and the direction have to agree; without that
// the bore-side rings, which the marching stores with the opposite winding, would be traversed
// backwards and the serpentine would shear across the wall instead of crossing it.
std::vector<Point> resample_ring(const Polygon &ring, size_t S, const Point &anchor, bool reference_ccw)
{
    Polygon r = ring;
    if (r.is_counter_clockwise() != reference_ccw)
        r.reverse();
    const size_t m = r.points.size();
    std::vector<double> cum(m + 1, 0.0);
    for (size_t i = 0; i < m; ++i)
        cum[i + 1] = cum[i] + (r.points[(i + 1) % m] - r.points[i]).cast<double>().norm();
    const double total = cum[m];
    std::vector<Point> out;
    if (total < 1e-9) { out.assign(S, r.points.empty() ? anchor : r.points.front()); return out; }
    // Where along the ring the anchor sits, so all rings start from the same place on the wall.
    double best = std::numeric_limits<double>::max(), off = 0.0;
    for (size_t i = 0; i < m; ++i) {
        const double d = (r.points[i] - anchor).cast<double>().squaredNorm();
        if (d < best) { best = d; off = cum[i]; }
    }
    out.reserve(S);
    size_t seg = 0;
    for (size_t k = 0; k < S; ++k) {
        double t = off + total * double(k) / double(S);
        if (t >= total) t -= total;
        while (seg + 1 < m && cum[seg + 1] < t) ++seg;
        while (seg > 0 && cum[seg] > t) --seg;
        const double span = cum[seg + 1] - cum[seg];
        const double f = (span > 1e-9) ? (t - cum[seg]) / span : 0.0;
        const Vec2d a = r.points[seg].cast<double>();
        const Vec2d b = r.points[(seg + 1) % m].cast<double>();
        const Vec2d p = a + (b - a) * f;
        out.emplace_back(coord_t(std::llround(p.x())), coord_t(std::llround(p.y())));
    }
    return out;
}

// Widest gap the sparsifier may leave between kept rings, in bead slots.
//
// Chosen against the PLASTIC, not against the tie metric, because the two disagree and the
// metric lost. On the de-lid pair, the dense object holds together with its worst interface tied
// 25% while the sparse object parts down the middle at 28.8% - so tie fraction does not predict
// peel resistance. What separates them is how WIDE the void is at an untied azimuth: 2.0 slots
// dense against 4.0 sparse. A wider void is both a bigger crack starter and less material to
// tear. Capping at 3.0 narrows the sparse case to 3.0 slots for about 5% more material; 2.5
// would match the dense geometry exactly but abolishes the density setting for walls this thin.
//
// Note this LOWERS the tie fraction (28.8% -> 19.1% on the sparse object). That is accepted
// deliberately rather than overlooked.
static const double CS_MAX_GAP_SLOTS = 3.0;

bool apply_radial_serpentine(ContourChain &cc, const WovenInfill &woven, coord_t spacing,
                             size_t layer_index, std::string *degrade_note)
{
    const size_t sl = woven.shell_loops;
    const size_t n  = cc.loops.size();
    const bool annular = cc.has_bore;
    // The buried rings, exactly as the sparsifier defines them: a solid skin of shell_loops on
    // every free face, and for a solid section the medial terminus stays put as well.
    const size_t first = sl;
    const size_t last  = annular ? (n - sl - 1) : (n - 2);
    if (n <= 2 * sl || last < first + 1) {
        add_note(degrade_note, "serpentine declined: fewer than two buried rings to traverse");
        return false;
    }
    const size_t M = last - first;                  // rings the traverse sweeps across

    // Traverses per circuit, from the requested density. Each traverse lays about one band
    // thickness of bead and there are 2N of them, so the fill fraction is 2*N*spacing/perimeter
    // whatever the band depth. Checked against the geometric prototype: 12.3 predicted, 12
    // delivering the requested 25%.
    const double perim = cc.loops[first + M / 2].length() * 1.0;
    const double density = 1.0 / double(std::max<size_t>(1, woven.stride));
    const size_t N = std::max<size_t>(3, size_t(std::lround(density * perim / (2.0 * double(spacing)))));

    // Sample densely enough that the reversal at each face is a corner and not a chord.
    const size_t per_traverse = std::max<size_t>(8, size_t(double(M) * 4.0));
    const size_t S = std::max<size_t>(2 * N * per_traverse, 256);

    const bool ccw = cc.loops[first].is_counter_clockwise();
    const Point anchor = cc.loops[first].points.empty() ? Point(0, 0) : cc.loops[first].points.front();
    std::vector<std::vector<Point>> band;
    band.reserve(M + 1);
    for (size_t j = first; j <= last; ++j) {
        band.push_back(resample_ring(cc.loops[j], S, anchor, ccw));
        if (band.back().size() != S) {
            add_note(degrade_note, "serpentine declined: a buried ring would not resample");
            return false;
        }
    }

    // The traverse interpolates between CORRESPONDING points of successive rings, which is only
    // meaningful while the chain really is a linear stack of nested rings. It is not always: a
    // section with more than one bore has three faces, not two, and a band that pinches or splits
    // breaks the nesting. Those are decomposed or declined upstream, but "upstream handles it" is
    // exactly the kind of assumption this subsystem has been caught on before, so check it here
    // rather than trust it. Corresponding points on adjacent rings sit about one ring gap apart;
    // if any pair is far outside that, the correspondence is not a crossing of the wall and the
    // traverse would cut across the section instead. Decline, loudly.
    double worst_step = 0.0;
    for (size_t k = 0; k + 1 < band.size(); ++k)
        for (size_t i = 0; i < S; ++i)
            worst_step = std::max(worst_step, (band[k + 1][i] - band[k][i]).cast<double>().norm());
    if (worst_step > 4.0 * double(spacing)) {
        add_note(degrade_note, "serpentine declined: rings do not correspond (chain is not a "
                               "simple nested stack - multi-bore or pinched section)");
        return false;
    }

    const double phase = double(layer_index) * woven.dphase;   // successive layers cross
    Polygon serp;
    serp.points.reserve(S);
    for (size_t i = 0; i < S; ++i) {
        const double u = double(i) / double(S);
        // Triangle wave over the ring stack: constant crossing speed, so the bead is spread
        // evenly through the depth instead of dwelling at the faces the way a sine would.
        const double x   = double(N) * u + phase / (2.0 * M_PI);
        const double tri = 4.0 * std::abs(x - std::floor(x + 0.5)) - 1.0;      // [-1, 1]
        const double f   = 0.5 * (tri + 1.0) * double(M);                      // [0, M]
        const size_t k   = std::min(size_t(f), M - 1);
        const double t   = f - double(k);
        const Vec2d a = band[k][i].cast<double>();
        const Vec2d b = band[k + 1][i].cast<double>();
        const Vec2d p = a + (b - a) * t;
        serp.points.emplace_back(coord_t(std::llround(p.x())), coord_t(std::llround(p.y())));
    }

    // [outer shells | serpentine | inner shells or terminus]. The slot map and any per-ring
    // widths are compacted in lockstep, as everywhere else: the serpentine spans every buried
    // slot, so it takes the outermost buried slot's entry rather than inventing one.
    Polygons kept;
    std::vector<size_t> kept_slot;
    std::vector<float>  kept_w;
    auto take = [&](size_t j) {
        kept.push_back(std::move(cc.loops[j]));
        if (!cc.ring_slot.empty())   kept_slot.push_back(cc.ring_slot[j]);
        if (!cc.ring_widths.empty()) kept_w.push_back(cc.ring_widths[j]);
    };
    for (size_t j = 0; j < first; ++j) take(j);
    kept.push_back(std::move(serp));
    if (!cc.ring_slot.empty())   kept_slot.push_back(cc.ring_slot[first]);
    if (!cc.ring_widths.empty()) kept_w.push_back(cc.ring_widths[first]);
    for (size_t j = last + 1; j < n; ++j) take(j);
    cc.loops = std::move(kept);
    if (!cc.ring_slot.empty())   cc.ring_slot   = std::move(kept_slot);
    if (!cc.ring_widths.empty()) cc.ring_widths = std::move(kept_w);
    cc.woven_infill = true;
    return true;
}

// One definition of "this chain has buried rings to thin", used by the sparsifier and by the
// C-lite deficit correction that must stand aside for it. They were two separate expressions
// and had to agree; now they cannot drift.
// Deliberately counts a full outer shell even where the outer end is a PINCH and the sparsifier
// will hold no shell there (sl_out). Correcting the arithmetic was tried and reverted: it makes a
// three-ring pinch chain eligible, and such a chain cannot be thinned anyway - two infill rings
// and a terminus, with the distribution obliged to anchor both ends - so the only effect is to
// start weaving it. Measured on the eccentric object's deep pocket, that turned its outer ring
// from a 41 mm loop into an 80 mm one sweeping 4.2 mm of radius, because the ring it weaves
// against is half its length and the gap to it is not a gap across the wall at all. The pocket
// still printed solid, and the "too thin to sparsify" note stopped being emitted for it.
// How far ring A is from ring B, at the given percentile of A's length.
//
// Consecutive rings of a chain are one offset apart everywhere, so this reads one spacing -
// EXCEPT where the offsets are not parallel, which is the whole reason it exists.
//
// Two things it must not do, both of which were done and measured wrong first:
//
//  * Sample A at its VERTICES. A marched ring is a raw offset polygon and its vertices sit
//    wherever the geometry put them - a square has four and all four are corners, where the
//    distance to the next offset is sqrt(2) times the distance along a face. Sample by arc length.
//
//  * Measure to B's nearest VERTEX. That is what a point KD-tree gives, and on a coarse ring it
//    is not the distance to the ring at all: from the middle of a square's edge the nearest
//    vertex of the next square is its corner, half a side away. It reported a 16 mm square's
//    adjacent offsets 3.458 mm apart when they are 0.446, and the error shrank with the ring, so
//    it looked like a plausible ramp rather than a bug. Measure to B's SEGMENTS.
double ring_gap(const Polygon &a, const Polygon &b, double pct = 0.5)
{
    if (a.points.size() < 2 || b.points.size() < 2)
        return 0.0;
    const AABBTreeLines::LinesDistancer<Line> to_b{ b.lines() };
    const size_t m = a.points.size();
    std::vector<double> cum(m + 1, 0.0);
    for (size_t i = 0; i < m; ++i)
        cum[i + 1] = cum[i] + (a.points[(i + 1) % m] - a.points[i]).cast<double>().norm();
    const double total = cum[m];
    if (total <= 0.0)
        return 0.0;
    constexpr size_t NS = 64;
    std::vector<double> d;
    d.reserve(NS);
    size_t seg = 0;
    for (size_t k = 0; k < NS; ++k) {
        const double at = total * double(k) / double(NS);
        while (seg + 1 < m && cum[seg + 1] < at) ++seg;
        const double len = cum[seg + 1] - cum[seg];
        const double t   = len > 0.0 ? (at - cum[seg]) / len : 0.0;
        const Vec2d p0 = a.points[seg].cast<double>();
        const Vec2d p1 = a.points[(seg + 1) % m].cast<double>();
        const Vec2d q  = p0 + (p1 - p0) * t;
        d.push_back(std::abs(to_b.distance_from_lines<false>(
            Point(coord_t(std::llround(q.x())), coord_t(std::llround(q.y()))))));
    }
    if (d.empty())
        return 0.0;
    const size_t k = std::min(d.size() - 1, size_t(pct * double(d.size())));
    std::nth_element(d.begin(), d.begin() + k, d.end());
    return d[k];
}

bool chain_thick_enough(const ContourChain &cc, const WovenInfill &w)
{
    return cc.has_bore ? cc.loops.size() > 2 * w.shell_loops + 1
                       : cc.loops.size() > w.shell_loops + 1;
}

bool chain_will_decimate(const ContourChain &cc, const WovenInfill &w)
{
    return w.enabled && w.stride > 1 && cc.loops.size() > 2 && chain_thick_enough(cc, w);
}

// The plan indexes its per-angle width and position tables by SLOT. Nothing downstream can
// notice a slot map that has drifted out of step with `loops` - it would simply read a
// neighbour's widths - and ChainBeadPlan::interp indexes its tables unchecked, so an
// off-by-one is an out-of-bounds read rather than a wrong bead. Check it explicitly.
bool validate_chain(const ContourChain &cc, std::string *fail_reason)
{
    auto fail = [&](const char *why) {
        if (fail_reason)
            *fail_reason = std::string("Continuous shell adaptive beading: ") + why;
        return false;
    };
    if (!cc.ring_widths.empty() && cc.ring_widths.size() != cc.loops.size())
        return fail("per-ring width array is out of step with the rings.");
    if (cc.plan == nullptr)
        return cc.ring_slot.empty() ? true : fail("a chain without a plan carries a slot map.");
    if (cc.ring_slot.size() != cc.loops.size())
        return fail("plan ring/slot map is out of step with the rings.");
    if (cc.ring_slot.empty())
        return fail("plan chain carries no slot map.");
    for (size_t i = 0; i < cc.ring_slot.size(); ++i) {
        if (cc.ring_slot[i] >= cc.plan->n_slots)
            return fail("plan slot index out of range.");
        if (i > 0 && cc.ring_slot[i] <= cc.ring_slot[i - 1])
            return fail("plan slot map is not strictly increasing.");
    }
    // Both surfaces are pinned and must survive any decimation: the ring form keeps the full
    // range, and the spiral form is truncated to exactly those two.
    if (cc.ring_slot.front() != 0 || cc.ring_slot.back() != cc.plan->n_slots - 1)
        return fail("plan chain lost one of its pinned surface rings.");
    return true;
}

// Decompose a single connected cross-section into the concentric loop sets that fill
// it as Fermat-style in/out spirals - ONE ContourChain per connected piece.
//
// For an annular wall the path must cover the whole wall thickness and leave the real
// bore hollow. We march inward one bead-spacing at a time. Each annular band has an
// outer contour and an inner hole contour; emitting all outer contours (outer edge ->
// medial axis) then all hole contours (medial axis -> bore) yields loops that spiral
// in, U-turn at the medial axis, and spiral back out to the bore. A solid cross-section
// (a cap) simply spirals inward to the centre. Where the region splits during marching
// (a dogbone waist, or a pinched non-uniform tube wall), each piece becomes its own
// chain so it is filled contiguously instead of interleaved with the others.
// Largest inset depth (in (0, hi]) at which `shape` still has material: bisection on
// offset_ex. Used to measure the residual half-thickness of a vanished band (the medial
// deficit) and the thickness of a sub-bead region. `hi` must be a known-empty depth.
coord_t vanish_depth(const ExPolygons &shape, coord_t hi)
{
    coord_t lo = 0;
    for (int it = 0; it < 4; ++it) {
        const coord_t mid = lo + (hi - lo) / 2;
        if (mid == lo)
            break;
        if (offset_ex(shape, -float(mid)).empty())
            hi = mid;
        else
            lo = mid;
    }
    return lo + (hi - lo) / 2;
}

// Sparse woven infill applied to ONE chain. Extracted from the marching loop so the
// PLAN path can call it too: before that, switching on continuous_shell_infill_density
// silently disabled adaptive beading for the whole object, because the plan branch
// returned its chain before ever reaching this code.
void apply_woven_infill(ContourChain &cc, const WovenInfill &woven, coord_t spacing,
                        size_t layer_index, std::string *degrade_note)
{
    if (!woven.enabled || cc.loops.size() <= 2)
        return;
    // The radial serpentine replaces the whole buried interior, so it is decided before any of
    // the ring-lattice machinery runs. It declines rather than degrades - a section without a
    // bore, or too thin to sweep, falls through to the ring lattice below with a note.
    if (woven.serpentine && woven.stride > 1 &&
        apply_radial_serpentine(cc, woven, spacing, layer_index, degrade_note))
        return;
    // The note is deferred to the OUTCOME, below: "no ring was dropped" is the thing that is
    // true of the plastic, where "chain_thick_enough said no" is a proxy for it. They part
    // company on a short chain whose span leaves the distribution no room to drop anything -
    // that chain prints solid and the proxy says nothing about it.
    const size_t rings_before = cc.loops.size();
    // Sparse WOVEN infill. Keep a solid shell of shell_loops rings on every free FACE and thin
    // the buried interior to the sparse stride, all carried by the SAME continuous bead. For a
    // SOLID chain the shell hugs the outer face and the fill thins toward the centre. For an
    // ANNULAR chain (a thick tube wall) the loop list runs outer-face -> medial -> bore-face, so
    // keep shell_loops at BOTH ends (a solid skin on the OD and the ID) and thin only the deep
    // middle - otherwise a very thick wall (e.g. OD100/ID10) would print 100% solid. With weave>0
    // each thinned ring is perturbed radially by A*sin(k*theta + phase), phase advancing per layer,
    // so a fixed angle bulges out on one layer and in on the next - diagonal inter-layer struts (a
    // gyroid-like 3D interlock, completed by the Z-wave in build_cross_section) instead of nested
    // tubes. Disabled unless infill density is set.
    if (woven.enabled && cc.loops.size() > 2) {
        const size_t sl = woven.shell_loops;
        // A shell exists to give a FACE a solid skin. Where this chain's outer end is not a face
        // but the point at which the previous chain pinched out, there is no surface to skin -
        // those rings are buried mid-wall, and skinning them printed a solid ring inside the
        // part. Measured on the eccentric object: 63% wall coverage against 44% for a plain
        // wall, with a run of one-slot gaps mid-wall where two chains' shells sat back to back.
        const size_t sl_out = cc.outer_pinch ? 0 : sl;
        const size_t K  = woven.stride;
        const size_t n  = cc.loops.size();
        const bool   annular = cc.has_bore;
        // is_infill on ORIGINAL indices: buried interior beyond the shell(s). Annular keeps a
        // shell at both ends; solid keeps the outer shell and the buried core (minus the terminus).
        // A ring that carries part of the solid mask is never thinned away, wherever it sits in
        // the stack. Granularity is a WHOLE RING, which is the natural unit here: the rings are
        // closed loops and the machinery downstream chains closed loops, so keeping half of one
        // would need open paths through the whole pipeline. It also happens to be the right unit
        // for the case that motivates this - a bore floor needs its INNER rings solid, and those
        // are whole rings.
        auto touches_mask = [&](size_t i) -> bool {
            if (woven.solid_mask == nullptr || woven.solid_mask->empty()) return false;
            const Points &pts = cc.loops[i].points;
            // Every eighth point: the mask is grown by a wall thickness before it gets here, so a
            // ring either sits well inside it or well outside, and sampling cannot straddle that.
            for (size_t k = 0; k < pts.size(); k += 8)
                for (const ExPolygon &m : *woven.solid_mask)
                    if (m.contains(pts[k])) return true;
            return false;
        };
        auto is_infill = [&](size_t i) -> bool {
            if (i < sl_out) return false;
            if (touches_mask(i)) return false;          // must stay solid here
            return annular ? (i + sl < n) : (i + 1 < n);
        };
        const bool thick = chain_thick_enough(cc, woven);   // one definition, shared
        // Perturbation centre = centroid of the outer boundary loop (shared by the
        // weave and the solid stagger below).
        Vec2d centre(0., 0.);
        if (!cc.loops.front().points.empty()) {
            for (const Point &p : cc.loops.front().points) centre += p.cast<double>();
            centre /= double(cc.loops.front().points.size());
        }
        // Angular-parametrization gate (see star_shaped_about): on a folded section
        // decline every perturbation that phases or clamps by angle-about-centre.
        // Checked over the FULL pre-decimation ring set (a superset of the kept
        // rings), so the kept rings are covered too. Decimation itself stays on.
        if (thick)
            for (const Polygon &ring : cc.loops)
                if (!star_shaped_about(ring, centre)) { cc.star_ok = false; break; }
        if (K > 1 && thick) {
            Polygons kept; std::vector<char> kept_infill;
            std::vector<size_t> kept_slot; std::vector<float> kept_w;
            kept.reserve(n); kept_infill.reserve(n);
            // WHICH buried rings survive. A fixed stride anchors the pattern at the outer shell
            // and dumps the whole remainder at the bore: printed on a 3.17 mm wall, the infill
            // tied the outer shell at 0.508 mm and stood 1.265 mm off the bore shell, so the
            // lattice joined one face and floated free of the other. Distribute the kept rings
            // across the buried span INCLUSIVE OF BOTH ENDS instead, so each shell gets an
            // adjacent ring and the leftover is shared through the middle.
            std::vector<char> keep(n, 0);
            // Rings that bound a MEDIAL POCKET are kept whatever the density asks, and pinned.
            //
            // A chain's ring list runs down one face to the medial and back up the other, so its
            // two middle rings are the last the marching could reach from either side. On a wall
            // of even thickness they meet one slot apart. On a wall of VARYING thickness they do
            // not: the marching closes on the thin azimuth first and hands the remainder on as a
            // separate chain, which fills the space between them. Those two rings are that
            // pocket's faces, and the chain filling it has nothing else to weld to.
            //
            // Nothing here could see that. The gap cap looks only WITHIN a chain, and so does the
            // weave, so the two chains' kept rings drifted four slots apart with a void between
            // them. Measured on the eccentric object: its middle pair is 7.41 slots apart and the
            // next chain's outermost ring sits exactly one slot inside EACH of them - yet the
            // stride stepped over both. The printed part came apart into three concentric pieces.
            //
            // The threshold is not chosen so much as read off: on every chain of every object in
            // the battery the consecutive gap is 0.70 to 1.11 slots, and the one pocket is 7.41.
            std::vector<char> pocket(n, 0);
            for (size_t i = 0; i + 1 < n; ++i)
                if (ring_gap(cc.loops[i], cc.loops[i + 1]) > 1.5 * double(spacing))
                    pocket[i] = pocket[i + 1] = 1;
            {
                // The span starts at whatever the outer end actually HOLDS: a face holds sl shell
                // rings, a pinch holds none. Starting at sl on a pinch chain spent that chain's
                // outermost rings - the ones meeting the chain that pinched out - as if they were
                // skin, and left the eccentric object's pocket three CONSECUTIVE rings, printed as
                // the solid crescent visible between its two walls.
                const size_t first = sl_out;
                const size_t last  = annular ? (n - sl - 1) : (n - 2);
                // The pocket's faces are forced into the keep set, ON TOP of one distribution
                // across the whole span. Splitting the span at the pocket and anchoring a run
                // either side is the tidier idea and was tried: it costs the same. A run bounded
                // by two anchors needs a free ring between them to do the weaving - both its ends
                // are pinned, and two pinned rings three slots apart weld over 0% of their length
                // - so the run has to carry three rings, which is the ring the forcing added. It
                // also came out worse where the arithmetic landed on two: chain 0 fell into four
                // detached pieces, welded only across the pocket.
                if (n > sl_out + sl && last >= first) {
                    const size_t span = last - first + 1;
                    // At least two where the span allows, so BOTH shells are tied; that floor is
                    // also what stops a thin interior degenerating into one ring down the middle.
                    //
                    // Pick the count that delivers the requested STRIDE, not the count that
                    // divides the span. Both ends are anchored, so m rings make m-1 gaps across
                    // span-1 slots and the delivered stride is (span-1)/(m-1) - which equals K
                    // only when the span happens to divide. `(span + K/2)/K` counts rings as if
                    // they were unanchored, and on the de-lid part shipped 4.67 slots against a
                    // requested 4: about 17% sparser than the density asked for, and wide enough
                    // that once the shell anchors are pinned there is not enough left to build a
                    // lattice. The span is small, so search it rather than juggle rounding.
                    //
                    // A gap also has a WIDTH LIMIT, independent of the requested density. Two
                    // rings tie by reaching for each other, so a wide gap asks each of them for a
                    // large excursion - and where the neighbour is a pinned anchor, it asks ONE
                    // ring for the whole thing. Past a point that swing dominates the ring's
                    // motion and starves its other interface. Measured on the de-lid part: 2.02
                    // slot gaps tie 25-52%, while a wall whose gaps came out 5.05/3.65/5.05 tied
                    // its middle interface at 0.0% - a clean parting plane down the wall, exactly
                    // where the printed part came apart. So the count is raised until no gap
                    // exceeds the limit, even where the requested density would allow wider.
                    // Denser walls are unaffected: their natural stride is already inside it.
                    size_t m = 1;
                    if (span > 1) {
                        auto max_gap = [&](size_t cand) {
                            size_t prev = 0, worst = 0;
                            for (size_t j = 0; j < cand; ++j) {
                                const size_t at = (cand == 1 ? (span - 1) / 2
                                                             : (j * (span - 1) + (cand - 1) / 2) / (cand - 1));
                                if (j) worst = std::max(worst, at - prev);
                                prev = at;
                            }
                            return worst;
                        };
                        double best = std::numeric_limits<double>::max();
                        for (size_t cand = 2; cand <= span; ++cand) {
                            if (double(max_gap(cand)) > CS_MAX_GAP_SLOTS && cand < span) continue;
                            const double err = std::abs(double(span - 1) / double(cand - 1) - double(K));
                            if (err < best - 1e-9) { best = err; m = cand; }   // ties keep the sparser
                        }
                        if (m == 1) m = span;          // nothing fits the limit: take the densest
                    }
                    for (size_t j = 0; j < m; ++j)
                        keep[first + (m == 1 ? (span - 1) / 2
                                             : (j * (span - 1) + (m - 1) / 2) / (m - 1))] = 1;
                }
            }
            for (size_t i = 0; i < n; ++i)
                if (pocket[i] != 0) keep[i] = 1;

            // Enforcing the gap cap in MILLIMETRES rather than slot index was tried here and
            // backed out on its own measurement. The reasoning was sound and the diagnosis it came
            // from still stands: a slot is a fixed distance only where consecutive offsets are
            // parallel, and on the eccentric object the real gap between kept rings runs to 8.4 mm
            // where every object that prints well tops out at 1.33. Inserting the midpoint ring
            // wherever the real gap exceeded the cap did narrow them - worst 90th-percentile gap
            // 6.96 mm to 4.50 - and cost 17% more material. It did not fix what it was aimed at:
            //
            //     eccentric        material  void/layer  through-void  moves/layer  channel
            //     clamp only        204.4       107.5      42.2 (39%)     17.2%     1.79 mm2
            //     + millimetre cap  238.7        74.1      32.8 (44%)     10.4%     1.42 mm2
            //
            // The void shrank because it was FILLED IN, not because it moved: as a fraction of the
            // void, the part that persists through ten layers went UP, and the lattice ended up
            // half as mobile again. That is the foam being traded for solid, which is the opposite
            // of the point. Wide gaps are therefore not the mechanism of the stacking - mobility
            // is - and the next thing to try is to move the voids rather than fill them.
            std::vector<char> kept_pocket;
            kept_pocket.reserve(n);
            for (size_t i = 0; i < n; ++i) {
                const bool inf = is_infill(i);
                if (!inf || keep[i]) {             // keep every shell/terminus ring + the chosen infill
                    kept_pocket.push_back(pocket[i]);
                    kept_infill.push_back(inf ? char(1) : char(0));
                    // The slot map and any per-ring widths are compacted in LOCKSTEP with the
                    // rings. Without this a surviving ring would read the widths of whatever
                    // slot happened to land at its new position - the reason adaptive beading
                    // and sparse infill could not previously coexist.
                    if (!cc.ring_slot.empty())   kept_slot.push_back(cc.ring_slot[i]);
                    if (!cc.ring_widths.empty()) kept_w.push_back(cc.ring_widths[i]);
                    kept.push_back(std::move(cc.loops[i]));
                }
            }
            cc.loops = std::move(kept);
            if (!cc.ring_slot.empty())   cc.ring_slot   = std::move(kept_slot);
            if (!cc.ring_widths.empty()) cc.ring_widths = std::move(kept_w);
            cc.woven_infill = true;   // spiralize_rings should tag [shell | infill | shell] accordingly

            if (woven.weave > 0.0) {
                // THE WEAVE, LOCAL TO EACH POINT.
                //
                // Two earlier formulations both failed, and for the same underlying reason - they
                // needed a global frame for something that is purely local.
                //
                // By ANGLE about the section centroid: a ring is only single-valued in angle if
                // the centroid sees all of it, so the whole weave had to be gated on `star_ok`.
                // Sections failing that got NO weave at all - thinned to the sparse stride but
                // never displaced, leaving plain nested tubes with nothing tying them together,
                // which is this subsystem's original failure mode. It sat undetected in the
                // non-star path for the entire campaign because every plate until PlateP was
                // star-shaped; PlateP's crescent was the first section to take it, and the one
                // object on that plate that peeled.
                //
                // By normalised ARC LENGTH along each ring: no centroid needed, but it assumes
                // equal fractions of arc length on adjacent rings lie across from each other.
                // On a wall of varying thickness they do not - the outer ring is much longer than
                // the inner one - so the correspondence check tripped and the ECCENTRIC object
                // stopped being woven too.
                //
                // What a point actually needs is its distance to the neighbouring ring and the
                // direction to it. Both are local: take the NEAREST point on that ring. No
                // parameterisation, no correspondence, no centroid, no star-shape assumption, and
                // nothing to gate on.
                //
                // It has to be the nearest point, though, and not a walk toward one. This used to
                // hill-climb a single index forwards and backwards from the previous sample's
                // answer, on the reasoning that adjacent points have adjacent neighbours. That
                // holds until the section is concave: a ring wraps AROUND a notch, so a point
                // facing it has one candidate across the gap and another round the far side of
                // the notch, and the distance function has two basins. Once the index settles in
                // the wrong one every later sample inherits it - the walk cannot cross the ridge
                // between them, and the step guard cannot either.
                //
                // Measured on the crescent, every ring, over its concave stretch: the walk was
                // wrong on 100% of samples, reporting a median gap of 5.206 mm where the true one
                // is 1.189. The gap feeds the amplitude law directly (cap = weave * share *
                // (gap - 0.8 * spacing)), so it permitted an excursion of 7.458 mm against a true
                // 0.873 - and the crescent's concave rings moved 7.384. The same index also picks
                // the DIRECTION of the excursion, and carries the phase field from ring to ring,
                // so one wrong answer spoiled the amplitude, the direction and the phase at once.
                // That is the crescent's concave side, which four measurements and five
                // hypotheses had failed to account for.
                //
                // Convex geometry is untouched: there the walk already agreed with the true
                // nearest on over 98% of samples, and the remainder differed by under 0.02 mm.
                auto nearest_on = [](const Points &from, const Points &to) {
                    std::vector<size_t> idx(from.size(), 0);
                    if (from.empty() || to.empty()) return idx;
                    auto coord = [&to](size_t i, size_t dim) -> double {
                        return dim == 0 ? double(to[i].x()) : double(to[i].y());
                    };
                    KDTreeIndirect<2, double, decltype(coord)> tree(coord, to.size());
                    for (size_t i = 0; i < from.size(); ++i) {
                        const std::array<double, 2> q{ double(from[i].x()), double(from[i].y()) };
                        const size_t k = find_closest_point(tree, q);
                        idx[i] = (k == decltype(tree)::npos) ? 0 : k;
                    }
                    return idx;
                };

                const size_t NR = cc.loops.size();
                // PIN the ring that hugs a shell: it stays concentric and welded to that shell the
                // whole way round, anchoring the lattice, and its neighbour covers the whole of
                // their gap instead of half. Without it the anchor's swing toward the interior
                // carried it off the shell for half the circumference, and the printed part came
                // apart into concentric sections along exactly that surface.
                std::vector<char> pinned(NR, 0);
                for (size_t j = 0; j < NR; ++j) {
                    if (!kept_infill[j]) continue;
                    if ((j > 0 && !kept_infill[j - 1]) || (j + 1 < NR && !kept_infill[j + 1]))
                        pinned[j] = 1;
                    // A pocket face anchors ANOTHER chain, and the reasoning that pins a ring
                    // against a shell applies unchanged: if it weaves it swings away from that
                    // chain's ring over half the circumference, which is the failure the shell
                    // pin was introduced to fix.
                    if (j < kept_pocket.size() && kept_pocket[j] != 0)
                        pinned[j] = 1;
                    // Pinning the ring on the OTHER side of that boundary as well - the pinch
                    // chain's own outermost ring, which hugs foreign material exactly as an
                    // anchor hugs a shell - was tried and reverted. It raised the boundary weld
                    // from 69% to 88% and then detached the pocket chain internally: its second
                    // ring is already pinned against its terminus, so pinning the first left two
                    // pinned rings three slots apart with nothing free between them to reach.
                    // The pin needs a free neighbour to do the reaching, and a three-ring chain
                    // has none to spare.
                }

                const double overshoot   = 4.0;
                const double drift_slots = 2.0;
                size_t n_infill = 0;
                for (const char c : kept_infill) n_infill += (c != 0) ? 1 : 0;
                const double drift = std::min(2.0 * M_PI / double(std::max<size_t>(2, n_infill)),
                                              drift_slots * 2.0 * std::asin(1.0 / overshoot));
                const double layer_phase = double(layer_index) * woven.dphase;

                // The rings are read from a snapshot so a ring already displaced this layer does
                // not shift the gaps its neighbour is about to measure.
                const Polygons before = cc.loops;

                // PHASE IS A FIELD ACROSS THE WALL, not a per-ring function of its own length.
                //
                // Phasing each ring by its OWN normalised arc length gives every ring the same
                // lobe COUNT. On an annulus that is harmless - all the rings are nearly the same
                // length, so anti-phase holds the whole way round. On a section whose rings differ
                // in length it is not: the inner rings of a crescent are far shorter than its
                // outer ones, so the same count means different wavelengths, and adjacent rings
                // BEAT against each other - in opposition in some places and in phase in others.
                // Two rings only tie where they oppose, so the ties collapse wherever the beat
                // puts them in phase. That is why the crescent still peeled on its concave side
                // after the weave was made to apply there at all: it was applying, but detuned.
                //
                // So the phase is carried ACROSS the wall instead. Ring 0 takes its phase from its
                // own arc length; every ring after it inherits the phase of its NEAREST point on
                // the ring outside it. Adjacent rings are then in exact anti-phase everywhere by
                // construction, whatever their lengths, and the lobes line up along the wall
                // rather than around each ring independently. On a circle the inherited phase is
                // the same angular position, so this reduces exactly to what annuli already do -
                // which is why their validated behaviour is unchanged.
                std::vector<std::vector<double>> basephase(NR);
                for (size_t j = 0; j < NR; ++j) {
                    const Points &pts = before[j].points;
                    basephase[j].assign(pts.size(), 0.0);
                    if (pts.empty()) continue;
                    if (j == 0) {
                        double run = 0.0;
                        std::vector<double> cum(pts.size() + 1, 0.0);
                        for (size_t i = 0; i < pts.size(); ++i)
                            cum[i + 1] = cum[i] + (pts[(i + 1) % pts.size()] - pts[i]).cast<double>().norm();
                        const double tot = cum[pts.size()];
                        for (size_t i = 0; i < pts.size(); ++i)
                            basephase[j][i] = (tot > 1e-9)
                                ? woven.lobes * 2.0 * M_PI * (cum[i] / tot) : 0.0;
                        (void)run;
                    } else {
                        const std::vector<size_t> up = nearest_on(pts, before[j - 1].points);
                        for (size_t i = 0; i < pts.size(); ++i)
                            basephase[j][i] = basephase[j - 1][up[i]];
                    }
                }

                size_t infill_ring = 0;
                for (size_t j = 0; j < NR; ++j) {
                    if (!kept_infill[j] || pinned[j]) continue;
                    const double ring_phase = layer_phase + double(infill_ring) * (M_PI + drift);
                    ++infill_ring;
                    const Points &me = before[j].points;
                    if (me.size() < 8) continue;
                    const std::vector<size_t> io = (j > 0)      ? nearest_on(me, before[j - 1].points)
                                                                : std::vector<size_t>();
                    const std::vector<size_t> ii = (j + 1 < NR) ? nearest_on(me, before[j + 1].points)
                                                                : std::vector<size_t>();
                    // A gap wider than a wall crossing CLAMPS the excursion. It does not cancel
                    // it, and the difference between those two is the whole point.
                    //
                    // The nearest point on the next ring is the right thing to weave toward while
                    // the two rings FACE each other. Where the next ring ends, or is much shorter
                    // than this one, the nearest point on it is round its end - along the wall,
                    // not across it - and the law then sizes an excursion from a distance that
                    // crosses no material at all. Worse, another chain's rings can be sitting in
                    // that space: band marching splits a wall into chains, and a chain measures
                    // gaps only within itself, so it reads a void where its neighbour chain's
                    // material is.
                    //
                    // Measured on the eccentric object, whose wall marches into two chains: the
                    // first chain's rings 3 and 4 sit a median 3.716 mm apart with the SECOND
                    // chain's three rings in between them. The law read that as the gap to close,
                    // permitted 5.063 mm of excursion, and flung both rings 2.53 mm into the
                    // pocket the other chain had already filled. That is the solid chunk that
                    // survived the pinch-boundary fix.
                    //
                    // Rejecting such a neighbour was the first fix and it was the wrong
                    // instrument. The weave does two jobs: in plane it welds rings to each other,
                    // and in Z it is the ONLY thing that makes one layer's lattice differ from the
                    // next, so that voids precess instead of stacking into channels. A reject
                    // gives up both - no weld AND no motion - which is the worst of the three
                    // outcomes. Measured across the battery as the fraction of ring length the
                    // weave actually moves:
                    //
                    //   P6_solid 44.0%  P1_sparse 39.1%  P4_lobed 38.6%    all came out well
                    //   P5_crescent 30.5%  P7_taper 24.7%
                    //   P3_ecc 17.5%                                      the one that did not
                    //
                    // The eccentric is the least mobile lattice in the suite and the only one
                    // whose voids stack: 47% of its void area is present on all ten consecutive
                    // layers against the crescent's 23%, and it changes 14.4% per layer against
                    // 26.1%. Its topology already freezes seven rings outright (two faces, a bore,
                    // a pocket's two sides, a terminus), so its few free rings have to move MORE
                    // than another object's, not less - and the reject had zeroed the amplitude
                    // over 47% of the one free ring in its pocket chain.
                    //
                    // Clamping keeps the motion and still bounds the excursion at what a legal gap
                    // would earn. The bound is what the sparsifier guarantees on the other side: it
                    // raises the kept-ring count until no gap exceeds CS_MAX_GAP_SLOTS, plus the
                    // half-slot of margin the reject carried, so every point that weaves today is
                    // byte-identical afterwards and the change only gives motion back where there
                    // was none. Per POINT, not per ring: a ring can face its neighbour over most of
                    // its length and run past its end for the rest.
                    //
                    // This is a PART of the fix, not all of it: it takes the eccentric from 17.5%
                    // to 20.9%, still short of the objects that print well. The rest is that the
                    // gap cap counts SLOTS, and a slot is a fixed distance only where consecutive
                    // offsets are parallel - see the handoff.
                    const double reach = (CS_MAX_GAP_SLOTS + 0.5) * double(spacing);
                    Points moved(me.size());
                    for (size_t i = 0; i < me.size(); ++i) {
                        const Vec2d here = me[i].cast<double>();
                        Vec2d out(0., 0.);
                        double g_out = 0.0, g_in = 0.0;
                        if (!io.empty()) {
                            const Vec2d d = before[j - 1].points[io[i]].cast<double>() - here;
                            const double g = d.norm();
                            if (g > 1e-9) { g_out = std::min(g, reach); out += d / g; }
                        }
                        if (!ii.empty()) {
                            const Vec2d d = here - before[j + 1].points[ii[i]].cast<double>();
                            const double g = d.norm();
                            if (g > 1e-9) { g_in = std::min(g, reach); out += d / g; }
                        }
                        const double len = out.norm();
                        if (len < 1e-9) { moved[i] = me[i]; continue; }
                        out /= len;
                        const double sh_out = (j > 0 && kept_infill[j - 1] && !pinned[j - 1]) ? 0.5 : 1.0;
                        const double sh_in  = (j + 1 < NR && kept_infill[j + 1] && !pinned[j + 1]) ? 0.5 : 1.0;
                        const double cap_out = std::max(0.0, woven.weave * sh_out * (g_out - 0.8 * double(spacing)));
                        const double cap_in  = std::max(0.0, woven.weave * sh_in  * (g_in  - 0.8 * double(spacing)));
                        const double e = std::clamp(
                            overshoot * std::sin(basephase[j][i] + ring_phase), -1.0, 1.0);
                        const Vec2d q = here + out * (e * (e >= 0.0 ? cap_out : cap_in));
                        moved[i] = Point(coord_t(std::llround(q.x())), coord_t(std::llround(q.y())));
                    }
                    cc.loops[j].points = std::move(moved);
                }
            }
        } else if (woven.solid_stagger > 0.0 && thick && !cc.star_ok && cc.plan == nullptr) {
            // The brick stagger is still parametrized by angle about the centroid, so unlike the
            // weave it is genuinely undefined on a folded section. Declining is right; declining
            // in silence is not - the interior rings of a solid cap then stack into an aligned
            // radial plane with nothing saying so.
            add_note(degrade_note, "solid-cap brick stagger declined: section is not star-shaped, "
                                   "so its interior rings stack aligned between layers");
        } else if (woven.solid_stagger > 0.0 && thick && cc.star_ok && cc.plan == nullptr) {
            // Solid cap/shell brick stagger (see build_cross_section). Shift only the interior rings
            // inward by solid_stagger this layer; the outer wall rings (< sl) and, for a tube, the
            // bore shell stay put. Bounded so a ring never invades the bore. The fill stays fully
            // solid (no decimation, no infill tag) - this only nudges interior geometry per layer.
            AngularProfile inner_prof;                       // per-angle bore-shell radius (shape-aware)
            if (annular)
                inner_prof = AngularProfile::build(cc.loops[n - sl], centre);
            for (size_t j = 0; j < n; ++j) {
                if (!is_infill(j)) continue;                 // shift interior only; shells stay put
                for (Point &p : cc.loops[j].points) {
                    const Vec2d d = p.cast<double>() - centre;
                    const double r = d.norm();
                    if (r < 1.0) continue;
                    const double th = std::atan2(d.y(), d.x());
                    double rr = r - woven.solid_stagger;
                    if (annular)
                        rr = std::max(rr, inner_prof.sample_nearest(th));
                    rr = std::max(rr, 1.0);
                    p = Point(coord_t(std::llround(centre.x() + rr * std::cos(th))),
                              coord_t(std::llround(centre.y() + rr * std::sin(th))));
                }
            }
        }
    }
    if (woven.stride > 1 && rings_before > 2 && cc.loops.size() == rings_before)
        add_note(degrade_note, "continuous shell: wall too thin to sparsify - printed solid");
}


std::vector<ContourChain> contour_spiral(const ExPolygon &region, const Flow &flow, size_t layer_index, const WovenInfill &woven,
                                         const ContinuousShellWidthLimits &limits,
                                         bool adaptive = false, double layer_height = 0.,
                                         bool spiral_wanted = false, double junction_ramp = 1., double seam_theta = 0.,
                                         std::string *fail_reason = nullptr, std::string *degrade_note = nullptr)
{
    const coord_t spacing = flow.scaled_spacing();
    const coord_t width   = flow.scaled_width();
    if (spacing <= 0 || width <= 0)
        return {};

    // Printable-width bounds; zeros derive the Arachne defaults from the nozzle.
    const double min_bead_width_mm   = limits.min_bead_width_mm   > 0. ? limits.min_bead_width_mm   : 0.85 * flow.nozzle_diameter();
    const double min_feature_size_mm = limits.min_feature_size_mm > 0. ? limits.min_feature_size_mm : 0.25 * flow.nozzle_diameter();

    // Adaptive beading (C-full): the plan is built from the REGION directly, before and
    // independent of band marching - the plan measures wall thickness itself and needs
    // nothing from the bands. (An earlier gate required the marching to produce exactly
    // one chain, which switched the plan OFF precisely where marching-based C-lite is
    // weakest: a lumpy wall that pinches during marching fragments into chains whose
    // deficit correction is skipped by design, underfilling the wall.) Out-of-scope
    // sections decline inside build_chain_plan and fall through to marching + C-lite;
    // a plan that contradicts its own invariants aborts the layer loudly.
    // Every route that skips the plan degrades this layer to the marching C-lite path, which
    // does NOT correct its medial deficit on a chain that pinches - a wall of varying thickness
    // can lose a third of its material that way. Each of those routes reports itself; a silent
    // fallback here is how a 31% underfill on eccentric bores went unnoticed (W0 fence).
    // `adaptive` already folds in force_solid, which this function cannot see, so the caller
    // pre-sets the note for the !adaptive cases and this only refines it.
    if (spiral_wanted && woven.enabled) {
        // Not belt-and-braces: a spiral plan chain carries only its two pinned surfaces, and
        // the sparsifier needs more than two rings, so without this sparse infill is a SILENT
        // NO-OP on spiral layers. They are alternatives for the same buried rings anyway - the
        // spiral reorganises exactly what sparse fill decimates - so say so and stop.
        if (fail_reason)
            *fail_reason = "Continuous shell: the interior spiral and sparse infill both reorganise "
                           "the buried interior - enable only one (set continuous_shell_infill_density "
                           "or continuous_shell_interior_spiral to 0).";
        return {};
    }
    if (adaptive && region.holes.size() != 1) {
        add_note(degrade_note, region.holes.empty()
                                   ? "adaptive beading out of scope: solid (no bore) section"
                                   : "adaptive beading out of scope: multi-bore section");
    }
    // Sparse infill no longer excludes the plan. Adaptive beading decides bead COUNTS and
    // WIDTHS; decimation decides which buried rings survive afterwards. They are orthogonal,
    // and coupling them cost every wall in a sparse part its bead count (a 2.035 mm wall got
    // four stretched beads instead of five, with a 0.078 mm void up the medial - printed).
    if (adaptive && region.holes.size() == 1) {
        std::string note, hard;
        std::unique_ptr<ChainBeadPlan> plan =
            build_chain_plan(region, flow, limits, layer_height, spiral_wanted, junction_ramp, seam_theta, &note, &hard);
        if (!hard.empty()) {
            if (fail_reason)
                *fail_reason = hard;
            return {};
        }
        if (degrade_note && !note.empty())
            // A plan that was BUILT can still carry a note (the spiral declining is not the plan
            // declining), so say which happened rather than reporting every note as a fallback.
            add_note(degrade_note, (plan ? "interior spiral declined: "
                                         : "adaptive beading declined: ") + note);
        if (plan) {
            ContourChain cc;
            Polygons pr = plan_rings(*plan);
            cc.has_bore = true;
            const size_t n_slots = plan->n_slots;
            if (plan->spiral && pr.size() >= 2) {
                cc.loops     = Polygons{ std::move(pr.front()), std::move(pr.back()) };
                cc.ring_slot = { 0, n_slots - 1 };          // the two pinned surfaces
            } else {
                cc.loops = std::move(pr);
                cc.ring_slot.resize(cc.loops.size());
                for (size_t i = 0; i < cc.ring_slot.size(); ++i)
                    cc.ring_slot[i] = i;                    // identity until decimation compacts it
            }
            cc.plan = std::move(plan);
            // Sparse infill on a PLANNED chain: the decimation runs here, on the plan's rings,
            // rather than only on marched ones. This is the whole point of carrying ring_slot -
            // before it, enabling continuous_shell_infill_density silently disabled adaptive
            // beading for the entire object.
            apply_woven_infill(cc, woven, spacing, layer_index, degrade_note);
            std::string bad;
            if (!validate_chain(cc, &bad)) {
                if (fail_reason)
                    *fail_reason = bad;
                return {};
            }
            std::vector<ContourChain> out;
            out.push_back(std::move(cc));
            return out;
        }
    }

    std::vector<std::vector<ExPolygon>> chains;
    std::vector<char> vanished;
    std::vector<char> outer_pinch;
    size_t point_count = 0;
    bool   multi_hole  = false;

    ExPolygons current = offset_ex(ExPolygons{ region }, -float(width) / 2.f);
    for (const ExPolygon &seed : current)
        if (!march_bands(seed, spacing, !seed.holes.empty(), chains, vanished, outer_pinch,
                         point_count, multi_hole))
            return {};
    if (multi_hole)
        add_note(degrade_note, "marching stopped early: a band split into several bores - the "
                               "material inward of that point is NOT filled");

    std::vector<ContourChain> out;

    if (chains.empty()) {
        // Sub-bead-width region: fall back to a medial-axis centreline pass so the wall
        // stays connected. Adaptive bead width is deliberately future work. Only a medial
        // that is GENUINELY (near-)closed may become a ring: force-closing an OPEN
        // centreline (a thin blade, an open horseshoe) would fabricate a chord from its
        // end straight back to its start - through air - and every downstream consumer
        // assumes rings. Decline those (empty result -> build_cross_section fails loudly)
        // rather than silently printing a phantom chord. All kept loops stay in ONE chain:
        // the medial of a thin closed wall can come back as several near-closed fragments,
        // and only the in-chain positional re-anchor (spiralize_rings) can hop between
        // them along the wall - the inter-chain stitcher has no rails at sub-bead width
        // (offset_ex by -1.5 beads is empty) and would chord straight across the bore.
        //
        // Adaptive width (C-lite): fit the bead to the wall instead of over-extruding at
        // nominal width. Measure the true thickness (twice the depth at which the raw
        // region vanishes; the w/2 inset above is known empty, so bisect below it),
        // decline outright below the minimum feature size (WideningBeadingStrategy
        // semantics), and otherwise print the centreline at the measured thickness
        // clamped up to the narrowest printable bead.
        const coord_t t_scaled = 2 * vanish_depth(ExPolygons{ region }, width / 2);
        const double  t_mm     = double(t_scaled) * SCALING_FACTOR;
        if (t_mm < min_feature_size_mm)
            return {};            // thinner than the minimum feature size -> decline loudly
        const float bead_mm = float(std::clamp(t_mm, min_bead_width_mm, double(flow.width())));

        Polygons loops;
        Polylines centerlines = region.medial_axis(double(width) / 5., double(2 * width));
        for (const Polyline &pl : centerlines) {
            if (pl.points.size() < 3)
                continue;
            Polygon loop(pl.points);
            if (loop.points.front() == loop.points.back())
                loop.points.pop_back();
            else if ((loop.points.front() - loop.points.back()).cast<double>().norm() > 2.0 * width)
                return {};        // an open centreline: no honest closed-bead path exists
            if (loop.points.size() >= 3)
                loops.push_back(std::move(loop));
        }
        if (loops.empty())
            return {};
        ContourChain cc;
        cc.loops = std::move(loops);
        cc.ring_widths.assign(cc.loops.size(), bead_mm);
        out.push_back(std::move(cc));
        return out;
    }

    for (size_t ci = 0; ci < chains.size(); ++ci) {
        const auto &chain = chains[ci];
        ContourChain cc;
        // Whether this chain's outer end is a real face or the point where a previous chain
        // pinched out. Only a real face gets a solid shell; see ContourChain::outer_pinch.
        cc.outer_pinch = ci < outer_pinch.size() && outer_pinch[ci] != 0;
        cc.loops.reserve(chain.size() * 2 + 1);
        // Outer contours: outer edge -> medial axis.
        for (const ExPolygon &band : chain)
            cc.loops.push_back(band.contour);
        const size_t contours_n = cc.loops.size();
        // Inner (hole) contours: medial axis -> bore (innermost first), for an annular
        // chain only. A piece pinched off solid contributes none.
        for (auto it = chain.rbegin(); it != chain.rend(); ++it)
            if (!it->holes.empty()) {
                cc.loops.push_back(it->holes.front());
                cc.has_bore = true;
            }
        // Medial deficit correction (adaptive bead width, C-lite). Marching stops when one
        // more band no longer fits, leaving a residual gap between the innermost beads'
        // edges. Measure it exactly: the last band vanished somewhere in (0, spacing], so
        // its residual half-thickness d_res is recoverable by bisection, and the void
        // between the innermost bead EDGES is v = 2*d_res - width. Because another band
        // would have fitted otherwise, v <= 2*spacing - width (~0.8 beads) ALWAYS - so the
        // gap never needs new medial rings (the old strip patch): it is closed by WIDENING
        // the innermost ring pair by v/2 each and shifting their centrelines v/4 toward the
        // medial, which extends each bead's inner flank by exactly v/2 (they meet at the
        // medial) while leaving its outer flank untouched. Skipped for chains about to be
        // sparsified (their interior is deliberately open) and for chains that ended by a
        // pinch or split (material continues elsewhere; there is no closed deficit).
        const bool will_decimate = chain_will_decimate(cc, woven);
        if (!cc.plan && !will_decimate && ci < vanished.size() && vanished[ci] && !chain.empty()) {
            const coord_t d_res = vanish_depth(ExPolygons{ chain.back() }, spacing);
            const coord_t v     = 2 * d_res - width;
            if (v > width / 10) {
                const float widened = flow.width() + float(double(v) * SCALING_FACTOR) / 2.f;
                cc.ring_widths.assign(cc.loops.size(), flow.width());
                // Recentre the innermost ring(s) v/4 toward the medial. v/4 < d_res always,
                // so the offset stays strictly inside the last band; if it still fragments
                // (numeric slivers), fall back to width-only (mild outward overlap instead
                // of a void - still strictly better than the gap).
                ExPolygons inner = offset_ex(ExPolygons{ chain.back() }, -float(v) / 4.f);
                const bool ok = inner.size() == 1 &&
                    (cc.has_bore ? inner.front().holes.size() == 1 : inner.front().holes.empty());
                if (ok) {
                    cc.loops[contours_n - 1] = inner.front().contour;
                    if (cc.has_bore)
                        cc.loops[contours_n] = inner.front().holes.front();
                }
                cc.ring_widths[contours_n - 1] = widened;
                if (cc.has_bore && contours_n < cc.loops.size())
                    cc.ring_widths[contours_n] = widened;
            }
        }
        apply_woven_infill(cc, woven, spacing, layer_index, degrade_note);
        if (!cc.loops.empty())
            out.push_back(std::move(cc));
    }

    // Drop negligible sliver pieces (the last necking bits where a non-uniform wall
    // pinches off). Each is a tiny leftover that adds no real fill but forces a long
    // bridge across the bore to reach it. Keep any piece at least a small fraction of
    // the largest piece's footprint, so genuine split lobes (a dogbone, a '+') stay.
    if (out.size() > 1) {
        double max_foot = 0.;
        for (const ContourChain &cc : out)
            max_foot = std::max(max_foot, std::abs(cc.loops.front().area()));
        const double min_foot = 0.05 * max_foot;
        out.erase(std::remove_if(out.begin(), out.end(),
                                 [&](const ContourChain &cc) { return std::abs(cc.loops.front().area()) < min_foot; }),
                  out.end());
    }
    return out;
}

// ---- Serpentine (connected rectilinear) fill for SOLID cap/shell planes -------------------------
// A concentric spiral is only ONE possible single-bead flat fill. For the base/top/solid planes a
// connected rectilinear pass - parallel scan lines joined end-to-end into one boustrophedon path,
// its direction rotated 45/-45 per layer - gives the classic cross-hatch inter-layer interlock and,
// with monotone ordering, the cleanest surface. The outer wall_loops stay concentric so the printed
// edge is smooth; the interior is the serpentine, all chained into one bead. Driven by surface pattern.
bool cs_serpentine_pattern(InfillPattern p)
{
    return p == ipRectilinear || p == ipAlignedRectilinear || p == ipMonotonic || p == ipMonotonicLine;
}

Polyline serpentine_fill(const ExPolygon &region, coord_t spacing, double angle_rad, const Point &center, bool *join_outside = nullptr)
{
    if (spacing <= 0)
        return {};
    ExPolygon rot = region;
    rot.rotate(-angle_rad, center);
    const BoundingBox bb = rot.contour.bounding_box();
    Polylines lines;
    for (coord_t y = bb.min.y() + spacing / 2; y <= bb.max.y(); y += spacing) {
        Polyline l;
        l.points.emplace_back(bb.min.x() - spacing, y);
        l.points.emplace_back(bb.max.x() + spacing, y);
        lines.push_back(std::move(l));
    }
    Polylines segs = intersection_pl(lines, rot);
    if (segs.empty())
        return {};
    for (Polyline &s : segs)                         // normalise every clipped run left -> right
        if (s.points.size() >= 2 && s.points.front().x() > s.points.back().x())
            std::reverse(s.points.begin(), s.points.end());
    std::sort(segs.begin(), segs.end(), [](const Polyline &a, const Polyline &b) {
        if (a.points.front().y() != b.points.front().y()) return a.points.front().y() < b.points.front().y();
        return a.points.front().x() < b.points.front().x();
    });
    // Boustrophedon: walk scan lines bottom -> top, flipping direction each row so the bead never lifts.
    // Every segment-to-segment connection is an EXTRUDING chord (a continuous bead cannot travel), so
    // a concave section is unfillable this way: joining split-row segments or hopping rows across a
    // notch would lay the chord through empty air. Collect every join and test its MIDPOINT against a
    // slightly inflated region (joins hug the contour, where exact containment fails numerically; a
    // real notch chord sits far outside the quarter-bead inflation).
    Polyline serp;
    Lines joins;
    coord_t last_y = std::numeric_limits<coord_t>::min();
    long row = -1;
    std::vector<Polyline> row_segs;
    auto flush_row = [&]() {
        if (row_segs.empty()) return;
        if (row % 2 == 1) {                          // odd rows run right -> left
            std::reverse(row_segs.begin(), row_segs.end());
            for (Polyline &s : row_segs) std::reverse(s.points.begin(), s.points.end());
        }
        for (const Polyline &s : row_segs) {
            if (!serp.points.empty())
                joins.emplace_back(serp.points.back(), s.points.front());
            for (const Point &p : s.points) serp.points.push_back(p);
        }
        row_segs.clear();
    };
    for (Polyline &s : segs) {
        if (s.points.empty()) continue;
        if (s.points.front().y() != last_y) { flush_row(); ++row; last_y = s.points.front().y(); }
        row_segs.push_back(s);
    }
    flush_row();
    const ExPolygons rot_grown = offset_ex(ExPolygons{rot}, float(spacing) / 4.f);
    for (const Line &join : joins) {
        const Point mid((join.a + join.b) / 2);
        bool inside = false;
        for (const ExPolygon &ep : rot_grown)
            if (ep.contains(mid)) { inside = true; break; }
        if (!inside) {
            if (join_outside != nullptr)
                *join_outside = true;
            return {};
        }
    }
    const double c = std::cos(angle_rad), sn = std::sin(angle_rad);   // rotate the path back to world
    for (Point &p : serp.points) {
        const double dx = double(p.x() - center.x()), dy = double(p.y() - center.y());
        p = Point(center.x() + coord_t(std::llround(dx * c - dy * sn)),
                  center.y() + coord_t(std::llround(dx * sn + dy * c)));
    }
    return serp;
}

// Fill a SOLID (no-bore) cap plane as wall_loops concentric rings + a per-layer-rotated serpentine
// interior, chained into one bead. Returns false (leaving `entities` untouched) if the region is not
// a clean single piece, so the caller falls back to the concentric spiral.
bool build_serpentine_cap(const ExPolygon &region, const Flow &flow, double layer_height,
                          const PrintRegionConfig *config, size_t layer_index,
                          ExtrusionEntityCollection &entities, std::string *reason)
{
    const coord_t spacing = flow.scaled_spacing();
    const coord_t width   = flow.scaled_width();
    if (spacing <= 0 || width <= 0) { if (reason) *reason = "Continuous shell serpentine: bad flow."; return false; }
    const size_t wl     = std::max<size_t>(1, size_t(config->wall_loops.value));
    const Point  center = region.contour.bounding_box().center();
    const double angle  = M_PI / 4.0 + double(layer_index % 2) * (M_PI / 2.0);   // 45 / 135 alternating

    // Collect the concentric wall rings (outer -> in) and the interior region they enclose.
    std::vector<Polygon> walls;
    ExPolygons cur = offset_ex(ExPolygons{ region }, -float(width) / 2.f);
    for (size_t i = 0; i < wl; ++i) {
        if (cur.size() != 1) { if (reason) *reason = "Continuous shell serpentine: wall split."; return false; }
        walls.push_back(cur.front().contour);
        cur = offset_ex(cur, -float(spacing));
        if (cur.empty()) break;                                                  // thin cap: walls only
    }
    if (walls.empty()) { if (reason) *reason = "Continuous shell serpentine: no walls."; return false; }
    Polyline serp;
    bool join_outside = false;
    if (cur.size() == 1)
        serp = serpentine_fill(cur.front(), spacing, angle, center, &join_outside);
    else if (!cur.empty()) { if (reason) *reason = "Continuous shell serpentine: interior split."; return false; }
    // A concave interior would need an extruding chord across the notch to stay one
    // bead; decline so the caller falls back to the concentric spiral, which follows
    // the region outline.
    if (join_outside) { if (reason) *reason = "Continuous shell serpentine: concave interior."; return false; }
    const bool  has_serp = serp.points.size() >= 2;

    // Anchor the WALL seam to the serpentine ENTRY point, so the whole part - all wall loops then the
    // serpentine - starts/ends at the same spot. Otherwise the wall finishes at an arbitrary angle and
    // the bead makes one long diagonal hop across the disc to reach the serpentine start, drawing a
    // connector straight over the zig-zag (the overlap seen in preview). With the seams aligned, every
    // wall->wall and wall->serpentine hop is a short local step.
    const Point anchor = has_serp ? serp.points.front() : walls.back().points.front();
    Polyline path;
    for (Polygon &ring : walls) {
        if (ring.points.empty()) continue;
        size_t best = 0; double bd = std::numeric_limits<double>::max();
        for (size_t k = 0; k < ring.points.size(); ++k) {
            const double d = (ring.points[k] - anchor).cast<double>().squaredNorm();
            if (d < bd) { bd = d; best = k; }
        }
        std::rotate(ring.points.begin(), ring.points.begin() + best, ring.points.end());
        for (const Point &p : ring.points) path.points.push_back(p);
        path.points.push_back(ring.points.front());                              // close the ring (back at the anchor)
    }
    if (has_serp)
        for (const Point &p : serp.points) path.points.push_back(p);
    if (path.points.size() < 4) { if (reason) *reason = "Continuous shell serpentine: empty fill."; return false; }
    ExtrusionMultiPath mp;
    mp.paths.emplace_back(make_path(Polyline3(path), flow, layer_height));
    entities.clear();
    entities.no_sort = false;
    entities.append(std::move(mp));
    return true;
}

} // namespace

bool ContinuousShellPathBuilder::build_cross_section(const Surfaces &surfaces, const Flow &flow, double layer_height,
                                                     const PrintRegionConfig *config, ExtrusionEntityCollection &entities,
                                                     std::string *reason, bool reverse_path, size_t layer_index, bool force_solid,
                                                     const ContinuousShellWidthLimits &width_limits,
                                                     const ExPolygons *solid_mask,
                                                     std::string *degrade_note)
{
    ExPolygons regions;
    regions.reserve(surfaces.size());
    for (const Surface &surface : surfaces) {
        if (surface.empty())
            continue;
        if (surface.is_bridge() || surface.surface_type == stInternalVoid || surface.surface_type == stPerimeter) {
            if (reason != nullptr)
                *reason = "Continuous shell mode does not support bridge, void, or perimeter-only surfaces.";
            return false;
        }
        regions.push_back(surface.expolygon);
    }

    if (regions.empty()) {
        if (reason != nullptr)
            *reason = "Continuous shell mode requires a non-empty layer cross-section.";
        return false;
    }

    const ExPolygons merged = union_ex(regions);
    if (merged.size() != 1) {
        if (reason != nullptr)
            *reason = "Continuous shell mode requires one connected layer cross-section.";
        return false;
    }

    // A solid cap/shell with a rectilinear-family surface pattern (top surface pattern): fill the
    // interior with a connected serpentine rotated 45/-45 per layer (cross-hatch inter-layer interlock,
    // clean surface) instead of concentric rings. SOLID (no-bore) planes only; on a split or annular
    // plane build_serpentine_cap declines and we fall through to the concentric spiral below.
    if (force_solid && config != nullptr && merged.front().holes.empty()
        && cs_serpentine_pattern(config->top_surface_pattern.value)
        && build_serpentine_cap(merged.front(), flow, layer_height, config, layer_index, entities, reason)) {
        return true;
    }

    WovenInfill woven;
    if (force_solid) {
        // Solid fill (no sparsification). If this is a multi-layer solid cap/shell inside a woven-
        // infill part, brick-interlock its interior rings between layers: on odd layers shift the
        // interior inward by half a bead so a course's beads nest between the courses above and below,
        // instead of the concentric rings stacking into an aligned radial weak plane. The outer wall
        // rings stay put, so the printed surface is unaffected.
        if (config != nullptr && config->continuous_shell_infill_density.value > 0.0) {
            woven.enabled       = true;
            woven.stride        = 1;                                            // keep every ring (fully solid)
            woven.shell_loops   = std::max<size_t>(1, size_t(config->wall_loops.value));
            woven.solid_stagger = (layer_index % 2 == 1) ? 0.5 * double(flow.scaled_spacing()) : 0.0;
        }
    } else {
        woven = woven_infill_params(config);
        // Only part of the section needs to be solid (an internal surface such as a bore floor).
        // The rings carrying it are kept whole; the rest of the wall still sparsifies.
        if (solid_mask != nullptr && !solid_mask->empty())
            woven.solid_mask = solid_mask;
    }
    // Adaptive beading (C-full) applies to forced-solid layers too.
    //
    // It used to be switched off for them, so a top/bottom shell - or any shell over air - fell
    // back to C-lite marching for its bead COUNTS while the layers above and below were planned.
    // The two disagree on how many beads a wall of given thickness takes, so the count stepped at
    // the boundary between a solid layer and its neighbours, and the step is visible: on the S5
    // cone it printed as a doubled lip at the bore, where the solid cap laid a different number
    // of beads than the wall it sat on.
    //
    // Beading decides bead counts and widths; force_solid decides whether the interior is
    // sparsified. Those are orthogonal, and coupling them is the same mistake W3 fixed for
    // sparse infill - which cost a 2.035 mm wall four stretched beads instead of five.
    //
    // The solid caps keep their serpentine and brick-stagger machinery untouched. The interior
    // SPIRAL stays gated on force_solid: it reorganises the buried interior into one migrating
    // arm, which is not what a solid cap is for.
    const bool adaptive = config != nullptr && config->continuous_shell_adaptive_beading.value;
    const bool   spiral_wanted = adaptive && !force_solid && config->continuous_shell_interior_spiral.value;
    const double junction_ramp = (config != nullptr)
                               ? std::clamp(config->continuous_shell_junction_ramp.value / 100.0, 0.05, 1.0) : 1.0;
    double cs_seam_theta = 0.0;
    if (config != nullptr) {
        const double step_rad   = 0.02 * (config->continuous_shell_seam_step.value / 100.0);
        const double jitter_rad = 0.10 * (config->continuous_shell_seam_jitter.value / 100.0);
        cs_seam_theta = double(layer_index) * step_rad + jitter_rad * seam_jitter_walk(layer_index);
    }
    // A forced-solid layer is no longer a reason for the plan to be off, so the only remaining
    // one at this point is the option itself.
    if (config != nullptr && ! adaptive)
        add_note(degrade_note, "adaptive beading disabled by continuous_shell_adaptive_beading");
    std::string plan_fail;
    const std::vector<ContourChain> chains = contour_spiral(merged.front(), flow, layer_index, woven, width_limits,
                                                            adaptive, layer_height, spiral_wanted, junction_ramp, cs_seam_theta,
                                                            &plan_fail, degrade_note);
    if (chains.empty()) {
        if (reason != nullptr)
            *reason = !plan_fail.empty() ? plan_fail : "Continuous shell mode could not construct a contour spiral.";
        return false;
    }
    // Was any chain actually thinned to sparse infill this layer? Only then does the Z-wave apply
    // (a wall too thin to sparsify keeps its plain buried-ring `internal` tag and must stay flat).
    bool any_woven = false;
    for (const ContourChain &c : chains) any_woven |= c.woven_infill;

    double step_ratio   = 1.0;
    double jitter_ratio = 0.0;
    if (config != nullptr) {
        step_ratio   = config->continuous_shell_seam_step.value / 100.0;
        jitter_ratio = config->continuous_shell_seam_jitter.value / 100.0;
    }
    // Seam control (both per-layer ANGLES, radius-independent):
    //  * seam_step   - a deterministic precession that turns the residual junction into a
    //    gentle HELIX. base_seam_step_rad ~ 1.15 deg/layer at 100%.
    //  * seam_jitter - the per-step size of an ACCUMULATED random walk (applied in
    //    spiralize_rings). At 0 the seam is a clean helix; turned up it wanders around the
    //    part over the whole print, SCATTERING the seam. The per-layer step stays small so
    //    inter-layer joins remain short. base_seam_jitter_rad ~ 5.7 deg/layer per unit
    //    step at 100%, which random-walks to a wide spread over a tall print.
    constexpr double base_seam_step_rad   = 0.02;
    constexpr double base_seam_jitter_rad = 0.10;
    const double seam_step_rad   = base_seam_step_rad * step_ratio;
    const double seam_jitter_rad = base_seam_jitter_rad * jitter_ratio;
    const double ring_spacing    = double(flow.scaled_spacing());

    // Brick-interlock amplitude (scaled). Adjacent buried rings alternate +amp/-amp, so
    // they differ by 2*amp; a "full" 100% interlock makes that half a bead height. See
    // spiralize_rings for how surfaces and the join rings are pinned to stay continuous.
    const double interlock_ratio = (config != nullptr) ? config->continuous_shell_interlock.value / 100.0 : 0.0;
    const coord_t interlock_amp  = coord_t(std::llround(scale_(interlock_ratio * layer_height / 4.0)));

    // Layer height for the interior arm's intrinsic Z-stagger (it takes a quarter of
    // this as its amplitude, and wraps its ramp-tilt compensation modulo the whole
    // layer). Half a layer of stagger was needed while the interior was a pure helical
    // thread, whose only defence against radial cleavage was interdigitation. Once the
    // arm holds ring-like columns (a junction ramp below 100) the topology carries that
    // load, and half a layer starts costing what it was buying: adjacent columns then
    // overlap only half their height, which reads in the section as reduced
    // column-to-column bonding. A quarter layer - the brick interlock's print-proven
    // envelope, and what the first peel-free spiral used - keeps the interdigitation
    // while restoring three quarters of the overlap.
    const coord_t cs_layer_h_scaled = coord_t(std::llround(scale_(layer_height)));

    // A2 - connector spread: distribute the per-ring surface<->interior handoffs around the
    // circumference instead of stacking them into one deep seam groove (passed to
    // spiralize_rings). A3 - seam coast: thin the flow across those handoffs so the residual
    // mark is softer (applied in the flow assignment below). Both default 0 (off).
    const double connector_spread = (config != nullptr) ? config->continuous_shell_connector_spread.value / 100.0 : 0.0;
    const double coast_pct        = (config != nullptr) ? config->continuous_shell_seam_coast.value : 0.0;

    // In-plane scarf (option A): taper the bead width across every ring-to-ring junction
    // (see spiralize_rings). The floor is transitional, not structural - a scarf bead
    // overlaps the junction zone - so it binds to a quarter nozzle, not min_bead_width.
    const coord_t scarf_len   = (config != nullptr) ? coord_t(std::llround(scale_(std::max(0., config->continuous_shell_scarf_length.value)))) : 0;
    const float   scarf_floor = 0.25f * flow.nozzle_diameter();

    // Fill each connected piece (chain) as its own Fermat spiral, then stitch the pieces
    // together greedily - each next piece is appended at whichever end lands nearest the
    // running tail - so a split cross-section's pieces are joined by the shortest hop
    // (along the wall) instead of long chords across the bore. A simple cross-section has
    // exactly one chain, so this is identical to a single spiral there.
    struct Sub { Polyline pl; std::vector<char> flags; std::vector<coord_t> zoff; std::vector<float> widths; bool zwave = false; };
    std::vector<Sub> subs;
    // Per-layer reversal exists to give SINGLE-chain solid layers short inter-layer
    // joins (core end meets core start). A stitched MULTI-chain tour instead closes
    // back to its own seam-anchored start below, so every layer already starts and
    // ends at the (precessing) outer anchor - and reversing would break exactly that
    // pairing: odd layers would START at the buried inner ring while the previous
    // layer ended at the outer anchor, chording across any concave notch between.
    // Ring dump, for analysis with TRUE ring identity. Off unless CS_DUMPRINGS names a layer.
    //
    // Kept rather than deleted because of what it cost to do without. Every other instrument this
    // subsystem has - ray crossings, angular profiles, nominal-band assignment, arc-length
    // correspondence - needs a centroid, a star-shaped outline, or rings of comparable length,
    // and a concave section has none of those. Five separate attempts to measure the crescent
    // produced numbers that were provably wrong (coincident-point gaps of 0.000mm, an inflated
    // section area, a tie profile built on stale nominal rings). Dumping the geometry itself is
    // the only measurement here that is correct by construction on any topology.
    static const char *dump_layer = std::getenv("CS_DUMPRINGS");
    if (dump_layer) {
        if (layer_index == size_t(std::atoi(dump_layer))) {
            if (FILE *f = std::fopen("cs_rings.txt", "a")) {
                for (size_t ci = 0; ci < chains.size(); ++ci) {
                    const ContourChain &c = chains[ci];
                    for (size_t j = 0; j < c.loops.size(); ++j) {
                        std::fprintf(f, "CHAIN %zu RING %zu BORE %d PINCH %d N %zu",
                                     ci, j, c.has_bore ? 1 : 0, c.outer_pinch ? 1 : 0,
                                     c.loops[j].points.size());
                        for (const Point &p : c.loops[j].points)
                            std::fprintf(f, " %.4f %.4f", unscale<double>(p.x()), unscale<double>(p.y()));
                        std::fputc(10, f);
                    }
                }
                std::fclose(f);
            }
        }
    }
    const bool reverse_chains = chains.size() > 1 ? false : reverse_path;
    for (const ContourChain &chain : chains) {
        std::vector<char> flags;
        std::vector<coord_t> zoff;
        std::vector<float> widths;
        Polyline sub = spiralize_rings(chain.loops, reverse_chains, seam_step_rad, seam_jitter_rad,
                                       ring_spacing, layer_index, chain.has_bore, &flags, interlock_amp, &zoff,
                                       connector_spread, chain.woven_infill ? woven.shell_loops : 0,
                                       &merged.front(), &chain.ring_widths, &widths, flow.width(),
                                       scarf_len, scarf_floor, chain.plan.get(), cs_layer_h_scaled,
                                       chain.ring_slot.empty() ? nullptr : &chain.ring_slot);
        if (sub.points.size() >= 2 && sub.points.size() == flags.size() && sub.points.size() == zoff.size()
            && sub.points.size() == widths.size())
            // A chain carries the Z-wave only when its interior was actually thinned AND its
            // rings passed the star-shape gate (angular phase is meaningless on a folded chain).
            subs.push_back({ std::move(sub), std::move(flags), std::move(zoff), std::move(widths),
                             chain.woven_infill && chain.star_ok });
    }
    if (subs.empty()) {
        if (reason != nullptr)
            *reason = "Continuous shell mode could not construct a continuous spiral path.";
        return false;
    }

    Polyline spiral;
    std::vector<char> internal;
    std::vector<char> connector;        // 1 on inter-chain connector points (thinned flow)
    std::vector<char> weavable;         // 1 where the owning chain may carry the Z-wave (star gate)
    std::vector<coord_t> interlock_z;   // per-point Z offset (brick interlock), parallel to spiral.points
    std::vector<float> pt_width;        // per-point extrusion width (adaptive bead width)
    auto append_sub = [&](Sub &s, bool rev) {
        if (rev) {
            std::reverse(s.pl.points.begin(), s.pl.points.end());
            std::reverse(s.flags.begin(), s.flags.end());
            std::reverse(s.zoff.begin(), s.zoff.end());
            std::reverse(s.widths.begin(), s.widths.end());
        }
        spiral.points.insert(spiral.points.end(), s.pl.points.begin(), s.pl.points.end());
        internal.insert(internal.end(), s.flags.begin(), s.flags.end());
        interlock_z.insert(interlock_z.end(), s.zoff.begin(), s.zoff.end());
        connector.insert(connector.end(), s.pl.points.size(), char(0));
        weavable.insert(weavable.end(), s.pl.points.size(), char(s.zwave ? 1 : 0));
        pt_width.insert(pt_width.end(), s.widths.begin(), s.widths.end());
    };

    // Inter-chain connector routing. A split cross-section's sub-chains can END deep
    // inside a lobe (a solid chain finishes at its own centre), so the shortest
    // endpoint hop may cut straight across the concave notch between lobes - through
    // air (the dogbone chords). Route those connectors through material instead:
    // straight when the hop stays inside the region, otherwise out to a "rail" one
    // bead inside the wall, along the wall (through the waist), and back in. Rails
    // sit at 1.5 bead-widths so the connector rides on buried beads, not the surface.
    Polygons rails;
    if (subs.size() > 1 && !merged.empty()) {
        // Densify the rails (Clipper offsets keep only corner vertices): the connector
        // enters/leaves at the NEAREST VERTEX, and with corner-only rails that approach
        // segment can cut across a nearby notch corner. At ~1 mm spacing the nearest
        // vertex is essentially the nearest point, so approaches run straight into the
        // wall through material.
        for (const ExPolygon &ep : offset_ex(ExPolygons{ merged.front() }, -1.5f * float(flow.scaled_width()))) {
            auto densify = [](const Polygon &poly) -> Polygon {
                const double step = scale_(1.0);
                const size_t n = std::clamp<size_t>(size_t(poly.length() / step), 16, 4096);
                std::vector<Point> pts = resample_closed(poly, n);
                return pts.empty() ? poly : Polygon(Points(pts.begin(), pts.end()));
            };
            rails.push_back(densify(ep.contour));
            for (const Polygon &h : ep.holes)
                rails.push_back(densify(h));
        }
    }
    auto route_connector = [&](const Point &from, const Point &to) -> Polyline {
        if ((to - from).cast<double>().norm() <= 2.0 * ring_spacing)
            return {};                                    // short local hop: leave straight
        Polyline direct;
        direct.points = { from, to };
        if (merged.front().contains(direct))
            return {};                                    // straight hop stays in material
        if (rails.empty())
            return {};
        // Pick the rail with the closest approach to both endpoints.
        size_t best_rail = 0, ia = 0, ib = 0;
        double best_sum = std::numeric_limits<double>::max();
        for (size_t r = 0; r < rails.size(); ++r) {
            const Points &pts = rails[r].points;
            if (pts.size() < 3) continue;
            size_t ka = 0, kb = 0; double da = 1e300, db = 1e300;
            for (size_t k = 0; k < pts.size(); ++k) {
                const double d1 = (pts[k] - from).cast<double>().squaredNorm();
                const double d2 = (pts[k] - to).cast<double>().squaredNorm();
                if (d1 < da) { da = d1; ka = k; }
                if (d2 < db) { db = d2; kb = k; }
            }
            const double sum = std::sqrt(da) + std::sqrt(db);
            if (sum < best_sum) { best_sum = sum; best_rail = r; ia = ka; ib = kb; }
        }
        if (best_sum >= 1e300 || ia == ib)
            return {};
        // Walk the shorter arc ia -> ib along the rail.
        const Points &pts = rails[best_rail].points;
        const size_t n = pts.size();
        double fwd = 0., bwd = 0.;
        for (size_t k = ia; k != ib; k = (k + 1) % n) fwd += (pts[(k + 1) % n] - pts[k]).cast<double>().norm();
        for (size_t k = ia; k != ib; k = (k + n - 1) % n) bwd += (pts[(k + n - 1) % n] - pts[k]).cast<double>().norm();
        Polyline arc;
        if (fwd <= bwd)
            for (size_t k = ia; ; k = (k + 1) % n) { arc.points.push_back(pts[k]); if (k == ib) break; }
        else
            for (size_t k = ia; ; k = (k + n - 1) % n) { arc.points.push_back(pts[k]); if (k == ib) break; }
        return arc;
    };

    std::vector<char> used(subs.size(), 0);
    append_sub(subs[0], false);   // start on the outermost chain (layer starts at a surface)
    used[0] = 1;
    for (size_t n = 1; n < subs.size(); ++n) {
        const Point tail = spiral.points.back();
        size_t best = 0; bool best_rev = false; double best_d = std::numeric_limits<double>::max();
        for (size_t i = 0; i < subs.size(); ++i) {
            if (used[i]) continue;
            const double ds = (subs[i].pl.points.front() - tail).cast<double>().squaredNorm();
            const double de = (subs[i].pl.points.back()  - tail).cast<double>().squaredNorm();
            if (ds < best_d) { best_d = ds; best = i; best_rev = false; }
            if (de < best_d) { best_d = de; best = i; best_rev = true; }
        }
        const Point head = best_rev ? subs[best].pl.points.back() : subs[best].pl.points.front();
        const Polyline arc = route_connector(tail, head);
        for (const Point &rp : arc.points) {
            spiral.points.push_back(rp);
            internal.push_back(1);          // rails are buried by construction
            interlock_z.push_back(0);
            connector.push_back(1);
            weavable.push_back(0);
            pt_width.push_back(flow.width());
        }
        const size_t head_idx = spiral.points.size();
        append_sub(subs[best], best_rev);
        if (head_idx < connector.size())
            connector[head_idx] = 1;        // thin the final hop segment onto the chain head
        used[best] = 1;
    }
    if (subs.size() > 1 && spiral.points.size() >= 2) {
        // Close the stitched tour back to its start. Single-chain layers keep their
        // inter-LAYER join short via reverse_path parity, but a stitched tour ends
        // deep inside its last lobe - with travels stripped, the NEXT layer's first
        // extrusion would then chord across the notch through air. Returning to the
        // (seam-anchored) start point makes the layer join a short surface hop again.
        const Point tail = spiral.points.back();
        const Point head = spiral.points.front();
        // NOTE: the arc must be materialized - ranging over `route_connector(...).points`
        // directly binds to a member of a destroyed temporary (dangling range-for).
        const Polyline back_arc = route_connector(tail, head);
        for (const Point &rp : back_arc.points) {
            spiral.points.push_back(rp);
            internal.push_back(1);
            interlock_z.push_back(0);
            connector.push_back(1);
            weavable.push_back(0);
            pt_width.push_back(flow.width());
        }
        spiral.points.push_back(head);
        internal.push_back(1);
        interlock_z.push_back(0);
        connector.push_back(1);
        weavable.push_back(0);
        pt_width.push_back(flow.width());
    }
    if (spiral.points.size() < 2) {
        if (reason != nullptr)
            *reason = "Continuous shell mode could not construct a continuous spiral path.";
        return false;
    }

    // Per-point flow multiplier combining two independent effects, both default 1.0 (off):
    //  * internal over-extrusion densifies provably-buried beads to close voids, and
    //  * seam coast (A3) thins the flow across the surface<->interior handoffs (the seam
    //    crossings) so the material deposited at the surface there is reduced, softening the
    //    once-per-layer seam groove.
    // Each point keeps its brick-interlock Z offset (interlock_z); apply_continuous_z later
    // ADDS the arc-length ramp on top, so smooth spiral + radial stagger compose into one bead.
    const bool   have_interlock = (interlock_z.size() == spiral.points.size());
    const bool   have_flags     = (internal.size() == spiral.points.size());
    const double over_pct       = (config != nullptr) ? config->continuous_shell_internal_overextrusion.value : 0.0;

    // Option D Z-wave. The radial weave (in contour_spiral) already migrates the buried infill
    // beads in/out per layer. Add a matching Z undulation so they also weave VERTICALLY: each
    // buried point gets a Z offset A_z*sin(k*theta + phase), same angular lobes and per-layer phase
    // as the radial wave but in QUADRATURE (zphase ~ pi/2), so a point traces a small circle in the
    // (r,z) plane as it goes around - a minimal-surface-like helical weave. With the per-layer
    // phase, adjacent layers' undulations are offset and nest into each other (vertical interlock),
    // completing the 3D gyroid-like structure. Composes onto interlock_z; apply_continuous_z adds
    // the ramp on top. Surfaces (shell) are left flat.
    if (have_interlock && have_flags && weavable.size() == spiral.points.size() && !merged.empty()
        && any_woven && woven.zamp > 0.0) {
        Vec2d centre(0., 0.);
        const Polygon &oc = merged.front().contour;
        if (!oc.points.empty()) {
            for (const Point &p : oc.points) centre += p.cast<double>();
            centre /= double(oc.points.size());
        }
        const double zamp  = std::clamp(woven.zamp, 0.0, 0.5) * scale_(layer_height);
        const double phase = double(layer_index) * woven.dphase + woven.zphase;
        for (size_t i = 0; i < spiral.points.size(); ++i) {
            // Buried (woven) infill only - not shells, not connectors, and not chains whose
            // rings failed the star-shape gate (their angular phase would be meaningless).
            if (!internal[i] || connector[i] || !weavable[i]) continue;
            const Vec2d d = spiral.points[i].cast<double>() - centre;
            const double th = std::atan2(d.y(), d.x());
            interlock_z[i] += coord_t(std::llround(zamp * std::sin(woven.lobes * th + phase)));
        }
    }

    // Volumetric cap for internal over-extrusion. Adding a fixed % of extra flow to EVERY
    // buried bead (the old behaviour) dumps far more material than the cross-section can
    // absorb, so the excess bulges out to the surface and wrecks the finish. Instead bound the
    // TOTAL extra volume to over_pct% of the residual VOID between the solid model cross-section
    // and what the nominal beads already deposit - so the grand total extruded can never exceed
    // the model's own solid volume (no bulge), with the budget spread across the buried beads:
    //   V_model   = wall cross-section area * layer height
    //   V_nominal = nominal volumetric flow (mm^3/mm) * total bead length
    //   void      = max(0, V_model - V_nominal);   budget = over_pct% * void
    //   over_mult = 1 + budget / (flow * buried bead length), plus a hard per-bead safety cap
    // Per-point flow rate: widened beads (adaptive width) deposit more per mm, so both the
    // nominal-volume estimate below and the emitted sub-paths must use the point's own
    // width, not the nominal flow. Widths take only a handful of distinct values per layer,
    // so a one-entry memo covers the lookups.
    const bool have_widths = (pt_width.size() == spiral.points.size());
    float  memo_w = flow.width();
    double memo_q = flow.mm3_per_mm();
    auto mm3_of = [&](float w) -> double {
        if (w != memo_w) { memo_w = w; memo_q = double(flow.with_width(w).mm3_per_mm()); }
        return memo_q;
    };

    double over_mult = 1.0;
    if (over_pct > 0.0 && have_flags && !merged.empty()) {
        double v_nom = 0.0, v_buried = 0.0;
        for (size_t i = 1; i < spiral.points.size(); ++i) {
            const double seg = (spiral.points[i] - spiral.points[i - 1]).cast<double>().norm() * SCALING_FACTOR;
            const double q   = have_widths ? mm3_of(pt_width[i]) : flow.mm3_per_mm();
            v_nom += q * seg;
            if (internal[i]) v_buried += q * seg;   // over-extrusion only lands on buried beads
        }
        const double v_model  = std::abs(merged.front().area()) * SCALING_FACTOR * SCALING_FACTOR * layer_height;
        const double budget   = (over_pct / 100.0) * std::max(0.0, v_model - v_nom);
        if (v_buried > 1e-9)
            over_mult = std::min(1.0 + budget / v_buried, 1.5); // safety: never more than +50% on a single bead
    }
    const double coast_mult     = 1.0 - std::clamp(coast_pct, 0.0, 95.0) / 100.0;

    std::vector<double> fmul(spiral.points.size(), 1.0);
    if (have_flags) {
        if (over_pct > 0.0)
            for (size_t i = 0; i < fmul.size(); ++i)
                if (internal[i] && !connector[i]) fmul[i] *= over_mult;
        if (coast_pct > 0.0) {
            // Thin a small window around every surface<->interior transition (the seam crossing).
            constexpr size_t win = 2;
            std::vector<char> at_seam(fmul.size(), 0);
            for (size_t i = 0; i + 1 < internal.size(); ++i)
                if (internal[i] != internal[i + 1]) {
                    const size_t lo = (i >= win) ? i - win : 0;
                    const size_t hi = std::min(fmul.size() - 1, i + 1 + win);
                    for (size_t j = lo; j <= hi; ++j) at_seam[j] = 1;
                }
            for (size_t i = 0; i < fmul.size(); ++i)
                if (at_seam[i]) fmul[i] *= coast_mult;
        }
    }
    // Inter-chain connectors ride along already-printed buried beads (the rail one
    // bead inside the wall), so they deposit at half flow: enough to stay fused and
    // continuous, without double-filling the wall they travel over.
    if (connector.size() == fmul.size()) {
        constexpr double connector_mult = 0.5;
        for (size_t i = 0; i < fmul.size(); ++i)
            if (connector[i]) fmul[i] *= connector_mult;
    }

    // Emit maximal runs of equal (flow multiplier, bead width) as sub-paths (keyed on the
    // segment's second point, matching the original over-extrusion convention), carrying
    // interlock Z. Each run is built with ITS OWN Flow so path.width and mm3_per_mm stay
    // mutually consistent - path.width is what the G-code ;WIDTH: tag and the preview
    // render, and it feeds real logic (volumetric speed caps), so it must not lie about
    // a widened bead.
    ExtrusionMultiPath path;
    size_t start = 0;
    while (start + 1 < spiral.points.size()) {
        const double cur   = fmul[start + 1];
        const float  cur_w = have_widths ? pt_width[start + 1] : flow.width();
        size_t end = start + 1;
        while (end + 1 < spiral.points.size() && fmul[end + 1] == cur
               && (!have_widths || pt_width[end + 1] == cur_w))
            ++end;
        Polyline sub;
        sub.points.assign(spiral.points.begin() + start, spiral.points.begin() + end + 1);
        Polyline3 p3(sub);
        if (have_interlock)
            for (size_t i = 0; i < p3.points.size(); ++i)
                p3.points[i].z() = interlock_z[start + i];
        const Flow run_flow = (cur_w == flow.width()) ? flow : flow.with_width(cur_w);
        ExtrusionPath ep = make_path(p3, run_flow, layer_height);
        ep.mm3_per_mm *= cur;
        path.paths.emplace_back(std::move(ep));
        start = end;
    }
    if (path.paths.empty()) {
        Polyline3 p3(spiral);
        if (have_interlock)
            for (size_t i = 0; i < p3.points.size(); ++i)
                p3.points[i].z() = interlock_z[i];
        path.paths.emplace_back(make_path(p3, flow, layer_height));
    }

    entities.clear();
    entities.no_sort = false;
    entities.append(std::move(path));
    return true;
}

void ContinuousShellPathBuilder::collect_paths(ExtrusionEntity *entity, std::vector<ExtrusionPath *> &out)
{
    if (entity == nullptr)
        return;
    if (auto *coll = dynamic_cast<ExtrusionEntityCollection *>(entity)) {
        for (ExtrusionEntity *e : coll->entities)
            collect_paths(e, out);
    } else if (auto *multi = dynamic_cast<ExtrusionMultiPath *>(entity)) {
        for (ExtrusionPath &p : multi->paths)
            out.push_back(&p);
    } else if (auto *loop = dynamic_cast<ExtrusionLoop *>(entity)) {
        for (ExtrusionPath &p : loop->paths)
            out.push_back(&p);
    } else if (auto *path = dynamic_cast<ExtrusionPath *>(entity)) {
        out.push_back(path);
    }
}

void ContinuousShellPathBuilder::apply_continuous_z(ExtrusionEntity *entity, double layer_height, double min_z_rel)
{
    std::vector<ExtrusionPath *> paths;
    collect_paths(entity, paths);
    if (paths.empty() || layer_height <= 0.)
        return;

    // Total XY arc length across the concatenated bead (paths share endpoints where
    // build_cross_section split them for over-extrusion, which contributes zero length).
    auto xy = [](const Point3 &p) { return Vec2d(double(p.x()), double(p.y())); };
    double total = 0.;
    const Point3 *prev = nullptr;
    for (const ExtrusionPath *p : paths) {
        for (const Point3 &pt : p->polyline.points) {
            if (prev != nullptr)
                total += (xy(pt) - xy(*prev)).norm();
            prev = &pt;
        }
    }
    if (total <= 0.)
        return;

    // Walk again, assigning each point a Z offset from the layer's nominal Z. The bead
    // rises one layer_height across its whole XY arc, so it starts a layer below nominal
    // and ends exactly at nominal. Arc-length (not depth) keeps the reverse_path join
    // continuous: consecutive layers meet at opposite arc-fractions whose z_diff differ
    // by exactly one layer, cancelling the nominal-Z step. Any brick-interlock offset that
    // build_cross_section pre-loaded into pt.z() is PRESERVED (added to, not overwritten),
    // so the smooth ramp and the radial stagger compose into one continuous bead.
    double acc = 0.;
    prev = nullptr;
    // Physical floor: the ramp plus any pre-loaded weave/interlock dip must never push a
    // point below the caller-supplied minimum (e.g. just above the bed on layer 1, where
    // a woven Z-wave could otherwise dive under a thin first layer).
    const coord_t min_z_scaled = coord_t(std::llround(scale_(min_z_rel)));
    for (ExtrusionPath *p : paths) {
        for (Point3 &pt : p->polyline.points) {
            if (prev != nullptr)
                acc += (xy(pt) - xy(*prev)).norm();
            const double factor = acc / total;
            const double z_diff = -(1.0 - factor) * layer_height;
            pt.z() += coord_t(std::llround(scale_(z_diff)));
            if (pt.z() < min_z_scaled)
                pt.z() = min_z_scaled;
            prev = &pt;
        }
        p->z_contoured = true;
    }
}

} // namespace Slic3r
