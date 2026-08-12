#include "ContinuousShell.hpp"

#include "ClipperUtils.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
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
                         const ExPolygon *region = nullptr)
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
        for (size_t k = 0; k < grid[0].size(); ++k) { spiral.points.push_back(grid[0][k]); flag(false); zoff(0); }
        spiral.points.push_back(grid[0][0]); flag(false); zoff(0);
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
    for (const size_t idx : order) {
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
        for (size_t k = 0; k < c.size(); ++k) { spiral.points.push_back(c[k]); flag(internal_ring); zoff(zo); }
    }

    return spiral;
}

// Collect concentric offset bands by marching inward from `band`, appending the nested
// sequence as ONE chain to `chains`. When the region splits during inward offset - a
// SOLID waist (e.g. a dogbone) or a pinched/branching annular wall (a non-uniform tube
// wall that necks off, common at fine spacing) - it ends the current chain and recurses
// so each piece becomes its OWN chain (filled contiguously, joined by a single bridge).
// Returns false only for topologies we still do not model (a solid band that grows a
// hole) or on a band/point budget overflow, so the caller can hard-stop.
bool march_bands(ExPolygon band, coord_t spacing, bool has_hole,
                 std::vector<std::vector<ExPolygon>> &chains, size_t &point_count)
{
    constexpr size_t max_chains = 4096;
    constexpr size_t max_points = 400000;
    std::vector<ExPolygon> chain;
    for (;;) {
        // A wall band must keep exactly one hole; when it pinches we hit the medial.
        if (has_hole && band.holes.size() != 1)
            break;
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
        if (next.empty())
            break;                           // vanished at the medial axis
        if (next.size() == 1) {              // still one piece -> keep marching
            band = std::move(next.front());
            continue;
        }
        // The chain split into several pieces (a solid waist, or a pinched annular
        // wall). End this chain and recurse into each piece as its own chain, using
        // each piece's OWN hole status (a piece may keep the bore or close up solid).
        if (!chain.empty())
            chains.push_back(std::move(chain));
        for (const ExPolygon &comp : next)
            if (!march_bands(comp, spacing, !comp.holes.empty(), chains, point_count))
                return false;
        return true;
    }
    if (!chain.empty())
        chains.push_back(std::move(chain));
    return true;
}

// One connected piece of a (possibly split) cross-section: its concentric loops in
// radial order (outer contours outer->medial, then, if annular, hole contours
// medial->bore) plus whether that piece is annular (so spiralize_rings knows its
// innermost loop is a bore surface, not a buried centre).
struct ContourChain {
    Polygons loops;
    bool     has_bore = false;
    bool     woven_infill = false;   // true once this chain's interior was thinned to sparse infill
    // False when the chain's rings are not star-shaped about the perturbation centre, so
    // every angular-parametrized perturbation (radial weave, solid stagger, and this
    // chain's share of the Z-wave) was/must be declined. Sparse DECIMATION is unaffected -
    // dropping rings needs no parametrization - so a declined chain still prints plain
    // concentric sparse rings.
    bool     star_ok = true;
};

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
std::vector<ContourChain> contour_spiral(const ExPolygon &region, const Flow &flow, size_t layer_index, const WovenInfill &woven)
{
    const coord_t spacing = flow.scaled_spacing();
    const coord_t width   = flow.scaled_width();
    if (spacing <= 0 || width <= 0)
        return {};

    std::vector<std::vector<ExPolygon>> chains;
    size_t point_count = 0;

    ExPolygons current = offset_ex(ExPolygons{ region }, -float(width) / 2.f);
    for (const ExPolygon &seed : current)
        if (!march_bands(seed, spacing, !seed.holes.empty(), chains, point_count))
            return {};

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
        out.push_back({ std::move(loops), false });
        return out;
    }

    for (const auto &chain : chains) {
        ContourChain cc;
        cc.loops.reserve(chain.size() * 2 + 1);
        // Outer contours: outer edge -> medial axis.
        for (const ExPolygon &band : chain)
            cc.loops.push_back(band.contour);
        // Medial gap fill. Marching inward from both surfaces stops when the annular band
        // can no longer be offset, but the last band's two fronts (its innermost outer
        // contour and innermost bore contour) can be more than one bead-width apart -
        // leaving an unfilled medial strip that shows up as a vertical void ring between
        // the inner and outer walls. Offsetting the last band in by width/2 isolates that
        // uncovered strip (empty where the wall already tiles, so this self-gates); its
        // outer and bore contours sit just inside the two fronts and bracket the gap.
        // Emitted here (between the outer and bore contour runs) so the path still spirals
        // in, crosses the medial, and spirals back out. Only ONE offset is used - offsetting
        // the already-thin strip again makes it near-degenerate and fragments into many
        // slivers (O(n^2) blow-up), so we take the strip's own contours directly.
        // The concentric loops emitted above are consumed by spiralize_rings as a strictly
        // radially-nested stack: the Fermat interleave connects loop i to loop i+2 assuming
        // each is one bead further in. A medial loop is only safe to insert if it, too, is a
        // clean concentric annulus that brackets the SAME bore. On a uniform wall the last
        // band's medial gap is exactly that. On a thick/non-uniform base (a flared foot) the
        // last band is lopsided: offsetting it by width/2 fragments into off-centre blobs or
        // opens the hole into the contour. Inserted into the radial order those blobs get
        // bridged straight across the hollow bore (bird's-nest chords). So gate strictly:
        // add the medial ring ONLY when the inward offset stays a single one-holed annulus.
        if (!chain.empty() && chain.back().holes.size() == 1) {
            // Coarsen the band before offsetting: the medial fill does not need full
            // slice-contour precision, and offsetting the detailed contour every layer is
            // the dominant cost (offset_ex is ~O(n^1.5) in contour points).
            ExPolygons simplified = chain.back().simplify(double(width) / 4.);
            if (simplified.size() == 1) {
                ExPolygons strips = offset_ex(ExPolygons{ simplified.front() }, -float(width) / 2.f);
                if (strips.size() == 1 && strips.front().holes.size() == 1) {
                    const ExPolygon &strip = strips.front();
                    if (strip.contour.points.size() >= 3 && strip.contour.length() > 4 * width)
                        cc.loops.push_back(strip.contour);
                    const Polygon &hole = strip.holes.front();
                    if (hole.points.size() >= 3 && hole.length() > 4 * width)
                        cc.loops.push_back(hole);
                }
            }
        }
        // Inner (hole) contours: medial axis -> bore (innermost first), for an annular
        // chain only. A piece pinched off solid contributes none.
        for (auto it = chain.rbegin(); it != chain.rend(); ++it)
            if (!it->holes.empty()) {
                cc.loops.push_back(it->holes.front());
                cc.has_bore = true;
            }
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
            const size_t K  = woven.stride;
            const size_t n  = cc.loops.size();
            const bool   annular = cc.has_bore;
            // is_infill on ORIGINAL indices: buried interior beyond the shell(s). Annular keeps a
            // shell at both ends; solid keeps the outer shell and the buried core (minus the terminus).
            auto is_infill = [&](size_t i) -> bool {
                if (i < sl) return false;
                return annular ? (i + sl < n) : (i + 1 < n);
            };
            const bool thick = annular ? (n > 2 * sl + 1) : (n > sl + 1);
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
                kept.reserve(n); kept_infill.reserve(n);
                for (size_t i = 0; i < n; ++i) {
                    const bool inf = is_infill(i);
                    if (!inf || (i - sl) % K == 0) {   // keep every shell/terminus ring + every K-th infill ring
                        kept_infill.push_back(inf ? char(1) : char(0));
                        kept.push_back(std::move(cc.loops[i]));
                    }
                }
                cc.loops = std::move(kept);
                cc.woven_infill = true;   // spiralize_rings should tag [shell | infill | shell] accordingly

                if (woven.weave > 0.0 && cc.star_ok) {
                    // Radial amplitude sized so that ANTI-PHASE neighbours just touch. Sparse rings
                    // sit K*spacing apart; giving consecutive infill rings a pi phase offset makes ring
                    // j bulge OUT where ring j+1 bulges IN, so their closest approach is
                    // K*spacing - 2*amp. Targeting a one-bead centre gap (amp = weave*(K-1)/2*spacing)
                    // brings the bead EDGES into contact at the lobe peaks - the rings knit into one
                    // connected in-plane lattice ("lobes along a slice") instead of separate wavy tubes,
                    // while the valleys stay open for the sparse voids. Capped below half the spacing so
                    // they never cross.
                    const double amp = std::min(woven.weave * 0.5 * double(K - 1), 0.49 * double(K)) * double(spacing);
                    // The lobes must stay INSIDE the infill band: an outward bulge that crossed the
                    // solid outer shell would push nubs past the object surface (a rough wall); an
                    // inward one would invade the bore shell. The shell rings are un-perturbed, so
                    // clamp every infill point to the shell radius AT ITS OWN ANGLE - shape-aware, so
                    // it works for non-round cross-sections too (a mean radius would flatten corners).
                    auto build_profile = [&](const Polygon &ring) {
                        const int NB = 360;
                        std::vector<double> prof(NB, 0.0); std::vector<int> cnt(NB, 0);
                        for (const Point &q : ring.points) {
                            const Vec2d d = q.cast<double>() - centre;
                            int b = int((std::atan2(d.y(), d.x()) + M_PI) / (2.0 * M_PI) * NB);
                            b = ((b % NB) + NB) % NB;
                            prof[b] += d.norm(); ++cnt[b];
                        }
                        for (int b = 0; b < NB; ++b) if (cnt[b]) prof[b] /= cnt[b];
                        double last = 0.0;                                  // fill empty bins from the last seen
                        for (int b = 0; b < 2 * NB; ++b) { int i = b % NB; if (prof[i] > 0.0) last = prof[i]; else if (last > 0.0) prof[i] = last; }
                        return prof;
                    };
                    auto sample = [](const std::vector<double> &prof, double th) {
                        const int NB = int(prof.size());
                        int b = int((th + M_PI) / (2.0 * M_PI) * NB); b = ((b % NB) + NB) % NB;
                        return prof[b];
                    };
                    const std::vector<double> outer_prof = build_profile(cc.loops[sl - 1]);                        // innermost outer-shell ring
                    const std::vector<double> inner_prof = annular ? build_profile(cc.loops[cc.loops.size() - sl]) // outermost bore-shell ring
                                                                   : std::vector<double>();
                    // Per-ring anti-phase: each successive infill ring flips by pi so neighbours interlock.
                    // The per-LAYER advance (dphase) still rotates the whole lattice so stacked layers weave.
                    const double layer_phase = double(layer_index) * woven.dphase;
                    size_t infill_ring = 0;
                    for (size_t j = 0; j < cc.loops.size(); ++j) {
                        if (!kept_infill[j]) continue;   // weave only the thinned infill, not the solid shells
                        const double phase = layer_phase + double(infill_ring) * M_PI;
                        ++infill_ring;
                        for (Point &p : cc.loops[j].points) {
                            const Vec2d d = p.cast<double>() - centre;
                            const double r = d.norm();
                            if (r < 1.0) continue;
                            const double th = std::atan2(d.y(), d.x());
                            double rr = r + amp * std::sin(woven.lobes * th + phase);
                            rr = std::min(rr, sample(outer_prof, th));                 // never cross the outer shell
                            if (annular) rr = std::max(rr, sample(inner_prof, th));    // never cross the bore shell
                            rr = std::max(rr, 1.0);
                            p = Point(coord_t(std::llround(centre.x() + rr * std::cos(th))),
                                      coord_t(std::llround(centre.y() + rr * std::sin(th))));
                        }
                    }
                }
            } else if (woven.solid_stagger > 0.0 && thick && cc.star_ok) {
                // Solid cap/shell brick stagger (see build_cross_section). Shift only the interior rings
                // inward by solid_stagger this layer; the outer wall rings (< sl) and, for a tube, the
                // bore shell stay put. Bounded so a ring never invades the bore. The fill stays fully
                // solid (no decimation, no infill tag) - this only nudges interior geometry per layer.
                const int NB = 360;
                std::vector<double> inner_prof;
                if (annular) {                                   // per-angle bore-shell radius (shape-aware)
                    inner_prof.assign(NB, 0.0); std::vector<int> cnt(NB, 0);
                    for (const Point &q : cc.loops[n - sl].points) {
                        const Vec2d d = q.cast<double>() - centre;
                        int b = int((std::atan2(d.y(), d.x()) + M_PI) / (2.0 * M_PI) * NB); b = ((b % NB) + NB) % NB;
                        inner_prof[b] += d.norm(); ++cnt[b];
                    }
                    for (int b = 0; b < NB; ++b) if (cnt[b]) inner_prof[b] /= cnt[b];
                    double last = 0.0; for (int b = 0; b < 2 * NB; ++b) { int i = b % NB; if (inner_prof[i] > 0.0) last = inner_prof[i]; else if (last > 0.0) inner_prof[i] = last; }
                }
                for (size_t j = 0; j < n; ++j) {
                    if (!is_infill(j)) continue;                 // shift interior only; shells stay put
                    for (Point &p : cc.loops[j].points) {
                        const Vec2d d = p.cast<double>() - centre;
                        const double r = d.norm();
                        if (r < 1.0) continue;
                        const double th = std::atan2(d.y(), d.x());
                        double rr = r - woven.solid_stagger;
                        if (annular) {
                            int b = int((th + M_PI) / (2.0 * M_PI) * NB); b = ((b % NB) + NB) % NB;
                            rr = std::max(rr, inner_prof[b]);
                        }
                        rr = std::max(rr, 1.0);
                        p = Point(coord_t(std::llround(centre.x() + rr * std::cos(th))),
                                  coord_t(std::llround(centre.y() + rr * std::sin(th))));
                    }
                }
            }
        }
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
                                                     std::string *reason, bool reverse_path, size_t layer_index, bool force_solid)
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
    }
    const std::vector<ContourChain> chains = contour_spiral(merged.front(), flow, layer_index, woven);
    if (chains.empty()) {
        if (reason != nullptr)
            *reason = "Continuous shell mode could not construct a contour spiral.";
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

    // A2 - connector spread: distribute the per-ring surface<->interior handoffs around the
    // circumference instead of stacking them into one deep seam groove (passed to
    // spiralize_rings). A3 - seam coast: thin the flow across those handoffs so the residual
    // mark is softer (applied in the flow assignment below). Both default 0 (off).
    const double connector_spread = (config != nullptr) ? config->continuous_shell_connector_spread.value / 100.0 : 0.0;
    const double coast_pct        = (config != nullptr) ? config->continuous_shell_seam_coast.value : 0.0;

    // Fill each connected piece (chain) as its own Fermat spiral, then stitch the pieces
    // together greedily - each next piece is appended at whichever end lands nearest the
    // running tail - so a split cross-section's pieces are joined by the shortest hop
    // (along the wall) instead of long chords across the bore. A simple cross-section has
    // exactly one chain, so this is identical to a single spiral there.
    struct Sub { Polyline pl; std::vector<char> flags; std::vector<coord_t> zoff; bool zwave = false; };
    std::vector<Sub> subs;
    // Per-layer reversal exists to give SINGLE-chain solid layers short inter-layer
    // joins (core end meets core start). A stitched MULTI-chain tour instead closes
    // back to its own seam-anchored start below, so every layer already starts and
    // ends at the (precessing) outer anchor - and reversing would break exactly that
    // pairing: odd layers would START at the buried inner ring while the previous
    // layer ended at the outer anchor, chording across any concave notch between.
    const bool reverse_chains = chains.size() > 1 ? false : reverse_path;
    for (const ContourChain &chain : chains) {
        std::vector<char> flags;
        std::vector<coord_t> zoff;
        Polyline sub = spiralize_rings(chain.loops, reverse_chains, seam_step_rad, seam_jitter_rad,
                                       ring_spacing, layer_index, chain.has_bore, &flags, interlock_amp, &zoff,
                                       connector_spread, chain.woven_infill ? woven.shell_loops : 0,
                                       &merged.front());
        if (sub.points.size() >= 2 && sub.points.size() == flags.size() && sub.points.size() == zoff.size())
            // A chain carries the Z-wave only when its interior was actually thinned AND its
            // rings passed the star-shape gate (angular phase is meaningless on a folded chain).
            subs.push_back({ std::move(sub), std::move(flags), std::move(zoff),
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
    auto append_sub = [&](Sub &s, bool rev) {
        if (rev) {
            std::reverse(s.pl.points.begin(), s.pl.points.end());
            std::reverse(s.flags.begin(), s.flags.end());
            std::reverse(s.zoff.begin(), s.zoff.end());
        }
        spiral.points.insert(spiral.points.end(), s.pl.points.begin(), s.pl.points.end());
        internal.insert(internal.end(), s.flags.begin(), s.flags.end());
        interlock_z.insert(interlock_z.end(), s.zoff.begin(), s.zoff.end());
        connector.insert(connector.end(), s.pl.points.size(), char(0));
        weavable.insert(weavable.end(), s.pl.points.size(), char(s.zwave ? 1 : 0));
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
        }
        spiral.points.push_back(head);
        internal.push_back(1);
        interlock_z.push_back(0);
        connector.push_back(1);
        weavable.push_back(0);
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
    double over_mult = 1.0;
    if (over_pct > 0.0 && have_flags && !merged.empty()) {
        double total_len = 0.0, buried_len = 0.0;
        for (size_t i = 1; i < spiral.points.size(); ++i) {
            const double seg = (spiral.points[i] - spiral.points[i - 1]).cast<double>().norm() * SCALING_FACTOR;
            total_len += seg;
            if (internal[i]) buried_len += seg;   // over-extrusion only lands on buried beads
        }
        const double q        = flow.mm3_per_mm();
        const double v_model  = std::abs(merged.front().area()) * SCALING_FACTOR * SCALING_FACTOR * layer_height;
        const double v_nom    = q * total_len;
        const double v_buried = q * buried_len;
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

    // Emit maximal runs of equal flow multiplier as sub-paths (keyed on the segment's
    // second point, matching the original over-extrusion convention), carrying interlock Z.
    ExtrusionMultiPath path;
    size_t start = 0;
    while (start + 1 < spiral.points.size()) {
        const double cur = fmul[start + 1];
        size_t end = start + 1;
        while (end + 1 < spiral.points.size() && fmul[end + 1] == cur)
            ++end;
        Polyline sub;
        sub.points.assign(spiral.points.begin() + start, spiral.points.begin() + end + 1);
        Polyline3 p3(sub);
        if (have_interlock)
            for (size_t i = 0; i < p3.points.size(); ++i)
                p3.points[i].z() = interlock_z[start + i];
        ExtrusionPath ep = make_path(p3, flow, layer_height);
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
