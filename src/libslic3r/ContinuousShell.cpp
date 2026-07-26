#include "ContinuousShell.hpp"

#include "ClipperUtils.hpp"
#include "Geometry/ConvexHull.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <mutex>
#include <vector>

namespace Slic3r {

namespace {

// Sparse WOVEN infill for a SOLID continuous-shell cross-section. Past the solid outer shell
// (wall_loops rings) the fill loops are thinned to a sparse spacing and (optionally) perturbed
// into a 3D interlocking weave - all carried by the SAME continuous bead. Driven by config:
//   continuous_shell_infill_density  0 => solid (disabled); else target interior density
//   continuous_shell_infill_weave    0 => plain concentric; 100 => full radial+Z weave
//   wall_loops                       => solid shell ring count
// Boolean env gate: ON only for a non-empty value other than "0". Some shells/harnesses
// "unset" a variable by passing it with an EMPTY value, which getenv() reports as set -
// a bare null check would silently flip experimental modes on.
bool cs_env_on(const char *name)
{
    const char *v = std::getenv(name);
    return v != nullptr && *v != '\0' && *v != '0';
}

struct WovenInfill {
    bool   enabled     = false;
    size_t shell_loops = 1;
    size_t stride      = 1;         // keep every Nth interior loop (1 = solid)
    double lobes       = 6.0;       // angular lobe count of the weave
    double weave       = 0.0;       // 0 = plain concentric, 1 = full weave (scales amp/zamp)
    double dphase      = 0.9;       // per-layer phase advance (rad)
    double zamp        = 0.0;       // Z amplitude, fraction of layer_height
    double zphase      = 1.5708;    // quadrature offset (rad) between the radial and Z waves
    bool   dyadic      = false;     // per-ring power-of-2 lobe hierarchy (bubble mode only)
    bool   cell_split  = false;     // pre-split deeply concave sections into near-star cells
    bool   bore_obstacles = false;  // SMALL bores are obstacles the global lattice swoops around,
                                    // not partition sources (the clean-puck / streamline mode)
    bool   bubble      = false;     // 3D bubble lattice (phase-locked + Z-carrier); config or CS_BUBBLE
    double bubble_h    = 0.0;       // bubble height (mm); 0 = auto (env CS_BUBBLE_H overrides)
    double bubble_pz   = 0.0;       // carrier period (layers); <4 = auto (env CS_BUBBLE_PZ overrides)
    double bubble_w    = 0.0;       // dyadic cell pitch (mm); 0 = auto (env CS_BUBBLE_W overrides)
    double cell_size   = 0.0;       // virtual-site cell target (mm); 0 = auto (4*stride beads)
    bool   adaptive_stride = true;  // small cells step DOWN to the largest hostable stride
                                    // (false: they print solid outright, the die-5 v7 behaviour)
    double solid_stagger = 0.0;     // radial inward shift (scaled) of a SOLID cap's interior rings this layer
};
WovenInfill woven_infill_params(const PrintRegionConfig *config)
{
    WovenInfill w;
    if (config == nullptr)
        return w;
    // Concave-cell split applies to SOLID parts too (it also prevents mid-march island
    // pinching), so it is read before the density gate below.
    w.cell_split = config->continuous_shell_cell_split.value || cs_env_on("CS_CELL_SPLIT");
    w.cell_size  = config->continuous_shell_cell_size.value;
    w.adaptive_stride = config->continuous_shell_adaptive_stride.value;
    w.bore_obstacles  = config->continuous_shell_bore_obstacles.value || cs_env_on("CS_BORE_OBSTACLES");
    // Bubble knobs: real config options; the env vars remain as global overrides so the
    // CLI A/B workflow keeps working (env value wins when set).
    w.bubble    = config->continuous_shell_bubble.value || cs_env_on("CS_BUBBLE");
    w.bubble_h  = config->continuous_shell_bubble_height.value;
    w.bubble_pz = double(config->continuous_shell_bubble_period.value);
    w.bubble_w  = config->continuous_shell_bubble_pitch.value;
    if (const char *ph = std::getenv("CS_BUBBLE_H"); ph != nullptr && std::atof(ph) > 0.0)
        w.bubble_h = std::atof(ph);
    if (const char *pz = std::getenv("CS_BUBBLE_PZ"); pz != nullptr && std::atof(pz) >= 4.0)
        w.bubble_pz = std::atof(pz);
    if (const char *pw = std::getenv("CS_BUBBLE_W"); pw != nullptr && std::atof(pw) > 0.1)
        w.bubble_w = std::atof(pw);
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
        // CS_BUBBLE (prototype): the bubble lattice needs the OPPOSITE of the per-layer rotation. Rotating
        // the lobes every layer means a crest never sits on the crest below it, so no vertical cell wall
        // can form - the print reads as decoherent wavy rings (measured: the m=lobes phase advances
        // exactly dphase per layer in sliced g-code). Bubbles want crests STACKED into welded columns,
        // with the vertical structure supplied by a slow Z-carrier on the amplitude (see contour_spiral),
        // so bubble mode pins the lattice phase across layers regardless of the configured shift.
        if (w.bubble)
            w.dphase = 0.0;
        // Dyadic cell hierarchy: a real (per-object overridable) config knob, with the env var
        // kept as a global force-on for quick experiments.
        w.dyadic = config->continuous_shell_infill_dyadic.value || cs_env_on("CS_BUBBLE_DYADIC");
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

// Resample a closed polygon to exactly N points spaced evenly by arc length.
// This gives adjacent concentric rings a common parametrization so we can morph
// smoothly between them.
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
// deep_concave: the caller determined the region has a genuine notch (deeper than
// ~2 bead widths past the convex hull). Rings then anchor by ARC FRACTION with
// gap-skipping joins instead of by absolute angle - see the anchor block below.
Polyline spiralize_rings(const Polygons &rings_in, bool reverse, double seam_step_rad, double seam_jitter_rad,
                         double ring_spacing, size_t layer_index, bool has_bore, std::vector<char> *internal,
                         coord_t interlock_amp = 0, std::vector<coord_t> *z_off = nullptr,
                         double connector_spread = 0.0, size_t infill_shell_loops = 0,
                         const ExPolygon *region = nullptr, std::vector<float> *fill_mul = nullptr,
                         bool deep_concave = false, bool thin_long_junctions = false,
                         const Polygons *junction_rails = nullptr)
{
    Polyline spiral;
    if (rings_in.empty())
        return spiral;
    const Polygons &rings = rings_in;    // original radial order: outer (0) -> inner

    size_t N = 0;
    for (const Polygon &r : rings)
        N = std::max(N, r.points.size());
    N = std::clamp<size_t>(N, size_t(64), size_t(256));

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

    // On a CONCAVE cross-section the angle->point map is discontinuous: a tiny seam
    // step can teleport the nearest-in-angle anchor across a notch, and since the
    // layer join lands on ring 0's anchor, consecutive LAYERS then join with a chord
    // through air (unfixable by in-layer containment - it happens between layers).
    // For concave regions anchor ring 0 by ARC FRACTION instead: map the seam angle
    // to a position along the ring's own perimeter, measured from a geometrically
    // stable origin (the x-extremal vertex - an extreme point is always on the convex
    // hull, so it moves continuously as the part tapers). Arc position is continuous
    // in the seam angle for ANY shape; convex regions keep the pure angle anchor
    // (identical behaviour, byte-stable output).
    std::vector<std::vector<Point>> grid;
    grid.reserve(rings.size());
    // Per-ring anchor twist expressed as a point-index offset (rings share N points,
    // evenly spaced by arc), used by the arc-anchored path to keep the Fermat
    // connectors spiraling inward exactly like the angle-anchored path does.
    const long twist_pts = std::lround(double(N) * ring_twist / (2.0 * M_PI));
    for (size_t ri = 0; ri < rings.size(); ++ri) {
        // Every ring must sweep the same rotational direction, otherwise morphing a
        // CCW ring (outer contour) to a CW ring (a hole contour, per Slic3r's
        // ExPolygon winding convention) maps each point to its antipode and the
        // interpolated chords cut straight across the bore (the "leaf" artifact).
        Polygon oriented = rings[ri];
        oriented.make_counter_clockwise();
        std::vector<Point> pts = resample_closed(oriented, N);
        if (pts.empty())
            return {};
        if (!deep_concave) {
            // Convex/gently-concave: the original absolute-angle anchor (byte-stable).
            anchor_to_angle(pts, ring_angle[ri]);
        } else if (ri == 0) {
            // Deep notch: the angle->point map is discontinuous (a small seam step can
            // teleport the anchor across the notch BETWEEN LAYERS, chording through
            // air). Anchor ring 0 by ARC FRACTION instead - position along the ring's
            // own perimeter from a stable origin (the x-extremal vertex; an extreme
            // point is always on the convex hull, so it moves continuously as the part
            // tapers). Arc position is continuous in the seam angle for ANY shape.
            size_t origin = 0;
            for (size_t k = 1; k < pts.size(); ++k)
                if (pts[k].x() > pts[origin].x() ||
                    (pts[k].x() == pts[origin].x() && pts[k].y() > pts[origin].y()))
                    origin = k;
            auto frac_of = [&](size_t layer) -> double {
                const double theta = double(layer) * seam_step_rad
                                   + seam_jitter_rad * seam_jitter_walk(layer);
                double f = std::fmod(theta, 2.0 * M_PI) / (2.0 * M_PI);
                return f < 0.0 ? f + 1.0 : f;
            };
            auto idx_of = [&](double frac) -> size_t {
                return (origin + size_t(std::llround(frac * double(N)))) % N;
            };
            // GAP SKIP: as the seam precesses across a slit / finger gap, the raw
            // inter-layer join (previous layer's anchor -> this layer's anchor) would
            // chord the gap through air. Advance this layer's start until the join
            // segment stays in material, so the seam JUMPS past the gap in one layer
            // instead of dragging through it. The previous anchor is recomputed with
            // the same rule (one level of lookback - the seam angle is a pure function
            // of layer index, so this stays deterministic and layer-local).
            auto adjusted_idx = [&](size_t prev_idx, size_t raw_idx) -> size_t {
                size_t idx = raw_idx;
                for (size_t step = 0; step < N; ++step) {
                    if (region->contains(Line(pts[prev_idx], pts[idx])))
                        return idx;
                    idx = (idx + 1) % N;
                }
                return raw_idx;   // nothing contained (degenerate) - keep raw
            };
            size_t start = idx_of(frac_of(layer_index));
            if (layer_index > 0) {
                size_t prev = idx_of(frac_of(layer_index - 1));
                if (layer_index > 1)
                    prev = adjusted_idx(idx_of(frac_of(layer_index - 2)), prev);
                start = adjusted_idx(prev, start);
            }
            if (start > 0)
                std::rotate(pts.begin(), pts.begin() + start, pts.end());
        } else {
            // Deep notch, inner rings: COUPLE each ring to the previous ring's anchor
            // positionally (nearest point is radially adjacent - rings are one bead
            // apart everywhere), then advance by the ring twist as an arc offset. This
            // keeps the inward connector spiral of the angle-anchored design without
            // re-introducing the discontinuous angle map (whose independent per-ring
            // anchors caused the gear/spiral regressions in the first experiment).
            const Point &prev_anchor = grid[ri - 1].front();
            size_t bestk = 0; double bd = std::numeric_limits<double>::max();
            for (size_t k = 0; k < pts.size(); ++k) {
                const double d = (pts[k] - prev_anchor).cast<double>().squaredNorm();
                if (d < bd) { bd = d; bestk = k; }
            }
            long s = (long(bestk) + twist_pts) % long(N);
            if (s < 0) s += long(N);
            if (s > 0)
                std::rotate(pts.begin(), pts.begin() + s, pts.end());
        }
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
    auto fmul = [&](float m) { if (fill_mul) fill_mul->push_back(m); };

    // Wrap-window fill (seam-valley compensation, ANNULAR chains). The arc-length ramp
    // spans the whole tour, so each ring climbs only its share of one layer height per
    // revolution. Between consecutive layers this leaves, once per wrap at each ring's
    // own start meridian, a slot (1 + share)*h tall instead of h - the bead there is
    // systematically underfilled by its ring's ramp share. On the surfaces this prints
    // as the helical seam valley; on buried rings it is a stack of aligned micro-voids
    // through the wall (a weak line). Compensate exactly: the points laid within the
    // first seam-step of arc after each ring's start (the ones under the tall slot) get
    // flow * (1 + share_r). Derived, not tuned - and zero when the tour has no ramp
    // discontinuity to heal (single ring).
    std::vector<double> ring_share(n_rings, 0.0);
    size_t wrap_pts = 0;
    if (has_bore && n_rings > 1) {
        double total_len = 0.0;
        std::vector<double> len(n_rings, 0.0);
        for (size_t ri = 0; ri < n_rings; ++ri) {
            const std::vector<Point> &c = grid[ri];
            for (size_t k = 0; k < c.size(); ++k)
                len[ri] += (c[(k + 1) % c.size()] - c[k]).cast<double>().norm();
            total_len += len[ri];
        }
        if (total_len > 0.0)
            for (size_t ri = 0; ri < n_rings; ++ri)
                ring_share[ri] = len[ri] / total_len;
        // Points within the first seam-step of arc; at least 2 because the emitter keys
        // a segment's flow on its SECOND point - a window must cover at least one full
        // segment past the ring start to deposit anything (a zero/tiny step still heals
        // the point-like discontinuity at the join itself).
        const double step = std::max(seam_step_rad, 0.0);
        wrap_pts = std::max<size_t>(2, size_t(std::llround(double(N) * step / (2.0 * M_PI))));
    }

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
        for (size_t k = 0; k < N; ++k) { spiral.points.push_back(grid[0][k]); flag(false); zoff(0); fmul(1.f); }
        spiral.points.push_back(grid[0][0]); flag(false); zoff(0); fmul(1.f);
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

    spiral.points.reserve(n_rings * N);
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
                // Candidates by distance, and take the NEAREST WHOSE JUNCTION STAYS IN
                // MATERIAL - the plain nearest point is not enough when an obstacle sits
                // between the rings: two swooped rings hugging opposite flanks of a bore
                // have their nearest points straight across the hole (the swoop-mode ring
                // junction chorded the bore). Falls back to the nearest point when nothing
                // within the scan is contained (degenerate slivers) - the old behaviour.
                std::vector<size_t> ks(c.size());
                for (size_t k = 0; k < c.size(); ++k) ks[k] = k;
                std::sort(ks.begin(), ks.end(), [&](size_t a2, size_t b2) {
                    return (c[a2] - prev).cast<double>().squaredNorm() < (c[b2] - prev).cast<double>().squaredNorm();
                });
                size_t bestk = ks.front();
                bool   found = false;
                for (size_t t2 = 0; t2 < ks.size(); ++t2)     // full scan: the nearest stretch of
                    if (region->contains(Line(prev, c[ks[t2]]))) { bestk = ks[t2]; found = true; break; }
                                                              // candidates can ALL sit across a bore
                if (!found && junction_rails != nullptr && !junction_rails->empty()) {
                    // No straight junction stays in material - both rings hug the same bore and
                    // every chord grazes the hole. Ride AROUND it on the fan-corridor rail
                    // (mirrors the inter-chain connector router): nearest rail vertices to the
                    // two ends, shorter arc between them, deposited at half flow like any ride.
                    const Point tgt = c[bestk];
                    size_t br = 0, ia = 0, ib = 0; double bs = std::numeric_limits<double>::max();
                    for (size_t r2 = 0; r2 < junction_rails->size(); ++r2) {
                        const Points &rp = (*junction_rails)[r2].points;
                        if (rp.size() < 3) continue;
                        size_t ka = 0, kb = 0; double da = 1e300, db2 = 1e300;
                        for (size_t k = 0; k < rp.size(); ++k) {
                            const double d1 = (rp[k] - prev).cast<double>().squaredNorm();
                            const double d2 = (rp[k] - tgt).cast<double>().squaredNorm();
                            if (d1 < da) { da = d1; ka = k; }
                            if (d2 < db2) { db2 = d2; kb = k; }
                        }
                        const double sum = std::sqrt(da) + std::sqrt(db2);
                        if (sum < bs) { bs = sum; br = r2; ia = ka; ib = kb; }
                    }
                    if (bs < 1e300 && ia != ib) {
                        const Points &rp = (*junction_rails)[br].points;
                        const size_t rn = rp.size();
                        double fwd = 0., bwd = 0.;
                        for (size_t k = ia; k != ib; k = (k + 1) % rn) fwd += (rp[(k + 1) % rn] - rp[k]).cast<double>().norm();
                        for (size_t k = ia; k != ib; k = (k + rn - 1) % rn) bwd += (rp[(k + rn - 1) % rn] - rp[k]).cast<double>().norm();
                        if (fwd <= bwd)
                            for (size_t k = ia; ; k = (k + 1) % rn) {
                                spiral.points.push_back(rp[k]); flag(true); zoff(0); fmul(0.5f);
                                if (k == ib) break;
                            }
                        else
                            for (size_t k = ia; ; k = (k + rn - 1) % rn) {
                                spiral.points.push_back(rp[k]); flag(true); zoff(0); fmul(0.5f);
                                if (k == ib) break;
                            }
                    }
                }
                if (bestk > 0)
                    std::rotate(c.begin(), c.begin() + bestk, c.end());
            }
        }
        // Long ring-to-ring junctions ride over already-printed beads (a swooped chain's
        // rerouted junctions hop the fan corridor beside a bore, every layer, at the same
        // spots) - deposit them at HALF flow like the inter-chain connectors, keyed on the
        // junction segment's second point (ring point 0). Opt-in per chain (swoop mode):
        // classic chains keep their full-flow Fermat connectors (the gear's radial spokes).
        const bool long_junction = thin_long_junctions && !spiral.points.empty() &&
                                   (c.front() - spiral.points.back()).cast<double>().norm() > 2.0 * ring_spacing;
        for (size_t k = 0; k < N; ++k) {
            spiral.points.push_back(c[k]);
            flag(internal_ring);
            zoff(zo);
            fmul(k == 0 && long_junction ? 0.5f
                                         : (k < wrap_pts ? float(1.0 + ring_share[idx]) : 1.f));
        }
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
        // A wall band must keep exactly one hole. TWO+ holes is unmodeled topology (the
        // upstream partitions own it) - end the chain. ZERO holes means the annulus
        // BREACHED: an eccentric bore ate through the thin side of its cell and the wall
        // opened into a solid crescent - keep marching it as SOLID. Breaking there (the
        // old behaviour, written for centred annuli whose pinch IS the medial finish)
        // silently discarded the whole crescent: the counterbored-disc gauntlet's solid
        // cap layers lost ~60% of their material to exactly this.
        if (has_hole && band.holes.size() != 1) {
            if (!band.holes.empty())
                break;
            has_hole = false;
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

// ---- Multi-bore partition (task #4 prototype, gated on CS_MULTIBORE) -----------------------
// march_bands models a SINGLE-hole annular wall; a cross-section with >1 hole trips its
// holes.size()!=1 guard and degrades to the medial-axis fallback (which smears fill across the
// bores). Split such a section into one single-bore cell per hole by a Voronoi partition of the
// hole centroids: each cell is the region closer to its own hole than to any other, so it holds
// exactly one bore and marches as an ordinary annulus. The cut lines are the perpendicular
// bisectors between hole centroids - they depend only on the hole positions, which drift smoothly
// up the part, so the partition is stable layer-to-layer (no per-frame flicker).

// Sutherland-Hodgman clip of `poly` to the halfplane of points at least as close to `keep` as to
// `other` (i.e. the `keep` side of the two sites' perpendicular bisector).
static Polygon clip_to_halfplane(const Polygon &poly, const Point &keep, const Point &other)
{
    if (poly.points.size() < 3)
        return Polygon{};
    const Vec2d  n   = (other - keep).cast<double>();        // bisector normal, away from `keep`
    const Vec2d  mid = 0.5 * (keep + other).cast<double>();
    const double c   = n.dot(mid);                           // halfplane: n.x <= c is closer to keep
    Polygon out;
    const size_t m = poly.points.size();
    for (size_t i = 0; i < m; ++i) {
        const Vec2d a = poly.points[i].cast<double>();
        const Vec2d b = poly.points[(i + 1) % m].cast<double>();
        const double da = n.dot(a) - c, db = n.dot(b) - c;
        const bool ai = da <= 0.0, bi = db <= 0.0;
        if (ai)
            out.points.emplace_back(coord_t(std::lround(a.x())), coord_t(std::lround(a.y())));
        if (ai != bi) {                                      // edge crosses the bisector: add cut point
            const double t = da / (da - db);
            const Vec2d  x = a + t * (b - a);
            out.points.emplace_back(coord_t(std::lround(x.x())), coord_t(std::lround(x.y())));
        }
    }
    return out;
}

// One concave-split or bore-bisector cut, recorded so the weave can treat the interface like
// intra-cell neighbours (see the interface field in contour_spiral). Declared here because
// partition_bores also records candidates when the modern cell pipeline is active.
struct CutLine;

// Partition a multi-hole seed into single-bore cells. Any cell that still carries >1 hole (two
// bores so close they share a Voronoi cell) is recursed on until each piece holds one bore.
// With `lines` set (modern cell pipeline), every site-pair bisector is recorded as a CANDIDATE
// cut line - extents zeroed; modernize_bore_cuts later measures the real shared edges, drops
// non-neighbour pairs, recesses the cells and registers the survivors for the interface weave.
static void partition_bores(const ExPolygon &seed, ExPolygons &out, int depth, std::vector<CutLine> *lines);
static void partition_bores(const ExPolygon &seed, ExPolygons &out, int depth = 0)
{
    partition_bores(seed, out, depth, nullptr);
}

// Concave-cell split (the woven/bubble generalisation to arbitrary outlines). The polar weave
// needs each chain to be near-STAR-SHAPED about its own centroid - lobes are laid along r(theta)
// and the anti-phase contacts assume a single-valued radius - and a deeply concave outline also
// pinches into disconnected islands mid-march (the TopPlate bird-nest). Rather than weave across
// a bad shape, split it FIRST along geometry-derived lines: find the contour's deepest concavity
// (max distance from a contour vertex to the convex-hull boundary), cut the cell with the line
// through that notch apex whose normal is the hull-bridge direction (the classic resolve-the-notch
// ACD heuristic: a dogbone gets vertical neck cuts that free the lobes, an L is cut along the
// diagonal through its inner corner), and recurse until every piece's residual concavity is
// shallow. The cut depends only on the outline, which drifts smoothly up the part, so cells are
// stable layer-to-layer - which the bubble Z-carrier needs (cells must stack vertically). Each
// piece then marches as its own chain: its rings hug the cut like a wall (the zero-width internal
// wall), it gets its own centroid/rosette/morph, and the connector router joins the cells into
// the one continuous bead exactly as it already does for multi-bore cells.

struct ContourNotch { Vec2d apex = Vec2d(0., 0.); Vec2d bridge = Vec2d(1., 0.); double depth = 0.0; };

// One concave-split cut, recorded so the weave can treat the interface like intra-cell
// neighbours instead of a dead wall: origin/tangent/normal of the cut line plus the parent
// cell's extent along it (the field tapers off beyond the real shared boundary, so the cut's
// ends - which land on the object surface - never wave the skin).
struct CutLine { Vec2d o = Vec2d(0., 0.); Vec2d t = Vec2d(1., 0.); Vec2d n = Vec2d(0., 1.); double smin = 0.0; double smax = 0.0; };

static void partition_bores(const ExPolygon &seed, ExPolygons &out, int depth, std::vector<CutLine> *lines)
{
    if (seed.holes.size() <= 1 || depth > 6) {               // already single-bore (or give up)
        out.push_back(seed);
        return;
    }
    std::vector<Point> sites;
    sites.reserve(seed.holes.size());
    for (const Polygon &h : seed.holes)
        sites.push_back(h.centroid());

    if (lines != nullptr)
        for (size_t i = 0; i < sites.size(); ++i)
            for (size_t j = i + 1; j < sites.size(); ++j) {
                const Vec2d a = sites[i].cast<double>(), b = sites[j].cast<double>();
                const Vec2d d = b - a;
                const double len = d.norm();
                if (len < 1.0) continue;
                CutLine cl;
                cl.o = 0.5 * (a + b);
                cl.n = d / len;
                cl.t = Vec2d(-cl.n.y(), cl.n.x());
                lines->push_back(cl);                        // candidate; extents measured later
            }

    const BoundingBox bb  = seed.contour.bounding_box();
    const coord_t     pad = scale_(5.0);
    const Point mn(bb.min.x() - pad, bb.min.y() - pad), mx(bb.max.x() + pad, bb.max.y() + pad);
    Polygon box;
    box.points = { mn, Point(mx.x(), mn.y()), mx, Point(mn.x(), mx.y()) };

    for (size_t i = 0; i < sites.size(); ++i) {
        Polygon cell = box;
        for (size_t j = 0; j < sites.size() && cell.points.size() >= 3; ++j)
            if (j != i)
                cell = clip_to_halfplane(cell, sites[i], sites[j]);
        if (cell.points.size() < 3)
            continue;
        for (ExPolygon &piece : intersection_ex(seed, Polygons{ cell }))
            partition_bores(piece, out, depth + 1, lines);   // a piece may still hold >1 bore
    }
}

// Efficient cell placement for MULTI-BORE sections (die-5 lesson). Voronoi over the bore
// centroids alone forces the bore cells to annex every far-flung region; the leftovers are
// slivers and corners that no lattice fits - they printed as stranded voids (v4) or solid
// slabs (v6). Instead AUGMENT the site set: farthest-point sampling adds a VIRTUAL site
// wherever material lies farther than ~4 sparse strides from every site, so open regions
// become their own compact SOLID cells with healthy rosettes and every bore keeps a tight
// annulus. Deterministic grid scan order keeps placement stable layer-to-layer. Single-bore
// seeds are untouched (a uniform annulus - the gear - wants no extra cells; the eccentric
// case is handled by the eccentricity split).
static void partition_cells(const ExPolygon &seed, ExPolygons &out, coord_t spacing, double target,
                            std::vector<CutLine> *lines)
{
    std::vector<Point> sites;
    sites.reserve(seed.holes.size() + 8);
    for (const Polygon &h : seed.holes)
        sites.push_back(h.centroid());

    const BoundingBox bb = seed.contour.bounding_box();
    const coord_t     gs = coord_t(2 * spacing);
    std::vector<Point> grid;
    for (coord_t gx = bb.min.x() + gs / 2; gx <= bb.max.x(); gx += gs)
        for (coord_t gy = bb.min.y() + gs / 2; gy <= bb.max.y(); gy += gs) {
            const Point g(gx, gy);
            if (seed.contains(g))
                grid.push_back(g);
        }
    // The cap is a runaway guard, not a coverage budget: it must sit ABOVE what any sane part
    // needs, or the last farthest gaps silently go uncovered (die-5: the rim consumed ~18 sites
    // and one of the four bore-diagonal gaps hit the old cap of 24 - its cell printed with a
    // 2.5mm void wedge while the other three diagonals got their rosette cells).
    for (int add = 0; add < 64; ++add) {
        double best_d = 0.0;
        Point  best;
        bool   have = false;
        for (const Point &g : grid) {
            double dmin = std::numeric_limits<double>::max();
            for (const Point &s2 : sites)
                dmin = std::min(dmin, (g - s2).cast<double>().norm());
            if (dmin > best_d) { best_d = dmin; best = g; have = true; }
        }
        if (!have || best_d <= target)
            break;
        sites.push_back(best);
    }
    if (sites.size() <= 1) {
        out.push_back(seed);
        return;
    }

    if (lines != nullptr)
        for (size_t i = 0; i < sites.size(); ++i)
            for (size_t j = i + 1; j < sites.size(); ++j) {
                const Vec2d a = sites[i].cast<double>(), b = sites[j].cast<double>();
                const Vec2d d = b - a;
                const double len = d.norm();
                if (len < 1.0) continue;
                CutLine cl;
                cl.o = 0.5 * (a + b);
                cl.n = d / len;
                cl.t = Vec2d(-cl.n.y(), cl.n.x());
                lines->push_back(cl);
            }

    const coord_t pad = scale_(5.0);
    const Point mn(bb.min.x() - pad, bb.min.y() - pad), mx(bb.max.x() + pad, bb.max.y() + pad);
    Polygon box;
    box.points = { mn, Point(mx.x(), mn.y()), mx, Point(mn.x(), mx.y()) };

    for (size_t i = 0; i < sites.size(); ++i) {
        Polygon cell = box;
        for (size_t j = 0; j < sites.size() && cell.points.size() >= 3; ++j)
            if (j != i)
                cell = clip_to_halfplane(cell, sites[i], sites[j]);
        if (cell.points.size() < 3)
            continue;
        for (ExPolygon &piece : intersection_ex(seed, Polygons{ cell })) {
            if (piece.holes.size() > 1)
                partition_bores(piece, out, 1, lines);   // two bores in one cell: classic recursion
            else
                out.push_back(std::move(piece));
        }
    }
}

// Shave the half-bead tapered recess strip off `piece` on side `sgn` of the cut (full standoff
// mid-span, bevelled to zero at the extent ends so skins/junctions stay continuous). Shared by
// the concave split and the bore-bisector modernization.
static ExPolygons shave_cut_strip(const ExPolygon &piece, const CutLine &cl, double sgn, coord_t spacing)
{
    const double h2    = 0.5 * double(spacing);
    const double bevel = std::min(4.0 * double(spacing), (cl.smax - cl.smin) / 3.0);
    auto pt = [&](double sv, double dv) {
        const Vec2d p = cl.o + sv * cl.t + dv * cl.n;
        return Point(coord_t(std::llround(p.x())), coord_t(std::llround(p.y())));
    };
    Polygon strip;
    if (sgn > 0)
        strip.points = { pt(cl.smin, 0.), pt(cl.smax, 0.), pt(cl.smax - bevel, h2), pt(cl.smin + bevel, h2) };
    else
        strip.points = { pt(cl.smin, 0.), pt(cl.smin + bevel, -h2), pt(cl.smax - bevel, -h2), pt(cl.smax, 0.) };
    ExPolygons out;
    for (ExPolygon &piece2 : diff_ex(ExPolygons{ piece }, Polygons{ strip })) {
        if (std::abs(piece2.contour.area()) < 4.0 * double(spacing) * double(spacing))
            continue;
        out.push_back(std::move(piece2));
    }
    return out;
}

// Strip sub-bead sliver holes from partitioned cells. Cut shaving and halfplane clipping can
// leave a hair-thin closed slot inside a cell (effective width 2A/P under a bead); the march
// would keep it hollow as a fake bore and lattice a void wedge around unprintable air. Real
// bores are round - their effective width equals their radius, always at least a bead.
static void strip_sliver_holes(ExPolygons &cells, coord_t spacing)
{
    const double s = double(spacing);
    for (ExPolygon &cell : cells)
        cell.holes.erase(std::remove_if(cell.holes.begin(), cell.holes.end(),
                                        [&](const Polygon &h) {
                                            const double a   = std::abs(h.area());
                                            const double per = h.length();
                                            const double aw  = per > 1.0 ? 2.0 * a / per : 0.0;
                                            return aw < 0.75 * s || a < 4.0 * s * s;
                                        }),
                         cell.holes.end());
}

// Modernize the bore-bisector cuts (cell_split pipeline): measure each candidate line's REAL
// shared edge from the cells' on-line boundary points, drop non-neighbour pairs, merge collinear
// duplicates, recess every adjacent cell by the tapered strip (kills the coincident 200%-flow
// bisector wall line), and register the survivors so the interface weave fields them exactly
// like concave cuts - the bore-cell walls then kiss/fuse/breathe instead of standing dead flat.
static void modernize_bore_cuts(ExPolygons &cells, const std::vector<CutLine> &candidates,
                                coord_t spacing, std::vector<CutLine> &out_cuts)
{
    const double on_tol = 0.35 * double(spacing);
    // A cell EDGE (adjacent vertex pair) lying on the line - isolated on-line VERTICES are not
    // enough: Voronoi vertices of NON-neighbour site pairs also land on those pairs' bisectors
    // (e.g. the centre cell's opposite corners on a corner-to-corner diagonal), and treating
    // them as a shared edge shaves channels straight through cell interiors.
    auto online_edge_extent = [&](const ExPolygon &cell, const CutLine &cl, double &smin, double &smax) {
        bool any = false;
        const Points &pts = cell.contour.points;
        for (size_t i = 0; i < pts.size(); ++i) {
            const Vec2d a = pts[i].cast<double>() - cl.o;
            const Vec2d b = pts[(i + 1) % pts.size()].cast<double>() - cl.o;
            if (std::abs(a.dot(cl.n)) >= on_tol || std::abs(b.dot(cl.n)) >= on_tol)
                continue;
            const double sa = a.dot(cl.t), sb = b.dot(cl.t);
            if (std::abs(sa - sb) < 2.0 * double(spacing))
                continue;                                   // a sliver touch, not a real edge
            smin = std::min({ smin, sa, sb });
            smax = std::max({ smax, sa, sb });
            any = true;
        }
        return any;
    };
    for (const CutLine &cand : candidates) {
        CutLine cl = cand;
        double smin = std::numeric_limits<double>::max(), smax = std::numeric_limits<double>::lowest();
        bool any = false;
        for (const ExPolygon &cell : cells)
            any = online_edge_extent(cell, cl, smin, smax) || any;
        if (!any || smax - smin < 4.0 * double(spacing))
            continue;                                       // not a real shared edge
        cl.smin = smin; cl.smax = smax;
        bool dup = false;
        for (CutLine &e : out_cuts)
            if (std::abs(e.n.dot(cl.n)) > 0.999 && std::abs(e.n.dot(cl.o - e.o)) < 0.5 * double(spacing)) {
                const double s0 = cl.smin * cl.t.dot(e.t) + (cl.o - e.o).dot(e.t);
                const double s1 = cl.smax * cl.t.dot(e.t) + (cl.o - e.o).dot(e.t);
                e.smin = std::min({ e.smin, s0, s1 });
                e.smax = std::max({ e.smax, s0, s1 });
                dup = true; break;
            }
        if (dup)
            continue;
        ExPolygons next;
        next.reserve(cells.size() + 2);
        for (ExPolygon &cell : cells) {
            double ca = std::numeric_limits<double>::max(), cb = std::numeric_limits<double>::lowest();
            if (!online_edge_extent(cell, cl, ca, cb)) {       // shave only cells sharing a real edge
                next.push_back(std::move(cell));
                continue;
            }
            // Shave with the CELL'S OWN on-line extent, not the candidate's union across all
            // cells: grid-aligned virtual sites make DIFFERENT pairs' bisectors collinear, so
            // the union can span far past this cell's shared edge - and where the boundary
            // turns off the line, the over-long strip dives through the cell INTERIOR, cutting
            // a sliver slot that closes into a phantom hole the march then keeps hollow as a
            // fake bore (die-5: 0.2-0.8mm2 sliver holes, one per over-shaved cell).
            CutLine ccl = cl;
            ccl.smin = std::max(cl.smin, ca);
            ccl.smax = std::min(cl.smax, cb);
            if (ccl.smax - ccl.smin < 2.0 * double(spacing)) {
                next.push_back(std::move(cell));
                continue;
            }
            const double sgn = (cell.contour.centroid().cast<double>() - cl.o).dot(cl.n) >= 0.0 ? 1.0 : -1.0;
            for (ExPolygon &piece : shave_cut_strip(cell, ccl, sgn, spacing))
                next.push_back(std::move(piece));
        }
        cells = std::move(next);
        out_cuts.push_back(cl);
    }
}

// Concavities of `contour` at least `min_depth` inside the convex-hull boundary, deepest first:
// each is a vertex, its depth, and the direction of the hull edge it hides behind (the pocket
// mouth). Capped to the 8 deepest - the split only ever consumes a couple per recursion.
static std::vector<ContourNotch> deep_notches(const Polygon &contour, double min_depth)
{
    std::vector<ContourNotch> out;
    if (contour.points.size() < 4)
        return out;
    const Polygon hull = Geometry::convex_hull(contour.points);
    const size_t  hn   = hull.points.size();
    if (hn < 3)
        return out;
    for (const Point &p : contour.points) {
        const Vec2d v = p.cast<double>();
        double d_min = std::numeric_limits<double>::max();
        Vec2d  e_dir(1., 0.);
        for (size_t i = 0; i < hn; ++i) {
            const Vec2d a  = hull.points[i].cast<double>();
            const Vec2d ab = hull.points[(i + 1) % hn].cast<double>() - a;
            const double len2 = ab.squaredNorm();
            if (len2 < 1.0)
                continue;
            const double t = std::clamp((v - a).dot(ab) / len2, 0.0, 1.0);
            const double d = (v - (a + t * ab)).norm();
            if (d < d_min) { d_min = d; e_dir = ab / std::sqrt(len2); }
        }
        if (d_min < std::numeric_limits<double>::max() && d_min >= min_depth)
            out.push_back({ v, e_dir, d_min });
    }
    std::sort(out.begin(), out.end(), [](const ContourNotch &a, const ContourNotch &b) { return a.depth > b.depth; });
    if (out.size() > 8)
        out.resize(8);
    return out;
}

// Split `cell` at its deepest notch until every piece is near-star-shaped. The threshold is
// generous (8 beads / 2.5mm) so tooth-scale concavities - gear flanks, keyways - never split:
// those already fill correctly via the pinch/merge path, and needless cuts cost connectors.
// The cut reuses clip_to_halfplane by synthesising two sites straddling the wanted line. A cell
// carrying a bore keeps it whole: if the cut would cross the bore, try the orthogonal line
// through the same apex, else leave the cell unsplit.
static void partition_concave(const ExPolygon &cell, coord_t spacing, ExPolygons &out, int depth = 0,
                              std::vector<CutLine> *cuts = nullptr)
{
    const double min_depth = std::max(8.0 * double(spacing), double(scale_(2.5)));
    std::vector<ContourNotch> notches = deep_notches(cell.contour, min_depth);
    // ECCENTRIC-ANNULUS trigger (die-5 lesson): a single-bore cell whose annulus width varies
    // wildly - the bore hugging one side of a wide cell - marches badly: the bore-side front
    // stops when the NARROW side pinches, stranding the bore with no nearby rings on its far
    // side (measured: a ~2mm void annulus bridged only by connector strands). Split such a
    // cell perpendicular to the bore->far direction, placed to keep a compact annulus around
    // the bore; the far chunk becomes its own solid cell with a healthy rosette. The gear's
    // tooth-scale variation (~1.5mm) and a centred bore stay far below the ~10-bead threshold.
    if (std::getenv("CS_DBG_ECC") != nullptr)
        std::fprintf(stderr, "[ecc] depth=%d holes=%zu notches=%zu area=%.1f\n",
                     depth, cell.holes.size(), notches.size(),
                     std::abs(cell.contour.area()) * 1e-12);
    if (notches.empty() && depth <= 4 && cell.holes.size() == 1) {
        const Vec2d b = cell.holes.front().centroid().cast<double>();
        double hr = 0.0;
        for (const Point &q : cell.holes.front().points)
            hr = std::max(hr, (q.cast<double>() - b).norm());
        double dmin = std::numeric_limits<double>::max(), dmax = 0.0;
        Vec2d far = b;
        for (const Point &q : cell.contour.points) {
            const double dd = (q.cast<double>() - b).norm();
            dmin = std::min(dmin, dd);
            if (dd > dmax) { dmax = dd; far = q.cast<double>(); }
        }
        if (std::getenv("CS_DBG_ECC") != nullptr)
            std::fprintf(stderr, "[ecc]   hr=%.2f dmin=%.2f dmax=%.2f thr=%.2f\n",
                         hr * 1e-6, dmin * 1e-6, dmax * 1e-6, 10.0 * double(spacing) * 1e-6);
        if (dmax - dmin > 10.0 * double(spacing)) {
            const Vec2d dir = (far - b).normalized();
            const double lo = hr + 8.0 * double(spacing);
            const double hi = dmax - 6.0 * double(spacing);
            if (lo < hi) {
                ContourNotch ecc;
                ecc.apex   = b + dir * std::clamp(0.5 * (hr + dmax), lo, hi);
                ecc.bridge = dir;                        // downstream: cut normal = bridge => cut ⊥ bore->far
                ecc.depth  = dmax - dmin;
                notches.assign(1, ecc);
            }
        }
    }
    if (notches.empty() || depth > 4) {
        out.push_back(cell);
        return;
    }
    const ContourNotch &n = notches.front();
    // Prefer PAIRING two opposing notches: a dogbone's two neck indentations resolve into the one
    // vertical neck cut, and a plus's opposing inner corners into a diagonal through the centre.
    // Without pairing, the second notch of a symmetric pair ends up resolved against a DIAGONAL
    // hull-bridge edge created by the first cut - valid cells, but skewed walls that cost
    // perimeter. Pair when the second apex is well separated and the connecting segment stays
    // inside the material (sampled); the cut line then passes through BOTH apexes.
    Vec2d nrm = n.bridge;                       // default: cut normal = bridge dir (across the pocket)
    {
        // The connecting segment must cross real INTERIOR - two notches of the same pocket
        // connect along its floor, and on-boundary points pass a plain contains() - so sample
        // against the cell eroded by one bead. Among passing pairs take the SHORTEST segment:
        // that is the thinnest waist (a dogbone's neck), not a valid-but-skewed diagonal.
        const ExPolygons interior = offset_ex(ExPolygons{ cell }, -float(spacing));
        double best_len = std::numeric_limits<double>::max();
        for (size_t k = 1; k < notches.size() && !interior.empty(); ++k) {
            const Vec2d seg = notches[k].apex - n.apex;
            const double len = seg.norm();
            if (len < 4.0 * double(spacing) || len >= best_len)
                continue;
            bool inside = true;
            for (double t : { 0.2, 0.4, 0.5, 0.6, 0.8 }) {
                const Vec2d q = n.apex + t * seg;
                const Point qq(coord_t(std::llround(q.x())), coord_t(std::llround(q.y())));
                bool in = false;
                for (const ExPolygon &ip : interior)
                    if (ip.contains(qq)) { in = true; break; }
                if (!in) { inside = false; break; }
            }
            if (inside) {
                best_len = len;
                nrm = Vec2d(-seg.y() / len, seg.x() / len);  // line through both apexes
            }
        }
    }
    if (!cell.holes.empty()) {
        const Vec2d hc = cell.holes.front().centroid().cast<double>();
        double hr = 0.0;
        for (const Point &q : cell.holes.front().points)
            hr = std::max(hr, (q.cast<double>() - hc).norm());
        auto clears = [&](const Vec2d &nn) { return std::abs(nn.dot(hc - n.apex)) > hr + double(spacing); };
        if (!clears(nrm)) {
            const Vec2d ortho(-nrm.y(), nrm.x());
            if (!clears(ortho)) { out.push_back(cell); return; }
            nrm = ortho;
        }
    }
    // The cut as a CutLine (extent = the cell's span along the line), recorded for the
    // interface weave when requested and used for the tapered recess either way.
    CutLine cl;
    cl.o = n.apex; cl.n = nrm; cl.t = Vec2d(-nrm.y(), nrm.x());
    {
        double smin = std::numeric_limits<double>::max(), smax = std::numeric_limits<double>::lowest();
        for (const Point &q : cell.contour.points) {
            const double sv = (q.cast<double>() - n.apex).dot(cl.t);
            smin = std::min(smin, sv); smax = std::max(smax, sv);
        }
        cl.smin = smin; cl.smax = smax;
    }
    if (cuts != nullptr)
        cuts->push_back(cl);

    // Recess each side HALF A BEAD from the cut via the shared tapered strip (see
    // shave_cut_strip): full standoff along the interior span - facing rings land one bead
    // apart, edges tiling the line, no 200%-flow double extrusion - bevelled to zero toward
    // the ends, so the two skins converge and meet exactly at the mouth and the perimeter
    // stays continuous.
    const double d = double(scale_(1.0));
    auto site = [&](double off) {
        return Point(coord_t(std::llround(n.apex.x() + off * nrm.x())), coord_t(std::llround(n.apex.y() + off * nrm.y())));
    };

    const BoundingBox bb  = cell.contour.bounding_box();
    const coord_t     pad = scale_(5.0);
    const Point mn(bb.min.x() - pad, bb.min.y() - pad), mx(bb.max.x() + pad, bb.max.y() + pad);
    Polygon box;
    box.points = { mn, Point(mx.x(), mn.y()), mx, Point(mn.x(), mx.y()) };

    ExPolygons pieces;
    for (int side = 0; side < 2; ++side) {
        const double sgn  = side == 0 ? -1.0 : 1.0;
        const Polygon plane = clip_to_halfplane(box, side == 0 ? site(-d) : site(+d), side == 0 ? site(+d) : site(-d));
        if (plane.points.size() < 3)
            continue;
        for (ExPolygon &clipped : intersection_ex(cell, Polygons{ plane }))
            for (ExPolygon &piece : shave_cut_strip(clipped, cl, sgn, spacing))
                pieces.push_back(std::move(piece));
    }
    if (pieces.size() < 2) {                    // degenerate cut: everything landed on one side
        if (cuts != nullptr && !cuts->empty())
            cuts->pop_back();                   // no split happened: drop the recorded cut
        out.push_back(cell);
        return;
    }
    for (ExPolygon &piece : pieces)
        partition_concave(piece, spacing, out, depth + 1, cuts);
}

// One connected piece of a (possibly split) cross-section: its concentric loops in
// radial order (outer contours outer->medial, then, if annular, hole contours
// medial->bore) plus whether that piece is annular (so spiralize_rings knows its
// innermost loop is a bore surface, not a buried centre).
struct ContourChain {
    Polygons loops;
    bool     has_bore = false;
    bool     woven_infill = false;   // true once this chain's interior was thinned to sparse infill
    int      peel_depth   = 0;       // 0 = no peel; 1 = ext ring peeled (cell keeps 1 outer shell);
                                     // 2 = full wall stack peeled (cell runs SHELL-FREE outside:
                                     //     its outermost structure is the bucket-0 woven ring)
    bool     from_solid   = false;   // marched piece had NO real holes: any band holes are the
                                     // convergence core, not a bore (degenerate rings droppable)
    bool     peel_band    = false;   // the global peeled-skin chain (never sliver-filtered)
    Polygons obstacles;              // bore keep-out outlines (CCW) the lattice must swoop around
                                     // (bore-obstacle mode; polygons, not circles, so non-round
                                     // bores work unchanged). Empty on every other chain.
    Polygons ride_rails;             // preferred connector-ride lanes (swoop mode): rings in the
                                     // reserved corridor between each bore's outermost fan arc and
                                     // the crest-flatten boundary. The default hole rails sit ON
                                     // the bore wall crown - rides there mounded ~a full extra
                                     // half-flow lap per layer onto the wall (die5_swoop print).
    double   footprint    = 0.;      // MATERIAL area of the marched bands (scaled^2). The sliver
                                     // filter must compare material, not contour extent: the peel
                                     // band's contour spans the whole part while its material is
                                     // a thin annulus, and compact virtual-site cells are small by
                                     // design - contour-area ratios mis-drop both.
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

    // Multi-bore partition is ON by default (a multi-hole section otherwise degrades to the
    // medial fallback, which smears fill across the bores). CS_NO_MULTIBORE opts out for A/B.
    const bool multibore = !cs_env_on("CS_NO_MULTIBORE");
    // A single-bore seed whose wall PINCHES and splits mid-march (sharp teeth, an eccentric keyed
    // bore, a herringbone base) fragments into an outer full-thickness band + a separately-marched
    // interior, with the tooth tips pinched off as their own tiny SOLID fragments. Each annular
    // fragment then holds only PART of the wall (the outer band keeps the outermost AND the
    // innermost/bore ring, the middle stranded elsewhere), so its loop list has a radial GAP that
    // the Fermat interleave chords straight across - two full-width surface->bore beads per layer
    // that pile onto the seam and bulge the exposed face (the gear top-face artifact). Since the
    // piece has exactly ONE bore, ALL its annular fragments are radially-nested pieces of that one
    // wall: merge them into one chain (concatenate every band, sort outer->in by |contour area|) so
    // the converter's outer-run/bore-run grouping yields a contiguous monotone stack again - what an
    // unsplit march would emit. Solid fragments (pinched-off tooth tips) stay their own chains. An
    // unsplit piece has one annular fragment already sorted outer->in, so this is byte-identical
    // there. Applied to a plain seed AND to each single-bore Voronoi cell of a multi-bore section.
    std::vector<char> chain_thin;                  // parallel to `chains`: peel depth behind the cell (0/1/2)
    std::vector<char> chain_solid;                 // parallel to `chains`: chain came from hole-less material
    std::vector<char> chain_band;                  // parallel to `chains`: 1 = the global peeled-skin chain
    auto march_and_merge = [&](const ExPolygon &piece, bool single_bore, char thin) -> bool {
        std::vector<std::vector<ExPolygon>> frags;
        if (!march_bands(piece, spacing, !piece.holes.empty(), frags, point_count)) {
            if (cs_env_on("CS_DBG_CELLS"))
                fprintf(stderr, "[march-fail] L%zu piece holes=%zu area=%.1f pts=%zu\n",
                        layer_index, piece.holes.size(),
                        unscale<double>(unscale<double>(std::abs(piece.contour.area()))),
                        piece.contour.points.size());
            return false;
        }
        if (single_bore) {
            std::vector<ExPolygon> annular_bands;
            for (std::vector<ExPolygon> &frag : frags) {
                if (!frag.empty() && frag.front().holes.size() == 1)
                    for (ExPolygon &b : frag) annular_bands.push_back(std::move(b));
                else {
                    chains.push_back(std::move(frag));      // solid pinched-off piece: keep as-is
                    chain_thin.push_back(thin);
                    chain_band.push_back(char(false));
                    chain_solid.push_back(char(true));      // a pinch-off holds no bore: its degenerate
                                                            // slivers (wedge-cell medial hooks) are droppable
                }
            }
            if (annular_bands.size() > 1)
                std::sort(annular_bands.begin(), annular_bands.end(),
                          [](const ExPolygon &a, const ExPolygon &b) {
                              return std::abs(a.contour.area()) > std::abs(b.contour.area());
                          });
            if (!annular_bands.empty()) {
                chains.push_back(std::move(annular_bands));
                chain_thin.push_back(thin);
                chain_band.push_back(char(false));
                chain_solid.push_back(char(false));         // the real bore wall
            }
        } else {
            for (std::vector<ExPolygon> &frag : frags) {
                chains.push_back(std::move(frag));
                chain_thin.push_back(thin);
                chain_band.push_back(char(false));
                chain_solid.push_back(char(piece.holes.empty()));
            }
        }
        return true;
    };

    ExPolygons current = offset_ex(ExPolygons{ region }, -float(width) / 2.f);
    std::vector<CutLine> region_cuts;              // all modern cuts (concave + bore bisectors), for the interface weave
    std::vector<Polygons> chain_obs;               // parallel to `chains` (may run short): swoop keep-outs
    for (const ExPolygon &seed : current) {
        // BORE-OBSTACLE ("swoop") mode: small bores are not partition sources - the lattice is
        // generated as if the puck were SOLID (one global ring field anchored at the outer skin
        // only, uniform radial buckets, unbroken hoop continuity) and the kept rings are later
        // rerouted around each bore's keep-out outline (the swoop, applied in the woven block).
        // Each bore keeps a dedicated shell annulus as its own mini-chain for its visible skin.
        // This exists because partitioning around SMALL bores reorganizes the whole fill for a
        // local feature (the die-5 arc: every cell scheme traded voids for density); rerouting
        // treats them as the local perturbation they are.
        // MIXED MODE: holes are classified PER HOLE. Every swoopable hole (small, clear of the
        // skin and of other holes) becomes a keep-out + shell annulus; any hole that fails -
        // oversized (the counterbored disc's taper), skin-adjacent, crowded - STAYS a hole of
        // the EFFECTIVE seed, and the cell/classic partition below runs on that reduced seed.
        // The cells' rings then swoop around the small bores exactly like the pure-puck path
        // (their chains carry the keep-outs). Previously one oversized bore demoted the WHOLE
        // part to cells (the cb_disc printed near-solid above the taper's gate crossing).
        // stride > 1: swoop is a SPARSE-lattice mode. force_solid layers arrive with
        // woven.enabled=true but stride 1 (the brick-stagger carrier) - they must use the
        // classic solid multibore fill (the swoop reroute lives inside the K>1 gate and
        // can never run for them; taking this branch anyway collapsed the solid cap/floor
        // layers of the counterbored-disc gauntlet to a fifth of their material).
        ExPolygon eff = seed;                      // effective partition seed (mixed mode strips swoopable holes)
        Polygons  seed_keepouts;
        std::vector<ExPolygon> seed_annuli;
        const size_t seed_ci0 = chains.size();
        auto finish_mixed = [&]() {
            if (seed_keepouts.empty())
                return;
            chain_obs.resize(chains.size());
            for (size_t ci2 = seed_ci0; ci2 < chains.size(); ++ci2)
                chain_obs[ci2] = seed_keepouts;
            for (const ExPolygon &ann : seed_annuli) {
                chains.push_back({ ann });
                chain_thin.push_back(char(0));
                chain_band.push_back(char(true));  // essential bore skin: never sliver-filtered
                chain_solid.push_back(char(false));
            }
        };
        if (woven.enabled && woven.stride > 1 && woven.bore_obstacles && !seed.holes.empty()) {
            const size_t slb  = woven.shell_loops;
            Polygons keepouts;
            std::vector<ExPolygon> annuli;
            Polygons big_holes;
            for (const Polygon &h : seed.holes) {
                bool ok = true;
                Polygon hc = h;
                hc.reverse();                                     // holes are CW; offsets want CCW
                // Wall = slb BEADS: the hole curve is already the innermost bead centreline
                // (the seed is centreline-offset), so the outer wall bead sits (slb-1) beads
                // out. The keep-out (= the rank-0 fan arc centreline) sits ONE bead beyond
                // the outer wall bead - fan lanes weld directly onto the wall.
                const size_t wall_off = std::max<size_t>(1, slb - 1);
                Polygons shell_out = offset(hc, float(double(wall_off) * double(spacing)));
                Polygons keep      = offset(hc, float(double(wall_off + 1) * double(spacing)));
                if (shell_out.size() != 1 || keep.size() != 1)
                    ok = false;
                if (ok) {
                    const BoundingBox kb = keep.front().bounding_box();
                    const double extent = std::max(double(kb.max.x() - kb.min.x()), double(kb.max.y() - kb.min.y()));
                    if (extent > 8.0 * double(woven.stride) * double(spacing))
                        ok = false;
                }
                if (ok) {
                    // Skin clearance: the keep-out plus the full fan allowance must stay inside
                    // the contour - a bore that close to the skin breaks the wall.
                    Polygons margin = offset(keep.front(), float(6.0 * double(spacing)));
                    if (margin.size() != 1)
                        ok = false;
                    else
                        for (const Point &mp : margin.front().points)
                            if (!seed.contour.contains(mp)) { ok = false; break; }
                }
                if (ok) {
                    // Hole clearance is TIGHTER (keep-out + one bead): the only real hazard is
                    // this bore's keep-out reaching another hole's opening. The full fan margin
                    // here was over-strict - it demoted the cb_disc's slant bore to the cell
                    // set merely for sharing the upper layers with the fat taper, so mixed
                    // mode never engaged there.
                    Polygons hm = offset(keep.front(), float(1.0 * double(spacing)));
                    if (hm.size() != 1)
                        ok = false;
                    else
                        for (const Point &mp : hm.front().points) {
                            for (const Polygon &h2 : seed.holes)
                                if (&h2 != &h && h2.contains(mp)) { ok = false; break; }
                            if (!ok) break;
                        }
                }
                if (ok) {
                    ExPolygon ann;
                    ann.contour = shell_out.front();
                    ann.holes.push_back(h);
                    annuli.push_back(std::move(ann));
                    keepouts.push_back(keep.front());
                } else {
                    big_holes.push_back(h);
                }
            }
            const bool all_swoopable = big_holes.empty() && !keepouts.empty() &&
                                       deep_notches(seed.contour, std::max(8.0 * double(spacing), scale_(2.5))).empty();
            if (all_swoopable) {
                // The PUCK chain goes FIRST: the tour starts (and, closed, ends) on chain 0,
                // and the whole seam recipe - join placement, coast flow trim, precession -
                // belongs on the puck's OUTER SKIN ring. With an annulus first, the layer
                // join landed on a bore wall: the coast trimmed a near-full half-flow lap
                // off it every layer and the join scar precessed around the bore interior.
                const size_t puck_ci = chains.size();
                if (!march_and_merge(ExPolygon(seed.contour), false, char(0)))
                    return {};
                chain_obs.resize(chains.size());
                for (size_t ci2 = puck_ci; ci2 < chains.size(); ++ci2)
                    chain_obs[ci2] = keepouts;
                // Bore skins: one PRE-BUILT single-band chain each. Never march these - a
                // 2-bead annulus pinches immediately and shatters into dozens of sliver
                // fragments (die5_swoop: ~100 chains/layer), each needing two connector
                // hops that lapped the bore wall at half flow every layer (the mounds).
                // The assembly loop's medial gap fill still covers any mid-wall gap.
                for (const ExPolygon &ann : annuli) {
                    chains.push_back({ ann });
                    chain_thin.push_back(char(0));
                    chain_band.push_back(char(true));             // essential skin: never sliver-filtered
                    chain_solid.push_back(char(false));
                }
                continue;
            }
            if (!keepouts.empty()) {
                // MIXED: strip the swoopable holes from the partition seed. The cell/classic
                // paths below run on the REDUCED seed (oversized holes only - or none, for a
                // concave outline whose bores all qualify) and finish_mixed() attaches the
                // keep-outs to every chain they produce plus the bore-skin annuli.
                eff.holes     = big_holes;
                seed_keepouts = std::move(keepouts);
                seed_annuli   = std::move(annuli);
            }
            // no swoopable holes at all: fall through with the seed unchanged
        }
        if (woven.cell_split) {
            // TRIAL pass: would ANY modern cut split this seed - a bore-bisector (multi-hole
            // section) or a concave notch? The trial's geometry is discarded; it only decides
            // whether to peel. A seed with neither marches exactly like the classic path.
            const double cell_target = woven.cell_size > 0.1
                                     ? scale_(woven.cell_size)
                                     : (woven.stride > 1 ? 4.0 * double(woven.stride) : 16.0) * double(spacing);
            ExPolygons tb; std::vector<CutLine> tlines;
            if (multibore && eff.holes.size() > 1)
                partition_cells(eff, tb, spacing, cell_target, &tlines);
            else
                tb.push_back(eff);
            ExPolygons tcells; std::vector<CutLine> tcuts;
            for (const ExPolygon &cell : tb)
                partition_concave(cell, spacing, tcells, 0, &tcuts);
            if (!tlines.empty() || !tcuts.empty()) {
                // PEEL THE SKIN. Splitting the seed directly hands each cell its own stretch of
                // the object's exterior, so the visible perimeter is stitched from per-cell arcs
                // meeting at the cut mouths. The one-continuous-bead premise only forces a single
                // skin<->interior handoff per layer, so instead: emit the seed's outer contour as
                // its own SINGLE-LOOP chain FIRST - the whole exterior prints as one contiguous
                // bead at the start (end, on reversed layers) of the tour, exactly like an unsplit
                // part, seam recipe and all - and partition only the INTERIOR one bead behind it.
                // Cut mouths then terminate against the skin's inner face instead of interrupting
                // it. Cells behind the peel run one shell ring thinner (the peel ring is wall #1),
                // which also thins the cut interfaces from 2+2 to 1+1 walls. Bore bisectors go
                // through modernize_bore_cuts: recessed, extent-measured, and registered for the
                // interface weave exactly like concave cuts.
                // Peel the FULL wall stack when possible: ext perimeter + wall #2 as ONE global
                // annular chain. Interior cells then never back the skin at all, so they run with
                // ZERO outer shells (thin_shell cells: sl_outer = 0) - their outermost structure
                // is the bucket-0 woven ring itself, and cell interfaces become lattice kissing
                // lattice with no solid walls between. Falls back to the single-ring peel (cells
                // one bead behind, one shell) when the second ring vanishes (very thin parts).
                Polygons wall2 = offset(eff.contour, -float(spacing));
                const bool deep_peel = wall2.size() == 1 && wall2.front().points.size() >= 3;
                const ExPolygons inner = intersection_ex(
                    eff, offset(eff.contour, deep_peel ? -2.f * float(spacing) : -float(spacing)));
                if (!inner.empty()) {
                    if (deep_peel) {
                        ExPolygon band;
                        band.contour = eff.contour;
                        band.holes.push_back(std::move(wall2.front()));
                        band.holes.back().reverse();          // ExPolygon holes are CW
                        chains.push_back({ std::move(band) });
                    } else {
                        chains.push_back({ ExPolygon(eff.contour) });
                    }
                    chain_thin.push_back(char(false));
                    chain_band.push_back(char(true));     // the peel skin is never sliver-filtered
                    chain_solid.push_back(char(false));
                    ExPolygons cells2;
                    for (const ExPolygon &ip : inner) {
                        ExPolygons bc; std::vector<CutLine> blines;
                        if (multibore && ip.holes.size() > 1)
                            partition_cells(ip, bc, spacing, cell_target, &blines);
                        else
                            bc.push_back(ip);
                        const size_t bc_pre = bc.size(), bl_n = blines.size();
                        modernize_bore_cuts(bc, blines, spacing, region_cuts);
                        if (cs_env_on("CS_DBG_CELLS"))
                            fprintf(stderr, "[cells] L%zu ip holes=%zu bc pre=%zu lines=%zu post=%zu cuts=%zu\n",
                                    layer_index, ip.holes.size(), bc_pre, bl_n, bc.size(), region_cuts.size());
                        for (const ExPolygon &cell : bc)
                            partition_concave(cell, spacing, cells2, 0, &region_cuts);
                    }
                    strip_sliver_holes(cells2, spacing);
                    if (cs_env_on("CS_DBG_CELLS"))
                        fprintf(stderr, "[cells] L%zu cells2=%zu\n", layer_index, cells2.size());
                    if (const char *dump = std::getenv("CS_DUMP_CELLS"); dump != nullptr && dump[0] != '\0') {
                        static std::mutex cell_dump_mutex;
                        std::lock_guard<std::mutex> lock(cell_dump_mutex);
                        if (FILE *f = std::fopen(dump, "a")) {
                            for (const ExPolygon &cell : cells2) {
                                fprintf(f, "L %zu H %zu C", layer_index, cell.holes.size());
                                for (const Point &p : cell.contour.points)
                                    fprintf(f, " %lld %lld", (long long) p.x(), (long long) p.y());
                                fprintf(f, "\n");
                                for (const Polygon &h : cell.holes) {
                                    fprintf(f, "L %zu HOLE C", layer_index);
                                    for (const Point &p : h.points)
                                        fprintf(f, " %lld %lld", (long long) p.x(), (long long) p.y());
                                    fprintf(f, "\n");
                                }
                            }
                            std::fclose(f);
                        }
                    }
                    for (const ExPolygon &cell : cells2)
                        if (!march_and_merge(cell, cell.holes.size() == 1, char(deep_peel ? 2 : 1)))
                            return {};
                    if (cs_env_on("CS_DBG_CELLS"))
                        fprintf(stderr, "[cells] L%zu chains=%zu\n", layer_index, chains.size());
                    finish_mixed();
                    continue;                       // seed fully handled by the peel path
                }
                // Interior vanished (skin-only part): march the trial cells unpeeled.
                region_cuts.insert(region_cuts.end(), tcuts.begin(), tcuts.end());
                strip_sliver_holes(tcells, spacing);
                for (const ExPolygon &cell : tcells)
                    if (!march_and_merge(cell, cell.holes.size() == 1, char(0)))
                        return {};
                finish_mixed();
                continue;
            }
            for (const ExPolygon &cell : tcells)    // no cuts: tcells is the plain seed
                if (!march_and_merge(cell, cell.holes.size() == 1, char(0)))
                    return {};
            finish_mixed();
            continue;
        }
        // Classic path (cell_split off): the multi-bore Voronoi partition only, unchanged.
        ExPolygons cells;
        if (multibore && eff.holes.size() > 1)
            partition_bores(eff, cells);
        else
            cells.push_back(eff);
        for (const ExPolygon &cell : cells)
            if (!march_and_merge(cell, cell.holes.size() == 1, char(0)))
                return {};
        finish_mixed();
    }

    std::vector<ContourChain> out;

    if (chains.empty()) {
        if (cs_env_on("CS_DBG_CELLS"))
            fprintf(stderr, "[medial-fallback] L%zu holes=%zu\n", layer_index, region.holes.size());
        // Sub-bead-width region: fall back to a single medial-axis centreline pass so
        // the wall stays connected. Adaptive bead width is deliberately future work.
        Polygons loops;
        Polylines centerlines = region.medial_axis(double(width) / 5., double(2 * width));
        for (const Polyline &pl : centerlines) {
            if (pl.points.size() < 3)
                continue;
            Polygon loop(pl.points);
            if (loop.points.size() > 1 && loop.points.front() == loop.points.back())
                loop.points.pop_back();
            if (loop.points.size() >= 3)
                loops.push_back(std::move(loop));
        }
        if (loops.empty())
            return {};
        out.push_back({ std::move(loops), false });
        return out;
    }

    for (size_t ci = 0; ci < chains.size(); ++ci) {
        const auto &chain = chains[ci];
        ContourChain cc;
        cc.peel_depth = ci < chain_thin.size() ? int(chain_thin[ci]) : 0;
        cc.from_solid = ci < chain_solid.size() && chain_solid[ci] != 0;
        cc.peel_band  = ci < chain_band.size() && chain_band[ci] != 0;
        if (ci < chain_obs.size())
            cc.obstacles = chain_obs[ci];
        cc.loops.reserve(chain.size() * 2 + 1);
        // A loop is admissible only if its CURVE stays clear of the REGION's holes - a
        // genuine ring (any offset of the wall) never enters a bore; a PINCHED band's
        // residual void pocket or a stray medial strip does (sampled per edge; corner-
        // sparse loops have multi-mm straight edges whose vertices all sit in material).
        auto clear_of_holes = [&](const Polygon &loop) {
            if (region.holes.empty())
                return true;
            const Points &lp = loop.points;
            for (size_t k2 = 0; k2 < lp.size(); ++k2) {
                const Point &a2 = lp[k2];
                const Point &b2 = lp[(k2 + 1) % lp.size()];
                const int ns2 = std::max(1, int((b2 - a2).cast<double>().norm() / scale_(0.3)));
                for (int s3 = 0; s3 <= ns2; ++s3) {
                    const Point q2 = a2 + (b2 - a2) * s3 / ns2;
                    for (const Polygon &h : region.holes)
                        if (h.contains(q2))
                            return false;
                }
            }
            return true;
        };
        // Outer contours: outer edge -> medial axis.
        for (const ExPolygon &band : chain) {
            cc.footprint += std::abs(band.area());
            cc.loops.push_back(band.contour);
        }
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
                    // The strip must also stay clear of the REGION's holes (see
                    // clear_of_holes above): a merged PINCHED cell can end on a stranded
                    // middle fragment whose offset strip cuts across the real bore.
                    const ExPolygon &strip = strips.front();
                    if (strip.contour.points.size() >= 3 && strip.contour.length() > 4 * width &&
                        clear_of_holes(strip.contour))
                        cc.loops.push_back(strip.contour);
                    const Polygon &hole = strip.holes.front();
                    if (hole.points.size() >= 3 && hole.length() > 4 * width &&
                        clear_of_holes(hole))
                        cc.loops.push_back(hole);
                }
            }
        }
        // Inner (hole) contours: medial axis -> bore (innermost first), for an annular
        // chain only. A piece pinched off solid contributes none. Gate each on
        // clear_of_holes: a genuine bore ring (any grown offset of the bore) never enters
        // the bore, while a PINCHED band's residual VOID POCKET - a scrap of grown-bore
        // space sealed off by the pinch - also lands in this run sitting ON TOP of the raw
        // bore, and its boundary then chords straight across it (the pip chords every
        // die-5 print shared: two pocket loops per corner cell at the solid layers).
        // Convergence cores on hole-less regions pass trivially (the from_solid
        // degenerate filter owns those).
        for (auto it = chain.rbegin(); it != chain.rend(); ++it)
            if (!it->holes.empty()) {
                const Polygon &hl = it->holes.front();
                if (!clear_of_holes(hl))
                    continue;
                cc.loops.push_back(hl);
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
            // Shell counts split by END. Behind a SHALLOW peel (depth 1) the global perimeter
            // chain is wall #1, so the cell keeps one outer shell fewer. Behind a DEEP peel
            // (depth 2: the full wall stack is global) the cell runs SHELL-FREE outside - its
            // outermost structure is the bucket-0 woven ring, so cell interfaces are lattice
            // kissing lattice with no solid wall between. Bore-side shells are always kept:
            // the bore surface is visible skin.
            const size_t sl_outer = cc.peel_depth >= 2 ? 0
                                  : cc.peel_depth == 1 ? std::max<size_t>(1, woven.shell_loops - 1)
                                                       : woven.shell_loops;
            const size_t sl_bore  = cc.peel_depth > 0 ? std::max<size_t>(1, woven.shell_loops - 1)
                                                      : woven.shell_loops;
            const size_t K  = woven.stride;
            const size_t n  = cc.loops.size();
            const bool   annular = cc.has_bore;
            // is_infill on ORIGINAL indices: buried interior beyond the shell(s). Annular keeps
            // shells at both ends; solid keeps the outer shells and the buried core (minus the
            // terminus).
            auto is_infill = [&](size_t i) -> bool {
                if (i < sl_outer) return false;
                return annular ? (i + sl_bore < n) : (i + 1 < n);
            };
            // Sparsify only when the cell can actually HOST a lattice: at stride K' the kiss
            // geometry needs at least ~two kept rings plus their voids (~2K' raw rings beyond
            // the shells). ADAPTIVE PER-CELL STRIDE (die-5 v7 rim mosaic): a cell too small for
            // the configured K used to print SOLID outright; instead step DOWN to the largest
            // stride Keff in [2, K] the cell can host, so compact cells still get a lattice -
            // just a tighter one, density landing between the configured value and 100% - and
            // only truly tiny cells (under ~5 rings beyond the shells) print SOLID: a dense
            // rosette beats one stranded ring in empty space (die-5 v4 print: sparsified caps
            // came out as bare corners crossed by connector strands). The cell-LOCAL geometry
            // (bucket pitch, amplitude, dyadic cell-pitch default) scales with Keff; the
            // CROSS-cell machinery keeps the configured K (the bubble Z-carrier period and the
            // interface field's wavelength/reach must agree on both sides of every cut).
            const size_t shells = sl_outer + (annular ? sl_bore : 0);
            size_t Keff = K;
            if (woven.adaptive_stride)
                while (Keff > 2 && n <= shells + 2 * Keff)
                    --Keff;
            const bool thick = n > shells + 2 * Keff;
            // The interface behaviour applies to EVERY cell with modern cuts nearby - including
            // small solid ones - so the carrier and field run before the sparsify gate.
            const bool cs_bubble = woven.bubble || cs_env_on("CS_BUBBLE");
            double zcarrier = 1.0;
            if (cs_bubble) {
                double P = std::ceil(M_PI * double(K + 1) / 0.75);
                if (woven.bubble_pz >= 4.0)
                    P = woven.bubble_pz;
                if (woven.bubble_h > 0.0 && double(flow.height()) > 1e-3)
                    P = std::max(6.0, 2.0 * woven.bubble_h / double(flow.height()));
                zcarrier = std::cos(2.0 * M_PI * double(layer_index) / P);
            }
            if (woven.weave > 0.0 && !region_cuts.empty()) {
                // INTERFACE WEAVE (see the block comment at the original site, now hoisted so
                // small SOLID cells' boundaries also kiss their neighbours across the cuts).
                const double A_int  = 0.5 * double(spacing);
                const double lambda = 3.0 * double(K) * double(spacing);
                const double reach  = double(K) * double(spacing);
                for (size_t j = 0; j < cc.loops.size(); ++j) {
                    Points &pts = cc.loops[j].points;
                    if (pts.size() < 3) continue;
                    {   // densify: shell rings of a faceted cell are corner-points only
                        const Polygon poly(pts);
                        const size_t nd = std::clamp<size_t>(size_t(poly.length() / (0.5 * double(spacing))), 48, 1024);
                        std::vector<Point> rs = resample_closed(poly, nd);
                        if (rs.size() >= 8) pts.assign(rs.begin(), rs.end());
                    }
                    for (Point &q : pts) {
                        const Vec2d p = q.cast<double>();
                        Vec2d delta(0., 0.);
                        for (const CutLine &cl : region_cuts) {
                            const Vec2d  rel = p - cl.o;
                            const double dd  = rel.dot(cl.n);
                            const double ss  = rel.dot(cl.t);
                            const double ad  = std::abs(dd);
                            if (ad > reach) continue;
                            const double over = std::max({ 0.0, cl.smin + reach - ss, ss - (cl.smax - reach) });
                            const double w    = (1.0 - ad / reach) * std::clamp(1.0 - over / reach, 0.0, 1.0);
                            if (w <= 0.0) continue;
                            const double sigma = dd >= 0.0 ? 1.0 : -1.0;
                            delta += cl.n * (sigma * A_int * zcarrier * std::sin(2.0 * M_PI * ss / lambda) * w);
                        }
                        if (delta.squaredNorm() > 1.0)
                            q = Point(coord_t(std::llround(p.x() + delta.x())),
                                      coord_t(std::llround(p.y() + delta.y())));
                    }
                }
            }
            if (K > 1 && thick) {
                Polygons kept; std::vector<char> kept_infill;
                kept.reserve(n); kept_infill.reserve(n);
                // RADIAL-BUCKET sparsification. K-striding the loop LIST assumes list position ==
                // uniform radial pitch, which breaks on a non-uniform annulus: the list is
                // outer-run + medial + bore-run, so kept rings land at irregular radial gaps, the
                // K*s kiss geometry fails, and the "lattice" degenerates into freestanding rings
                // bridged by sagging connector strands (the die-5 wedge voids). Instead bucket
                // every infill ring by its MEASURED depth from the chain's outer face and keep one
                // ring per K*s bucket (nearest to the bucket centre wins): uniform radial pitch by
                // construction, list structure irrelevant. The bucket index also supplies the
                // weave parity below, so anti-phase is radially correct across the medial too.
                const Polygon &outer_ref = cc.loops.front();
                auto ring_depth = [&](const Polygon &ring) {
                    double sum = 0.0; size_t cnt = 0;
                    const size_t step = 1 + ring.points.size() / 48;
                    for (size_t k = 0; k < ring.points.size(); k += step) {
                        const Vec2d p = ring.points[k].cast<double>();
                        double best = std::numeric_limits<double>::max();
                        const Points &op = outer_ref.points;
                        for (size_t e2 = 0; e2 < op.size(); ++e2) {
                            const Vec2d a  = op[e2].cast<double>();
                            const Vec2d ab = op[(e2 + 1) % op.size()].cast<double>() - a;
                            const double l2 = std::max(1.0, ab.squaredNorm());
                            const double t  = std::clamp((p - a).dot(ab) / l2, 0.0, 1.0);
                            best = std::min(best, (p - (a + t * ab)).squaredNorm());
                        }
                        sum += std::sqrt(best); ++cnt;
                    }
                    return cnt > 0 ? sum / double(cnt) : 0.0;
                };
                std::vector<double> depth_of(n, -1.0);
                double d0 = std::numeric_limits<double>::max();      // shallowest infill ring anchors bucket 0
                for (size_t i = 0; i < n; ++i)
                    if (is_infill(i)) {
                        depth_of[i] = ring_depth(cc.loops[i]);
                        d0 = std::min(d0, depth_of[i]);
                    }
                std::vector<long>   bucket_of(n, -1);
                std::map<long, size_t> pick;                 // bucket -> ring index nearest its centre
                std::map<long, double> pick_err;
                for (size_t i = 0; i < n; ++i) {
                    if (depth_of[i] < 0.0) continue;
                    const double rel = depth_of[i] - d0;
                    const long   b   = std::lround(rel / (double(Keff) * double(spacing)));
                    const double err = std::abs(rel - double(b) * double(Keff) * double(spacing));
                    bucket_of[i] = b;
                    auto it = pick.find(b);
                    if (it == pick.end() || err < pick_err[b]) { pick[b] = i; pick_err[b] = err; }
                }
                std::vector<int> kept_bucket;
                for (size_t i = 0; i < n; ++i) {
                    const bool inf = is_infill(i);
                    const bool keep = !inf || (bucket_of[i] >= 0 && pick[bucket_of[i]] == i);
                    if (keep) {
                        kept_infill.push_back(inf ? char(1) : char(0));
                        kept_bucket.push_back(inf ? int(bucket_of[i]) : -1);
                        kept.push_back(std::move(cc.loops[i]));
                    }
                }
                // Centre-convergence nub (task #9): the inward march converges at a genuinely
                // SOLID piece's core into DEGENERATE rings - a tiny loop, or on an elongated
                // cell a flat medial sliver traversed down-and-back at effectively double flow
                // (the elongated case also routes through the medial-strip/bore machinery, so
                // the degenerate loops can sit MID-LIST and even flag the chain annular) - and
                // the weave scribbles them into a knot of excess plastic. For from_solid chains
                // (the marched piece had no real holes - any apparent bore is the convergence
                // core) drop EVERY degenerate loop beyond the outer shells: effective bead width
                // 2*area/perimeter under one bead, or perimeter under ~8 beads. Real bores are
                // untouchable: their chains are not from_solid.
                if (cc.from_solid) {
                    Polygons          fkept;
                    std::vector<char> finf;
                    std::vector<int>  fbuck;
                    fkept.reserve(kept.size()); finf.reserve(kept.size()); fbuck.reserve(kept.size());
                    for (size_t i2 = 0; i2 < kept.size(); ++i2) {
                        bool degen = false;
                        if (i2 >= sl_outer) {
                            const double per = kept[i2].length();
                            const double aw  = per > 1.0 ? 2.0 * std::abs(kept[i2].area()) / per : 0.0;
                            degen = per < 8.0 * double(spacing) || aw < 1.0 * double(spacing);
                        }
                        if (!degen) {
                            fkept.push_back(std::move(kept[i2]));
                            finf.push_back(kept_infill[i2]);
                            fbuck.push_back(kept_bucket[i2]);
                        }
                    }
                    kept = std::move(fkept);
                    kept_infill = std::move(finf);
                    kept_bucket = std::move(fbuck);
                }
                cc.loops = std::move(kept);
                cc.woven_infill = true;   // spiralize_rings should tag [shell | infill | shell] accordingly

                // THE SWOOP (bore-obstacle mode): reroute every kept ring that crosses a bore
                // keep-out around it, like streamlines around a cylinder. Region booleans keep
                // the invariants the tour needs (closed, simple, nested): a ring whose enclosed
                // region contains the bore bulges around its FAR flank (union with the keep-out);
                // one passing on the near side dips around the NEAR flank (difference). Rings
                // crossing the same flank nest at one-bead pitch - keep-out grown by each ring's
                // same-side closeness rank - forming a dense passing-lane fan beside the bore
                // (free bore-margin reinforcement); the weave amp cap below flattens the fan.
                // Rings swallowed whole by a keep-out (the march's convergence core landing
                // inside a central bore) vanish, ending the chain hugging that bore.
                std::vector<char> kept_swooped(cc.loops.size(), 0);   // rings rerouted around a bore
                Polygons amp_keepouts;    // per-bore FAN outlines (base keep-out grown by the deepest
                                          // detour rank + the 0.85-bead flatten clearance) - the
                                          // amplitude cap must flatten crests BESIDE the outermost
                                          // fan arc, not at the base keep-out: crests that reach the
                                          // base land ON the solid fan lanes and mound excess plastic
                                          // in a glossy ring around every bore (die5_swoop print).
                if (!cc.obstacles.empty()) {
                    auto dist_to_poly = [](const Vec2d &p, const Polygon &poly) {
                        double best = std::numeric_limits<double>::max();
                        const Points &pp = poly.points;
                        for (size_t e2 = 0; e2 < pp.size(); ++e2) {
                            const Vec2d a  = pp[e2].cast<double>();
                            const Vec2d ab = pp[(e2 + 1) % pp.size()].cast<double>() - a;
                            const double l2 = std::max(1.0, ab.squaredNorm());
                            const double t  = std::clamp((p - a).dot(ab) / l2, 0.0, 1.0);
                            best = std::min(best, (p - (a + t * ab)).norm());
                        }
                        return best;
                    };
                    std::vector<char> dead(cc.loops.size(), 0);
                    for (const Polygon &keep : cc.obstacles) {
                        const Point kc = keep.centroid();
                        struct Cross { size_t j; double pen; bool far_side; };
                        std::vector<Cross> crossers;
                        for (size_t j = 0; j < cc.loops.size(); ++j) {
                            if (dead[j] || j < sl_outer) continue;   // every non-shell ring, core included
                            const Polygon &ring = cc.loops[j];
                            if (ring.points.size() < 3) continue;
                            const bool inside = ring.contains(kc);
                            // Crossing = the ring curve enters the keep-out (or the keep-out
                            // swallowed the whole ring - deep convergence rings near a central
                            // bore). PENETRATION (deepest incursion of the curve into the
                            // keep-out) is the nesting key: on BOTH flanks the deepest-cutting
                            // ring hugs the base keep-out and shallower ones nest outward, which
                            // provably preserves ring nesting (union side: penetration falls
                            // with ring size; difference side: it falls as rings shrink away
                            // from the bore). Distance-to-boundary is ~0 for every crosser and
                            // ranks them arbitrarily.
                            double pen = 0.0;
                            for (const Point &q : ring.points)
                                if (keep.contains(q))
                                    pen = std::max(pen, dist_to_poly(q.cast<double>(), keep));
                            if (pen <= 0.0)
                                continue;                            // curve stays clear of this bore
                            crossers.push_back({ j, pen, inside });
                        }
                        size_t bore_max_rank = 0;
                        bool   bore_has_fan  = false;
                        for (int side = 0; side < 2; ++side) {
                            std::vector<Cross> ss;
                            for (const Cross &c : crossers)
                                if (c.far_side == (side == 1)) ss.push_back(c);
                            std::sort(ss.begin(), ss.end(), [](const Cross &a, const Cross &b) { return a.pen > b.pen; });
                            for (size_t rank = 0; rank < ss.size(); ++rank) {
                                bore_has_fan  = true;
                                bore_max_rank = std::max(bore_max_rank, rank);
                                const size_t j = ss[rank].j;
                                Polygons grown = rank == 0 ? Polygons{ keep }
                                                           : offset(keep, float(double(rank) * double(spacing)));
                                if (grown.size() != 1) continue;
                                ExPolygons rr = ss[rank].far_side
                                              ? union_ex(Polygons{ cc.loops[j], grown.front() })
                                              : diff_ex(Polygons{ cc.loops[j] }, Polygons{ grown.front() });
                                // keep the largest surviving piece; a vanished ring is dropped
                                size_t bi = 0; double ba = -1.0;
                                for (size_t pi2 = 0; pi2 < rr.size(); ++pi2) {
                                    const double a2 = std::abs(rr[pi2].contour.area());
                                    if (a2 > ba) { ba = a2; bi = pi2; }
                                }
                                if (rr.empty() || ba < 16.0 * double(spacing) * double(spacing))
                                    dead[j] = 1;
                                else {
                                    cc.loops[j] = rr[bi].contour;
                                    kept_swooped[j] = 1;
                                }
                            }
                        }
                        // The crest-flatten boundary for this bore: beside the OUTERMOST fan
                        // arc, with the same 0.85-bead clearance the shells use. No fan: the
                        // base keep-out already carries the shell clearance.
                        const double fan_grow = bore_has_fan
                                              ? (double(bore_max_rank) + 0.85) * double(spacing)
                                              : 0.0;
                        Polygons gk = fan_grow > 0.0 ? offset(keep, float(fan_grow)) : Polygons{ keep };
                        if (gk.size() == 1)
                            amp_keepouts.push_back(std::move(gk.front()));
                        else
                            amp_keepouts.push_back(keep);
                        // Connector-ride lane: mid-corridor between the outermost fan arc and
                        // the crest-flatten boundary - reserved empty by the cap above, so a
                        // half-flow ride lands in open space beside the fan instead of mounding
                        // onto the bore wall crown.
                        if (bore_has_fan) {
                            Polygons rl = offset(keep, float((double(bore_max_rank) + 0.4) * double(spacing)));
                            if (rl.size() == 1)
                                cc.ride_rails.push_back(std::move(rl.front()));
                        }
                    }
                    Polygons          fkept;
                    std::vector<char> finf, fsw;
                    std::vector<int>  fbuck;
                    for (size_t j = 0; j < cc.loops.size(); ++j)
                        if (!dead[j]) {
                            fkept.push_back(std::move(cc.loops[j]));
                            finf.push_back(kept_infill[j]);
                            fbuck.push_back(kept_bucket[j]);
                            fsw.push_back(kept_swooped[j]);
                        }
                    cc.loops     = std::move(fkept);
                    kept_infill  = std::move(finf);
                    kept_bucket  = std::move(fbuck);
                    kept_swooped = std::move(fsw);
                }

                if (woven.weave > 0.0) {
                    // Perturbation centre = centroid of the outer boundary loop.
                    Vec2d centre(0., 0.);
                    const Polygon &outer = cc.loops.front();
                    if (!outer.points.empty()) {
                        for (const Point &p : outer.points) centre += p.cast<double>();
                        centre /= double(outer.points.size());
                    }
                    // Radial amplitude sized so that ANTI-PHASE neighbours just touch. Sparse rings
                    // sit K*spacing apart; giving consecutive infill rings a pi phase offset makes ring
                    // j bulge OUT where ring j+1 bulges IN, so their closest approach is
                    // K*spacing - 2*amp. Targeting a one-bead centre gap (amp = weave*(K-1)/2*spacing)
                    // brings the bead EDGES into contact at the lobe peaks - the rings knit into one
                    // connected in-plane lattice ("lobes along a slice") instead of separate wavy tubes,
                    // while the valleys stay open for the sparse voids. Capped below half the spacing so
                    // they never cross.
                    // CS_BUBBLE (prototype): drive the crossbonds to genuine FUSION instead of kisses.
                    // The kiss weave sizes amp so anti-phase neighbours just touch (a marginal ~0.1mm bond
                    // that does not stiffen the part); bubble mode overlaps them by ~one bead at the lobe
                    // peaks - real welds that tie the concentric rings into one network - leaving big open
                    // cells between. closest approach = K*spacing - 2*amp; one-bead overlap => amp=(K+1)/2.
                    const double amp = (cs_bubble ? 0.5 * (double(Keff) + 1.0)
                                                  : std::min(woven.weave * 0.5 * double(Keff - 1), 0.49 * double(Keff)))
                                       * double(spacing);
                    // Bubble Z-carrier: multiply the radial wave by a slow cosine of the LAYER index. With
                    // the lattice phase pinned (dphase=0 in bubble mode) the crest welds form vertical
                    // columns; the carrier makes those columns wax, vanish, and reappear SIGN-FLIPPED half
                    // a period later - crests become valleys - so the contact columns alternate between the
                    // two crest-angle sets. Each void cell is then pinched on four sides: left/right by this
                    // half-period's welds, top/bottom by the flipped welds above and below - a 3D checker-
                    // board of ellipsoidal bubbles (open-cell: throats stay open at the carrier's zero
                    // crossings, where the rings relax to plain concentric). Period from printability: the
                    // crest migrates radially amp*2*pi/P per layer, kept under ~3/4 bead so each bead still
                    // lands mostly on the one below => P = ceil(pi*(K+1)/0.75) (~21 layers at K=4, giving
                    // ~2.1mm bubbles vs the 1.68mm radial cell - near-round). CS_BUBBLE_PZ overrides P.
                    // (cs_bubble / zcarrier / the interface field are computed and applied at the
                    // top of the woven gate, before the sparsify decision, so small SOLID cells
                    // get the interface behaviour too.)
                    // The lobes must stay INSIDE the infill band: an outward bulge that crossed the
                    // solid outer shell would push nubs past the object surface (a rough wall); an
                    // inward one would invade the bore shell. The shell rings are un-perturbed, so
                    // clamp every infill point to the shell radius AT ITS OWN ANGLE - shape-aware, so
                    // it works for non-round cross-sections too (a mean radius would flatten corners).
                    auto build_profile = [&](const Polygon &ring) {
                        const int NB = 360;
                        // Densify the source ring first: a faceted shell ring is little more than its
                        // corner points, which leaves most angular bins EMPTY - and the fill-forward
                        // below then smears a corner's radius across whole faces, so the clamp lets
                        // crests poke past flat walls (found as ~0.1mm nubs outside a split cell's
                        // outer face; the cube's shell had the same smear all along).
                        Points dense = ring.points;
                        if (dense.size() < 3 * NB) {
                            std::vector<Point> rs = resample_closed(ring, size_t(3 * NB));
                            if (rs.size() >= 8) dense.assign(rs.begin(), rs.end());
                        }
                        std::vector<double> prof(NB, 0.0); std::vector<int> cnt(NB, 0);
                        for (const Point &q : dense) {
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
                    // Conservative clamp samples: the outward limit takes the MIN of the bin and its
                    // neighbours (and the bore limit the MAX), so bin quantisation near corners - where
                    // the radius changes fastest - can never let a crest overshoot a wall.
                    auto sample_min = [](const std::vector<double> &prof, double th) {
                        const int NB = int(prof.size());
                        int b = int((th + M_PI) / (2.0 * M_PI) * NB); b = ((b % NB) + NB) % NB;
                        return std::min({ prof[(b + NB - 1) % NB], prof[b], prof[(b + 1) % NB] });
                    };
                    auto sample_max = [](const std::vector<double> &prof, double th) {
                        const int NB = int(prof.size());
                        int b = int((th + M_PI) / (2.0 * M_PI) * NB); b = ((b % NB) + NB) % NB;
                        return std::max({ prof[(b + NB - 1) % NB], prof[b], prof[(b + 1) % NB] });
                    };
                    // Outer clamp reference: the innermost outer-shell ring - or, SHELL-FREE
                    // (sl_outer == 0), the bucket-0 ring's own current (field-waved, pre-lobe)
                    // geometry, i.e. the cell boundary path: outward crests then flatten exactly
                    // at the boundary, where the neighbour cell's lattice meets them a bead away.
                    const std::vector<double> outer_prof = build_profile(cc.loops[sl_outer > 0 ? sl_outer - 1 : 0]);
                    const std::vector<double> inner_prof = annular ? build_profile(cc.loops[cc.loops.size() - sl_bore]) // outermost bore-shell ring
                                                                   : std::vector<double>();
                    double outer_mean_r = 1.0;                      // mean radius of the outer shell (bubble morph reference)
                    if (cs_bubble && !outer_prof.empty()) {
                        double s = 0.0; for (double v : outer_prof) s += v;
                        outer_mean_r = std::max(1.0, s / double(outer_prof.size()));
                    }
                    // Per-ring anti-phase: each successive infill ring flips by pi so neighbours
                    // interlock. The per-LAYER advance (dphase) rotates the whole lattice so stacked
                    // layers weave. Perturb along ARC LENGTH and the local NORMAL - NOT by centroid
                    // angle/radius. A raw sin(lobes*theta) radial push only works on a circle: on a
                    // square/gear the centroid angle is non-uniform (corners subtend little angle but
                    // lots of perimeter), so the lobes bunch at the corners, and a radial push is not
                    // perpendicular to a slanted wall, so bumps shove diagonally - together that
                    // scrambles the lattice. Arc length spaces the lobes evenly along the perimeter and
                    // the local normal keeps each bump perpendicular to the wall, giving clean, radially
                    // aligned, interlocking lobes on round AND faceted cross-sections alike. On a circle
                    // arc-fraction is proportional to angle and the normal is radial, so this reduces to
                    // the old behaviour exactly.
                    const double layer_phase = double(layer_index) * woven.dphase;
                    // PASS 1 - densify + (bubble) morph every kept infill ring, recording its loop
                    // index, perimeter, capped amplitude, and the classic parity phase. The dyadic
                    // assignment below needs every ring's perimeter BEFORE any ring is perturbed.
                    struct WeaveRing { size_t loop; double perim; double amp; double lobes; double phase; };
                    std::vector<WeaveRing> wrings;
                    {
                        size_t infill_ring = 0;
                        for (size_t j = 0; j < cc.loops.size(); ++j) {
                            if (!kept_infill[j]) continue;   // weave only the thinned infill, not the solid shells
                            // Parity from the ring's RADIAL bucket (not its list position): anti-phase
                            // then interlocks correctly even across the medial junction of a
                            // non-uniform annulus, and is stable against list reshuffles.
                            const long parity = (j < kept_bucket.size() && kept_bucket[j] >= 0)
                                                    ? long(kept_bucket[j]) : long(infill_ring);
                            const double parity_phase = layer_phase + double(parity) * M_PI;
                            ++infill_ring;
                            Points &pts = cc.loops[j].points;
                            if (pts.size() < 3) continue;
                            // Densify first: march-offset rings of a POLYGON are sparse - a square is just its
                            // 4 corner points, and no lobe wave (nor a fusion bump) can be drawn on 4 points, so
                            // the perturbation below would be a near-no-op on any faceted cross-section. Resample
                            // to ~2 points per bead. Round rings arrive dense from the mesh, so this barely
                            // changes them (arc-fraction/normal already reduce to the old behaviour there).
                            // (Skipped when the interface field ran: it already densified every kept ring, and
                            // resampling a field-waved ring would smooth the wave back out.)
                            if (region_cuts.empty()) {
                                const Polygon poly(pts);
                                const size_t cap = cs_bubble ? 1024 : 512;   // dyadic outer rings carry more lobes
                                const size_t nd = std::clamp<size_t>(size_t(poly.length() / (0.5 * double(spacing))), 48, cap);
                                std::vector<Point> rs = resample_closed(poly, nd);
                                if (rs.size() >= 8) pts.assign(rs.begin(), rs.end());
                            }
                            double amp_ring = amp;
                            // (Self-collision is handled PER-POINT in the perturb pass: a whole-ring
                            // width cap muted elongated rings entirely - die-5 v4's corner-cell rings
                            // barely reached their neighbours - when only their narrow ENDS need it.)
                            if (cs_bubble) {
                                double rmean = 0.0;
                                for (const Point &q : pts) rmean += (q.cast<double>() - centre).norm();
                                rmean /= double(pts.size());
                                // Morph the ring toward a CIRCLE of its own mean radius, more the deeper it
                                // sits (circle at the core, full contour shape at the wall) - the polar-Fourier
                                // "circle-in / shape-out" base that gives smooth cells and natural contact zones.
                                // SWOOPED rings are exempt: circularizing one drags its bore detour straight
                                // back through the bore, and the keep-out projection then snaps consecutive
                                // points to OPPOSITE flanks - a full-width chord across the hole.
                                const bool swooped_ring = j < kept_swooped.size() && kept_swooped[j] != 0;
                                const double beta = swooped_ring ? 0.0 : std::clamp(1.0 - rmean / outer_mean_r, 0.0, 1.0);
                                if (beta > 0.01)
                                    for (Point &q : pts) {
                                        const Vec2d d = q.cast<double>() - centre;
                                        const double r = d.norm(); if (r < 1e-6) continue;
                                        const double rr = (1.0 - beta) * r + beta * rmean;
                                        q = Point(coord_t(std::llround(centre.x() + rr * d.x() / r)),
                                                  coord_t(std::llround(centre.y() + rr * d.y() / r)));
                                    }
                                amp_ring = std::min(amp_ring, 0.5 * rmean);   // keep a deep ring's bump from reaching the centre
                            }
                            wrings.push_back({j, Polygon(pts).length(), amp_ring, woven.lobes, parity_phase});
                        }
                    }
                    // CS_BUBBLE_DYADIC (prototype) - per-ring lobe counts in a power-of-2 hierarchy.
                    // A single global count stretches cells with perimeter (a 30mm cube's outer cells are
                    // ~20mm wide while the core's are ~2mm). Instead give each ring the lobes*2^m nearest
                    // perimeter / CS_BUBBLE_W (target cell pitch, mm; default 3*K*spacing), monotone
                    // non-increasing inward. Rings sharing a count interlock anti-phase as usual; at a
                    // q = 2^dm boundary the phase rule
                    //     phi_outer = q*phi_inner - (q+1)*pi/2      (q=1 reduces to the classic +pi flip)
                    // lands EVERY inner crest on an outer crest, while the outer ring's remaining crests
                    // stay free - smaller cells that bond only within their own band, the accepted trade
                    // for near-constant cell size. On carrier-negative half-cycles crests swap with
                    // valleys, which would un-bond the boundaries; flipping every DOUBLED (m>=1) ring by
                    // pi restores them (the exact correction is (2^m - 1)*pi ≡ pi mod 2pi for all m>=1),
                    // and the flip lands at the carrier zero crossing where the amplitude vanishes, so
                    // the geometry stays continuous in Z.
                    if (cs_bubble && woven.dyadic && !wrings.empty()) {
                        const double w_target = woven.bubble_w > 0.1 ? scale_(woven.bubble_w)
                                                                     : 3.0 * double(Keff) * double(spacing);
                        std::vector<int> ms(wrings.size(), 0);
                        int prev_m = 5;
                        for (size_t t = 0; t < wrings.size(); ++t) {          // outer -> inner
                            const double want = wrings[t].perim / std::max(1.0, w_target);
                            int best_m = 0; double best_d = std::abs(woven.lobes - want);
                            for (int mm = 1; mm <= 5; ++mm) {
                                const double d = std::abs(woven.lobes * double(1 << mm) - want);
                                if (d < best_d) { best_d = d; best_m = mm; }
                            }
                            best_m = std::min(best_m, prev_m);                 // perimeters shrink inward: never double inward
                            prev_m = best_m; ms[t] = best_m;
                            wrings[t].lobes = woven.lobes * double(1 << best_m);
                        }
                        double phi = 0.0;                                      // march innermost -> outward
                        wrings.back().phase = phi;
                        for (size_t t = wrings.size() - 1; t-- > 0;) {
                            const double q = wrings[t].lobes / wrings[t + 1].lobes;
                            phi = q * phi - (q + 1.0) * M_PI / 2.0;
                            wrings[t].phase = phi;
                        }
                        if (zcarrier < 0.0)
                            for (size_t t = 0; t < wrings.size(); ++t)
                                if (ms[t] > 0) wrings[t].phase += M_PI;
                    }
                    // PASS 2 - perturb along arc length / local normal. The weld flat lands BESIDE the
                    // shell, never on top: clamping a crest exactly to the shell ring's own path (the old
                    // behaviour) double-deposits a bead along the wall, and the displaced melt telegraphs
                    // through the 2-bead shell as alternating surface bumps tracking the internal weld
                    // pads. Keep ~0.85 bead of centreline clearance instead: the flattened crest bonds
                    // sideways with ~15% edge overlap - the same weld every print uses between infill and
                    // walls - and the wall stays flat.
                    // Shell-free cells (sl_outer == 0) clamp AT the boundary (gap 0): the flattened
                    // crest IS the interface contact surface. Shelled cells keep the clearance so
                    // crests bond beside the wall instead of riding it.
                    const double shell_gap = (cs_bubble && sl_outer > 0) ? 0.85 * double(spacing) : 0.0;
                    for (const WeaveRing &wr : wrings) {
                        Points &pts = cc.loops[wr.loop].points;
                        const size_t m = pts.size();
                        // Cumulative arc length + perimeter; origin anchored to the +x-most point so the
                        // fraction is radially aligned across the concentric rings (their anti-phase
                        // lobes then line up radially and actually interlock).
                        std::vector<double> arc(m, 0.0);
                        for (size_t k = 1; k < m; ++k)
                            arc[k] = arc[k - 1] + (pts[k] - pts[k - 1]).cast<double>().norm();
                        const double perim = arc[m - 1] + (pts[0] - pts[m - 1]).cast<double>().norm();
                        if (perim < 1.0) continue;
                        size_t origin = 0; double best = std::numeric_limits<double>::max();
                        for (size_t k = 0; k < m; ++k) {
                            const Vec2d d = pts[k].cast<double>() - centre;
                            const double a = std::abs(std::atan2(d.y(), d.x()));       // nearest angle 0 (+x)
                            if (a < best) { best = a; origin = k; }
                        }
                        const double s0 = arc[origin];
                        // PER-POINT self-collision cap: local amplitude limited by the ring's own
                        // clearance to its NON-neighbouring segments (subsampled). An elongated ring
                        // waves at full amplitude along its long sides - so it still reaches the next
                        // ring out - and flattens only at its narrow ends where its own sides close in.
                        std::vector<double> local_amp(m, wr.amp);
                        {
                            const size_t astep = std::max<size_t>(1, m / 96);
                            const size_t sstep = std::max<size_t>(1, m / 128);
                            for (size_t k = 0; k < m; k += astep) {
                                const Vec2d p = pts[k].cast<double>();
                                double dmin = std::numeric_limits<double>::max();
                                for (size_t e3 = 0; e3 < m; e3 += sstep) {
                                    const size_t sep = e3 > k ? e3 - k : k - e3;
                                    if (std::min(sep, m - sep) < m / 8)
                                        continue;                    // own neighbourhood along the loop
                                    const Vec2d a  = pts[e3].cast<double>();
                                    const Vec2d ab = pts[(e3 + sstep) % m].cast<double>() - a;
                                    const double l2 = std::max(1.0, ab.squaredNorm());
                                    const double t2 = std::clamp((p - a).dot(ab) / l2, 0.0, 1.0);
                                    dmin = std::min(dmin, (p - (a + t2 * ab)).norm());
                                }
                                double la = wr.amp;
                                if (dmin < std::numeric_limits<double>::max())
                                    la = std::min(wr.amp, std::max(0.0, 0.5 * (dmin - 0.6 * double(spacing))));
                                // Bore-obstacle cap: a crest may reach the FAN outline (outermost
                                // detour arc + flatten clearance) but never past it - capping at
                                // the base keep-out let crests land ON the solid fan lanes, which
                                // mounded a glossy over-extruded ring around every bore (die5_swoop
                                // print). (Outlines are polygons - non-round bores cap identically.)
                                for (const Polygon &keep : (amp_keepouts.empty() ? cc.obstacles : amp_keepouts)) {
                                    const Points &kp = keep.points;
                                    double db = std::numeric_limits<double>::max();
                                    for (size_t e4 = 0; e4 < kp.size(); ++e4) {
                                        const Vec2d a2  = kp[e4].cast<double>();
                                        const Vec2d ab2 = kp[(e4 + 1) % kp.size()].cast<double>() - a2;
                                        const double l22 = std::max(1.0, ab2.squaredNorm());
                                        const double t3  = std::clamp((p - a2).dot(ab2) / l22, 0.0, 1.0);
                                        db = std::min(db, (p - (a2 + t3 * ab2)).norm());
                                    }
                                    la = std::min(la, std::max(0.0, db));
                                }
                                for (size_t f2 = k; f2 < std::min(m, k + astep); ++f2)
                                    local_amp[f2] = la;
                            }
                        }
                        Points out; out.reserve(m);
                        for (size_t k = 0; k < m; ++k) {
                            const Vec2d p = pts[k].cast<double>();
                            const double f = std::fmod(arc[k] - s0 + perim, perim) / perim;   // arc fraction [0,1)
                            const double delta = local_amp[k] * zcarrier * std::sin(2.0 * M_PI * wr.lobes * f + wr.phase);
                            // Local outward normal (central difference), flipped to point AWAY from the
                            // centre so it is correct regardless of the ring's winding.
                            Vec2d nrm((pts[(k + 1) % m] - pts[(k + m - 1) % m]).y(),
                                     -(pts[(k + 1) % m] - pts[(k + m - 1) % m]).x());
                            if (nrm.dot(p - centre) < 0.0) nrm = -nrm;
                            const double nl = nrm.norm();
                            if (nl < 1e-6) { out.push_back(pts[k]); continue; }
                            const Vec2d q = p + (delta / nl) * nrm;
                            // Poke-through guard: keep the perturbed point inside the solid shell (and,
                            // for an annulus, outside the bore shell) via the per-angle shell profiles,
                            // less the flatten clearance above.
                            const Vec2d dd = q - centre;
                            double rr = dd.norm();
                            const double th = std::atan2(dd.y(), dd.x());
                            rr = std::min(rr, sample_min(outer_prof, th) - shell_gap);
                            if (annular) rr = std::max(rr, sample_max(inner_prof, th) + shell_gap);
                            rr = std::max(rr, 1.0);
                            Vec2d fin(centre.x() + rr * std::cos(th), centre.y() + rr * std::sin(th));
                            // Bore-obstacle safety net: whatever the morph/weave/clamp composed,
                            // no point may land inside a keep-out - snap escapees to the nearest
                            // boundary point (the amp cap above keeps this a rare last resort).
                            for (const Polygon &keep : cc.obstacles) {
                                if (!keep.contains(Point(coord_t(std::llround(fin.x())), coord_t(std::llround(fin.y())))))
                                    continue;
                                const Points &kp = keep.points;
                                double db = std::numeric_limits<double>::max();
                                Vec2d  nearest = fin;
                                for (size_t e4 = 0; e4 < kp.size(); ++e4) {
                                    const Vec2d a2  = kp[e4].cast<double>();
                                    const Vec2d ab2 = kp[(e4 + 1) % kp.size()].cast<double>() - a2;
                                    const double l22 = std::max(1.0, ab2.squaredNorm());
                                    const double t3  = std::clamp((fin - a2).dot(ab2) / l22, 0.0, 1.0);
                                    const Vec2d cand = a2 + t3 * ab2;
                                    const double d3  = (fin - cand).norm();
                                    if (d3 < db) { db = d3; nearest = cand; }
                                }
                                fin = nearest;
                            }
                            out.push_back(Point(coord_t(std::llround(fin.x())), coord_t(std::llround(fin.y()))));
                        }
                        pts = std::move(out);
                    }
                }
            } else if (woven.solid_stagger > 0.0 && thick
                       && (!cc.has_bore ||
                           (cc.loops.size() >= sl_bore + 1 &&
                            cc.loops[cc.loops.size() - sl_bore].contains(cc.loops.front().centroid())))) {
                // Solid cap/shell brick stagger (see build_cross_section). Shift only the interior rings
                // inward by solid_stagger this layer; the outer wall rings (< sl) and, for a tube, the
                // bore shell stay put. Bounded so a ring never invades the bore. The fill stays fully
                // solid (no decimation, no infill tag) - this only nudges interior geometry per layer.
                // GATED to chains whose bore ENCLOSES the chain centre: the bore clamp below is a
                // per-angle POLAR profile about the centre, which is only single-valued for a
                // centred annulus (the tube walls this was built for). On a multi-bore Voronoi
                // cell the bore is ECCENTRIC - its polar map is two-valued, the bin average
                // collapses the ring's near and far sides into a mid-bore radius, and the
                // max-clamp then shoved interior ring points INTO the bore (the odd-layer chords
                // visible through every die-5 print's pips). Eccentric-bore cells simply skip
                // the stagger - a solid-layer nicety they can live without.
                Vec2d centre(0., 0.);
                const Polygon &outer0 = cc.loops.front();
                if (!outer0.points.empty()) {
                    for (const Point &p : outer0.points) centre += p.cast<double>();
                    centre /= double(outer0.points.size());
                }
                const int NB = 360;
                std::vector<double> inner_prof;
                if (annular) {                                   // per-angle bore-shell radius (shape-aware)
                    inner_prof.assign(NB, 0.0); std::vector<int> cnt(NB, 0);
                    for (const Point &q : cc.loops[n - sl_bore].points) {
                        const Vec2d d = q.cast<double>() - centre;
                        int b = int((std::atan2(d.y(), d.x()) + M_PI) / (2.0 * M_PI) * NB); b = ((b % NB) + NB) % NB;
                        inner_prof[b] += d.norm(); ++cnt[b];
                    }
                    for (int b = 0; b < NB; ++b) if (cnt[b]) inner_prof[b] /= cnt[b];
                    double last = 0.0; for (int b = 0; b < 2 * NB; ++b) { int i = b % NB; if (inner_prof[i] > 0.0) last = inner_prof[i]; else if (last > 0.0) inner_prof[i] = last; }
                }
                for (size_t j = 0; j < n; ++j) {
                    if (!is_infill(j)) continue;                 // shift interior only; shells stay put
                    // Only NESTED rings (enclosing the chain centre) may be clamped against the
                    // per-angle bore profile. A pinched band's residual void POCKET sits mid-list
                    // without enclosing the centre; max-projecting its points onto the bore circle
                    // scattered them around it and the edges between them chorded straight
                    // THROUGH the bore (the die-5 pip chords' final mechanism).
                    if (annular && !cc.loops[j].contains(Point(coord_t(std::llround(centre.x())),
                                                               coord_t(std::llround(centre.y())))))
                        continue;
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

    if (cs_env_on("CS_DBG_CELLS")) {
        for (size_t ci2 = 0; ci2 < out.size(); ++ci2)
            for (size_t li2 = 0; li2 < out[ci2].loops.size(); ++li2) {
                const Points &lp = out[ci2].loops[li2].points;
                size_t inh = 0;
                for (size_t k2 = 0; k2 < lp.size(); ++k2) {
                    const Point &a2 = lp[k2];
                    const Point &b2 = lp[(k2 + 1) % lp.size()];
                    const double ll = (b2 - a2).cast<double>().norm();
                    const int ns2 = std::max(1, int(ll / scale_(0.3)));
                    for (int s3 = 0; s3 <= ns2 && inh == 0; ++s3) {
                        const Point q2 = a2 + (b2 - a2) * s3 / std::max(1, ns2);
                        for (const Polygon &h : region.holes)
                            if (h.contains(q2)) { ++inh; break; }
                    }
                    if (inh > 0) {
                        fprintf(stderr, "[ring-in-hole] L%zu chain=%zu loop=%zu/%zu pts=%zu edge=%zu (%.2f,%.2f)->(%.2f,%.2f) len=%.2f bore=%d solid=%d\n",
                                layer_index, ci2, li2, out[ci2].loops.size(), lp.size(), k2,
                                unscale<double>(a2.x()), unscale<double>(a2.y()),
                                unscale<double>(b2.x()), unscale<double>(b2.y()),
                                unscale<double>(ll), int(out[ci2].has_bore), int(out[ci2].from_solid));
                        for (size_t l3 = 0; l3 < out[ci2].loops.size(); ++l3) {
                            const Polygon &pl3 = out[ci2].loops[l3];
                            const Point c3 = pl3.centroid();
                            fprintf(stderr, "[loops]   chain=%zu loop=%zu pts=%zu area=%.1f cen=(%.2f,%.2f) ccw=%d\n",
                                    ci2, l3, pl3.points.size(),
                                    unscale<double>(unscale<double>(std::abs(pl3.area()))),
                                    unscale<double>(c3.x()), unscale<double>(c3.y()),
                                    int(pl3.is_counter_clockwise()));
                        }
                        break;
                    }
                }
            }
    }

    // Drop negligible sliver pieces (the last necking bits where a non-uniform wall
    // pinches off). Each is a tiny leftover that adds no real fill but forces a long
    // bridge across the bore to reach it. Keep any piece at least a small fraction of
    // the largest piece's footprint, so genuine split lobes (a dogbone, a '+') stay.
    // Compare MATERIAL footprints (marched band area), never contour extents: the peel
    // band's contour spans the whole part while its material is a thin annulus, and the
    // virtual-site cells are compact by design - contour-area ratios dropped every
    // interior cell of the die-5 behind its peel (skin-only gcode). The peel band itself
    // is the object's skin and is exempt outright.
    if (out.size() > 1) {
        double max_foot = 0.;
        for (const ContourChain &cc : out)
            if (!cc.peel_band)
                max_foot = std::max(max_foot, cc.footprint);
        const double min_foot = 0.05 * max_foot;
        out.erase(std::remove_if(out.begin(), out.end(),
                                 [&](const ContourChain &cc) { return !cc.peel_band && cc.footprint < min_foot; }),
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
                                                     const ExPolygons *bridge_over)
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

    ExPolygons merged = union_ex(regions);
    if (merged.size() != 1) {
        if (reason != nullptr)
            *reason = "Continuous shell mode requires one connected layer cross-section.";
        return false;
    }

    // BRIDGE PATCH (a blind bore's cap layer): each interior island sitting over AIR is
    // reintroduced as a HOLE for the ring march - the rings then wrap the opening exactly
    // like the bore below did, instead of sagging full-flow loops across it - and a
    // serpentine patch is appended to the tour (after the sanctity pass, which would
    // otherwise drop it) to span the opening with taut anchored lines. The NEXT solid
    // layer closes over the patch normally: classic two-stage bridging.
    ExPolygons bridge_islands;
    if (bridge_over != nullptr && !bridge_over->empty()) {
        const double min_area = 4.0 * double(flow.scaled_spacing()) * 4.0 * double(flow.scaled_spacing());
        for (const ExPolygon &isl : *bridge_over)
            if (std::abs(isl.contour.area()) >= min_area)
                bridge_islands.push_back(isl);
        if (!bridge_islands.empty()) {
            ExPolygons cut = diff_ex(merged, bridge_islands);
            if (cut.size() == 1)
                merged = std::move(cut);
            else
                bridge_islands.clear();          // odd geometry: print as before rather than split the region
        }
    }

    // PATHOLOGY GUARD: refuse outlines the boundary-driven fill cannot sanely cover, instead
    // of producing them badly. CS marches concentric offsets of the outline; a shape whose
    // inward offsets SHATTER (a fractal silhouette, filigree, engraved text) explodes into
    // thousands of chains whose rings overlap across fragments - the Mandeltree test case
    // sliced to ~1.4x its own SOLID volume in "sparse" fill (221g, 73m, 6.5h of soup).
    // Probe: offset inward by 2 and 6 beads and count surviving pieces. Legitimate CS parts
    // stay low (tube 1, die-5 <=6, herringbone gear ~25 pinched tooth tips); pathological
    // outlines fragment into hundreds. The threshold leaves generous headroom over the gear.
    {
        size_t worst = 0;
        for (float d : { 2.f, 6.f })
            worst = std::max(worst, offset_ex(merged, -d * float(flow.scaled_spacing())).size());
        if (worst > 64) {
            if (reason != nullptr)
                *reason = "Continuous shell mode cannot fill this outline: it fragments into " +
                          std::to_string(worst) + " pieces within a few bead widths of the surface "
                          "(detail finer than the bead). Simplify or smooth the outline, thicken "
                          "features to at least two bead widths, or slice without continuous shell.";
            return false;
        }
    }

    // Debug region dump (CS_DUMP_REGION=<path>): append each layer's merged contour so
    // divergences between unit-test polygons and real sliced contours can be replayed
    // exactly (GUI-subsystem builds cannot print to a console - HANDOFF file-logger rule).
    if (const char *dump = std::getenv("CS_DUMP_REGION")) {
        static std::mutex dump_mutex;
        std::lock_guard<std::mutex> lock(dump_mutex);
        if (FILE *f = std::fopen(dump, "a")) {
            std::fprintf(f, "L %zu W %.4f H %.4f S %lld R %d F %d C %zu", layer_index,
                         double(flow.width()), double(flow.height()),
                         (long long) flow.scaled_spacing(), int(reverse_path), int(force_solid),
                         merged.front().contour.points.size());
            for (const Point &p : merged.front().contour.points)
                std::fprintf(f, " %lld %lld", (long long) p.x(), (long long) p.y());
            std::fprintf(f, " H %zu", merged.front().holes.size());
            for (const Polygon &hole : merged.front().holes) {
                std::fprintf(f, " %zu", hole.points.size());
                for (const Point &p : hole.points)
                    std::fprintf(f, " %lld %lld", (long long) p.x(), (long long) p.y());
            }
            std::fprintf(f, "\n");
            std::fclose(f);
        }
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

    // DEEP-concavity gate (computed once per cross-section): a point of the noise-
    // simplified contour lying deeper inside the convex hull than two bead widths.
    // Round/convex parts stay miles under the threshold and keep the byte-stable
    // absolute-angle anchors; notched parts (gear teeth, slit collars, combs, spiral
    // bands, dogbones) switch to arc-fraction anchoring with gap-skipping joins,
    // un-reversed closed tours, and rail-routed connectors.
    bool deep_concave = false;
    {
        Polygons simplified = merged.front().contour.simplify(0.25 * double(flow.scaled_spacing()));
        if (!simplified.empty() && simplified.front().points.size() >= 4) {
            const Polygon &sc = simplified.front();
            const Polygon hull = Geometry::convex_hull(sc.points);
            const double  threshold = 2.0 * double(flow.scaled_spacing());
            const Lines   hull_edges = hull.lines();
            for (const Point &p : sc.points) {
                double d = std::numeric_limits<double>::max();
                for (const Line &e : hull_edges)
                    d = std::min(d, e.distance_to(p));
                if (d > threshold) { deep_concave = true; break; }
            }
        }
    }

    // Fill each connected piece (chain) as its own Fermat spiral, then stitch the pieces
    // together greedily - each next piece is appended at whichever end lands nearest the
    // running tail - so a split cross-section's pieces are joined by the shortest hop
    // (along the wall) instead of long chords across the bore. A simple cross-section has
    // exactly one chain, so this is identical to a single spiral there.
    struct Sub { Polyline pl; std::vector<char> flags; std::vector<coord_t> zoff; std::vector<float> fill; };
    std::vector<Sub> subs;
    // Per-layer reversal exists to give SINGLE-chain solid layers short inter-layer
    // joins (core end meets core start). A stitched MULTI-chain tour instead closes
    // back to its own seam-anchored start below, so every layer already starts and
    // ends at the (precessing) outer anchor - and reversing would break exactly that
    // pairing: odd layers would START at the buried inner ring while the previous
    // layer ended at the outer anchor, chording across any concave notch between.
    // DEEP-CONCAVE solids get the same treatment as multi-chain (no reversal + tour
    // closure below): reversal makes odd layers start at the buried core, where the
    // arc-anchor gap-skipping cannot protect the inter-layer join (observed as slit
    // hops on the clamp-collar shape).
    const bool closed_tour     = chains.size() > 1 || deep_concave;
    const bool reverse_chains  = closed_tour ? false : reverse_path;
    for (const ContourChain &chain : chains) {
        std::vector<char> flags;
        std::vector<coord_t> zoff;
        std::vector<float> fill;
        Polyline sub = spiralize_rings(chain.loops, reverse_chains, seam_step_rad, seam_jitter_rad,
                                       ring_spacing, layer_index, chain.has_bore, &flags, interlock_amp, &zoff,
                                       connector_spread, chain.woven_infill ? woven.shell_loops : 0,
                                       &merged.front(), &fill, deep_concave, !chain.obstacles.empty(),
                                       chain.ride_rails.empty() ? nullptr : &chain.ride_rails);
        if (sub.points.size() >= 2 && sub.points.size() == flags.size() && sub.points.size() == zoff.size() &&
            sub.points.size() == fill.size())
            subs.push_back({ std::move(sub), std::move(flags), std::move(zoff), std::move(fill) });
        else if (cs_env_on("CS_DBG_CELLS"))
            fprintf(stderr, "[spiral-drop] L%zu loops=%zu pts=%zu flags=%zu zoff=%zu fill=%zu bore=%d woven=%d peel=%d\n",
                    layer_index, chain.loops.size(), sub.points.size(), flags.size(), zoff.size(), fill.size(),
                    int(chain.has_bore), int(chain.woven_infill), chain.peel_depth);
    }
    if (subs.empty()) {
        if (reason != nullptr)
            *reason = "Continuous shell mode could not construct a continuous spiral path.";
        return false;
    }

    Polyline spiral;
    std::vector<char> internal;
    std::vector<char> connector;        // 1 on inter-chain connector points (thinned flow)
    std::vector<coord_t> interlock_z;   // per-point Z offset (brick interlock), parallel to spiral.points
    std::vector<float> wrap_fill;       // seam-valley compensation multiplier per point
    auto append_sub = [&](Sub &s, bool rev) {
        if (rev) {
            std::reverse(s.pl.points.begin(), s.pl.points.end());
            std::reverse(s.flags.begin(), s.flags.end());
            std::reverse(s.zoff.begin(), s.zoff.end());
            std::reverse(s.fill.begin(), s.fill.end());
        }
        spiral.points.insert(spiral.points.end(), s.pl.points.begin(), s.pl.points.end());
        internal.insert(internal.end(), s.flags.begin(), s.flags.end());
        interlock_z.insert(interlock_z.end(), s.zoff.begin(), s.zoff.end());
        wrap_fill.insert(wrap_fill.end(), s.fill.begin(), s.fill.end());
        connector.insert(connector.end(), s.pl.points.size(), char(0));
    };

    // Inter-chain connector routing. A split cross-section's sub-chains can END deep
    // inside a lobe (a solid chain finishes at its own centre), so the shortest
    // endpoint hop may cut straight across the concave notch between lobes - through
    // air (the dogbone chords). Route those connectors through material instead:
    // straight when the hop stays inside the region, otherwise out to a "rail" one
    // bead inside the wall, along the wall (through the waist), and back in. Rails
    // sit at 1.5 bead-widths so the connector rides on buried beads, not the surface.
    Polygons rails;
    if ((subs.size() > 1 || closed_tour) && !merged.empty()) {
        // Densify the rails (Clipper offsets keep only corner vertices): the connector
        // enters/leaves at the NEAREST VERTEX, and with corner-only rails that approach
        // segment can cut across a nearby notch corner. At ~1 mm spacing the nearest
        // vertex is essentially the nearest point, so approaches run straight into the
        // wall through material.
        auto densify = [](const Polygon &poly) -> Polygon {
            const double step = scale_(1.0);
            const size_t n = std::clamp<size_t>(size_t(poly.length() / step), 16, 4096);
            std::vector<Point> pts = resample_closed(poly, n);
            return pts.empty() ? poly : Polygon(Points(pts.begin(), pts.end()));
        };
        // Swoop mode publishes dedicated ride lanes (mid-fan-corridor rings, reserved
        // empty by the crest cap) - use those INSTEAD of the hole-derived rails, which
        // sit on the bore wall crown and mounded every ride onto it.
        Polygons swoop_rails;
        for (const ContourChain &c2 : chains)
            for (const Polygon &rr2 : c2.ride_rails)
                swoop_rails.push_back(rr2);
        if (cs_env_on("CS_DBG_CELLS"))
            fprintf(stderr, "[rails] L%zu chains=%zu swoop_rails=%zu\n", layer_index, chains.size(), swoop_rails.size());
        for (const ExPolygon &ep : offset_ex(ExPolygons{ merged.front() }, -1.5f * float(flow.scaled_width()))) {
            rails.push_back(densify(ep.contour));
            if (swoop_rails.empty())
                for (const Polygon &h : ep.holes)
                    rails.push_back(densify(h));
        }
        for (const Polygon &sr : swoop_rails)
            rails.push_back(densify(sr));
    }
    auto route_connector = [&](const Point &from, const Point &to) -> Polyline {
        const bool dbg = std::getenv("CS_DEBUG_ROUTE") != nullptr;
        const double hop = (to - from).cast<double>().norm();
        if (hop <= 2.0 * ring_spacing)
            return {};                                    // short local hop: leave straight
        Polyline direct;
        direct.points = { from, to };
        // A straight hop that stays in material is only harmless while it is SHORT. A LONG
        // contained hop - the multi-bore inter-cell join - lays a bead straight across the
        // already-filled interior, overlapping the cell fills into a raised ridge. Route long
        // contained hops along a rail too, so the connector rides on buried beads / hugs a
        // bore rail (near the shared cell edge) instead of ridging across the surface fill.
        const bool long_hop = hop > scale_(4.0);
        if (merged.front().contains(direct) && !long_hop) {
            if (dbg) std::printf("ROUTE (%.2f,%.2f)->(%.2f,%.2f): contained+short, straight\n",
                                 from.x() / 1e6, from.y() / 1e6, to.x() / 1e6, to.y() / 1e6);
            return {};                                    // short straight hop stays in material
        }
        if (rails.empty()) {
            if (dbg) std::printf("ROUTE (%.2f,%.2f)->(%.2f,%.2f): NO RAILS\n",
                                 from.x() / 1e6, from.y() / 1e6, to.x() / 1e6, to.y() / 1e6);
            return {};
        }
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
        const Points &pts = rails[best_rail].points;
        const size_t n = pts.size();
        // The APPROACH segments (from -> rail entry, rail exit -> to) must also stay
        // in material: entering at the blindly-nearest vertex can clip a notch corner
        // obliquely. Slide the entry/exit along the rail (nearest first, expanding)
        // until the approach segment is contained; rails are dense (~1mm), so a few
        // steps suffice. Falls back to the nearest vertex if nothing within reach
        // is contained (degenerate slivers) - still no worse than before.
        auto contained_approach = [&](const Point &p, size_t k) {
            Polyline seg; seg.points = { p, pts[k] };
            return merged.front().contains(seg);
        };
        auto slide_to_contained = [&](const Point &p, size_t k0) -> size_t {
            if (contained_approach(p, k0)) return k0;
            for (size_t off = 1; off <= std::min<size_t>(n / 2, 64); ++off) {
                const size_t kf = (k0 + off) % n, kb = (k0 + n - off) % n;
                if (contained_approach(p, kf)) return kf;
                if (contained_approach(p, kb)) return kb;
            }
            return k0;
        };
        ia = slide_to_contained(from, ia);
        ib = slide_to_contained(to, ib);
        if (ia == ib) {
            Polyline arc; arc.points = { pts[ia] };
            return arc;
        }
        // Walk the shorter arc ia -> ib along the rail.
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

    // Visit the sub-chains in a nearest-neighbour order on their CENTROIDS so consecutive
    // cells are spatial NEIGHBOURS (short inter-cell hops). The old per-step greedy on nearest
    // ENDPOINT could strand a far cell for last -> a long chord clear across a multi-bore part
    // (observed: a 26 mm jump on the 7-bore carrier). A single-chain part has one sub, so this
    // is identical to a plain spiral there.
    std::vector<Vec2d> cen(subs.size(), Vec2d(0.0, 0.0));
    for (size_t i = 0; i < subs.size(); ++i) {
        for (const Point &p : subs[i].pl.points) cen[i] += p.cast<double>();
        cen[i] /= double(std::max<size_t>(1, subs[i].pl.points.size()));
    }
    std::vector<size_t> order; order.reserve(subs.size());
    std::vector<char> visited(subs.size(), 0);
    order.push_back(0); visited[0] = 1;             // start on the first (outer) chain
    for (size_t k = 1; k < subs.size(); ++k) {
        size_t best = 0; double bd = std::numeric_limits<double>::max();
        for (size_t j = 0; j < subs.size(); ++j) {
            if (visited[j]) continue;
            const double d = (cen[j] - cen[order.back()]).squaredNorm();
            if (d < bd) { bd = d; best = j; }
        }
        order.push_back(best); visited[best] = 1;
    }
    if (cs_env_on("CS_DBG_CELLS"))
        fprintf(stderr, "[tour] L%zu sub=%zu(chain %zu) start=0 size=%zu\n",
                layer_index, size_t(0), order[0], subs[order[0]].pl.points.size());
    append_sub(subs[order[0]], false);
    for (size_t oi = 1; oi < order.size(); ++oi) {
        const size_t i = order[oi];
        const Point tail = spiral.points.back();
        const double ds = (subs[i].pl.points.front() - tail).cast<double>().squaredNorm();
        const double de = (subs[i].pl.points.back()  - tail).cast<double>().squaredNorm();
        const bool rev = de < ds;                    // enter at whichever end lands nearer the tail
        const Point head = rev ? subs[i].pl.points.back() : subs[i].pl.points.front();
        const Polyline arc = route_connector(tail, head);
        if (std::getenv("CS_DEBUG_ROUTE") != nullptr)
            std::printf("STITCH tail(%.2f,%.2f) -> head(%.2f,%.2f) rev=%d arc=%zu\n",
                        tail.x() / 1e6, tail.y() / 1e6, head.x() / 1e6, head.y() / 1e6,
                        int(rev), arc.points.size());
        for (const Point &rp : arc.points) {
            spiral.points.push_back(rp);
            internal.push_back(1);          // rails are buried by construction
            interlock_z.push_back(0);
            wrap_fill.push_back(1.f);
            connector.push_back(1);
        }
        const size_t head_idx = spiral.points.size();
        if (cs_env_on("CS_DBG_CELLS"))
            fprintf(stderr, "[tour] L%zu sub=%zu(chain %zu) start=%zu size=%zu rev=%d\n",
                    layer_index, oi, i, head_idx, subs[i].pl.points.size(), int(rev));
        append_sub(subs[i], rev);
        if (head_idx < connector.size())
            connector[head_idx] = 1;        // thin the final hop segment onto the chain head
    }
    if (closed_tour && spiral.points.size() >= 2) {
        // Close the tour back to its start. Convex single-chain layers keep their
        // inter-LAYER join short via reverse_path parity, but a stitched multi-chain
        // tour (or a deep-concave solid, which prints un-reversed) ends deep inside
        // its last lobe / at its buried core - with travels stripped, the NEXT layer's
        // first extrusion would then chord across the notch through air. Returning to
        // the (gap-skipped) start point makes the layer join a short surface hop again.
        const Point tail = spiral.points.back();
        const Point head = spiral.points.front();
        // NOTE: the arc must be materialized - ranging over `route_connector(...).points`
        // directly binds to a member of a destroyed temporary (dangling range-for).
        const Polyline back_arc = route_connector(tail, head);
        for (const Point &rp : back_arc.points) {
            spiral.points.push_back(rp);
            internal.push_back(1);
            interlock_z.push_back(0);
            wrap_fill.push_back(1.f);
            connector.push_back(1);
        }
        spiral.points.push_back(head);
        internal.push_back(1);
        interlock_z.push_back(0);
        wrap_fill.push_back(1.f);
        connector.push_back(1);
    }
    if (spiral.points.size() < 2) {
        if (reason != nullptr)
            *reason = "Continuous shell mode could not construct a continuous spiral path.";
        return false;
    }

    // BORE SANCTITY pass: no segment of the finished tour may cross a HOLE of the region.
    // Every routing mechanism above tries to avoid it, but each has fallback branches that
    // can emit a straight chord when nothing better is found (observed on every die-5 print:
    // the solid bottom shells' inter-cell connectors chorded two bores at one layer, visible
    // straight through the pips). Catch them ALL here: any hole-crossing segment is rerouted
    // around that hole on a densified ring one bead outside it, at connector half flow.
    if (!merged.empty() && !merged.front().holes.empty()) {
        struct HoleRail { Polygon hole_ccw; BoundingBox bb; Polygon rail; };
        std::vector<HoleRail> hrails;
        for (const Polygon &h : merged.front().holes) {
            Polygon hc = h;
            hc.reverse();
            Polygons r2 = offset(hc, 1.2f * float(flow.scaled_spacing()));
            if (r2.size() != 1) continue;
            const double step = scale_(0.6);
            const size_t nn = std::clamp<size_t>(size_t(r2.front().length() / step), 16, 2048);
            std::vector<Point> rs = resample_closed(r2.front(), nn);
            if (rs.size() < 8) continue;
            BoundingBox bb = hc.bounding_box();
            bb.offset(coord_t(scale_(0.5)));
            hrails.push_back({ std::move(hc), bb, Polygon(Points(rs.begin(), rs.end())) });
        }
        if (!hrails.empty()) {
            auto seg_crosses = [](const Point &a, const Point &b, const HoleRail &hr) {
                if (!hr.bb.overlap(BoundingBox(Points{ a, b })))
                    return false;
                // crossing = midpoint inside, an endpoint inside, or the segment intersects
                // the hole outline (sampled - segments here are short post-resample)
                if (hr.hole_ccw.contains(Point((a + b) / 2))) return true;
                if (hr.hole_ccw.contains(a) || hr.hole_ccw.contains(b)) return true;
                const double len = (b - a).cast<double>().norm();
                const int ns = std::min(16, std::max(2, int(len / scale_(0.3))));
                for (int s2 = 1; s2 < ns; ++s2)
                    if (hr.hole_ccw.contains(a + (b - a) * s2 / ns)) return true;
                return false;
            };
            auto in_hole = [&](const Point &p) -> const HoleRail * {
                for (const HoleRail &hr : hrails)
                    if (hr.bb.contains(p) && hr.hole_ccw.contains(p))
                        return &hr;
                return nullptr;
            };
            auto detour = [&](const Point &a, const Point &b, const HoleRail &hr,
                              Points &np, std::vector<char> &ni, std::vector<char> &nc,
                              std::vector<coord_t> &nz, std::vector<float> &nf) {
                const Points &rp = hr.rail.points;
                const size_t rn = rp.size();
                size_t ia = 0, ib = 0; double da = 1e300, db2 = 1e300;
                for (size_t k = 0; k < rn; ++k) {
                    const double d1 = (rp[k] - a).cast<double>().squaredNorm();
                    const double d2 = (rp[k] - b).cast<double>().squaredNorm();
                    if (d1 < da) { da = d1; ia = k; }
                    if (d2 < db2) { db2 = d2; ib = k; }
                }
                if (ia == ib) return;
                double fwd = 0., bwd = 0.;
                for (size_t k = ia; k != ib; k = (k + 1) % rn) fwd += (rp[(k + 1) % rn] - rp[k]).cast<double>().norm();
                for (size_t k = ia; k != ib; k = (k + rn - 1) % rn) bwd += (rp[(k + rn - 1) % rn] - rp[k]).cast<double>().norm();
                auto emit = [&](const Point &p) {
                    np.push_back(p); ni.push_back(1); nc.push_back(1); nz.push_back(0); nf.push_back(1.f);
                };
                if (fwd <= bwd)
                    for (size_t k = ia; ; k = (k + 1) % rn) { emit(rp[k]); if (k == ib) break; }
                else
                    for (size_t k = ia; ; k = (k + rn - 1) % rn) { emit(rp[k]); if (k == ib) break; }
            };
            Points            np;
            std::vector<char> ni, nc;
            std::vector<coord_t> nz;
            std::vector<float>   nf;
            const size_t n0 = spiral.points.size();
            np.reserve(n0 + 64); ni.reserve(n0 + 64); nc.reserve(n0 + 64); nz.reserve(n0 + 64); nf.reserve(n0 + 64);
            const HoleRail *gap_hole = nullptr;   // a run of in-hole points was DROPPED
            for (size_t i = 0; i < n0; ++i) {
                const Point &b = spiral.points[i];
                if (const HoleRail *hr = in_hole(b)) {
                    if (cs_env_on("CS_DBG_CELLS"))
                        fprintf(stderr, "[sanctity] L%zu drop i=%zu/%zu (%.2f,%.2f) internal=%d connector=%d wrap=%.2f\n",
                                layer_index, i, n0, unscale<double>(b.x()), unscale<double>(b.y()),
                                int(internal[i]), int(connector[i]), double(wrap_fill[i]));
                    gap_hole = hr;                // drop the point outright - nothing may
                    continue;                     // deposit inside a bore, ever
                }
                if (!np.empty()) {
                    const HoleRail *cross = gap_hole;
                    if (cross == nullptr)
                        for (const HoleRail &hr : hrails)
                            if (seg_crosses(np.back(), b, hr)) { cross = &hr; break; }
                    if (cross != nullptr)
                        detour(np.back(), b, *cross, np, ni, nc, nz, nf);
                }
                gap_hole = nullptr;
                np.push_back(b); ni.push_back(internal[i]); nc.push_back(connector[i]);
                nz.push_back(interlock_z[i]); nf.push_back(wrap_fill[i]);
            }
            if (np.size() >= 2) {
                spiral.points = std::move(np);
                internal      = std::move(ni);
                connector     = std::move(nc);
                interlock_z   = std::move(nz);
                wrap_fill     = std::move(nf);
            }
        }
    }

    // Lay the bridge patches (appended AFTER the sanctity pass: they deliberately cross the
    // reintroduced cap holes, which the pass would otherwise strip). Each patch is a solid
    // serpentine over the island grown ~one bead, so the line ends anchor on the material
    // ringing the opening; the approach rides the rails at connector half-flow.
    for (const ExPolygon &isl : bridge_islands) {
        ExPolygons grown = offset_ex(ExPolygons{ isl }, 1.2f * float(flow.scaled_spacing()));
        if (grown.size() != 1)
            continue;
        bool jo = false;
        // Alternate the scan direction per layer so a stacked patch (a bore cap's shell
        // stack) cross-hatches like any solid shell instead of laying parallel courses.
        const double patch_angle = (layer_index % 2 != 0) ? M_PI / 4.0 : -M_PI / 4.0;
        Polyline bp = serpentine_fill(grown.front(), flow.scaled_spacing(), patch_angle, isl.contour.centroid(), &jo);
        if (bp.points.size() < 2)
            continue;
        const Point tail = spiral.points.back();
        const Polyline arc = route_connector(tail, bp.points.front());
        for (const Point &rp : arc.points) {
            spiral.points.push_back(rp);
            internal.push_back(1); interlock_z.push_back(0); wrap_fill.push_back(1.f); connector.push_back(1);
        }
        for (size_t bi2 = 0; bi2 < bp.points.size(); ++bi2) {
            spiral.points.push_back(bp.points[bi2]);
            internal.push_back(1); interlock_z.push_back(0); wrap_fill.push_back(1.f);
            connector.push_back(bi2 == 0 ? char(1) : char(0));   // thin only the hop onto the patch
        }
    }

    // Per-point flow multiplier combining three effects:
    //  * internal over-extrusion densifies provably-buried beads (knob, default off),
    //  * wrap-window seam-valley fill compensates the once-per-wrap tall slot each ring's
    //    ramp share creates at its start meridian (derived from geometry, always on), and
    //  * seam pressure trim thins both sides of the layer join to fight the bulge
    //    (seam_coast knob, symmetric pre+post-join de-pressurization, default off).
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
    if (have_interlock && have_flags && !merged.empty() && any_woven && woven.zamp > 0.0) {
        Vec2d centre(0., 0.);
        const Polygon &oc = merged.front().contour;
        if (!oc.points.empty()) {
            for (const Point &p : oc.points) centre += p.cast<double>();
            centre /= double(oc.points.size());
        }
        const double zamp  = std::clamp(woven.zamp, 0.0, 0.5) * scale_(layer_height);
        const double phase = double(layer_index) * woven.dphase + woven.zphase;
        for (size_t i = 0; i < spiral.points.size(); ++i) {
            if (!internal[i] || connector[i]) continue;   // buried (woven) infill only - not shells, not connectors
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
    std::vector<double> fmul(spiral.points.size(), 1.0);
    if (have_flags && over_pct > 0.0)
        for (size_t i = 0; i < fmul.size(); ++i)
            if (internal[i] && !connector[i]) fmul[i] *= over_mult;

    // Wrap-window seam-valley fill (derived, on by default): points under the
    // once-per-wrap tall slot get their ring's ramp-share of extra flow. See
    // spiralize_rings. CS_WRAP_FILL=<percent> scales the dose for print A/B tests
    // (0 disables, 100 = full derived dose) - print evidence 2026-07-16: full dose
    // is a lateral move on seam appearance (wider band, no ridge, slight valley
    // remains), kept for its buried-void-closing strength rationale pending the
    // strength-sweep verdict.
    double fill_scale = 1.0;
    if (const char *fs = std::getenv("CS_WRAP_FILL"))
        fill_scale = std::max(0.0, std::atof(fs) / 100.0);
    if (wrap_fill.size() == fmul.size())
        for (size_t i = 0; i < fmul.size(); ++i)
            fmul[i] *= 1.0 + (double(wrap_fill[i]) - 1.0) * fill_scale;

    // Thin-wall flow taper. Where the part is locally thinner than TWO beads, opposing
    // wall passes share the same space - the skin loop's two sides down a tapering tail,
    // or any sub-2-bead wall stretch - and constant flow over-stuffs it (approaching 200%
    // at a cusp), venting as dots of excess plastic on the surface. Taper each NON-BURIED
    // point's flow by the local part thickness: full flow at >= 2 beads, fair sharing
    // below (f = t / 2w, floored at 0.4 for extrusion stability). Thickness is probed by
    // stepping along the local normal and testing region containment; a single deep
    // inward probe accepts the common thick case, so the ladder only runs in thin zones.
    if (have_flags && !merged.empty()) {
        const double bw = double(flow.scaled_width());
        const ExPolygon &reg = merged.front();
        auto inside = [&](const Vec2d &q) {
            return reg.contains(Point(coord_t(std::llround(q.x())), coord_t(std::llround(q.y()))));
        };
        const size_t np = spiral.points.size();
        for (size_t i = 0; i < np; ++i) {
            if (internal[i] || connector[i]) continue;
            const Vec2d p = spiral.points[i].cast<double>();
            Vec2d t = (spiral.points[std::min(i + 1, np - 1)] - spiral.points[i == 0 ? 0 : i - 1]).cast<double>();
            if (t.squaredNorm() < 1.0 && i + 1 < np)
                t = (spiral.points[i + 1] - spiral.points[i]).cast<double>();   // U-turn apex: one-sided
            if (t.squaredNorm() < 1.0 && i > 0)
                t = (spiral.points[i] - spiral.points[i - 1]).cast<double>();
            const double tl = t.norm();
            if (tl < 1.0) continue;
            const Vec2d n(-t.y() / tl, t.x() / tl);
            // SYMMETRIC probing - a bead has its own lane whenever material extends ~2 beads
            // on EITHER side (a bore-side shell has its bore right behind it but plenty of
            // wall in front; picking a single "inward" direction mis-read that as thin).
            if (inside(p + n * (1.6 * bw)) || inside(p - n * (1.6 * bw)))
                continue;                                       // full flow
            // genuinely thin: ladder both directions for the local total thickness
            auto ladder = [&](const Vec2d &dir) {
                double best = 0.0;
                for (double s : { 0.3, 0.6, 0.9, 1.2, 1.45 })
                    if (inside(p + dir * (s * bw))) best = s * bw; else break;
                return best;
            };
            const double thickness = ladder(n) + ladder(-n) + 0.15 * bw;  // + probe resolution allowance
            fmul[i] *= std::clamp(thickness / (2.0 * bw), 0.4, 1.0);
        }
    }

    // Seam pressure-transient trim (seam_coast). SYMMETRIC (2026-07-18): the physical
    // seam bulge is melt oozing at the layer-join corner. Reducing flow AFTER the join
    // (the original post-join trim) only lowered the region's baseline, not the bulge
    // PEAK - the ooze comes from pressure stored on APPROACH, which firmware pressure
    // advance does not cancel on these z-contoured moves. So trim BOTH sides of the
    // join with the same exponential decay: the PRE-join half (backward from the tour
    // end) de-pressurizes the nozzle before the corner so there is less to ooze; the
    // POST-join half (forward from the tour start) trims whatever still surges out. The
    // tour end and start are the two sides of the once-per-layer join. coast% = trim
    // depth at the join; decays over ~16mm; 0 = off (default).
    if (coast_pct > 0.0 && fmul.size() >= 2) {
        const double c0  = std::clamp(coast_pct, 0.0, 95.0) / 100.0;
        const double tau = scale_(5.0);            // decay length ~5mm -> dead by ~16mm
        double s = 0.0;                            // POST-join: forward from tour start
        for (size_t i = 1; i < fmul.size(); ++i) {
            s += (spiral.points[i] - spiral.points[i - 1]).cast<double>().norm();
            if (s > 4.0 * tau) break;
            fmul[i] *= 1.0 - c0 * std::exp(-s / tau);
        }
        s = 0.0;                                   // PRE-join: backward from tour end
        for (size_t i = fmul.size() - 1; i > 0; --i) {
            s += (spiral.points[i] - spiral.points[i - 1]).cast<double>().norm();
            if (s > 4.0 * tau) break;
            fmul[i] *= 1.0 - c0 * std::exp(-s / tau);
        }
    }
    // Quantize to 2% steps so the run-length path emitter below produces a handful of
    // sub-paths per layer instead of one per point along the smooth trim/fill ramps.
    for (double &m : fmul)
        m = std::round(m * 50.0) / 50.0;
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

    // Debug output dump (CS_DUMP_PATH=<path>): the final planned polyline, for
    // bisecting artifacts between the path builder and the gcode pipeline.
    if (const char *dump = std::getenv("CS_DUMP_PATH")) {
        static std::mutex dump_mutex;
        std::lock_guard<std::mutex> lock(dump_mutex);
        if (FILE *f = std::fopen(dump, "a")) {
            std::fprintf(f, "L %zu", layer_index);
            for (const ExtrusionPath &ep : path.paths)
                for (const Point3 &p : ep.polyline.points)
                    std::fprintf(f, " %lld %lld", (long long) p.x(), (long long) p.y());
            std::fprintf(f, "\n");
            std::fclose(f);
        }
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
