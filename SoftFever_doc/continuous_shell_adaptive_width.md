# Continuous Shell: Adaptive Bead Width — Design

Status: C-lite, A (scarf; print-tested NEGATIVE for TPU, default stays 0 — see
`test_models/tpu_suite/RESULTS.md`), and C-full (`continuous_shell_adaptive_beading`,
default OFF pending print validation) implemented. Sequencing agreed 2026-08-12:
**C-lite → A (scarf) → C-full → B (spiral)**; B remains.

## 1. What the helical seam actually is

Three distinct artifacts, treated differently by each option:

1. **Per-ring join diagonal.** Every ring is a complete closed loop; the bead then hops one
   spacing inward (two on the Fermat arms) over a short arc. The diagonal crosses a region the
   neighbouring rings already cover → a local ~50% overfill ridge. `continuous_shell_seam_coast`
   thins flow there; `continuous_shell_connector_spread` fans the ridges apart. Neither removes them.
2. **Surface dive** — the once-per-layer point where the bead leaves ring 0. Same mechanism, visible.
3. **Layer join** — solved by the arc-length Z ramp (`apply_continuous_z`).

Topological invariant: at a given radius the bead covers the circumference exactly once, so the
one-spacing-per-revolution descent must happen *somewhere*. Every design option is a choice of
**how long an arc the descent spreads over, per ring class**:

- **Option A (scarf)** = short transition windows at every join.
- **Option B (interior spiral)** = full-revolution transitions for buried rings; surfaces stay
  closed loops with compact scarfed joins.
- **Option C (BeadingStrategy)** = the machinery deciding bead **count and widths** — which is
  what makes B's boundary geometry close without voids (see §4). They compose; they don't compete.

## 2. Platform constraints (bind all options)

Established by code survey 2026-08-12; line numbers as of commit after the review fixes.

| Constraint | Where |
|---|---|
| `ExtrusionPath` has no per-point width — scalar `width`/`mm3_per_mm`/`height` only. Variable width = many constant-width sub-paths (Arachne quantizes width into 0.05 mm buckets: `VariableWidth.cpp:29-63`) or flow-only modulation (precedents: `z_contoured` ratio `GCode.cpp:7087-7106`, scarf `e_ratio` `GCode.cpp:7114-7123`) | `ExtrusionEntity.hpp:162-318` |
| Sub-paths of one multipath emit back-to-back, no travel/retract/seam between them | `GCode.cpp:6117-6142` |
| Preview width is sticky-forced from `;WIDTH:` tags (`m_forced_width`); changing `mm3_per_mm` without `path.width` makes the preview lie; changing width without flow doesn't change the print. Keep them consistent via `Flow::with_width` | `GCode.cpp:6867-6871`, `GCodeProcessor.cpp:3943` |
| `path.width` feeds real logic, not just preview: volumetric speed cap (`GCode.cpp:6629`), seam gap, `polygons_covered_by_width` | various |
| Hard flow floor: `Flow::mm3_per_mm` throws below width ≈ 0.215×height; practical floor `min_bead_width` (85% nozzle default) | `Flow.cpp:227` |
| Volumetric cap makes F ∝ 1/mm3_per_mm per sub-path → feedrate steps at width changes | `GCode.cpp:6689` |
| Resonance-avoidance bi-quantizes external-perimeter speed → large F jumps at transitions. **Gated off under continuous shell** (phase-1 patch) | `GCode.cpp:6691-6725` |
| Adaptive PA is absent from the spiral-vase pipeline (its `;PA_CHANGE` tags are never consumed) → flow transitions must be long/gradual | `GCode.cpp:3763-3770` |
| PressureEqualizer runs after SpiralVase, splits/ramps F at flow steps, preserves per-point Z — helpful | `PressureEqualizer.cpp:471-509` |
| SpiralVase XY smoothing is forced on in continuous shell and pulls points toward the previous layer by (1−factor) — smears deliberate centerline-shift geometry (previous layer's transition is at a precessed angle). Width/flow modulation is immune; escape hatch `spiral_mode_max_xy_smoothing = 0` | `SpiralVase.hpp:27`, `SpiralVase.cpp:248-277` |
| SpiralVase short-segment filter is lossy in E; threshold 2×`resolution` ≈ 0.02 mm | `SpiralVase.cpp:182,213` |
| Whole shell is `erExternalPerimeter`: `dont_slow_down_outer_wall` disables layer-time cooling slowdown for everything | `CoolingBuffer.cpp:432-437` |
| Per-sub-path Z-sync can emit a bare mid-bead `G1 Z`, silently swallowed by the SpiralVase filter — watch when debugging Z discontinuities | `GCode.cpp:6463-6469`, `SpiralVase.cpp:200-210` |

## 3. Arachne reuse boundary

The **entire `BeadingStrategy` stack is pure, stateless, standalone-callable** (proven in-tree:
`tests/libslic3r/test_arachne_walls.cpp:247-267`). `compute(thickness, bead_count)` with a count
that deliberately disagrees with `getOptimalBeadCount`, blended against `compute(t, count+1)` by
`SkeletalTrapezoidation::interpolate` (~20 lines, `SkeletalTrapezoidation.cpp:1712-1731`), is the
entire transition mechanism; the strategies guarantee both branches vary continuously in thickness.

Graph-bound and NOT reusable: transition anchor placement + arclength ramp walk
(`SkeletalTrapezoidation.cpp:746-1200`). The angular re-derivation for rings is ~100–150 lines:
root-find θ* where t(θ) crosses `getTransitionThickness(n)`, straddle with an arc window of
`transition_length` split by `getTransitionAnchorPos`, blend beadings inside the window.

Gotchas when reusing: feed **spacings** as the preferred widths; derive split/add thresholds from
extrusion **widths** (`WallToolPaths.cpp:516-522`); handle `LimitedBeadingStrategy`'s zero-width
marker walls; `transitioning_angle` is inert inside the strategies.

Config: the relevant options (`wall_transition_length/angle/filter_deviation`,
`wall_distribution_count`, `min_bead_width`, `min_feature_size`) are **object-scope percent of the
min nozzle diameter** — reuse them rather than minting new ones (GUI note: they may be hidden when
the wall generator isn't Arachne; toggle logic needs a continuous-shell clause when C-full lands).

## 4. Phases

### C-lite (implemented, phase 1)

Insight that collapsed the scope: marching only stops when another band cannot fit, so the medial
void `v = 2·d_res − w ≤ 2·spacing − w ≈ 0.8·w` **always** — no medial-bead-count planning needed.

- Measure the vanished band's residual half-thickness by bisection (`vanish_depth`); close the
  void by widening the innermost ring pair by v/2 each, centrelines shifted v/4 toward the medial
  (inner flank extends by exactly v/2 → beads meet at the medial; outer flank untouched).
  Replaces the heuristic medial-strip patch. Skipped for sparsified chains and pinch/split chains.
- Sub-bead fallback fits the bead to the measured wall (clamped to `min_bead_width`), declines
  below `min_feature_size` (WideningBeadingStrategy semantics).
- Emission: per-point width array → sub-paths split on (flow multiplier, width), each built with
  its own `Flow::with_width` so `path.width`/`mm3_per_mm` stay consistent (truthful preview).
- Prereq patches: resonance-avoidance gated off in continuous-shell mode; honest per-sub-path width.
- The strategy stack itself enters at C-full; C-lite implements its degenerate uniform case
  directly and validates all the emission plumbing.

### A — scarf at the joins (implemented, phase 2)

Per-point width taper over a window at each ring join (`continuous_shell_scarf_length`, mm,
default 0 = off): the outgoing ring's tail narrows linearly to a floor, the incoming ring's
head widens back, and the junction diagonal itself is laid at the floor. Width/flow only —
no centerline moves (SpiralVase smoothing smears those; flow passes through untouched).
Design decisions that emerged in implementation:

- **Floor is transitional, not structural**: `max(w/2, nozzle/4)` — a scarf bead overlaps the
  junction zone, so it is not bound by `min_bead_width` (which is for standalone beads).
- **Free-surface rings are never tapered** — narrowing a surface bead recedes its outer flank
  into a visible relief groove; only the diagonal and the buried partner soften those joins.
- **Window capped at a quarter of the ring's own arc** so aggressive scarf lengths never thin a
  small inner ring wall-to-wall (at most half of any ring is below full width).
- Tapers compute from the ring's ORIGINAL width (overlapping head/tail windows on short rings
  otherwise compound below the floor; overlaps keep min()).
- Widths quantized to Arachne's 0.05 mm buckets so a taper is a handful of sub-paths.
- Woven chains excluded: their inter-ring diagonals are structural struts across open space.
- `seam_coast` is superseded but kept functional (both compose; tooltip documents).

### C-full — per-angle beading (implemented + print-validated, phase 3, 2026-08-15)

Implemented as `continuous_shell_adaptive_beading` (region bool, **default ON** — the T4
A/B re-print eliminated both the visible transition band and the stretched-regime medial
void line; see `test_models/tpu_suite/RESULTS.md`). Scope: dense annular star single-chain sections; everything else keeps C-lite
bit-identically (byte-identity regression-tested). Thickness field t(θ) from 720-bin radial
profiles with **nearest-distance** thickness (radial difference overstates 1/cos α on oblique
walls); rings are iso-index curves placed by radial FRACTION (exact perpendicular depth on
straight walls; single-valued r(θ) ⇒ self-intersection-free on star sections). Per-point widths
are a function of θ evaluated at emission (`ChainBeadPlan::width_at`), so resampling can never
desynchronize a width array. Strategy stack per §3's recipe (`build_chain_plan`), with the
following decisions made in implementation — they refine the printed-T4/T5 brief:

- **Count law** (`squeeze_biased_count`): never stretch an inner bead past **+10%** of nominal
  spacing (T4 voided at +29%; the flip is NOT at the equidistant thickness), clamped to
  [optimal, optimal+1] — single-step transitions by construction. The appearing bead is
  TRANSITIONAL while narrow: its floor is w_term = max(0.30 mm, 0.6·min_bead_width) — the T5
  die-swell floor — not the structural `min_bead_width` (else the law could not flip before
  the stretch reached the void regime).
- **Uniform interior distribution**: the factory gets an unbounded distribution radius — the
  two surface beads stay pinned at nominal (Redistribute) and ALL interior beads share the
  squeeze uniformly, the consolidation pattern of T4's squeezed band. `wall_distribution_count`
  is deliberately not consumed.
- **`SkeletalTrapezoidation::interpolate` is NOT reused** (its index-wise pairing would ramp the
  bore bead and un-pin the bore surface, and it drops the appearing bead at equal thickness) —
  `blend_ring_beadings` pairs front/back with the unpaired middle slot as the appearing bead.
- **Width-birth transitions — no open arcs.** The plan's medial bead never terminates: through a
  thin sector (eccentric wall) or below full blend (taper flip band, 0.25 spacings of thickness
  ≈ 20+ layers of height) it thins to w_term as a FULL ring and its neighbours squeeze to absorb
  the difference (exact conservation). Open coverage arcs were implemented and rejected: a young
  arc's far tip always costs either a bore-crossing chord or a displaced layer join (measured
  9 mm), and the printed T5 verdict says bead terminations gap — so the topology changes only at
  the single flip layer, as a one-floor-width step. This supersedes the brief's "stagger the
  flip plane": there is no flip plane left to stagger.

B's boundary wedge remains the next consumer of the plan machinery.

### B — true interior spiral (implemented, phase 4; v4 Z-wave 2026-08-17)

Implemented as `continuous_shell_interior_spiral` (region bool, default OFF pending a print;
requires adaptive beading and ≥4 buried rings; everything else keeps ring mode bit-identically).

**The print campaign's most important lesson (revised across v1→v3): radial INTERDIGITATION —
not thermal adjacency — is the master strength variable.** v1 (double-arm) delivered the seam
goal (junction census ~k → ≤1/layer, fingernail-flat seam) but delaminated into helical sheets;
we first blamed the armA↔armB cold joints. v2 (single arm, warm one-revolution-apart radial
welds) still peeled — falsifying the thermal explanation — and v3 (fixed traversal, doubling
gone) peeled too, the unrolled peel front tracking the precessing seam helix. The real
discriminator: ring mode's Fermat ORDER scrambles the arc-length ramp fraction across radial
neighbours (adjacent rings sit up to ~(k−1)/k·h apart in Z, varying with angle — the graded
incommensurate nesting the interlock A/B tear test proved tear-proofs the wall), and it holds
even though its outer interfaces are laid up to k revolutions apart in time. A monotone spiral
orders path position radially, so radial neighbours differ by only h/k — every interface is a
smooth laminate cleavage plane. Thermal adjacency is demoted to a secondary weld-quality factor.

**v2: a SINGLE arm at pitch one slot per revolution.** Radially adjacent turns are consecutive
revolutions of the same arm — deposited exactly one revolution apart, warmer than ring mode's
own worst case. This also deletes the double arm's whole complexity budget: no interleave
constraint (the moiré analysis), no U-turn, no s→2s wedge, no π-per-layer phase. What remains:

- Comb has the SAME bead count and pitch as ring mode (no stretch); arm tips tangent-merge onto
  the closed surface rings at clearance widths floored at the 0.30 transitional width — no bead
  terminates in air; both junctions ride the (precessing) seam azimuth.
- Traversal is FIXED outer→bore every layer (v3, 2026-08-17): the v2 per-layer reversal printed
  a staggered bowtie-void lattice — reverse-join doubling, the codebase's own documented reason
  annular ring mode never reverses (a reversed layer puts each spatial bead at the opposite
  arc-ramp fraction, pairing consecutive layers' beads with a two-layer void between pairs).
  The inter-layer join is one radial in-material diagonal across the wall at the precessing
  seam azimuth — a single event per layer, total radial travel equal to ring mode's hops.
- Generated analytically in plan-slot space (fractional-slot interpolation over the C-full
  tables, `frac_at_f`/`width_at_f`); the historical "pairwise ring-morph parametrization" is
  unnecessary. The explicit interlock knob is N/A on the spiral interior. Width-birth count
  transitions flow through the shared tables.
- **v4 — the intrinsic Z-wave (2026-08-17): printed, NO PEEL — mechanism confirmed.** Buried
  arm points carried z = (h/4)·sin(2πg·v), g the golden-ratio conjugate — a Z-only change, and
  the tube held. But the section showed organized void rows: a sinusoid makes the
  interdigitation depth oscillate ALONG each interface, so every interface passes through the
  void-producing intermediate/coplanar configurations somewhere. Ring mode's interfaces hold
  one constant depth each (both neighbours' ramp fractions advance together).
- **v5 — the per-revolution staircase (2026-08-17):** each buried revolution holds one constant
  Z level from the golden low-discrepancy sequence in [−h/2, +h/2] (amplitude = ring mode's own
  bead-Z spread), stepped over a 0.08-rev smoothstep arc just past the seam azimuth — where
  ring mode's own order boundaries live. Adjacent-interface depths are the golden gaps, always
  0.38 h or 0.62 h: never coplanar, constant per interface, scattered across interfaces — ring
  mode's full property set inside a spiral. Not a knob (the spiral has no natural bricklaying
  to fall back on); tapers out over one revolution before each tangent merge; surfaces pinned;
  layer 0 flat. Prediction: void rows collapse to ring-mode pinpricks with peel still gone.

## 5. Seam-impact / risk summary

| Phase | Seam impact | Other wins | Risk |
|---|---|---|---|
| C-lite | indirect (void-free, width-continuous medial; fewer improvised joins) | kills medial gap, ring-pop voids, sub-bead over-extrusion | low |
| A | direct: de-lumps every join ridge + surface dive | subsumes seam_coast | low-med |
| C-full | enables B; per-angle density uniformity | non-uniform walls | med |
| B | interior joins eliminated; 2 joins/layer total | most uniform interior | high |
