# Continuous Shell — Handoff (2026-08-14)

For the session that picks up **once the TPU test prints are in hand** — likely on a
**different machine**, so this doc is self-contained: everything described here
(code, test suite, sliced `.gcode.3mf` files, tooling scripts, this doc) **is
committed to the repo**. Do not rely on agent memory files from the original
machine — they do not transfer; anything they held that matters is inlined below.
Companion: `continuous_shell_adaptive_width.md` (the full design analysis — read it
before touching geometry code).

## Where the work stands

The fork's continuous-shell mode (one unbroken bead per object, per-point-Z ramp,
travel-free) got, in order:

1. **Review fixes** — commit `73c26c11ad`: star-shape gate (declines
   weave/stagger/Z-wave on folded sections), per-ring adaptive resampling, hardened
   medial fallback, debug-crumb removal.
2. **C-lite (adaptive bead width, phase 1)**: medial deficit closed by
   widening+recentring the innermost ring pair (measured via vanish-depth
   bisection; `v = 2·d_res − w ≤ ~0.8 beads` always, so widening always suffices);
   sub-bead walls width-fitted with a `min_feature_size` decline; per-point width
   emission via `Flow::with_width` (truthful preview); resonance avoidance gated off
   in continuous-shell mode.
3. **Option A: junction scarf**: `continuous_shell_scarf_length` (region option, mm,
   default 0) tapers bead width to `max(w/2, nozzle/4)` across every ring junction;
   surface rings never tapered; window capped at ring_arc/4; 0.05 mm width buckets;
   woven chains excluded. Supersedes `seam_coast` (kept functional).
4. **TPU test suite** — `test_models/tpu_suite/`: STLs, generator-derived README with
   the full test matrix and weight targets, flattened config overlays, the **sliced,
   print-ready `.gcode.3mf` files** (thumbnails injected, right nozzle, generic TPU
   @230C), and the tooling scripts under `tools/`. The sliced `.gcode.3mf` files are
   stored via **Git LFS** (`.gitattributes` tracks `sliced/*.gcode.3mf`): run
   `git lfs install` once on the new machine BEFORE cloning/pulling, or you'll get
   3-line pointer files instead of printable 3MFs (fix after the fact with
   `git lfs pull`). Future re-slices dropped in `sliced/` pick up LFS automatically.

Items 2-4 are committed together with this doc. All 27 `[ContinuousShell]` unit
tests and the full 156-case libslic3r suite pass at that commit.

## The physical experiment now printing

Two jobs in `test_models/tpu_suite/sliced/` (H2D, right nozzle, generic TPU @230C,
engineering plate, 0.45/0.25 flow — spacing 0.3963):

- **PlateA_T1-T4** — four tubes, by-object. T1 control (1.693 wall, tiles exactly),
  T2 medial-deficit (1.993 wall → innermost pair widened to ~0.61), T3 sub-bead
  (0.35 wall → fitted bead), T4 taper (2.60→1.50 wall, ring-count transitions).
- **T5_trio_fmb** — three 60 mm seam tubes (2.485 wall, exact tiling), front/middle/back
  = scarf 0/6/12 mm via per-object overrides. Verified: 0 / 9.74 / 18.70 m of tapered bead.

### Reading the results (decision table)

| Observation | Interpretation | Action |
|---|---|---|
| T2 weight ≈ T1×1.19 (both near predicted; see README table) | medial correction works | none — lock in |
| T2 light / void ring visible in section | deficit not closing on real hardware | check bisection resolution (4 iters, s/16) and first-layer interplay; add debug output |
| T3 near 1.29 g, smooth wall | width-fitting works | none |
| T3 rippled/gappy | fitted width too low for TPU die swell | raise the fallback clamp floor (currently `min_bead_width`, 85% nozzle) |
| T4 inner surface smooth at transition heights | graceful ring-count flips | none |
| T4 gap/bulge ring at a height | transition still abrupt | this is the C-full motivation — prioritize per-angle beading |
| T5 seam softens 0→6→12 under raking light | scarf works | pick a default (likely 6-8), consider making it default-on |
| T5 scarf regions gappy/stringy | floor too low for TPU | raise `scarf_floor` from nozzle/4 (one line in `build_cross_section`) |
| T5 scarf no visible difference | seam is dominated by the surface dive, not junction ridges | revisit surface-ring exclusion — maybe taper ring 0's tail after all, or move to C-full/B |
| weights uniformly off by a few % | flow calibration / density, not geometry | rescale expectations, compare ratios only |

Weight targets and per-print predictions are in `test_models/tpu_suite/README.md`
(predictions at 1.24 g/cm³; the T2/T1 ratio is density-independent).

> **2026-08-14 — prints done, results in `test_models/tpu_suite/RESULTS.md`.**
> C-lite validated on hardware (weights within 1%; T2/T1 = 1.188 vs 1.192;
> medial closed in section). Scarf (Option A) verdict NEGATIVE in TPU: seam
> creasing increases with scarf length, scarf-0 most isotropic — keep default
> 0, seam workstream moves to C-full/B. RESULTS.md also carries the 5-point
> C-full design brief (ring parity double-flip, squeeze-vs-stretch voids).

> **2026-08-15 — C-full implemented AND print-validated** as
> `continuous_shell_adaptive_beading` (region bool, **default ON** after the T4 A/B
> re-print: no visible transition band, stretched-regime medial void line eliminated —
> `test_models/tpu_suite/RESULTS.md`). Per-angle bead plan for dense annular star
> sections via the standalone Arachne stack; +10%-max-stretch squeeze-biased count law;
> single-step ring transitions via WIDTH-BIRTH (the medial ring thins to the 0.30 mm
> transitional floor across the flip band's ~20 layers instead of terminating — open
> coverage arcs were built and rejected: tips gap per T5, and a young arc's far tip costs a
> bore chord or a 9 mm layer join). All out-of-scope sections keep C-lite bit-identically
> (regression-tested). Design details in `continuous_shell_adaptive_width.md` §4.

> **2026-08-15 (later) — phase B implemented** as `continuous_shell_interior_spiral`
> (region bool, default OFF pending a print; requires adaptive beading + ≥4 buried
> rings). Full-span same-winding double-arm spiral in plan-slot space: k-bead no-stretch
> comb, tangent merges onto the surface rings, π-per-layer phase (two-start helix seam),
> junction census ~k → ≤1 per layer. Geometry settled numerically after hand derivation
> failed twice — see the design doc §4-B. The roadmap is now fully implemented.

## Next phases (agreed order: C-lite ✓ → A ✓ → C-full ✓ → B ✓ (unprinted))

- **C-full**: per-angle thickness field → per-angle bead plan via Arachne's standalone
  `BeadingStrategy` stack (pure/stateless — see design doc §3 for the exact reuse recipe
  and the ~150-line angular transition re-derivation). Inherits the star-shape gate.
  Print feedback from T4 should set the transition-length default.
- **B**: interior double-arm Fermat spiral; needs C-full to close the spiral-to-loop
  boundary wedge; requires deliberate pairwise ring morphing (see design doc §4).

Deferred/open items: weave arc-length reparametrization (star-gated shapes),
`build_cross_section` variant extraction, serpentine-cap last-row width quantization,
GUI toggle wiring for `wall_transition_*` options under continuous shell (C-full),
`pick_1.png` object-ID encoding in headless 3MFs is a best-guess (skip-object UI may
mispick on multi-object plates).

## Environment / workflow crib (self-contained — memory files do NOT transfer)

**Building** (original machine was Windows-on-ARM with VS 2026; a new machine needs
its own deps build — build trees are not in git):
- Deps: `cd deps && mkdir build-<arch> && cmake ../ -G "<your VS generator>" -A <arch>
  -DCMAKE_BUILD_TYPE=Release && cmake --build . --config Release --target deps -- -m`.
  Set `CMAKE_POLICY_VERSION_MINIMUM=3.5` for CMake ≥ 4. **OpenSSL needs Perl on PATH**
  (exit 9009 from dep_OpenSSL = perl missing; Strawberry Perl works). Match `-A` to the
  host arch — mixed-arch configures fail deep inside Boost.
- App/tests: configure the top level with `-DBUILD_TESTS=ON` (default OFF — without it
  no test project exists), build targets `libslic3r_tests` (fast iteration) and
  `OrcaSlicer` + `OrcaSlicer_app_gui` (produces `orca-slicer.exe` for CLI slicing).
- Run tests: `<build>/tests/libslic3r/Release/libslic3r_tests.exe "[ContinuousShell]"
  --order rand --warn NoAssertions`.

**Headless slicing → Bambu Connect** (all encoded in `test_models/tpu_suite/tools/`):
- `gen_tpu_suite.py` — parametric STL generator; edit W/H at the top for a different flow.
- `slice_tpu_suite.py` — flattens BBL profile `inherits` chains into self-contained
  JSONs and drives the CLI. Non-negotiable settings it encodes: `filament_map=["2"]` +
  `filament_map_mode="Manual"` (**CLI defaults to the LEFT nozzle**; David's hardened
  0.4 is on the RIGHT of the H2D); `print_sequence="by object"` (required by continuous
  shell for multi-object plates); generic-TPU temps overridden to 230C (stock profile's
  240C exceeds the spool's 235 max); `nozzle_volume_type=["Standard","Standard"]`
  (**Bambu Connect refuses files whose G-code CONFIG_BLOCK declares a nozzle flow
  vector it can't resolve against the printer** — the option's 1-element default
  emits `; nozzle_volume_type = Standard` for a two-extruder printer, leaving the
  right nozzle undefined. Connect reads the CONFIG_BLOCK comments, not the 3MF
  metadata. David's H2D hotends are both 0.4 hardened steel STANDARD flow. The
  committed 3MFs were normalized via `patch_nozzle_flow.py`, which fixes metadata +
  the comment line and recomputes the gcode md5).
- `make_t5_trio.py` — multi-object plate with per-object region overrides via
  `<metadata>` on `<object>` in model_settings.config. Never edit `3D/3dmodel.model`
  build items (required `p:UUID` attrs). 60 mm objects pass sequential clearance in a
  front/middle/back column at the arranger's 128 mm pitch (100 mm fails).
- `inject_thumbnails.py` — **run after EVERY headless slice**: CLI output has no
  thumbnail PNGs and Bambu Connect refuses to load such files. Injects
  plate/small/no_light/top/pick PNGs rendered from `Metadata/plate_1.json`.

The mode's error philosophy is fail-loudly: no silent fallbacks. Keep it that way.

## Work plan from the 2026-08-23 audit

Ordering rationale: fix what is wrong in **already-validated territory** first (density 0,
default ON), then the prerequisites the decoupling would otherwise expose, then the
decoupling itself. Each workstream is one commit with its own verification.

| # | Workstream | Scope | Why here | Verified by |
|---|---|---|---|---|
| **W0** | Golden `CsSig` harness | Byte-identity table at density 0 over annulus / eccentric / star tube / taper flip / sub-bead / solid / the four force_solid sections / spiral-on | Fence for every later step. Any diff is a bug, never a threshold to retune | Green before W1 |
| **W1** | `floor_transitional` re-derives `toolpath_locations` | `:599-616`, consumed `:626` | **Default ON at density 0** - the only item that touches print-validated output today | Golden table moves ONLY on transition layers; P2 print |
| **W2** | Weave clamps | `:1680-1705`: clamp to shell radius -/+ one spacing; build `inner_prof` for the solid (no-bore) case | Prerequisite: the decoupling newly exposes *planned* surface rings to a bead-on-bead collision | G-code: no kept rib within a spacing of a shell; P4 section |
| **W3** | Decoupling: `ContourChain::ring_slot` | Steps 1-8 of the brief; sparsifier extracted to a free function and called from the plan branch | The headline regression | Desync detector test; bit-identity at density 0; P3 |
| **W4** | `force_solid` keeps adaptive | drop `&& !force_solid` at `:1970`, keep spiral gated on it | Owns the printed doubled lip; density>0 only, so inert at density 0 | P3 cone lip vs the S5 in hand |
| **W5** | Surface the decline notes | 7 sites write `note`, nothing reads it (`:1489-1491`) | W3/W4 widen the set of silent fallbacks | A declining part logs a warning naming the layer |
| **W6** | `force_solid` -> solid *mask* | accumulate `diff_ex` over shell layers, offset by `wall_loops*spacing` | Biggest behaviour change; wants its own print | Bit-identical on a flat lid; P5 |

**Not scheduled:** the `SlotSchedule` single-source-of-truth refactor. It is the better end
state and both judges said so, but it rewrites the print-validated sparsifier and carries an
unprinted rib-placement law. Revisit after P4 says whether weave x per-angle width composes.

## Print campaign

Every print states what it must show. A print that cannot falsify anything is not run.

| Plate | Objects | Question | Baseline | Must show |
|---|---|---|---|---|
| **P0** `PlateD_data` (ready) | D1_ladder (wall 8/7/6/5/4 beads), D2_core | Where the sparsify gate really flips, and what a wall looks like either side | density-0 reslice of the same objects | C-lite baseline for P3. Gate fires at 6+ beads only; the 5-bead step is 4 stretched beads |
| **P1** `S5_conesparse` (in flight) | S5 geometry, `bottom_shell_layers=0` | Does sparsifying the cone wedge help or hurt the overhanging lip? | S5_sparse (printed) | Lip better (less over-extrusion feeding the droop) or worse (less backing). Either answer is data |
| **P2** after W1 | T4_taper A/B, density 0 | Does re-deriving bead locations close the transition-band skin deficit? | T4_adaptive (printed) | Flip bands lose their void rows; skins gain ~49 um where k=3 |
| **P3** after W3+W4 | S5 geometry + D1_ladder reprint, density 25 | **The load-bearing question: is C-full-under-sparse actually stronger than C-lite-under-sparse?** | S5_sparse and P0 (both printed) | Ladder 5-bead step 4.0 -> 5.0 beads/ray; base-wall void gone; cone lip undoubled; **tear test vs the C-lite coupon** |
| **P4** after W2 | Weaved sparse annulus, density 20-30, sectioned | Does a kept rib at its slot's planned width, displaced by the weave, still overlap its neighbours? | P0 ladder 8-bead step | Ribs overlapping, NOT the tangent-with-no-overlap signature the T4 microscopy found |
| **P5** after W6 | Capped cone + flat-lid part | Does a solid *mask* beat the whole-layer boolean without hurting the lid? | S5_conesparse (P1) | Lid unchanged; wedge sparse only where nothing needs support |

**P3 is the gate on the whole workstream.** RESULTS.md has no sparse + C-full entry, so the
claim that honest bead counts matter once the interior is decimated is inference. If the tear
test cannot separate C-full-under-sparse from C-lite-under-sparse, W3 is correctness work with
no strength dividend and the effort belongs elsewhere.

Two effects ride along on P3 and are unmeasured today: Fermat parity and `spread_twist` are
computed from POST-decimation indices (`:819-821, :1207-1220`), so a decimated plan chain gets
different radial-neighbour Z offsets than the same wall at density 0 - a different non-periodic
interdigitation, which is exactly what north star #4 buys. Tear both orientations.

### Measured baselines already in hand (2026-08-23)

- `PlateD_data` 12.00 g / 57 min vs 17.34 g / 79 min at density 0. Ladder beads per ray,
  dense -> sparse: 8 -> 4.8, 7 -> 4.8, 6 -> 4.6 (stride); **5 -> 4.0 with the gate OFF**, so
  that bead is lost to C-lite stretch, not to the sparsifier; 4 -> 4.0 unchanged.
- D2_core does NOT demonstrate the `:1704` clamp gap: a solid centre is filled either way
  (density 0 reaches r=0.2412 mm, density 25 r=0.2563 mm). The code gap is real; detecting it
  needs a ring self-crossing test, not an `rmin` probe.

### W1b (new, 2026-08-23): the medial deficit on walls of varying thickness

The marching-pinch fix (c17d5e36fd) restored the BULK material on eccentric bores - 69% ->
101% at a 2 mm offset - but measuring the worst ray rather than the median shows one
oversized gap left per ray, and it sits at the MEDIAL, where the two marching fronts meet:

| tube | wall thin/thick | material | worst medial gap |
|---|---|---|---|
| E1_ecc  (offset 2.0, pinches)        | 4.0 / 8.0 mm | 101% | 0.82 mm (~1 bead) |
| E2_mild (offset 0.5, does NOT pinch) | 5.5 / 6.5 mm |  94% | 1.29 mm (~2 beads) |

Every other gap on both is 0.40-0.47 mm. E2 measured 94% before the pinch fix as well, so
its void is pre-existing and untouched by it - the two are the same defect at different
severities, not a regression introduced by the fix.

Root cause, and why the obvious patch is wrong: C-lite's deficit correction is SCALAR. It
takes one `d_res` from `vanish_depth` and widens the innermost ring pair uniformly. On a
wall whose thickness varies with azimuth the residual varies too, so no single number is
right. Simply marking the pinched chain `vanished = 1` makes it worse, not better:
`vanish_depth` bisects until the whole shape vanishes, so on a pinched chain it returns the
THICK side's half-thickness and would widen the innermost beads by ~0.6 mm each.

The per-angle machinery that can express this already exists - it is C-full - and it
declines on exactly these walls:

    adaptive beading declined: wall thickness crosses more than one bead-count threshold

So the correct fix is to widen C-full's scope to walls spanning more than one bead-count
step, which means a per-angle bead COUNT rather than a single `n_slots` with per-angle
widths. The width-birth mechanism already prints a transitional ring floored at `w_term`
with the deficit redistributed, so the ingredients are there, but the count law and
`blend_ring_beadings` pair only k and k-1 today. That is a design-sized change and should
get the same treatment C-full itself got, not an improvisation.

`sliced/PlateE_eccentric.gcode.3mf` (19.25 g, 87 min) is the print that scopes it: E1 and
E2, 30 mm tall so both can be cut and torn. Two questions it answers that G-code cannot -
(1) does the recovered material BOND at the pinch boundary, where a separate chain is
joined at a place that used to be a void, or has a void been traded for a weak plane, and
(2) how bad is a 0.8-1.3 mm medial void ring in plastic, which sets the priority of the
C-full scope work against the rest of the queue.

**W1b sharpened by the PlateE print (2026-08-23).** Both predictions held, and the risk the
print existed to test did not materialise: E1 - the tube the marching fix took from 69% to
101%, where recovered material arrives as a separate chain joined where a void used to be -
shows near-perfect bricklaying with extremely few voids on either side of the cut. No weak
plane at the new interface. E2 shows its 1.29 mm medial void exactly where predicted, and
is otherwise essentially solid.

So the fix left the EXTREME case cleaner than the mild one, and the fence pins where the
behaviour changes:

| bore offset | chains | coverage | path |
|---|---|---|---|
| 0.5 mm | 1 | 0.906 | never pinches - 9.4% short, visible medial void |
| 1.0 mm | 4 | 1.059 | pinches, marching fix applies |
| 2.0 mm | 4 | 1.056 | pinches, marching fix applies |

A 0.5 mm bore offset on a 6 mm wall is an utterly ordinary geometry, so this is now the
worst case in the family. W1b's first step is therefore NOT to widen C-full's scope (which
would also work, and is far bigger): it is to instrument march_bands and find why a
5.5/6.5 mm wall does not take the pinch path at all. By hand the thin side should close at
k~6.9 and the thick side at k~8.2, which ought to produce a hole-less crescent and route
into the already-validated pinch branch. It does not, and the single-chain result says the
band keeps a hole until the whole thing vanishes - most likely a numerical sliver hole
surviving the offset. If that is the cause, the fix is to treat a collapsed hole as a
pinch, which routes this geometry through the mechanism the print just validated rather
than inventing a new one.

### W1b DONE for drifting walls (2026-08-23) - and what it changes downstream

The ring-ladder probe showed the marching path cannot be rescued. It lays rings at fixed
spacing from each face, so the leftover at the medial is a function of azimuth, and after
the pinch fix the gap simply MOVED: at a 2 mm bore offset the thick side became perfect
(18 rings, 0.448) and the THIN side inherited a 0.821 mm gap. No scalar correction can
express an azimuth-varying residual.

The fix is on the C-full side, and it turned out to be a pre-filter that was measuring the
wrong thing. `n_max - n_min > 1` rejected any wall whose per-bin counts differed by two -
but the count law assigns per-bin counts, so a wall that merely DRIFTS gets counts two
apart even where one count serves the whole section inside the same +-10% stretch budget.
The check now tests that directly (two probes, at the thinnest and thickest bins) and, when
one count fits, uses it uniformly - which runs the section down the well-tested
constant-count path instead of asking the k/k-1 blend to span two steps, which it cannot
represent. When no count fits it still DECLINES rather than hard-fails: refusing to slice a
part that used to slice is a worse failure than the underfill.

Measured, real slicer, tubes of identical volume:

| bore offset | before W1b | after |
|---|---|---|
| 0.5 mm | 94%, 1.29 mm void | **100%, worst gap 0.551 mm** |
| 1.0 - 3.0 mm | 101-102% | 101-102% |

Fence: only the varying-wall rows move. bore@+0.5 coverage 0.906 -> 1.002, bore@+1
1.059 -> 0.997, and the lobed star tube 0.931 -> 1.028 (its 7% deficit was the same
defect). Uniform walls are byte-identical. Widths stay inside the envelope: 0.300 (the
transitional floor) to 0.550 (exactly the eps = 0.10 cap).

**Downstream phase updates warranted by this:**

- **W3 (decoupling) is worth MORE than costed.** Its value was always "adaptive beading
  survives sparse infill", but the plan used to decline on any wall crossing more than one
  bead-count step - i.e. on most real geometry. It now engages far more widely, so W3
  carries the benefit to correspondingly more parts.
- **W2 (weave clamps) matters MORE, and is more clearly a prerequisite.** The decoupling
  now exposes *planned* surface rings to the bead-on-bead weave collision on many more
  cross-sections than when only near-uniform tubes took the plan.
- **W1b-v2 (per-angle bead COUNTS) drops in priority but stays open.** Walls where no
  single count fits still decline to marching. Coverage there is now 101-102% because the
  pinch fix carries them, but the fence still shows a 0.821 mm thin-side ladder gap at a
  2 mm offset. It is no longer a material-loss bug, only a local void.
- **NEW RISK to watch in print: width-run explosion.** Planning a lobed wall produces 334
  separate width runs in one layer (vs 1 on a uniform tube, 4 on a taper). Each is a 0.05 mm
  bucket step so the variation is smooth, and T4 validated width ladders - but 334 flow
  changes per layer is far beyond anything printed so far. The next print with a non-uniform
  wall must be checked for surface banding; the mitigation if it shows is coarser bucketing.
- **PlateD_data still needs its second object replaced** (D2_core cannot detect the :1704
  clamp gap) before it is worth filament. Fold into W2.
