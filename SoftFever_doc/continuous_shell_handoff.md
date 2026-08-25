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

### PlateF built (2026-08-23), and two corrections to the record

`sliced/PlateF_varying.gcode.3mf` (20.02 g, 90 min): E2_mild re-sliced from the SAME STL as
the printed part, plus a 5-lobe tube. Both at density 0, so W2 cannot invalidate this plate -
the weave clamps live inside the decimation block, which only runs at density > 0. That is
why this print goes BEFORE W2 rather than after.

**Correction 1 - the per-object splitter was wrong, and it poisoned the first report.** It
split a sequential plate on a Z drop scanned from `G1` lines only; the inter-object move on
this plate lands on a `G0`, so both objects merged into one and every per-object number was
computed over the pair. It reported 791 width runs and a 0.710 mm gap, neither of which was a
property of any object. All three plate tools now split on the explicit `; printing object
NAME` markers Orca emits, which also name the object instead of inferring it from height.

**Correction 2 - reporting a maximum conflated a seam with a void ring.** One ray carrying a
big gap is a seam artefact; a void ring shows as most rays carrying one. Measured on E2_mild
at z=15, same STL, before and after W1b/W1:

| | median | p90 | max | rays with gap > 0.55 mm |
|---|---|---|---|---|
| before (PlateE, printed) | 0.811 | 1.271 | 1.631 | **48 / 72** |
| after (PlateF) | 0.399 | 0.435 | 0.690 | **1 / 72** |

So the void ring IS gone and the typical gap is now one bead spacing. Quoting the max alone
made the fix look half as effective as it is - a mistake in the opposite direction from the
single-layer 0.551 mm figure quoted earlier, which was equally unrepresentative. The tools
now report the whole distribution.

**Risk downgraded: the width-run explosion is not one at this profile.** The fence measured
334 runs in a layer on a lobed wall, and the first PlateF report said 791 - but that was the
merged-object bug. Measured per object at the campaign's flow: 19 runs per layer on the lobed
tube, 9 on the eccentric. The fence's 334 came from a different flow (0.5 mm nozzle) and
geometry. Still worth a look at the printed lobe flanks, but this is no longer a reason to
hold W3, and coarser bucketing is not needed on this evidence.

One thing to check on the print: F2_lobed carries beads down to 0.246 mm, below the 0.30 mm
transitional floor, on segments long enough not to be corner noise. Look at the lobe flanks
for thin or broken beads.

### PlateG (2026-08-24): the strength attribution control

`sliced/PlateG_control.gcode.3mf` (20.02 g, 90 min), density 0 - so W2 can proceed in
parallel without invalidating it.

- **G1_conc**, a plain concentric tube, is geometry that NONE of the wall-filling work
  touches: W1b only widens the plan's scope on walls of VARYING thickness, and W1 plus the
  floor guard only fire where a transitional bead is floored, which a uniform wall of this
  thickness never does. The density-0 fence confirms it - the concentric rows are
  byte-identical across every one of those commits, and this slice reads a flat 0.398 mm
  median gap, 0/72 rays over 0.55, no sub-floor extrusion, widths 0.394-0.400.

  It exists because PlateF came off "extremely difficult to cut even with a scalpel" and that
  result is UNATTRIBUTED. If this control cuts like the older campaign parts, the gain belongs
  to the wall-filling work. If it resists equally, the cause is something this campaign did
  not change and the strength claim has to be withdrawn.

- **G2_lobed** is PlateF's lobed STL re-sliced with the floor-affordability guard. Expect the
  void column REDUCED, not absent: sub-floor extrusion is 6.75% of the layer, one transitional
  ring's worth, and that residue belongs to the count law rather than the rebate.

Note: `PlateF_varying.gcode.3mf` in the repo is now the POST-guard slice, because it was
re-sliced to verify the guard. It no longer matches the printed PlateF part. PlateG supersedes
it as the canonical artifact for this geometry.

### W2 done (2026-08-24): weave clamps

Two audit findings, both density>0 only, both in the weave perturbation block:

- **#7** the clamps bounded infill points by the SHELL RING'S OWN RADIUS, so a lobe could
  close the pitch to a shell to zero and lay bead on bead. The amplitude law sizes lobes for
  infill-to-infill neighbours K*spacing apart, but the first kept infill ring's outward
  neighbour is an immovable shell one spacing away. Now clamped to one spacing clear of both
  shells, with a band narrower than two spacings left unwoven rather than forced.
- **#5** a SOLID (no-bore) chain had no inward clamp at all - `if (annular)` guarded it - so
  the innermost kept ring was bounded only by `max(rr, 1.0)`, a nanometre from the centroid.
  It now profiles the terminus ring; a degenerate terminus reads zero, which floors the bound
  at one spacing from the centre, also correct.

**Verification, stated honestly.** Full suite green (171 cases) and the density-0 fence is
byte-identical on all 20 rows, so the change is confined to density>0 as intended. Indirect
but attributable evidence that #7 bites: on the D1_ladder at density 25, with nothing else
changed, sparse beads per ray rose 4.8 -> 5.3 (8-bead wall), 4.8 -> 5.1 (7) and 4.6 -> 5.2
(6) - rings that used to be clamped onto a shell's centreline now resolve as separate
crossings a spacing away.

**What could NOT be verified from G-code, and why.** Four probes for #5 were built and all
failed for the same principled reason: a solid chain legitimately runs its terminus ring to
the centre, so neither the innermost extrusion radius (0.2412 mm at density 0 vs 0.2563 at
25 while the bug was live) nor the smallest adjacent-ring gap can separate a collapsed ring
from normal convergence. Nor can a shell-to-infill pitch be isolated by ray crossings, because
the weave deliberately brings anti-phase infill neighbours' bead EDGES into contact, so small
gaps are the design working. #5's fix rests on code reasoning plus containment, not on a
measurement - the honest status.

**PlateD_data's dud object is replaced.** D2_core (solid cylinder) becomes D2_weave, a thick
annulus (r_out 10 / bore 2.5, 7.5 mm wall). It has real boundaries at both faces, exercises
the full woven lattice, and is the geometry where a bad clamp shows physically as
over-extrusion against a shell. That is the W2 print question, and it is a question for the
printer rather than for a G-code script. Plate is 11.93 g / 56 min.

### W3 DONE (2026-08-24): adaptive beading and sparse infill compose

The headline regression is closed. `ContourChain` gains `ring_slot`, the plan slot each loop
came from. Position in `loops` used to BE the slot, which is exactly why the two could not
coexist: decimation drops rings, so every surviving ring past the first gap would have read a
neighbour's per-angle widths.

The design's load-bearing claim was verified before building on it: in ring mode the plan is
read by loop index in exactly ONE place, `point_w` (`:1003`), called from the ring emitter and
the outer-ring emitter. The spiral branch takes explicit slots and never calls it; `width_of`
at the scarf uses `ring_widths`, not the plan. So routing `point_w` through the map re-routes
every ring-mode width, and phase B is untouchable by the change.

Shape of it:

- the sparsifier is extracted to `apply_woven_infill(cc, woven, spacing, layer_index, note)`
  so the PLAN path can call it too - it previously returned its chain before ever reaching
  that code;
- decimation compacts `ring_slot` and `ring_widths` in lockstep with `loops`;
- `chain_thick_enough` / `chain_will_decimate` replace two separate expressions that had to
  agree and could drift;
- `validate_chain` checks the map is present, strictly increasing, in range, and still holds
  both pinned surfaces - `ChainBeadPlan::interp` indexes its tables unchecked, so a drifted
  map is an out-of-bounds read, not merely a wrong bead;
- the solid-cap brick stagger is skipped on planned chains: it is a rigid half-spacing
  translation of the interior block and the plan owns radial placement;
- spiral x sparse is now REJECTED rather than silently ignored (a spiral chain carries only
  its two pinned surfaces, and the sparsifier needs more than two rings, so sparse infill was
  a no-op on spiral layers). Mirrored in `Print::validate` so it is caught before slicing
  rather than at layer 60 of a long job.

**Verified.** Full suite green (172 cases) and the density-0 fence byte-identical on all 20
rows - the change is confined to density>0. End to end, the pre-registered prediction on
D1_ladder is met exactly: the 5-bead step reads **5.0 beads per ray, up from 4.0**, while the
4-bead step stays 4.0 and the 6/7/8 steps stay sparsified.

One existing test was INVERTED, which was expected: "woven sparse annulus keeps its lattice"
asserted that at density 30 the adaptive option changed nothing. That was the defect written
down as a requirement. It is replaced by a test that the plan engages, that decimation still
happens, and by the DESYNC DETECTOR - on a varying wall with the weave off, decimation never
moves a point, so any surviving point must keep exactly the width it had when nothing was
decimated. A slot desync corrupts a whole ring (hundreds of points) and cannot hide.

Building that detector took three passes and both dead ends are worth knowing: coordinates
shared by two rings cannot say which slot a later run meant (49 false hits), and the point ON
a width-run boundary belongs to both ExtrusionPaths and legitimately carries either width
(2 more). Excluding both leaves zero.

### W7 (2026-08-24): the weave's lateral tie is geometry-dependent

Printed section of a thick woven annulus showed concentric wavy cylinders joined at points -
top-down support, almost no lateral tie. Measured per adjacent ring pair: 13.0% of the
circumference with bead edges touching, 5.1% with real overlap, one mid pair at 0.1%.

**Cause (arithmetic, not a bug).** The drive amplitude equalled the contact bound, so the
closest approach is exactly one bead spacing at the instant the sine peaks - peak meets valley
at a tangent POINT. Contact needs `sin >= (K*spacing - width)/(2*amp)` ~ 0.98, about 6% of
every period, and it is scale-free: sweeping the lobe count over its whole range moved bonding
only 13.5% -> 9.6% (finer lobes are WORSE - they sweep through the approach faster). Neither
the amplitude nor the bound has changed since the original commit, so it never was a
regression, and no user setting escapes it.

**Fix (fed6bae35d):** the drive now overshoots the bound 1.6x and is CLIPPED to it, so each
meeting is a parallel run at one spacing rather than a graze. Contact holds while
`|sin| >= 1/overshoot`, predicting ~28% of each period.

**Why the first measurement missed that prediction, and what it means.** On the thick annulus
the change gave only 14.0% touching, against 28% predicted. The annulus's 7.5 mm wall holds
about FOUR kept interior rings, so most of them are shell-adjacent - and W2's shell clamp caps
their outward swing, which is exactly what suppresses the weave. On a SOLID disc (no bore,
r 15 mm, same settings) where shell-adjacent rings are a small minority, bonding is 26.7%,
matching the prediction, and it rises toward the centre:

| radius band | edges touching |
|---|---|
| 3-6 mm | 30.9% |
| 6-9 mm | 30.3% |
| 9-12 mm | 23.2% |
| 12-15 mm | 25.5% |

That radial trend is the remaining design question. Lobes are a fixed ANGULAR count, so the
cell aspect ratio - arc wavelength against ring pitch - drifts with radius: equiaxed cells and
a real net near the centre, long thin cells and tube-like behaviour far out. Making the lobe
count scale with radius (constant arc wavelength) would hold the cell shape across the part.
Each ring must close, so the count has to stay integral per ring, which means adjacent rings
would carry different counts and meet on a beat rather than everywhere in phase - worth
designing deliberately rather than improvising.

**Practical guidance meanwhile:** the woven infill earns its lattice on LARGE, SOLID sections.
On a thin or few-ring interior most rings are shell-adjacent and it degenerates toward
concentric tubes; there, higher density (more kept rings) helps more than any weave setting.

### W7 step 2 (DESIGNED, not implemented): dynamic density floor

David's idea, and it targets the failure mode exactly. A lattice needs enough buried rings to
BE a lattice: with too few, most kept rings are shell-adjacent, W2's shell clamp caps their
swing, and the weave degenerates into concentric tubes. So set the stride from the ring budget,
not from density alone:

    n_buried = rings between the shells
    K_eff    = clamp(n_buried / MIN_LATTICE_RINGS, 2, K_from_density)

The shape is right because it is FREE where the lattice already works. The solid disc has ~35
buried rings, so a floor of 6 needs K <= 5 and it keeps its requested K = 4 - no change, no
extra material. The thin annulus has 15 and densifies. Material is spent only where the lattice
would otherwise fail to form, which is also where there is least area to spend it on.

Measured on the 7.5 mm-wall annulus (tie fraction includes ties to the shells):

| floor | K | material | tie |
|---|---|---|---|
| none (today) | 4 | 3.94 g | 16.2% |
| 4 rings | 3 | 4.25 g (+8%) | ~19% |
| 6 rings | 2 | 5.10 g (+29%) | 22.3% |
| *solid disc reference* | *4* | *unchanged* | *26.7%* |

Suggested default: a floor of 4 - most of the benefit for a third of the cost - with the floor
exposed rather than hard-coded.

**Deliberately not implemented yet.** It overrides the density the user asked for (request 25%
on a thin wall, get 33-50% there), and the open question is whether 6 points of tie fraction is
visible in plastic at all. This whole W7 thread has been wrong-until-corrected-by-a-print, so
the density comparison goes on the bench BEFORE a policy is built around the number. When it is
implemented it must report through the W5 degrade_note ("interior only 15 rings, densified
stride 4 -> 3 so the weave can tie") rather than silently inflating a part.
