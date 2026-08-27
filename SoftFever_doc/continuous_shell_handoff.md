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

### W7 step 2: the ring-count floor is NOT worth building (2026-08-25)

Costed before implementing, and the numbers killed it. Kept infill rings m against a floor:

| geometry | buried span | natural m | floor 2 | floor 4 | floor 6 |
|---|---|---|---|---|---|
| D2_weave, 7.5 mm wall | 15 | 4 | 4 | **4** | 6 |
| D3_dense, same at density 50 | 15 | 8 | 8 | 8 | 8 |
| ladder 8-bead step | 4 | 1 | 2 | 4 (= solid) | 4 |
| solid disc r15 | 34 | 9 | 9 | 9 | 9 |

**A floor would not have changed D2 at all.** D3 was preferred because it is denser everywhere
(m = 8 against 4), not because a floor fired. The floor's real reach is a thin interior whose
natural stride leaves ONE lonely ring - and that case is already handled by the floor of two
that came with the bore-gap fix (7119045180). Floor 4 mostly turns 3-4 bead interiors solid;
floor 6 does change D2, but that is just "use a higher density", which is an existing knob.

So the mechanism is done at m >= 2 and the lever for thick interiors is density. Recorded
rather than built, to avoid machinery that does not do the thing it was wanted for.

### PlateL: the de-lid plate (2026-08-25)

`sliced/PlateL_delid.gcode.3mf` (5.18 g, 28 min). Everything the woven infill has been judged
on so far shows the wrong thing: a section shows in-plane structure, a top view shows one
layer. Neither shows whether consecutive layers LOCK - which is what the Z-wave and the
per-layer phase advance exist to do, and the one claim in this subsystem with no physical
evidence behind it.

`top_shell_layers = 0` leaves the lattice open from above, so several layers are visible
through each other. Two objects, same annulus (r_out 10 / bore 2.5, 12 mm):

    L1_open25   density 25, stride 4
    L2_open50   density 50, stride 2, via a per-object override

Density is the comparison because density is the axis known to matter (PlateD: denser
preferred by eye), and because the floor turns out not to be.

Checked before publishing: the top layer runs 0.62-0.66x the length of a mid layer, so there
really is no solid skin - a lid would run several times longer.

**Also fixed a bug in the publish guard itself.** Builders import each other for shared
helpers, and an imported module runs its own `resolve_outdir` at import time, which CONSUMED
`--publish` from argv before the importing tool's call saw it. `--publish` therefore produced
a silent verify-only run - the very failure the guard exists to prevent, one level up. The
decision is now memoized.

### PlateL reprint (2026-08-25): the amplitude law was wrong, not the phase

The staggered-tie plate came off "better, not quite perfect": the dense object was nearly a
foam, missing only a few connections between adjacent columns; the sparse object still had
obvious gaps and its columns were **not connected at all**.

The right statistic turned out not to be the one used to justify the drift. Tie fraction is a
per-PAIR measure - does ring j touch ring j+1 anywhere - and it says nothing about whether any
given azimuth is braced. The wall reads as columns when an azimuth has **no weld anywhere
across its depth**, so that is what is measured now (`weldmap.py`: welds per azimuth per layer,
plus the azimuthal width of the dead bands and how many layers they stack).

**What the drift actually did**, measured on the two printed files:

| | mean welds/azimuth | azimuths with none | dead band | tallest dead stack |
|---|---|---|---|---|
| print 1, L1 sparse | 0.85 | 61.4% | 11.9 deg | 5 layers |
| print 2, L1 sparse | 0.80 | **40.0%** | **4.0 deg** | 8 layers |
| print 1, L2 dense | 1.18 | 58.0% | 17.9 deg | 5 layers |
| print 2, L2 dense | 1.03 | **43.4%** | **7.3 deg** | 4 layers |

So the drift was a real improvement - it broke the wide aligned gaps into narrow scattered ones
- and it slightly reduced weld depth doing it. Both match the part in hand. But 40% of azimuths
still carried nothing, and a sweep of saturation overshoot x drift **through the slicer** (16
slices, scored on the same weld map) hit a hard floor near 37%: no combination of the two knobs
did better. That is the signal that the knobs were not the problem.

### The actual defect

Slicing with the weave OFF shows the kept rings where the sparsifier puts them. On L1:

```
r:    2.71 3.11 | 3.51        5.52   6.96        8.97 | 9.37 9.76
gap:   0.40 0.40  2.00        1.45   2.00        0.40   0.40      (mm)
```

The gaps are **5, 4, 5 slots** - and the amplitude was a single scalar sized from the nominal
stride K=4. Every gap therefore closed by exactly `2*cap` whatever its width: a 4-slot gap
closed to one spacing (a weld), a 5-slot gap closed to two, which is air. **No setting reaches
it** - the amplitude that would close the 5-slot gap crosses at the 4-slot one, and crossing
retraces ground already laid in the layer. The printed result is exactly what the geometry
predicts: rings welded into doublets that touch neither the shells nor each other.

This is not peculiar to this part. Sweeping thickness and density shows walls whose kept-ring
gaps exceed K are the common case, not the exception (d=25 t=9mm: gaps 5,5,5; d=50 t=9mm: a
3-slot gap where the stride is 2).

### Fix: size the amplitude per gap, per angle

Measure the gap that actually exists between each ring and each of its neighbours, at every
angle, and size the excursion to close *that*. Two weaving neighbours approach each other, so
each covers half the excess over one spacing; an immovable shell covers none, so its neighbour
covers all of it. Both halves come from the same measured gap, so opposed rings meet at exactly
one spacing by construction and cannot cross. Per-angle because adaptive beading makes the gaps
vary around the ring.

The two wave constants were then re-swept against the new geometry and set to overshoot 4.0 and
a drift of 2 slot-widths - chosen because that pair minimises the WORST of the two objects, and
the sparse one is the failure mode.

| | mean welds/azimuth | azimuths with none | dead band | tallest dead stack |
|---|---|---|---|---|
| L1 sparse | 0.80 -> **1.17** | 40.0% -> **28.3%** | 4.0 -> **2.7 deg** | 8 -> 5 layers |
| L2 dense | 1.03 -> **1.34** | 43.4% -> **25.7%** | 7.3 -> **3.2 deg** | 4 -> 5 layers |

Cost: 5.19 g -> 5.34 g (+2.9%), 28 -> 29 min. Density-0 fence byte-identical; full suite green
(173 cases).

### What is deliberately NOT claimed

- **No azimuth is welded all the way through, and none can be.** Anti-phase neighbours mean that
  where one pair closes, the next is at its widest. The wall is a truss of alternating doublets,
  not a solid; connectivity comes from the alternation around the circumference, not from
  simultaneous welds. Measured: 0.00% of 9000 rays fully welded, before and after. Any future
  metric that asks for a full radial bridge is asking for the wrong thing.
- **Welds are a G-code measure** - beads whose centrelines come within a bead width. Whether they
  bond in TPU still needs the part.
- **The two constants are tuned on one plate at two densities.** Treat them as the current best
  estimate, not as settled.

### Regression test, and the trap it nearly fell into

`Every woven neighbour pair reaches a weld, whatever its gap` asserts the per-pair property on a
wall that mixes 4- and 5-slot gaps, and asserts that mix first so it cannot go vacuous. It took
four attempts to make it honest, and the first version **passed against the bug** - the wall it
used happened to have a single buried gap of exactly K slots, which the old scalar law serves
perfectly. Two further versions failed for measurement reasons worth recording:

- a wedge sampler splits one saturated ring into two radius clusters and renumbers every ring
  above it;
- matching rings by nearest nominal radius breaks **exactly at a weld**, where both rings have
  moved to within one spacing of the midpoint and each is nearest the other's nominal.

Rings are now identified by the nominal band they lie in, which is sound only because every
buried pair on that wall is infill-to-infill. Verified in both directions: the test fails
(2.07 slots vs 1.10 allowed) when the per-gap law is swapped back for the scalar.

### Open, and worth a decision

The sparsifier anchors a ring on each shell and shares the remainder through the middle, so m
rings make m-1 gaps across span-1 slots and the **delivered stride is (span-1)/(m-1), not K**.
On L1 that is 4.67 against a requested 4 - about 17% sparser than the density asked for. The
per-gap amplitude now welds those wider gaps regardless, so this is no longer a strength bug,
but it does mean the density setting is not being honoured. Correcting it costs roughly one
extra ring per chain (more material) and changes default output for every sparse part, so it
should be a deliberate call rather than a quiet fix.

**Next:** reprint PlateL to see whether the columns actually merge in plastic. Then W4
(`force_solid` layers keeping adaptive beading - the printed doubled bore lip), then W6.

### The peel test (2026-08-25): the anchor ring is the peel plane

With the bottom shell removed, the dense de-lid object **peels apart into three concentric
sections**. That is a different kind of evidence from everything before it: a peel finds
separation surfaces at particular RADII, and every metric in this campaign so far has been
per-azimuth. The weld map cannot see a cylinder of poor bond, because such a cylinder is equally
unbraced at every azimuth.

Three measurements, in order:

1. **It is not a void.** Radial coverage (what fraction of azimuths and heights carry material at
   each radius) never falls below ~32% anywhere in the wall. There is no empty cylinder. So the
   separation is beads that touch without bonding.
2. **The infill-to-shell interface is an order of magnitude weaker than the interior.** In
   contact at 4.5% of samples against 32% for infill-to-infill, and that ratio held at every
   Z-wave amplitude tested (0, 0.20, 0.45).
3. **The cause: the anchor ring does not stay anchored.** The sparsifier deliberately keeps a
   ring adjacent to each shell so the lattice ties to both faces - but that ring is woven like
   any other, and its swing toward its inner neighbour carries it off the shell.

| L2 dense, anchor to shell | bonded | median gap | p90 gap |
|---|---|---|---|
| bore face | **43.8%** | 0.454 mm | 0.603 mm |
| outer face | **42.6%** | 0.450 mm | 0.603 mm |

One and a half spacings at the 90th percentile, at a FIXED radius, repeated every layer. That is
a peel surface, and there is exactly one at each face - hence three sections.

### Two corrections to earlier records

- **The weld map takes no account of Z.** Two beads binned into the same layer can be most of a
  layer apart vertically and still be counted as welded. Measured dZ at contact is 0.35-0.64 of
  a layer. Every weld-map figure quoted before this therefore OVERSTATES bonding. The
  version-to-version comparisons remain valid - zamp was constant across all of them - but the
  absolute percentages are not a measure of fused material.
- **The Z stagger between radially adjacent rings is inherent, not a bug.** At zamp = 0 the
  interior still shows 0.36-0.50 of a layer between neighbours, because a CONTINUOUS shell ramps
  in Z as it spirals inward. Radially adjacent beads never meet flat. This is the same effect the
  campaign has been calling bricklaying; it cannot be removed without giving up continuity, and
  a quadrature Z wave cannot fix it either - the simulation gives ~0.28 of a layer at contact
  even with zero per-ring drift, because W7 deliberately widened contact from a point to a run
  and the Z wave keeps sweeping throughout that run.

### The fix, designed and measured but NOT landed

Pin the anchor ring: a kept ring that neighbours a shell does not weave at all. It stays
concentric, welded to that shell the whole way round, and its inner neighbour covers the whole of
their gap rather than half - exactly as if the pinned ring were itself a shell. Also removes an
old `min(c, 0.49*gap)` clamp in `gap_cap`, which was not a safety net but a leak: `share*(gap -
spacing)` already lands a pair one spacing apart at full excursion, and the 0.49 clamp capped the
excursion below what contact needs whenever a share-1.0 gap exceeded about two slots.

Patch kept at `pin_anchor_rings.patch` (54 lines). Measured:

- **The radial peel plane is gone.** Anchor-to-shell gap becomes 0.400 mm median with a p90 of
  0.401 - one spacing, constant, all the way round, both faces, both densities.
- **But bonded only reaches 40-61%**, because the pinned ring is still tagged as infill and so
  still carries the Z wave, which lifts it off the shell vertically. The radial peel plane has
  been traded for a vertical one.
- **And the sparse object regresses badly**: 28.3% -> 48.9% of azimuths with no weld, welds per
  azimuth 1.17 -> 0.70. Pinning two of L1's four infill rings leaves only two that weave, across
  gaps of 5 slots.
- The weld regression test fails on its sampling guard, correctly: pinning makes the next ring
  cover a whole gap, so it crosses the nominal midpoint and breaks the band-assignment premise
  the test documents as its own precondition.

### What has to be settled before this lands

1. **Exclude the pinned ring from the Z wave.** It needs to be tagged as shell rather than
   infill, which is the `chain.woven_infill ? woven.shell_loops : 0` argument at the
   spiralize_rings call. Tagging it shell also changes its extrusion role and the over-extrusion
   cap, which is arguably correct - it IS functionally part of the shell now - but is a wider
   change than the Z wave alone and needs the fence and a surface check.
2. **The sparse regression is really the stride bug.** L1's kept rings sit 5 slots apart because
   the delivered stride is (span-1)/(m-1), not K. At that spacing there is not enough material to
   build a lattice once two rings are pinned. Fixing the placement so the stride matches the
   requested density costs about one extra ring per chain and changes default output for every
   sparse part - the deliberate call already flagged.
3. Re-derive the weld test's ring identity, which pinning invalidates on that wall.

### Pinning the anchor ring, and honouring the density (2026-08-25)

Both changes from the peel diagnosis are now in, together, because neither stands alone.

**Pin the anchor ring.** A kept ring that neighbours a shell no longer weaves. It stays
concentric, welded to that shell the whole way round, and its inner neighbour covers the whole of
their gap rather than half - as if the pinned ring were itself a shell. The `min(c, 0.49*gap)`
clamp in `gap_cap` goes with it: `share*(gap - spacing)` already lands a pair exactly one spacing
apart at full excursion, so the clamp was not a safety net but a leak, capping the excursion
below what contact needs whenever a share-1.0 gap exceeded about two slots.

**Honour the requested density.** The ring count was `(span + K/2)/K`, which counts rings as if
they were unanchored. Both ends ARE anchored, so m rings make m-1 gaps across span-1 slots and
the delivered stride is `(span-1)/(m-1)`. The count is now chosen by searching the span for the m
whose delivered stride is closest to K. Rings sit on integer slots, so an exact K is not always
reachable - on the sparse de-lid object the best available moves 4.59 to 3.44 against a requested
4.0 - but the dense object was already optimal at 1.97 against 2.0 and is untouched.

### Results

| anchor to shell | printed | now |
|---|---|---|
| L2 bore | 54% in contact, p90 gap 0.603 mm | **100%**, p90 0.401 mm |
| L2 outer | 51%, p90 0.603 mm | **100%**, p90 0.401 mm |
| L1 bore | 60%, p90 1.203 mm | **100%**, p90 0.401 mm |
| L1 outer | 92%, p90 0.456 mm | **100%**, p90 0.401 mm |

One spacing, constant, all the way round, both faces, both densities: the peel surface is gone
from the geometry. Welds per azimuth rise from 3.68 to 4.76 (sparse) and 4.59 to 5.83 (dense).
5.34 g -> 5.43 g (+1.7%), 29 min. Density-0 fence byte-identical; full suite green at 173 cases.

### A metric that was lying, and how it was caught

The earlier weld map counted a pair as welded when its centreline gap was at most the mean bead
width. Beads sit one spacing apart when adjacent - 0.3963 mm against a 0.40 mm width - so the
threshold sat exactly where the data lands, and which side a pair fell on was numerical noise.
Read at that threshold, pinning looked like a REGRESSION on the sparse object (28.3% -> 48.9% of
azimuths unwelded). Read at 1.15x, every variant has 0.0% unwelded and pinning is the clear
winner. The ordering at 1.00x was an artifact.

Consequences, stated plainly:

- **Any "% of azimuths with no weld" figure in this file is threshold-dominated** and should not
  be quoted. The stable statistic is welds per azimuth at a tolerant threshold.
- The per-gap amplitude work that preceded this still stands, but on its own evidence - it took
  a pair that could never meet at any setting from 2.07 slots to 1.0, which is a large,
  threshold-insensitive change - not on the dead-azimuth percentages quoted alongside it.
- A knife-edge threshold is not a small inaccuracy. It inverted a conclusion.

### Tests

`Every woven neighbour pair reaches a weld, whatever its gap` gained a second claim: the
shell-side tie must hold ALL THE WAY ROUND, asserted as the fraction of rays within contact
rather than the worst case, since one stray crossing off a junction hop should not decide it.
Verified to fail with the pin disabled - tied on 49.8% of rays, opening to 2.74 slots, which
matches the 43% measured on the printed part.

Ring identity in that test also had to change again, to RANK on rays showing the full set (715 of
720). Pinning is what broke the previous scheme: a ring facing a ring that cannot move covers the
whole of their gap, so for a two-slot gap it finishes exactly on the midpoint and lands in its
neighbour's band. That is the fourth identity scheme this one test has needed; the three that
failed are recorded above and in the test's own comments.

### Still open

- The pinned ring is still tagged infill, so it still carries the Z wave and is lifted off the
  shell vertically even though it is now radially perfect. Tagging it shell (the
  `chain.woven_infill ? woven.shell_loops : 0` argument at the spiralize_rings call) would fix
  that, but also changes its extrusion role and the over-extrusion cap - a wider change needing
  the fence and a surface check.
- **Next print: peel the reprinted part, both objects.** The geometry no longer has a peel
  surface; whether the plastic agrees is the question, and the sparse object should be peeled
  too, which the last round did not do.

### The peel again (2026-08-25): the pin worked and the peel did not move

Reprinted with the pin and the stride fix, both objects peel **virtually the same as before**.
The user's read of the dense one: the middle column of foam is somewhat joined to the innermost,
still entirely disconnected from the outside, and the "somewhat" is a very weak joint.

**The pin did what it claimed.** Tie fraction by interface, measured on the printed part,
running bore-outward (rank matching on rays showing the full ring set):

| L2 dense | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| tied | 99.4% | 99.7% | 73.8% | **0.8%** | 15.1% | 63.0% | 13.6% | 43.7% | 21.4% | 99.8% | 99.9% |

The shell interfaces (0, 1, 9, 10) are now tied 99.4-100%, against 43% before. The peel plane did
not survive there. It moved INTO the lattice, where interface 3 is tied over 0.8% of its
circumference - a very nearly free cylinder - and 4, 6, 8 are tied 13-21%.

Note the shape of that list: it alternates. High, low, high, low. That is not noise.

**The mechanism: the lattice sets into permanent doublets.** Adjacent rings are meant to be in
anti-phase, so each interface should close wherever its two rings oppose - about 42% of the
circumference at overshoot 4, and every interface equally, because a ring spends half its period
facing each way. Measuring the correlation of adjacent rings' radial deviations, which should be
about -1 for anti-phase:

| iface | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| correlation | +0.00 | **+0.25** | -0.53 | -0.93 | -0.68 | -0.60 | -0.02 |
| tied | 74% | 1.6% | 13% | 64% | 14% | 46% | 22% |

Rings 3 and 4 move TOGETHER. Only one pair in the wall (-0.93) is properly opposed. Each ring
commits to one neighbour and stays committed all the way round, so the wall is a stack of bonded
doublets with unbonded cylinders between them - which is precisely what a peel finds, and
precisely what the hand reports.

**The drift is not the cause.** It weakens opposition to cos(d) - 30% lost on the dense object -
and that was worth testing, but setting it to zero moves the worst interface only from 0.8% to
6.6% and lowers the mean (57.3% -> 54.6%). Both are catastrophic for peel. Default left at 2
slots; changing it would be noise-chasing.

### What this means

Every fix so far has moved the weak surface rather than removed it. The scalar amplitude left
wide gaps unclosable; per-gap sizing fixed that and the weak surface moved to the shell face; the
pin fixed the shell face and the weak surface moved one ring inward. That is not bad luck - it
follows from the geometry:

**A lattice of concentric rings that tie by tangency has no radial members.** Its only radial
connections are the places where two rings happen to touch, so between them there is always a
continuous cylindrical void, and a peel needs exactly one such surface. Making individual ties
better cannot remove the class of failure; it relocates it to whichever interface currently ties
least.

The consequence for the roadmap is that peel resistance is not reachable by tuning the weave.
What would reach it is a generator whose path traverses the wall RADIALLY - a serpentine that
sweeps across the thickness as it goes round, so every cell is bounded by continuous material in
both directions - rather than concentric rings that undulate. That is a new interior mode, not a
parameter, and it should be weighed against what the ring lattice is good at before it is built.

### Not to be repeated

Three metrics in a row have flattered this subsystem: tie fraction (per-pair, blind to whether
any azimuth is braced), welds per azimuth (blind to radius, so blind to peel surfaces), and both
at a knife-edge threshold. The statistic that matches the hand is **tie fraction per interface,
in radial order, worst-first**. A mean over interfaces hides a free cylinder; a peel does not.

### Radial-traverse interior (2026-08-25): the generator that has radial members

The ring lattice cannot be made peel-resistant, because its only radial connections are where two
rings happen to touch and a continuous cylindrical void always sits between them. The serpentine
replaces the buried interior with one path that sweeps ACROSS the wall as it goes round, so every
traverse crosses every ring and a crack at any depth must sever traverses rather than separate.

Geometric prototype first, before any C++, on the de-lid geometry at density 25:

| | ring lattice | serpentine |
|---|---|---|
| members a crack must sever, WORST depth | 0 to 1 | **24** |
| ...at every other depth | 0 to 74, wildly uneven | 24, uniform |
| path length for the same fill | 1.00x | 1.01x |

min equals mean by construction, and it costs nothing: both lay the same area of bead. The
density law falls out clean - **N = density x perimeter / (2 x spacing)** - and the prototype
confirmed it (12.3 predicted, 12 delivering the requested 25%).

### Built in ring-index space, not in angle

The first cut parameterised the traverse by angle about a centroid, which would have confined it
to star-shaped single-bore sections. That is a small share of what the slicer meets, so it was
thrown away and rebuilt as a function of NORMALISED ARC LENGTH along the rings the marching
already produces: the traverse interpolates between corresponding points of successive rings.
Every topology the marching can produce is then carried for free, and nothing depends on a
centroid - which matters, since a crescent's centroid lies outside its own material.

Tested on four topologies, all green:

| topology | result |
|---|---|
| annulus | 30 wall crossings against the ring lattice's 2 |
| solid, no bore | sweeps to the medial terminus; 39% of solid length, so it really sparsifies |
| crescent, centroid OUTSIDE the material | 0 emitted points outside the section |
| two bores | refused identically with the option on and off |

Plus a byte-identity test at densities 0, 25 and 50 with the option off, and the density-0 fence.
Full suite green at 175 cases.

### Multi-bore, and the Voronoi cells

A two-bore section has THREE faces, so it is not a linear stack of nested rings and the
serpentine's correspondence between successive rings is meaningless on it. Measured rather than
assumed: continuous shell **refuses such a section outright today**, and has since long before
this work - the marching gives up when a band carries more than one hole. The serpentine neither
worsens nor papers over that, which is what the test asserts.

"Upstream handles it" is the kind of assumption this subsystem has been caught on repeatedly, so
the generator also checks it directly: after resampling, corresponding points on adjacent rings
must sit within about one ring gap. If any pair is far outside that, the chain is not a nested
stack, the traverse would cut across the section rather than through the wall, and it declines
with a note instead of emitting nonsense.

**The Voronoi cell pipeline is the right partner, and the fit is exact.** On
`continuous-shell-multibore-and-infill`, `partition_bores` splits a multi-hole section by a
Voronoi partition of the hole centroids, recursing until every cell holds exactly one bore - and
a single-bore cell IS the linear stack the serpentine requires. So the serpentine needs no
knowledge of cells at all. That is asserted now, while the two pieces are still built apart,
rather than discovered at merge time: a test partitions a two-bore section by hand along the
bisector (which for symmetric bores IS the Voronoi boundary) and confirms the cell traverses
properly - 26 wall crossings, no point outside.

**The integration risk to carry into that merge**, written down while it is fresh: a cut line
between two cells is an ARTIFICIAL face, not a real surface. The serpentine turns at it, so two
adjacent cells' interiors do not connect across the cut - which recreates exactly the failure
this whole generator exists to remove, only at the cut instead of between rings. The branch
already registers cut-line survivors for an "interface weave"; that weave is what must bond
across cuts, and the merge is not done until a peel across a cut line has been printed.

### Status

Off by default (`continuous_shell_interior_serpentine`), fence clean, no existing profile moves.
Nothing here is print-validated yet - the numbers above are all geometry. The next physical step
is a plate carrying the annulus, the solid section and the crescent with the lid off, peeled the
same way the ring lattice was.

### PlateN (2026-08-25): the serpentine's first physical test, with a control

Everything known about the serpentine was geometry. This plate is the first evidence in plastic,
and it carries its own control: **N0 and N1 are the same STL**, the de-lid annulus, differing only
in whether the serpentine is on. PlateG already taught this campaign that an uncontrolled result
is unattributable - an "enormous strength improvement" that turned out to be present in the
untouched part too - so the comparison is on one plate, one filament, one session.

    N0_rings      annulus r10 / bore 2.5, ring lattice          <- the control
    N1_serp       the same STL with the serpentine on
    N2_solid      disc r8, NO bore: sweeps to the medial terminus and back
    N3_crescent   concave, centroid OUTSIDE the material
    N4_ecc        bore offset 2 mm: wall thickness varies all the way round
    N5_lobed      5-lobe outer wall against a round bore
    N6_thin       too thin to traverse: must decline and print sanely

14.23 g, 72 min, lids off so the interior is visible without cutting.

### The measure, and the one that did not work

The obvious statistic - count direction reversals in depth - does not separate these generators
at all, and reported them BACKWARDS: 121 reversals for the ring lattice against 61 for the
serpentine. The lattice undulates six times a revolution on every one of its rings, so it racks
up reversals without ever crossing the wall. What separates them is the SPAN of each monotonic
run: a woven ring turns back within about twice its weave amplitude, a traverse runs face to face.

| object | wall | deep runs (a monotonic run crossing more than half the wall) |
|---|---|---|
| N0_rings, the control | 7.07 mm | **0** |
| N1_serp, same STL | 7.09 mm | **52** |
| N2_solid | 7.61 mm | 39 |
| N3_crescent | 13.97 mm | 17 |
| N4_ecc | 9.87 mm | 44 |
| N5_lobed | 7.88 mm | 49 |

Zero against fifty-two on identical geometry. The ring lattice cannot cross its own wall, which
is the whole reason printed parts peeled into concentric sections, and every non-round section
traverses too.

### Meshes are validated, not assumed

Nothing here needs a general polygon-with-holes tessellator: a tube is capped by a strip between
its two loops on a common parameter, a solid section by a fan from an anchor. Both are checked -
every edge shared by exactly two triangles, no degenerate triangles, and for a fan, that the
anchor actually sees the whole boundary. That check earned its place immediately: it caught the
crescent's bite arc running the wrong way round its circle and carrying a wrong y-coordinate,
before the mesh reached the plate. An unnoticed non-manifold mesh is a whole wasted print.

### What this plate does NOT settle

Everything above is still G-code. The serpentine's claim is that a crack at any depth must sever
traverses rather than separate, and only the peel answers that. **Peel N0 and N1 the same way**,
after removing the bottom shell; the rest are there to show the generator survives shapes the
angular formulation could never have handled, not to be peeled.

N6 is the decline case: a 1 mm wall has nothing to traverse, so it should fall back to the ring
lattice and print normally. Its "deep runs" reading is not meaningful at that thickness.

### Correction (2026-08-25): the tie profile was a measurement artifact

PlateN printed cleanly on every topology - the serpentine's shape claims held - but the verdict
on it was that **the previous foam was better**, wanting only the tubes joined to each other. So
the line went back to the ring lattice, and re-examining it turned up something that invalidates
several numbers recorded above.

`rank_rows` matched rings by their ORDER along a ray and kept only rays showing the modal number
of crossings. On a saturated weave that filter discards more than half the circumference - a ray
grazing a steep radial transition can cross one ring twice or miss it entirely - and the rays it
keeps are not a fair sample. It was keeping **162 of 360 rays**.

Three facts about the same geometry did not fit what it reported:

- adjacent rings correlate at **-0.87 to -0.98**, which is the exact anti-phase the code asks for;
- **every** interface reaches a minimum gap of exactly one spacing, including the ones reported
  as barely tied;
- an idealised model of the very same cap law gives a **uniform 43%** at every interface.

Re-measured with rings assigned by DISJOINT INTERVAL instead of by rank - each ring's excursion
is bounded by its own caps, those bounds do not overlap on this geometry, and nothing is
discarded - over 6455 rays instead of 162:

| config | tie % per interface (bore -> out) | worst |
|---|---|---|
| pin + drift 2 (what was printed) | 44 19 16 37 18 25 33 | **16.3%** |
| pin + drift 0 | 37 29 24 23 25 16 33 | 15.6% |
| no pin + drift 0 | 0 21 24 29 15 29 1 | 0.0% |

**What is withdrawn.** The "0.8% tied" interface, the 74/0.8/15/63 alternation, and the
"permanent doublets" reading built on top of them. The real profile is 16-44% and far flatter.
The mechanism story - that each ring commits to one neighbour and never swaps - was overstated.

**What stands.** The peel is physical evidence and is not in question. The infill-to-shell
interface was genuinely weak and pinning genuinely fixed it: worst interface 0% -> 16%, which the
anchor-hug regression test corroborates independently of this instrument (49.8% tied with the pin
disabled, matching the 43% measured on the plastic).

**What it changes.** The remaining question is now much better posed and much smaller: the model
says a uniform 43%, the slicer delivers 16-44%. That is a bug-sized gap in a specific place, not
evidence that concentric rings cannot be made to connect.

### The per-ring lobe count, tested

The proposal was to scale each ring's lobe count to its own circumference - so a bore ring and an
outer ring get the same lobe PITCH rather than the same lobe COUNT - and to choose the offsets to
set how many connections form. Modelled on the L2 geometry:

| scheme | worst interface |
|---|---|
| exact anti-phase, one lobe count | **43%**, uniform |
| per-ring count from circumference, wavelength 4-10 mm | 0-20% |
| forced step, adjacent counts differing by 1 to 4 | ~19% |

Detuning makes it worse, and the reason is the reverse of the intuition: two rings at exact
anti-phase are in opposition at EVERY azimuth, so they touch wherever both saturate - 43% of the
circumference. Give them different counts and they are only in opposition at the beat nodes, so
contact collapses to a few narrow zones. The pitch argument is right - a bore ring at six lobes
has less than half an outer ring's lobe pitch - but it cannot be bought by detuning neighbours
against each other.

### Rule for this subsystem, learned the hard way three times

Before a measurement is used to justify a design change, check what it DISCARDS. Both instruments
that misled this campaign did so by silently dropping inconvenient samples: a weld threshold
sitting exactly where the data lands, and a ring matcher keeping 45% of rays. Both produced
confident, ordered, plausible tables.

### The 43-vs-16 hunt (2026-08-26): the geometry was right, the target was wrong

The model of the shipped cap law predicted a uniform 43% tie at every interface; the slicer
delivered 16-44%. That gap was the last open question, and chasing it went through four wrong
answers before the right one, each ruled out by measurement rather than argument.

**Ruled out.** Amplitude - the actual p5..p95 envelopes match the prediction ring by ring.
**Phase** - each ring's phase extracted by Fourier component steps by 3.76 to 4.18 radians,
averaging 3.93 against the 3.927 the code asks for. **Arc fitting** - the config has it on, and
the G-code parser reads only G0/G1, so it would have silently dropped every arc; the file
contains zero G2/G3 moves. **Path simplification** - point pitch is a uniform 0.198 mm with a p90
of 0.201, so nothing was decimated. Unfolded, ring 8.15 shows a clean six-lobe wave at the
predicted radii with the saturated duty cycle, two bins high and three low per sixty degrees.

**The geometry was correct all along. The design target was not.**

One spacing is where two adjacent beads sit, so aiming a tie there leaves NO margin: full closure
lands at 0.396 mm against a 0.456 mm contact threshold - a seventh of a bead. Every
sub-tenth-millimetre effect in a real toolpath then flips a tie on or off, which is exactly how a
law that models as 43% delivers 16%.

Aiming at 0.8 of a spacing instead makes the beads at a tie OVERLAP rather than merely touch,
which is what a weld wants anyway. Swept through the slicer and measured with the interval
instrument:

| closure target | tie % per interface | worst | mean | mass |
|---|---|---|---|---|
| 1.00 (as shipped) | 44 19 16 37 18 25 33 | 16.3% | 27.5% | 5.43 g |
| 0.90 | 44 25 22 47 24 31 31 | 21.7% | 31.9% | 5.46 g |
| **0.80** | 44 26 25 52 27 34 27 | **25.0%** | **33.5%** | 5.48 g |
| 0.70 | 42 25 24 49 32 46 8 | 8.0% | 32.4% | 5.49 g |

0.8 lifts the worst interface by half again for 0.9% more material and one minute. By 0.7 it
destabilises - one interface collapses to 8% as excursions start to overrun their neighbours -
so this is a shallow optimum with a cliff on one side, not a monotone knob to keep turning.

Fence byte-identical, suite green at 175 cases, PlateL republished at 5.48 g / 30 min.

### What is still unexplained, and should not be papered over

25% is not 43%. The model still predicts more than the slicer delivers, and the remaining factor
of 1.7 has no explanation yet. Everything ruled out above stays ruled out, so the residue lives
somewhere not yet looked at - the likeliest candidates are the junction and seam regions, which
this analysis has consistently treated as noise to be suppressed rather than as geometry to be
accounted for, and the layer-to-layer phase advance, which every measurement here has averaged
over rather than examined.

**And the ceiling matters more than the residue.** Even a perfect 43% is under half the
circumference. If the foam needs its tubes joined most of the way round, no amount of amplitude
or phase work reaches it - that needs explicit radial ties added to the lattice, keeping the cells
and adding real spokes between them. This closure fix is worth having and worth printing, but it
is an improvement to a design whose ceiling is already known.

### Gap width limit (2026-08-26): chosen against the plastic, not the metric

The closure fix printed: the DENSE object came out fixed entirely, and the sparse object almost -
one central void where it peels, with, in the user's words, "no free-standing foam, all of it
connected to either the outer or inner shell". That is a single parting plane rather than the
column bundle this started as.

Chasing that plane produced a third instrument failure and one genuinely useful finding.

**The failure.** An interface measured at 0.0% tied, apparently pinpointing the plane. It was an
artifact: the band assignment used `delid_noweave.3mf` as its nominal reference, and that file
was sliced BEFORE the stride and closure changes, so the bands were built from stale ring
positions. Re-measured against a matched reference, the sparse wall reads 30/29/32/29 - worst
28.8%, and uniform.

**The finding, which is the useful part.** With correct references, the dense object holds
together with its worst interface tied 25.0% while the sparse object parts at 28.8%. **Tie
fraction does not predict peel resistance.** What separates the two is how WIDE the void is at an
untied azimuth: 2.0 slots dense against 4.0 sparse. A wider void is a bigger crack starter and
leaves less material to tear.

So the sparsifier now caps how wide a gap it may leave. Gaps land on integer slots, so the choice
is effectively binary - 3 or 4:

| limit | sparse gaps (slots) | worst tie | mass |
|---|---|---|---|
| none, as shipped | 3.0 3.8 2.8 **4.0** | 28.8% | 5.48 g |
| **3.0** | 3.0 3.0 1.6 **3.0** | 19.1% | 5.78 g |
| 2.5 | all **2.0** (identical to dense) | 25.0% | 6.16 g |

**3.0 landed.** 2.5 would reproduce the known-good dense geometry exactly, but it abolishes the
density setting for walls this thin - density 25 would emit density 50. 3.0 is the untested
middle, and the print decides.

Two things stated plainly rather than buried. The limit LOWERS the tie fraction, 28.8% to 19.1%,
and that is accepted deliberately because the metric has been shown not to track the outcome. And
the material cost is geometry-dependent and can be much larger than the 5% seen here: the
composition test's 4 mm wall went from 67% of solid to 79%, because on a thin wall the shells
already dominate and added rings land on top of that. Its assertion bar moved from 0.7 to 0.85
with that reason recorded in the test.

**Next:** print PlateL and peel both objects again. If the sparse one holds, a wider battery of
geometries follows, then W4.

### VALIDATED (2026-08-26): the foam holds, both densities

PlateL reprinted with the closure target at 0.8 spacings and the gap cap at 3 slots. **Both
objects passed the peel test.** Sectioned, both show a regular honeycomb - closed cells, connected
through the wall depth - rather than the concentric tubes, columns, or parting planes that every
earlier version produced.

That closes the arc this campaign has been on since the first de-lid print:

| print | sparse | dense |
|---|---|---|
| de-lid, original | columns, not connected at all | columns |
| + tie-band stagger | still columns | almost a foam |
| + per-gap amplitude | 3 concentric sections | 3 concentric sections |
| + anchor pin, stride | peels at the shell faces | peels at the shell faces |
| + closure 0.8 | one central parting plane | **holds** |
| + gap cap 3 slots | **holds** | **holds** |

The two changes that mattered, both landed in the last day, were unglamorous: aim the ties at 0.8
of a spacing instead of exactly one, and refuse to leave a gap wider than three slots. Neither is
a new mechanism. Both are corrections to targets that had been set at the edge of what works.

**What the metrics got wrong, kept here because it was expensive.** Tie fraction - the statistic
this subsystem was steered by for a week - does not predict peel resistance: the dense object
holds at 25% and the sparse one parted at 28.8%. Void WIDTH does track it. And the gap cap that
fixed the sparse case LOWERED its tie fraction to 19.1%, so optimising the metric would have
rejected the change that worked.

### PlateP: the battery, with one pre-registered prediction (2026-08-26)

PlateL proved the foam holds on ONE geometry. Everything that made it hold is geometry-dependent:
the gap cap decides ring counts from the buried span, the per-gap amplitude sizes excursions from
measured gaps, and the anchor pin needs a kept ring beside a shell. A round tube exercises none of
that. So: eccentric bore, lobed wall, crescent, solid section, a taper whose ring count changes
partway up, and a wall too thin to sparsify - 19.73 g / 95 min, lids off.

Measuring the ring gaps before printing found a hole in the fix:

| object | rings | ring gaps (slots) | max buried |
|---|---|---|---|
| P1_sparse | 10 | 3.0 3.0 1.6 3.0 3.0 | 3.0 |
| P2_dense | 12 | 2.0 x6 | 2.0 |
| **P3_ecc** | 8 | 2.1 **4.4** 2.0 | **4.4 - over cap** |
| P4_lobed | 10 | 2.8 2.8 1.5 2.8 2.8 | 2.8 |
| P6_solid | 4 | 3.0 | 3.0 |
| P7_taper | 8 | 3.1 3.2 3.1 | 3.2 |
| P8_thin | 4 | sub-spacing | declined, as intended |

**The cap is expressed in SLOTS, but adaptive beading makes bead width vary with azimuth.** On an
eccentric wall the thick side carries wider beads, so a three-slot gap there is 4.4 spacings of
physical void. The cap bounds ring COUNTS; it does not bound the quantity that actually predicts
peeling, which is void WIDTH in millimetres.

That is left in deliberately, as a **pre-registered prediction**: if void width governs peel
resistance, P3 comes apart and the others hold. Fixing it first would have tested the fix;
printing it tests the mechanism, and the mechanism is what the rest of this design rests on. The
fix, if needed, is to convert the cap through the LOCAL bead width rather than the nominal one.

P5_crescent measures 15.4 slots, but that number is meaningless - it comes from rays cast from a
centroid lying outside the material. There is no valid pre-print reading for it; the print is the
measurement.

P7 is the one to examine along its HEIGHT rather than across: its ring count changes partway up,
so the cap re-decides mid-part, and the layers either side of that change still have to tie.

### PlateP result (2026-08-26): 7 of 8 hold, the prediction was wrong, and the crescent is a GATE

Seven of eight objects came out acceptable. Only the crescent is peelable. A few objects showed
what read as overly-solid columns, noted for later.

**The pre-registered prediction failed, and that is the most useful thing here.** P3_ecc was
predicted to peel because a three-slot gap on its thick side is 4.4 spacings of physical void.
**P3 held.** The object that failed was the crescent - the one object explicitly recorded as
having no valid pre-print measurement. So void width does NOT govern peel resistance on its own;
the story built on it is at best incomplete, and should stop being leaned on. What it did do
correctly was motivate the gap cap, which fixed the sparse annulus - but the mechanism behind
that is now open again.

**The crescent's cause is a gate, not a tuning failure.** The weave is parameterised by ANGLE
ABOUT A CENTROID, so it is guarded by `star_ok`:

    if (woven.weave > 0.0 && cc.star_ok) { ...the entire weave... }

A crescent's centroid sits at (-3.66, 0) and the section is not star-shaped about it, so the gate
fails and the weave is skipped WHOLESALE. The crescent gets sparse concentric rings with no
undulation whatever - zero ties by construction. That is the original failure mode of this whole
campaign, preserved intact in the non-star path, and it has been there the entire time: every
plate until now was star-shaped, so nothing ever exercised it.

**The fix is one already done once.** The radial serpentine hit exactly this and was rebuilt in
ring-index space - normalised arc length along the rings the marching produces, which needs no
centroid and no star assumption. The weave needs the same treatment: modulate along each ring's
own arc-length parameter, displace along the ring's outward normal rather than along a radius
from a centroid, and delete the `star_ok` gate on the weave. `resample_ring` already exists.

Worth stating plainly: this means every foam result in this ledger is a result about STAR-SHAPED
sections. The lobed and eccentric objects are star-shaped; concave outlines were never tested
until this plate, and they were never woven at all.

### The crescent fixed, and what the eccentric object's solidity is (2026-08-26)

### The weave is now local to each point

Two formulations failed before this one, both because they needed a global frame for something
purely local.

**By ANGLE about the centroid** - the shipped version. A ring is single-valued in angle only if
the centroid sees all of it, so the weave was gated on `star_ok`, and sections failing that gate
got NO weave at all: thinned to the sparse stride but never displaced, leaving plain nested tubes.
That is this subsystem's original failure mode, and it sat in the non-star path for the entire
campaign because every plate until PlateP was star-shaped.

**By normalised ARC LENGTH along each ring** - the first attempt at a fix. No centroid needed, but
it assumes equal fractions of arc length on adjacent rings lie across from each other. On a wall
of varying thickness they do not, so the correspondence check tripped and the ECCENTRIC object
stopped being woven too. Measured: crescent still +0.0%, eccentric fell +8.0% -> 0.0%. Caught
before it went anywhere near a plate.

**By NEAREST POINT on the neighbouring ring** - what landed. A point needs only its distance to
the next ring and the direction to it, and both are local. No parameterisation, no correspondence,
no centroid, no star assumption, nothing to gate on. The nearest point is walked rather than
searched, so it stays linear.

| object | weave off | angular (printed) | nearest-point |
|---|---|---|---|
| P1_sparse | 386 mm | +7.8% | +7.8% |
| P2_dense | 460 mm | +4.4% | +4.4% |
| P3_ecc | 548 mm | +8.0% | +7.5% |
| P4_lobed | 424 mm | +6.1% | +5.7% |
| **P5_crescent** | 389 mm | **+0.0%** | **+13.3%** |
| P6_solid | 258 mm | +13.9% | +14.0% |
| P7_taper | 365 mm | +6.2% | +6.3% |
| P8_thin | 128 mm | +0.0% | +0.0% (correctly declines) |

The peel-validated annulus is unchanged against a matched weave-off reference: L1 worst 19.1%
(identical), L2 24.9% against 25.0%. Fence byte-identical, suite green at 175.

**The test that asserted the old behaviour was inverted, not deleted.** It used to require the
weave be declined on a crescent, and its stated fear - that clamping to an angle-binned profile
could drag a bead "half the section width off its ring" - was legitimate for THAT formulation.
The safety claim is kept and made sharper: every emitted point must lie inside the section (0 of
1038 outside), and no point may sit further than 15% of the section span from where the unwoven
path put it (measured 0.57 mm of 16 mm).

### The eccentric object's solidity

P3_ecc covers 63% of its wall against P1_sparse's 44% at the same nominal density, and it tracks
thickness: 42% on the thin side, 73% on the thick. Reading the kept rings by azimuth shows why -
on the thick side the gaps run `1.1 1.1 2.1 3.1 1.0 1.0 1.0 1.1 1.0 1.0 1.0 1.0 2.0 2.0 1.0 1.0`.
That long run of one-slot gaps is SOLID: consecutive rings, nothing skipped.

This is NOT yet explained. The likeliest cause is that the marching splits an eccentric wall into
more than one chain - the thin side closes first, so a hole-less band continues as a separate
SOLID chain - and each chain then carries its own shell allowance, or fails `chain_thick_enough`
and prints solid. That is a hypothesis, not a measurement: it has not been confirmed, and the
sparsifier's behaviour on multi-chain sections has not been looked at. Recorded so it is not
mistaken for a settled account.

### PlateP reprint (2026-08-27): crescent much improved, eccentric solidity DIAGNOSED

The crescent now peels only on its **concavity side** - the rest of it holds - and no other object
regressed. So weaving the non-star path was the right fix and it worked; what remains on that
object is localised to the concave region rather than general.

### The eccentric object's solid ring: confirmed, and it is a chain boundary

Two probes failed before the right one, both worth recording because they looked plausible:

- **Counting separate extrusion runs per layer** gave P3_ecc 4 and P4_lobed 242, which would have
  meant severe fragmentation. Both were artifacts: continuous shell emits Z on nearly every move
  for the ramp, and the run-counter treated Z-only moves as breaks.
- **Counting retractions and travels** returned ZERO of each for every object - correct, and it
  disproves fragmentation outright, because continuous shell links everything into one path with
  connectors rather than travels. Which also means travels can never reveal chain structure.

Instrumenting the slicer directly settles it:

| object | chains per layer |
|---|---|
| P1_sparse | 1 |
| P4_lobed | 1 |
| P5_crescent | 1 |
| **P3_ecc** | **2** |

An eccentric wall pinches closed on its thin side before its thick side, so `march_bands` ends the
annular chain there and continues the remainder as a SECOND chain. Each chain then applies
`shell_loops` at both of its own ends - but the end where the two chains meet is a PINCH
ARTIFACT, not a real surface. Those interior shells are the solid ring, and they explain the
whole measurement: 63% wall coverage against 44% for a plain wall, tracking thickness, with a run
of one-slot gaps mid-wall where two chains' shells sit back to back.

Note this is unrelated to the gap cap, which was the earlier guess. It is also not a decline: the
slicer emits no degrade note for this object at all.

**The fix, specified but not built:** a chain needs to know whether each of its ends is a real
face or a pinch, and apply the solid shell allowance only at real faces. The information exists
at the moment of the split - `march_bands` continues when a band loses its hole, and that
continuation point is exactly the artificial end.

### Pinch-boundary shells, and phase as a field (2026-08-27)

### The eccentric object's solid ring is gone

A chain now knows whether its OUTER end is a real face or the point at which a previous chain
pinched out, and the solid shell allowance is applied only at real faces. A shell exists to skin
a surface; where the end is a pinch artifact there is no surface, and skinning it printed a solid
ring buried mid-wall.

    P3_ecc wall coverage   62.4% -> 53.4%
    P3_ecc path per layer  589 mm -> 466 mm   (-20.7%)
    plate                  19.99 g -> 19.26 g

Every other object is byte-identical in path length, which is the check that matters: the flag is
only ever set on a chain the marching created by continuation, and no other object on the plate
produces one.

### Phase is now a FIELD across the wall

The crescent still peeled on its concave side after the weave was made to apply there. The weave
phased each ring by its OWN normalised arc length, which gives every ring the same lobe COUNT. On
an annulus that is harmless - all rings are nearly the same length, so anti-phase holds the whole
way round. On a crescent the inner rings are far shorter than the outer ones, so the same count
means different wavelengths and adjacent rings BEAT: in opposition in some places, in phase in
others. Two rings tie only where they oppose, so the ties collapse wherever the beat puts them in
phase. The weave was applying, but detuned.

Phase is now carried ACROSS the wall. Ring 0 takes its phase from its own arc length; every ring
after it inherits the phase of its nearest point on the ring outside it. Adjacent rings are then
in exact anti-phase everywhere by construction whatever their lengths, and the lobes line up along
the wall rather than around each ring independently.

| object | before | after |
|---|---|---|
| P1_sparse, P2_dense | 416 / 481 mm | **unchanged** |
| P5_crescent | 440 mm | +2.1% |
| P4_lobed | 448 mm | +1.7% |
| P3_ecc | 466 mm | +0.2% |
| P6, P7, P8 | | unchanged |

The annuli being untouched is the point: on a circle the inherited phase IS the same angular
position, so this reduces exactly to what they already did. The change lands only on sections
whose rings differ in length - which is where the defect was.

Peel-validated annulus, matched reference: L1 worst 19.4% against 19.1%, L2 25.1% against 24.9%.
Fence byte-identical, suite green at 175.

**Stated plainly: the beat is inferred, not measured.** Every instrument in this ledger that can
measure tie quality needs rays, a centroid, or a ring correspondence, and none of them work on a
crescent - three separate attempts to measure that object returned numbers known to be wrong
(coincident-point gaps of 0.000 mm, an inflated section area). The mechanism fits the geometry and
the fix moves exactly the objects it should, but the crescent print is what will decide it.

### PlateQ (2026-08-27): incremental on all three, and the crescent measured at last

The three-object plate came back an incremental improvement across the board. The eccentric's
solid chunk is still present though smaller; the crescent is mildly less peelable but its
concave-side outer ring still has little interconnect with the foam; the lobed object is
essentially unchanged, which is what it was there to confirm.

### An instrument that works on any topology

Five attempts to measure the crescent had failed, each returning numbers that were provably
wrong - coincident-point gaps of 0.000 mm, an inflated section area, a tie profile built on stale
nominal rings. Every instrument this subsystem had needs a centroid, a star-shaped outline, or
rings of comparable length, and a concave section has none of those.

`CS_DUMPRINGS=<layer>` now writes the chain rings themselves, so analysis has exact ring identity
and the outline's shape is irrelevant (`tools/ringanalyse.py`). It is kept rather than deleted:
the cost of working without it, across this campaign, has been several wrong conclusions.

### What it says about the crescent, and what remains unexplained

| interface | tied CONVEX | tied CONCAVE |
|---|---|---|
| 2 (anchor to foam) | 45.8% | **5.1%** |
| 3 | 27.1% | **0.0%** |

The user's reading confirmed and quantified: the concave side is tied roughly nine times less.

What was then ruled OUT, each by measurement:

- **Gap geometry.** Unwoven gaps are identical between the regions - 1.189 mm at interface 2,
  0.793 at interface 3, on both sides.
- **The weave failing to apply.** Ring 3 is displaced 0.238 mm (concave) against 0.252 mm
  (convex): it is weaving there.
- **Reach.** Both regions achieve full closure at their peaks - minimum gap 0.317 mm convex,
  0.280 mm concave, against a 0.456 mm tie threshold.
- **Lobe starvation.** Ring 3 gets 3 lobes over 34.7 mm of convex arc and 2 over 22.9 mm of
  concave arc - the same pitch either side.

So the ring reaches the anchor at its peaks, with the same lobe density and the same gap to
close, and is nonetheless tied over nine times less of its concave length. **That is not
explained.** Five hypotheses have died on this object today - the star gate (real, fixed, not
sufficient), arc-length correspondence, the beat, phase-field flattening, and lobe starvation -
and the honest position is that the cause is still unidentified rather than that the next guess
is likely right.

The instrument now exists to settle it, which it did not before. The next step on this object is
to plot ring 3's displacement against arc length either side and look, rather than to reason
about it.

### Queue (2026-08-27)

**Parked, with everything needed to resume:**

1. **The eccentric object's remaining solid chunk.** The pinch-boundary fix cut its path 20.7% and
   the chunk visibly shrank, so a second and smaller contributor is still there. Start from the
   fact that P3_ecc is the only object on the battery that emits TWO chains per layer; the first
   contributor was each chain skinning its own pinch end. Check what else is per-chain: the
   `chain_thick_enough` gate, the shell allowance at the chain's INNER end (the terminus), and
   whether the second chain's own gap cap re-decides its ring count independently.

2. **The crescent's concave side.** Tied 5.1% against 45.8% convex at the anchor interface, 0.0%
   against 27.1% one ring in. Cause UNIDENTIFIED - four candidates ruled out by measurement (gap
   geometry, weave not applying, failure to reach, lobe starvation), five hypotheses dead. The
   instrument to settle it now exists and did not before: `CS_DUMPRINGS=<layer>` plus
   `tools/ringanalyse.py`. The next step is a desk exercise, not a print: plot ring 3's
   displacement against arc length either side of the concavity and look at it.

Both are real but neither is the common case, and the annulus, eccentric, lobed, solid, taper and
thin-wall geometries are all in acceptable shape. Proceeding to W4.

> **Both closed on 2026-08-27** - see "The weave was measuring gaps that cross no wall" at the end
> of this file. They were the same defect: the amplitude law was fed a distance that crosses no
> wall. Awaiting print validation, not further analysis.

### Queue (2026-08-27, revised)

1. **The axial core of a solid section.** Diagnosed (see the correction above): the weave adds a
   FIXED `lobes * 4 * amplitude` of path to every ring regardless of its radius, so near the axis,
   where a ring is short, it is a 130% surcharge and the ring folds into a tight zigzag. The fix
   is to bound the amplitude by the local wavelength `2 pi r / lobes` as well as by the gap. That
   changes geometry on every part with a small-radius feature, so it is print-gated. The ring
   selection is NOT the lever - moving the inner anchor onto the terminus was tried and the
   amplitude law spends the saving.

2. **A combined battery over the Voronoi subsectioning**, which allows multi-bore geometry. Blocked
   on `continuous-shell-multibore-and-infill` being brought forward. The sharp question: whether a
   peel runs along a CELL CUT LINE - an artificial face the generator will treat as real. Note
   that the reach bound landed above makes this more interesting, not less: a Voronoi cell
   boundary is exactly the situation where one chain's rings sit in another's measured gap.

### W4 (2026-08-27): a forced-solid layer keeps adaptive beading

`force_solid` used to switch adaptive beading off, so a top/bottom shell - or any shell over air -
fell back to C-lite marching for its bead COUNTS while the layers above and below were planned.
The two disagree on how many beads a wall of a given thickness takes, so the count STEPPED at the
boundary between a solid layer and its neighbours. On the S5 cone that printed as a doubled lip at
the bore, where the solid cap laid a different number of beads than the wall it sat on.

Beading decides bead counts and widths; force_solid decides whether the interior is sparsified.
They are orthogonal, and coupling them is the same mistake W3 fixed for sparse infill - which cost
a 2.035 mm wall four stretched beads instead of five.

The fix is one clause. What makes it convincing is the signature table, where the forced-solid row
used to be byte-identical to the adaptive-OFF row and is now byte-identical to the adaptive one:

| fixture row | before | after |
|---|---|---|
| `adaptive=0` (C-lite) | `803b06085a5e` | `803b06085a5e` |
| `adaptive=1` | `86d421e06ffb` | `86d421e06ffb` |
| **`force_solid`** | **`803b06085a5e`** | **`86d421e06ffb`** |

2758 points and 624.2 mm become 2561 and 580.2 - exactly its neighbours' layout. Coverage moves
from 1.056 to 0.982, which is the adaptive path's own figure on this fixture (its plain row reads
0.982 and its spiral row 0.968), so this is the solid layer adopting the plan's coverage rather
than a new deviation. Absolute error against a full layer improves, 5.6% over to 1.8% under.

**The fence moved, deliberately, and only where it should.** Four force_solid rows exist; the two
annulus rows changed and the two solid-disc rows did not, because the plan declines on a section
with no bore. No other row moved.

On the plates: bottom-shell layers shift by 0.24 mm at p99 and **the mid-wall is bit-identical**
(0.0000 mm) on every object. The eccentric and crescent objects do not move even in their bottom
shells - the first is multi-chain, the second has no bore.

The interior SPIRAL stays gated on force_solid and has its own test, because the two used to ride
the same flag and only one of them was meant to change.

**Not yet print-validated.** The defect this fixes was seen on the S5 cone, which is not on any
current plate. The geometry to reprint is the S5 cone; what to look for is the bore lip, which
should no longer be doubled.

### W5 (2026-08-27): the notes were single-slot, so half of them never appeared

The plumbing was already in place - `degrade_note` threaded out of `build_cross_section` and
logged by LayerRegion with the layer id. What was missing was coverage and, worse, a defect in the
mechanism itself.

**Every site guarded on the note being EMPTY**, so a layer that took two weaker paths reported
one, and which one depended on the order the code happened to run in. That is not a cosmetic
problem: it is the same blindness W5 exists to remove, one level up. Notes now accumulate through
a single `add_note`, and a multi-bore section correctly reports both of its declines:

    adaptive beading out of scope: multi-bore section; marching stopped early: a band split
    into several bores - the material inward of that point is NOT filled

**Two silent declines surfaced:**

- **The multi-hole break in `march_bands`.** A band that splits into several bores ends the chain,
  and the material inward of that point is not filled. This is precisely the family that cost an
  eccentric bore 31% of its material with nothing in the log.
- **The solid-cap brick stagger on a non-star section.** Unlike the weave, which was rebuilt to
  need no angle, the stagger is still parametrized by angle about the centroid and is genuinely
  undefined on a folded outline. Declining is right; declining silently means a solid cap's
  interior rings stack into an aligned radial plane with nothing recording it.

**The test caught itself being vacuous.** Its first version asked only that a note be non-empty,
and passed while proving nothing - it was seeing an unrelated pre-existing note about adaptive
beading being off, because of the empty-slot guard above. It now asserts the specific text. That
is the fourth time in this campaign a green test turned out to be asserting nothing; the rule
stands - a test is not evidence until you have watched it fail for the reason you intend.

Every fence signature is byte-identical: W5 adds information, never geometry. Suite green at 178.

### Queue addition (2026-08-27)

**A combined battery over the Voronoi subsectioning.** The geometry battery covers what a single
`ContourChain` can be - annulus, eccentric, lobed, crescent, solid, taper, thin. It does not cover
what the multi-bore pipeline on `continuous-shell-multibore-and-infill` unlocks, and that is where
the interesting cases are: two bores, bores of different sizes, a bore near a wall, cell CUT LINES
where two Voronoi cells meet.

That last one is the point. A cut line is an artificial face, and everything this campaign learned
about faces says the interior generator will treat it as one - the serpentine turns at it, the
sparsifier skins it (until told not to, as the pinch fix now does), and the weave anchors to it.
Whether a peel runs along a cut line is the single question that merges the two workstreams, and
no plate has ever asked it.

Blocked on that branch being brought forward; it has diverged a long way from main. Recorded here
so it is not rediscovered late.

### W6 (2026-08-27): the solid region is a mask, not a whole layer

`force_solid` was a whole-layer boolean. If any part of the section had air within the shell
count, EVERY ring in the layer was kept. For a genuine top or bottom shell that is right - the
whole face is exposed. For an INTERNAL surface it is not: the floor of a blind bore needs solid
over the bore's footprint, and the wall around it is buried. That layer printed solid edge to
edge.

`cs_solid_mask` now answers WHERE rather than whether: the union of the uncovered regions over the
shell layers, grown by a wall thickness so the solid patch anchors into the surrounding wall
instead of ending at the surface's edge. On a 20 mm disc with the middle masked:

| | path per layer |
|---|---|
| sparse | 390.9 mm |
| **masked middle** | **446.9 mm** |
| fully solid, as before | 698.1 mm |

36% less path on such a layer, with solid still where the surface needs it. A whole-layer boolean
can only produce the two endpoints; that the masked case lands between them is the whole claim,
and it is what the test asserts.

**Granularity is a WHOLE RING**, and that is a deliberate limit rather than an oversight. The
rings are closed loops and every stage downstream chains closed loops, so keeping half a ring
would need open paths through the entire pipeline. It also happens to be the right unit for the
case that motivates this: a bore floor needs its inner rings, and those are whole rings. A surface
that cuts ACROSS the rings - a shelf on one side of a tube, say - will still promote every ring it
touches, which is more solid than strictly necessary. That case is not addressed here.

**Real top and bottom shells are untouched.** A mask covering 95% or more of the section is
treated as the whole-layer case and stays on the old path, so the common case is byte-identical -
every fence signature is unchanged.

**Not print-validated.** The geometry that exercises it is a blind bore: a tube with a floor part
way up, where the floor layer should now be solid over the bore and sparse in the wall. That is
the object to slice.

### PlateR (2026-08-27): the geometry W4 and W6 actually need

Neither change was exercised by any existing plate - every one of them has its lid off and no
internal surface, which is precisely the geometry that avoids force_solid. 9.05 g, 45 min.

    R1_tube       plain tube, LID ON. Its top layers are forced solid AND annular, which is the
                  case W4 changes. Look down the bore at the top face: no step where the cap
                  meets the wall.
    R2_blindbore  a bore starting part way up, so its FLOOR is an internal surface - the case W6
                  is for. Cut it at the floor: solid disc in the middle, foam in the ring around.
                  The bore opens UPWARD on purpose; a bore closed at the top would put a ceiling
                  over air, and asking TPU to bridge 8 mm would test the bridge, not the mask.
    R3_control    the validated annulus, unchanged by either change.

**W6 measured end to end**, through the real LayerRegion mask rather than the unit test's
hand-built one, by slicing the same plate with the mask suppressed:

| layer | z | whole-layer solid | W6 mask | saved |
|---|---|---|---|---|
| 19 | 4.75 mm | 680 mm | 443 mm | 35% |
| 20 | 5.00 mm | 813 mm | 577 mm | 29% |
| 21 | 5.25 mm | 776 mm | 460 mm | 41% |
| 22 | 5.50 mm | 813 mm | 575 mm | 29% |
| 23 | 5.75 mm | 776 mm | 457 mm | 41% |

21% less over the shell band under the floor, and the layers away from it are untouched - 0% at
layers 16, 17 and 26, which is the containment that matters. Whole object 27896 -> 26525 mm.

**What each object can and cannot settle.** R2 shows W6 in plastic directly: the floor either has
foam in its outer ring or it does not. R1 is weaker evidence for W4 - there is no same-plate
control, because W4's old behaviour cannot be reproduced per object (it was a disagreement
between two beading paths, and switching adaptive off changes both). The comparison available is
against the S5 cone already in hand, which is a different geometry. What R1 can show on its own is
the absence of a step at the cap boundary, which is what the defect looked like.

### PlateR result (2026-08-27): W4 and W6 validated in plastic

The plate came out well. The blind-bore object shows a compact solid column under the bore with
properly foamed material around and below it - which is W6 doing exactly what it was built to do.
No object showed a regression.

**One reading corrected, because the direction matters.** The solid column looked smaller than the
bore by eye. Measured, it is slightly LARGER: solid extends to r = 4.6 mm on the masked layers
against a bore radius of 4.0 mm. That is the bore footprint grown by about a wall thickness
(`wall_loops * spacing` = 0.79 mm, trimmed a little by whole-ring granularity), and the direction
is the whole point - a solid region NARROWER than the bore would leave the floor's rim sitting on
foam, which is the failure the grow term exists to prevent.

What reads as "small" is its HEIGHT: the column is only as tall as the shell count, about five
layers or 1.25 mm, where before W6 those same layers were solid edge to edge out to r = 10 mm.
The layer below the shell band measures solid only to r = 2.3 mm, which is not the mask at all -
it is the naturally dense core near the axis of a solid section.

So both changes are now print-validated:

- **W6** directly and unambiguously: the floor is solid where the bore sits and foam elsewhere,
  and the G-code showed 29-41% less path on those layers with the layers away from them untouched.
- **W4** by the absence of the defect: no step where the cap meets the wall. Weaker evidence, as
  recorded when the plate was built - there is no same-plate control for it - but the signature
  table already showed the forced-solid layer adopting its neighbours' layout exactly.

### The axial core in a solid section (2026-08-27) - pre-existing, queued

The blind-bore part has a solid column running the full height BELOW its floor, coaxial with the
bore and narrower than it. That is not the W6 mask, and not the shell band under the floor - it is
present on layers nowhere near either, and the slicer emits no decline note for the object at all.

**It is the medial core of a solid section.** Below the floor R2 is a plain cylinder, and a
concentric fill of a solid section has every ring converging on the axis. The terminus ring and
its pinned anchor are always kept, so within roughly r = 2.5 mm several rings survive whatever the
density asks for. Path lying inside r = 2.5 mm, by layer:

    layers 3-17:   45, 17, 45, 17, 45, 17 ... mm   every layer occupied, alternating
    layers 24+:    0 mm                            above the floor, where the bore is

The core area inside r = 2.5 mm is 19.6 mm2, and 45 mm of bead at 0.4 mm wide covers about 18 of
them - so odd layers are nearly solid there and even layers about a third. Stacked, that reads as
a solid column. It stops dead at the floor because the bore replaces it, which is why it appears
aligned with the bore: same axis, not the same cause.

**Not a regression.** Those layers are not forced solid, so neither W4 nor W6 touches them. This
is how a solid section has always been filled; the blind-bore geometry simply made it visible by
putting a bore directly above it.

**The odd/even alternation is unexplained.** It is not the brick stagger, which only runs on
forced-solid layers. `reverse_path` alternates with layer parity and would move where the terminus
lands, which is a candidate and nothing more.

**Queued** alongside the eccentric solid chunk, which is plausibly the same family: both are
places where the fill is denser than the density setting asks, and in both the mechanism sits
near a chain end rather than in the weave.

### Correction: the odd/even alternation was the instrument (2026-08-27)

The alternating 45/17 mm figures above are wrong, and the "unexplained" verdict with them. They
came from bucketing G-code segments by `int(round(z / layer_height))`. The continuous shell writes
a per-POINT Z - the Z-wave, and the layer-join ramp - so that expression does not name the layer a
segment belongs to, and it scattered each layer's path across two buckets.

Bucketed by the G-code's own `; CHANGE_LAYER` markers, the profile is flat:

    layers  99-114:  ~31 mm per layer      the plain cylinder below the floor
    layers 115-119:  ~49 mm per layer      the five shell layers under the bore floor (W6)
    layers 120+:      0 mm                 above the floor, where the bore is

So there is no alternation to explain. The core is a steady ~31 mm of path per layer inside
r = 2.5 mm, which is 12.4 mm2 of bead in a 19.6 mm2 disc - **63% dense where the setting asks
25%**.

The ring dump says where it comes from. Kept rings on layer 8, by radius:

    shell    9.775  9.379  8.982                        two shell rings plus a pinned anchor
    infill   7.819  6.617  5.428  4.238  3.049  1.860   stride 3 slots, weaving
    anchor   0.658                                      kept infill, pinned, un-woven
    terminus 0.262                                      always kept

Two contributions, and the smaller one is the obvious one:

- **The terminus and its anchor sit ONE slot apart** at the axis - 5.8 mm of path between them.
  Moving the distribution's inner anchor onto the terminus itself removes that doublet; it was
  tried, and it does not help, because the amplitude law then grants the next ring a bigger cap
  and it spends the saving.

- **The weave adds a FIXED length to every ring**, not a proportional one. The drive saturates
  (`clamp(4 sin, -1, 1)`), so a ring is a square wave of `lobes` periods and amplitude A, and the
  radial travel is about `lobes * 4 * A` whatever the ring's radius. At r = 7.8 mm that is a 21%
  surcharge on the circumference; at r = 1.86 mm the same 15.6 mm is a 130% surcharge, and the
  ring folds into a zigzag of 1.94 mm wavelength and 1.3 mm peak to peak. 24 of the 30 mm inside
  r = 2.5 mm is that one ring.

So the axial core is the weave's fixed length surcharge landing in the smallest area of the
section. It cannot be fixed by choosing different rings; it needs the amplitude bounded by the
LOCAL wavelength (2 pi r / lobes) as well as by the gap - which changes geometry on every part
with a small-radius feature, and so is a print-gated change rather than a desk one.

## The crescent and the eccentric: the weave was measuring gaps that cross no wall (2026-08-27)

Both parked objects turned out to be the same defect wearing two costumes. The amplitude law is

    cap = weave * share * (gap - 0.8 * spacing)

and `gap` is the distance to the nearest point on the neighbouring ring. That is the right
quantity while the two rings FACE each other across the wall. Neither object's did.

### The crescent: the nearest point was walked, and a concavity has two basins

`nearest_walk` carried an index from one sample to the next and hill-climbed it forwards and
backwards while the distance improved - on the reasoning that adjacent points have adjacent
neighbours. That holds until the section is concave. A ring wraps AROUND a notch, so a point
facing it has one candidate across the gap and another round the far side, and the distance
function along the neighbouring ring has two basins with a ridge between them. The walk cannot
cross a ridge, so once the index settled in the wrong basin every later sample inherited it.

Transcribing the walk into Python and running it against a true nearest, on the slicer's own
dumped rings:

| ring pair | region | walk says | truth | wrong on |
|---|---|---|---|---|
| 3 -> 2 | convex | 1.189 mm | 1.189 mm | 1.5% |
| 3 -> 2 | concave | 5.206 mm | 1.189 mm | **100%** |
| 1 -> 0 | concave | 5.573 mm | 0.396 mm | **100%** |

Every ring, over the whole concave stretch. The law then permitted 7.458 mm of excursion against a
true 0.873, and the printed rings moved 7.384. The same index also chooses the DIRECTION of the
excursion, and carries the phase field from ring to ring, so one wrong answer spoiled the
amplitude, the direction and the phase together. Four of the five dead hypotheses were true
observations of a downstream symptom.

Searching for the nearest point instead of walking toward it, same object, same layer:

| | welded, CONVEX | welded, CONCAVE |
|---|---|---|
| before | 41.5% | **5.7%** |
| after | 41.5% | **38.1%** |

Ring 3's displacement over the concave stretch falls from a mean of 2.256 mm (max 7.384) to a mean
of 0.544 (max 0.872 - the cap, to three digits). The convex side does not move.

### The eccentric: a chain reads its neighbour chain's material as a void

Band marching splits a wall into chains where it pinches, and each chain measures gaps only within
itself. On P3_ecc the first chain's rings 3 and 4 lie a median 3.716 mm apart - with the SECOND
chain's three rings sitting between them. The law read that as the gap to close, permitted 5.063
mm, and flung both rings 2.53 mm into a pocket that was already full. That is the solid chunk that
survived the pinch-boundary fix.

Here the walk was innocent: walk and truth agree at 3.716 mm. The gap is real; it is just not a
wall crossing.

### The bound

A neighbour further off than a wall crossing is not a neighbour. The bound is the one the
sparsifier already guarantees on the other side - it raises the kept-ring count until no gap
exceeds `CS_MAX_GAP_SLOTS`, and consecutive band contours are one spacing apart, so a genuine
crossing between kept rings is at most that many spacings. Half a slot of margin keeps a gap that
sits exactly at the cap (the crescent's 1.189 mm) from being excluded by rounding. Applied per
POINT: a ring can face its neighbour over most of its length and run past its end for the rest.

### Blast radius

The excursion the law permits, on the same unwoven rings, old law against new:

| object | old | new | object | old | new |
|---|---|---|---|---|---|
| P1_sparse | 0.885 | 0.885 | P6_solid | 0.873 | 0.873 |
| P3_ecc | **5.063** | **1.064** | P7_taper | 0.875 | 0.875 |
| P4_lobed | 1.016 | 1.016 | P8_thin | 0.035 | 0.035 |
| P5_crescent | **9.837** | **1.041** | | | |

Exactly the two objects that failed in the plastic, and nothing else. On the test fixtures the
same split appears as topology - a wall that marches as ONE chain is byte-identical, and only
walls that pinch into chains lose the excess:

| eccentric annulus r10/bore3 | surplus path, old | new |
|---|---|---|
| bore offset 0 mm (one chain) | 6.47% | 6.47% |
| bore offset 1 mm (one chain) | 6.69% | 6.69% |
| bore offset 2 mm (pinches) | **9.59%** | **2.46%** |
| bore offset 3 mm (pinches) | 2.11% | 1.18% |

Density-0 signature fence byte-identical across all 20 fixtures. Suite 49 -> 50 cases.

### Two things tried and reverted, with the reason recorded in the code

- **`chain_thick_enough` double-counts an outer shell on a pinch chain**, where the sparsifier
  holds none (`sl_out`). The arithmetic is genuinely wrong, and correcting it makes things worse:
  a three-ring pinch chain cannot be thinned anyway - two infill rings and a terminus, with the
  distribution obliged to anchor both ends - so the only effect is to start weaving it. On the
  eccentric pocket that turned the outer ring from a 41 mm loop into an 80 mm one sweeping 4.2 mm
  of radius, and silenced a "printed solid" note that was true.

- **Moving the distribution's inner anchor onto the terminus** for a bore-less chain, to remove
  the one-slot doublet at the axis. It removes the doublet and saves nothing, and it broke the
  weave on thin chains by spending one of the distribution's slots on a ring kept anyway.

### The note now reports the outcome

"wall too thin to sparsify - printed solid" fired on `chain_thick_enough` saying no. It now fires
on no ring having been dropped, which is the thing that is true of the plastic. They part company
on a short chain whose span leaves the distribution no room: that chain printed solid and nothing
said so.

### One test replaced rather than retuned

The crescent test asserted that the weave adds at least 4% of path. That proxy reads the wrong way
round: extra path comes MORE readily from a bead going where it should not than from one weaving,
and on that fixture the defect scored 10.5% against the fix's 2.2% - the assertion was satisfied
twice as well by the bug as by the fix. It is now a two-sided window on the same law, since the
surplus has a ceiling (`lobes * 4 * cap`) as well as a floor.

Also recorded there: the existing "no bead is displaced further than the wall allows" bound did
NOT catch this, and the numbers are written down so nobody re-derives it. A bead displaced 7 mm
from its own ring lands on top of a NEIGHBOURING ring, so its distance to the nearest unwoven bead
stays small - 0.555 mm with the walk against 0.592 mm with the search. That is also why the
printed defect was stacked material and an unwelded interface rather than a bead in mid-air.
