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
