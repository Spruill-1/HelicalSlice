# TPU Continuous-Shell Test Suite

> **Printed and measured 2026-08-14 — see `RESULTS.md`.** C-lite validated
> (weights within 1%, T2/T1 = 1.188 vs 1.192); scarf verdict NEGATIVE (keep
> default 0); C-full design brief extracted.

Generated 2026-08-12 for flow **0.4 nozzle / 0.45 line width / 0.25 layer height**
(bead spacing 0.3963 mm). If the profile changes, regenerate — the wall thicknesses
are tuned to this spacing and test nothing at a different flow.

All tubes OD 32 mm. Print **one object at a time** (or "By object" sequencing —
continuous shell requires it for multiples). Suggested common settings:
`continuous_shell_mode` on, `wall_loops` 2, seam step 100% (default), seam jitter 0,
`continuous_shell_interlock` **0** for T1–T4 (isolates the width changes from Z stagger;
T5 can keep defaults). Everything else per your normal TPU profile.

Weights below assume TPU at 1.21 g/cm³ — scale linearly for your filament's density.
Weigh on a 0.1 g scale; weight vs target is the primary density/void metric.

| Tube | Wall (mm) | H | Target weight | Tests | What to measure |
|---|---|---|---|---|---|
| T1_tile_exact | 1.693 | 30 | **5.85 g** | Control: 4 rings tile exactly; the medial correction must NOT fire | Weight ≈ target; section a slice — no widened beads, no void ring |
| T2_deficit | 1.993 | 30 | **6.82 g** | Medial deficit v = 0.30 mm; innermost pair widens to ~0.60 mm | Weight ≈ target (the OLD code would print ~1.0 g light, ≈ −15%); cut a ring and look for the closed medial — no continuous void line between wall halves |
| T3_thin | 0.350 | 30 | **1.26 g** | Sub-bead fitted bead (prints at 0.35–0.34 width, not 0.45) | Weight ≈ target (old behavior over-extruded ~28% → ~1.6 g and rippled); wall smoothness, water-tightness |
| T4_taper | 2.60 → 1.50 | 40 | **9.32 g** | Ring-count transitions mid-print (crosses 6-ring and 4-ring bands) | Look for the transition bands on the INSIDE surface: should be smooth width blends, no gap ring / no bulge ring at any height |
| T5_seam | 2.485 | 60 | **16.73 g** | Scarf A/B/C: print **three copies** at `continuous_shell_scarf_length` = **0 / 6 / 12** | Seam: rake light across the outer wall — the helical groove should soften visibly at 6 and more at 12; weigh all three (12 mm scarf should cost ≲ 3% weight vs 0); flex each tube at the seam line — scarfed tubes should not preferentially crease there |

## Sliced and ready (2026-08-12)

`sliced/` contains all seven jobs as `.gcode.3mf`, ready for Bambu Connect:
H2D, **right** nozzle (filament_map=2, pinned manually — the CLI default was the left!),
0.4 hardened steel, **Generic TPU at 230C nozzle / 35C engineering plate** (stock generic
profile runs 240C, above this spool's 235 max — overridden in-spec, range set 190–235),
0.25 layer / 0.45 line width, process base `0.20mm Standard @BBL H2D`. Overlays in `config/`.

Slicer-predicted deposition (at the profile's 1.24 g/cm³ — rescale to your spool):

**PlateA_T1-T4.gcode.3mf** carries T1+T2+T3+T4 on one bed, sequenced **by object**
(each tube prints start-to-finish as its own continuous bead — verified in the G-code).

**T5_trio_fmb.gcode.3mf** carries all three scarf variants on ONE bed in a
front/middle/back column (x=175, y=32/160/288 — the sequential arranger's
clearance-safe pitch for 60 mm objects), scarf **0 = front, 6 = middle, 12 = back**
via per-object setting overrides. Verified per object block in the G-code:
0 / 9.74 / 18.70 m of tapered bead — identical to the three single prints, which
remain in `sliced/` as fallbacks. Print order is front→back.

All final `.gcode.3mf` files carry the injected thumbnail set
(plate/small/no_light/top/pick PNGs rendered from the plate layout) — headless CLI
output omits them and Bambu Connect refuses to load such files. (The four
single-tube fallbacks T1–T4 were never given thumbnails; run `inject_thumbnails.py`
on them first if you ever need one.)

**Nozzle flow declaration (2026-08-14):** Bambu Connect validates the declared
nozzle *flow type* against the printer ("The nozzle flow setting does not
match") — and it reads it from the **G-code CONFIG_BLOCK comments**
(`; nozzle_volume_type = ...`), not from the 3MF metadata: patching
slice_info.config alone, and even stripping the declaration entirely, was
still refused. The original slice emitted a **1-element vector**
(`; nozzle_volume_type = Standard`) for the two-extruder H2D, leaving the
right nozzle's flow type unresolvable, so Connect refused even though the
value matched the hardware (both of David's H2D hotends report 0.4 hardened
steel **Standard** flow; High Flow hotends are an optional upgrade). All nine
files were normalized in place with `tools/patch_nozzle_flow.py` to
`Standard,Standard`: 3MF metadata + the CONFIG_BLOCK comment line + recomputed
`.gcode.md5`. The diff against the originals is exactly that comment line plus
metadata — no motion commands changed. `slice_tpu_suite.py` now pins
`nozzle_volume_type = ["Standard", "Standard"]` so future re-slices declare it
correctly from the start.

Per-tube deposition predictions (at the profile's 1.24 g/cm³ — rescale to your spool;
weigh the PlateA tubes individually after cutting them apart):

| Print | Filament | Predicted weight | vs model volume |
|---|---|---|---|
| T1_tile_exact | 4.59 cm³ | 5.69 g | 95% |
| T2_deficit | 5.46 cm³ | 6.77 g | **97% — deficit closed (old code ≈85%)** |
| T3_thin | 1.04 cm³ | 1.29 g | 100% |
| T4_taper | 7.78 cm³ | 9.65 g | 101% |
| T5_seam_scarf0 | 13.34 cm³ | 16.54 g | 96% |
| T5_seam_scarf6 | 13.02 cm³ | 16.14 g | −2.4% vs scarf0 |
| T5_seam_scarf12 | 12.74 cm³ | 15.80 g | −4.5% vs scarf0 |
Verified in the G-code: per-point-Z continuous bead (~30-60k XYZE moves each), scarf
width ladder 0.225-0.4 present only in scarf6/12 (9.7 m vs 18.7 m of tapered bead — 2x,
as designed), T2's widened 0.61 medial pair, T3's 0.35 fitted bead, T4's graded
0.49-0.61 sweep. The 0.64-0.69 widths are the layer-0 medial fix at first-layer flow.

## C-full validation print (2026-08-15)

`sliced/T4_adaptive.gcode.3mf` re-slices T4_taper with
`continuous_shell_adaptive_beading` ON (same flow/temps/right nozzle as the
suite; nozzle flow declared Standard,Standard natively). Predicted 9.39 g /
~46 min. G-code verified: single-step ring flips 6->5 at z~16.2mm and 5->4 at
z~32.0mm with 0.45->0.30 width ramps across each ~3.5mm band — vs the printed
T4's abrupt 6->4 double flip at z~20.5mm. A/B against the printed T4: the two
transition bands and the thin section should lose the rough band and the
periodic medial voids (cut it at the same heights).

## Notes

- **T2 is the C-lite money shot**: same print settings as T1, only the wall differs.
  If both land near target weight, the deficit correction is doing its job.
- **T5's wall tiles exactly** (6 rings) so scarf effects are isolated from the medial
  correction. Label the three copies immediately — they look identical from the outside
  until you find the seam.
- The scarf floor is nozzle/4 = 0.10 mm bound by w/2, so tapers bottom out at ~0.225 mm
  width. If TPU die-swell makes the floor segments visibly gappy, that's a finding —
  the floor fraction is one line to change.
- Weight deltas < ±3% of target are within flow-calibration noise; the effects being
  tested are 15–28% where they appear.
- A wall below 0.10 mm (min feature size) now *refuses to slice* by design — not in the
  suite since it produces no print, but worth knowing if you experiment with thinner walls.
