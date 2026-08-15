# TPU Continuous-Shell Suite — Physical Results (2026-08-14)

Printed on the Bambu H2D, right nozzle (0.4 hardened steel, Standard flow),
generic 95A TPU @230C, engineering plate, from the committed `sliced/` jobs
(after the nozzle-flow declaration fix, commit `4ea6a94278`). Weights on a
0.1 g scale; predictions at the profile's 1.24 g/cm³. One specimen per
condition — treat every conclusion as n=1, flex results as qualitative.

## Weights: everything on target

| Tube | Measured | Predicted (G-code) | Delta |
|---|---|---|---|
| T1 control | 5.6 g | 5.63 g | −0.5% |
| T2 deficit | 6.65 g | 6.71 g | −0.9% |
| T3 thin (brim removed) | 1.1 g | 1.11 g | −0.9% |
| T4 taper | 9.55 g | 9.58 g | −0.3% |
| T5 scarf 0 | 16.5 g | 16.54 g | −0.2% |
| T5 scarf 6 | 16.0 g | 16.14 g | −0.9% |
| T5 scarf 12 | 15.6 g | 15.80 g | −1.3% |

- **T2/T1 = 1.188 vs target 1.192** (density-independent): the medial-deficit
  correction deposited the full corrected volume. Old code ≈ 15% light.
- T3 at 1.1 g with a smooth wall: sub-bead width fitting works (old behavior
  ≈ 1.6 g, rippled). T3's ~0.13 g auto-brim is NOT in the targets.
- Uniform −0.3…−0.9% across all tubes → flow/density at profile values.
- T5 ladder steps (−0.5 g, −0.9 g) match the designed 2× tapered-bead ratio.
  Scarfed tubes run ~0.5–1% leaner than predicted — consistent with TPU
  pressure lag during rapid width taper.

Predictions derived per object block from the plate G-code (net extrusion ×
filament cross-section; runs ~1% under slicer stats uniformly — use ratios).

## Microscopy

- **T2 cut ring, off-seam:** regular dense bead lattice, no continuous void
  channel — the widened innermost pair closed the medial. C-lite confirmed
  structurally as well as by weight.
- **T2 cut, on-seam (scarf 0, interlock 0):** disordered junction column with
  crossing interfaces. Looks alarming; per the T5 flex result below it is NOT
  a mechanical weak line in TPU.
- **T4 axial cut — squeeze consolidates, stretch voids:** the squeezed 6-ring
  band shows only pinprick corner voids; the stretched 4-ring band shows a
  PERIODIC VOID LINE along the medial (innermost-pair meeting closes
  imperfectly once per revolution at ~0.51 spacing) despite on-target weight —
  a packing defect, not underfill. The transition plane itself is a
  void-rich band. Beads visibly nest even at interlock 0 — the per-point-Z
  ramp gives radial neighbors angle-dependent height offsets for free.

## T4 G-code forensics: ring parity and the double flip

Ring counts only exist in PAIRS (the out-and-back radial path must close each
revolution): T1 and T2 both print 4 rings; T3's single bead is the lone odd
case. So the taper cannot do 6→5→4: it does a double flip **6→4 at
z≈20.5 mm**, exactly where wall = 5×spacing (1.98 mm, the equidistant point),
stepping per-ring spacing ~0.34 → ~0.51 within ~0.25 mm of height. Above the
flip the whole 4-ring band runs over-wide including the surface ring — the
visibly rougher thin section.

## T5 scarf verdict: NEGATIVE — keep default 0

- Backlit, the scarf is cosmetically backwards: scarf-0 is a thin sharp
  helical line; 6/12 show progressively wider dark textured bands (the
  stacked taper windows). Dotted texture in the bands + lean weights =
  taper floor `max(w/2, nozzle/4)` = 0.225 mm is below the 0.25 layer
  height — die swell, micro-gaps.
- Fingernail across the seam: all three similar; **scarf-0 slightly
  smoothest**. The ridge the scarf exists to remove is not a problem in TPU.
- Flex: **seam creasing increases 0→6→12; scarf-0 feels most isotropic.**
  The thinned taper band becomes the preferential crease line — the scarf
  manufactures the weak line it was meant to prevent.

Decision: `continuous_shell_scarf_length` stays **default 0**; not
recommended for TPU. If ever revisited (rigid materials, where a ridge is a
crack initiator), raise the taper floor to ≥0.30 mm first. The seam-quality
workstream moves to C-full/B.

## C-full A/B validation (2026-08-15) — PASSED, default flipped ON

`sliced/T4_adaptive.gcode.3mf` (T4 geometry, `continuous_shell_adaptive_beading` on)
printed at 9.2 g vs 9.39 predicted (−2%, within flow noise; slight elephant's foot
at the base — first-layer squish, unrelated to beading). Against the original T4:

- **Exterior:** no visible transition band anywhere (old tube: distinct rough band
  at the z≈20.5 double flip and rougher thin section). The user's verdict: "no
  noticeable seam between the thick/thin sections whatsoever."
- **Thin-section cut:** the periodic medial void line of the stretched regime is
  GONE — uniform pinprick corner-void array, same texture as the thick section.
- **Transition cut:** gradual packing change; no void-rich flip plane, no notch.

This confirms the mechanism chain end-to-end: squeeze-biased single-step flips +
width-birth transitions eliminate both failure modes the original campaign
documented. `continuous_shell_adaptive_beading` now defaults ON (C-lite remains
available by turning it off, and is still the automatic fallback for non-star /
sparse / solid sections).

## C-full design brief (from this campaign)

1. Retire ONE ring at a time via an angular wedge (locally odd counts) —
   halves the transition step the pair-parity forces today.
2. Pin the surface ring at nominal width; absorb all width variation in
   interior rings.
3. Spread transitions over real height — the current effective transition
   length is ~0.25 mm.
4. Bias ring-count flips toward SQUEEZE (hold more rings down to the
   min-width clamp) — squeezed walls consolidate, stretched walls void
   along the medial.
5. Spread/stagger the flip plane — today it concentrates voids at one height.

Strength prediction to test destructively: stretched-regime walls fail along
the medial void line; squeezed walls have no preferred plane.
