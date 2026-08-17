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

## Interlock A/B on adaptive walls (2026-08-15) — NEGATIVE, default flipped to 0

`T4_adaptive_interlock` (same slice as T4_adaptive, XY-identical byte-verified,
plus `continuous_shell_interlock` 100) vs the plain adaptive tube:

- **Microscopy:** the stagger is geometrically active — corner voids become
  3-corner TRIANGULAR (T-junction) rows, the brick signature — but misfit
  concentrates in organized rows at the pinned/staggered ring boundaries, and
  the exterior gains a slight regular bumpiness (offset beads overlap instead
  of nesting; the crowding telegraphs through the pinned surface beads in TPU).
- **Tear test:** a tear could be STARTED and propagated on the interlocked tube
  (ductile, crack deflected across courses — not catastrophic). On the plain
  adaptive tube a tear could NOT be started at all (scalpel required).

Interpretation: the arc-length Z ramp is itself a bricklayer — radially adjacent
rings sit at path fractions that differ around the circumference, so neighbors
are offset in Z by a smoothly varying, NON-PERIODIC amount (graded,
incommensurate nesting; visible in the very first T4 cuts at interlock 0). The
explicit periodic +-h/4 stagger replaces that with coherent, repeating interface
geometry whose boundary void rows act as tear initiation sites. Natural beats
periodic. `continuous_shell_interlock` now defaults 0 (knob retained; n=1
caveat, but microscopy, tear initiation, and surface all point the same way).

## Phase B spiral: v1 double-arm print FAILED, v2 single-arm (2026-08-17)

T5_spiral v1 (double-arm Fermat interior), two prints, pink TPU: the SEAM GOAL was
met — the helical seam significantly shallower, barely registering on the
fingernail — but both tubes delaminated into helical "sheets" peeling from the
bore. Mechanism: in the double-arm interleave every radially adjacent interior
pair is an armA/armB interface deposited up to half the interior apart in time —
every intra-layer radial weld is a COLD JOINT. Ring mode's ~k junctions had a
hidden virtue all along: ring-adjacent beads land within a few revolutions.
**Thermal adjacency (when the neighbouring bead was laid) is a first-class design
constraint, co-equal with geometry.** v2 rebuilds the interior as a SINGLE arm at
one slot per revolution (radial neighbours = consecutive revolutions, one rev
apart — warmer than ring mode's worst case) with per-layer traversal reversal for
short alternating joins. Same seam benefit (2 junctions/layer), none of the
double-arm machinery. `sliced/T5_spiral.gcode.3mf` regenerated as v2.

## Phase B v2 print: seam goal held, reverse-join doubling struck (2026-08-17)

The single-arm v2 tube printed with the same excellent shallow seam, but the
cross-section showed LARGE staggered bowtie voids in a brick lattice. Mechanism:
v2 alternated the traversal direction per layer for short joins - and under the
arc-length Z ramp a reversed layer puts the same spatial bead at the OPPOSITE
ramp fraction, so consecutive layers' beads pair up with a two-layer void
between pairs, the pairing phase shifting with radius (the staggered lattice).
This is the codebase's own documented "reverse-join doubling" - the reason
annular ring mode never reverses - reintroduced by the spiral's reversal. v3
fixes the traversal direction (outer->bore every layer); the inter-layer join
becomes one radial in-material diagonal across the wall at the precessing seam
azimuth (single event per layer; total radial travel equal to ring mode's
per-ring hops combined). `sliced/T5_spiral.gcode.3mf` regenerated as v3.

## Phase B v3 print: peel persists - mechanism revised to laminate cleavage (2026-08-17)

Two decisive results arrived together: the v2 tube DID have the bore peel after
all (reported late), and the v3 reprint (fixed traversal, verified in the
printed G-code: 236/237 layers run outer->bore) shows the peel AND a void
pattern too. So the delamination survived fixing both cold arm-to-arm joints
(v2) and reverse-join doubling (v3) - neither was the root cause. Unrolling the
peeled wall showed the peel organized along a diagonal matching the precessing
seam helix (0.02 rad/layer is ~36 deg from horizontal at this radius): the
cleavage runs along helical interfaces and initiates at the seam column.

REVISED MECHANISM - radial interdigitation is the master strength variable:

- Ring mode's Fermat ORDER (0,2,4,...,5,3,1) scrambles the arc-length ramp
  fraction across radial neighbours: adjacent rings sit up to ~(k-1)/k of a
  layer apart in Z, varying with angle - the graded incommensurate nesting
  visible in the very first T4 cuts, the thing the interlock A/B tear test
  proved tear-proofs the wall. Those interfaces hold DESPITE being laid up to
  k revolutions apart in time (ring mode's outer interfaces are its coldest) -
  which is what falsifies the thermal-adjacency story.
- Any monotone spiral orders path position radially, so radial neighbours
  differ by only h/k (~0.04 mm): no interdigitation anywhere. Every interface
  is a smooth laminate plane; a crack starting at the seam column propagates
  freely around the wall. All three spiral revisions share this flaw; v1's
  cold joints and v2's doubling were real defects stacked on top, not the root.
- Thermal adjacency is demoted to a secondary weld-quality factor.

v4 response: the spiral arm now carries an intrinsic graded Z-wave,
z = (h/4)*sin(2*pi*g*v) with g the golden-ratio conjugate and v the winding
coordinate - radial-neighbour offsets up to ~0.47 h, graded and incommensurate
around the circumference so no repeating interface geometry appears (the
tear-initiation failure of the explicit periodic interlock), tapering to zero
at the tangent merges, surfaces and layer 0 pinned. Single-variable change vs
v3: XY byte-identical, only buried Z moves. Prediction if the revision is
right: the peel is gone (or dramatically harder to start) and the aligned void
pattern breaks into ring-mode-like incommensurate nesting. If v4 still peels,
the interior spiral is retired.

A/B caveat found while regenerating v4: the committed v3 file (the tube printed
above) was sliced by a STALE binary predating its own commit - it carries
~63 mm/layer of ~0.5 mm-wide bead on the inner wall (radii 14.0-14.7 on this
tube) and -3.5% total extrusion vs what the committed source produces. The v4
file is sliced from the committed source and verified strictly Z-only against a
pure-HEAD reslice: XY streams identical, E within rounding ulps, max |dZ| equal
to the wave amplitude. So the v4-vs-v3 print comparison carries a small
width-plan delta in addition to the wave; the qualitative verdict (peel gone or
not) is robust to it. v4 predicts 17.34 g (the stale v3 build predicted 16.75 g).

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
