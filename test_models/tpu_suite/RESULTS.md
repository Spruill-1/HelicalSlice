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

## Phase B v4 print: NO PEEL - mechanism confirmed; sine trades peel for void rows (2026-08-17)

The v4 tube did not delaminate: the first spiral tube that holds together, on a
Z-only change. The interdigitation mechanism is confirmed by construction -
three failures and one fix all isolate the same variable. Seam stays minimal.

But the section shows MORE voids than ring mode, in organized rows, and less
visible inter-layer stagger. Diagnosis: a sinusoid makes the interdigitation
depth OSCILLATE along each interface (period 1/g rev), so every interface
sweeps through the void-producing intermediate and coplanar configurations
somewhere on its circumference. Ring mode never does this: both radial
neighbours' ramp fractions advance together, so each interface holds ONE
constant depth all the way around, a different depth per interface. Constant
along the interface + scattered across interfaces is the full property set;
v4 delivered only the second half.

v5 copies ring mode's structure exactly: each buried revolution holds one
constant Z level from the golden low-discrepancy sequence in [-h/2, +h/2]
(amplitude raised from h/4 to h/2 = ring mode's own bead-Z spread), stepped
over a 0.08-revolution smoothstep arc just past the seam azimuth - where ring
mode's own order boundaries live. Golden gaps make adjacent-interface depths
always 0.38 h or 0.62 h: never coplanar, mean ~ ring parity. Prediction: void
rows collapse to ring-mode pinpricks, no peel, seam unchanged. Still Z-only
(XY identical to v3/v4 slices).

## Phase B v5 print: mixed - pinpricks return, but XY-plane splits at weak interfaces (2026-08-17)

The staircase tube shows ring-mode pinprick packing in part of the section - the
constant-depth-per-interface principle is right - but the wall still SPLIT along
some radial interfaces (two sections separated in the cut, triangular void rows
on the split faces). Diagnosis is quantitative: the v5 levels were designed on
the staircase ALONE, but the physical offset between adjacent revolutions is
staircase PLUS the arc-length ramp's own tilt (~0.13 h per revolution - the
inner neighbour is one revolution later in path). Adding that tilt to the
designed -0.38h/+0.62h golden gaps gives -0.25h/+0.75h: EVERY interface at an
effective ~0.25 h from coplanar - uniformly T-junction grade, the explicit
interlock's documented tear geometry. Secondary defects: simultaneous level
steps at the seam let opposite-sign gap pairs sweep through zero mid-step (a
short coplanar ribbon at the seam column), and the envelope taper left the
surface-adjacent interfaces barely staggered.

v6 designs the TOTAL offset: golden targets, accumulated ramp tilt subtracted
(estimated per revolution from plan radii), wrapped mod h (offsets are
nest-equivalent mod h) - every physical gap a true golden gap. Envelope
tightened to a quarter revolution of the merges. Steps stay simultaneous at
the seam: scattering them was worked through and rejected (on a single
continuous arm, a revolution holds its predecessor's level until its own step,
so scattered steps make whole arcs of every interface exactly coplanar); the
residual ribbon is confined to the 0.08-rev window on the already-sacrificial
seam column, far milder than ring mode's k-junction column.

Verified in the emitted G-code: XY identical to the whole v3-v6 chain (Z-only
throughout), and 98% of physically-measured adjacent-interface gaps folded mod
h land at 0.096-0.098 mm = 0.385 h, the designed golden depth (the v5 file
measures ~0.065 there). Prediction: pinpricks everywhere, no splits, seam
unchanged. This is the staircase concept at full strength - if interfaces
still split, the deficit is not geometric stagger depth and the spiral's
remaining gap to ring mode is something else (deposition order/thermal).

## Phase B v6 print: mostly solid, residual splits - the limit is TOPOLOGY (2026-08-17)

The ramp-compensated staircase printed mostly solid with pinprick packing - the
best spiral wall of the campaign - but some radial interfaces still split. The
stagger is now provably at its designed depth (98% of gaps measured 0.385 h in
the emitted G-code), so the geometric variable is exhausted. It bought the
progression full peel -> no peel -> mixed -> mostly solid and then stopped.

WHY IT STOPS - the columns' topology, not their nesting:

- In ring mode each bead column is a CLOSED CIRCLE joined to its neighbours by
  ~k-1 short, steep radial hops (measured: 6 per layer on this 7-bead wall).
  A closed circle cannot be pulled out radially without breaking it, and the
  hops staple the interfaces on top of that.
- The spiral's interior is ONE HELICAL THREAD. Its "columns" are turns of that
  thread, and a helix CAN unwind: a crack along an interface separates a
  helical strip without breaking anything. That is exactly the "sheet peeling
  from the bore" of v1-v3, and the peel front bounded by the seam helix in the
  unrolled photo is the single ligament line (1 steep crossing per layer, at
  the layer join) holding the strip on.
- Z-stagger raises the interfacial bond and mechanical interlock - which is why
  v4-v6 improved so much - but it cannot remove a failure MODE that needs no
  material to break. No stagger depth closes a topological gap.

This unifies all six prints under one mechanism, including v1's two interleaved
sheets (the double arm's A/B interface had ~2 ligaments per layer, the fewest
of any revision, and failed the worst).

Also measured while looking for a ring-mode seam route: `connector_spread` is
nearly INERT on a thin tube. At 100% the junction azimuths disperse only 0.07
(circular dispersion; 0.00 at spread 0), because the radius-aware clamp must
stop each hop's chord from crossing the bore void. Dispersing junctions on an
annulus requires migrating radially ALONG the path, which is what the spiral
does - the two knobs are the ends of one continuum.

Standing recommendation: `continuous_shell_interior_spiral` stays default OFF
and is documented as a seam-only option that trades radial peel strength. The
open design between the ends of the continuum is a RAMPED junction - migrate
each ring-to-ring transition over an arc (spreading the overfill ridge and
keeping hops out of the bore) while each column still spans >180 degrees, so
radial extraction stays blocked. That is a new feature, not a spiral revision;
it is not implemented.

## Ramped junction implemented - the middle of the continuum (2026-08-18)

`continuous_shell_junction_ramp` (percent, default 100 = the pure spiral, so nothing
changes by default) sets how much of each revolution the interior arm spends
migrating between bead columns. Below 100 the arm HOLDS each column at a constant
radius for (1 - ramp) of a turn and concentrates the migration into one short ramp.
Measured on the T5 tube at ramp 20 against ramp 100:

| | pure spiral (100) | ramped (20) | ring mode |
|---|---|---|---|
| longest constant-radius run | 0.45 turn | **1.10 turn** | closed ring |
| migration concentrated in | 52% of arm | 42% of arm | n/a (hops) |
| crossings per ray | 6-7 (wedge) | 6 (13/16 rays) | 6 |
| worst radial gap | 0.55 mm | 0.54 mm | - |
| predicted weight | 17.34 g | 15.67 g | 16.55 g |

Two findings came out of building it, both measured, both non-obvious:

1. **A junction spread over an arc costs REDUNDANT material, linearly in the arc.**
   The transit arc lies in the interstice between two columns, where the neighbouring
   beads already reach - unlike ring mode's hop, which is ~1 pitch long. Uncompensated,
   ramp 20 came out at 18.87 g, +14% against ring mode. The fix is to pay for it in
   width: the arm is scaled (uniformly, preserving the plan's per-angle variation,
   floored at the transitional width) to the deposit the same plan would lay over those
   columns as closed rings. This is the same trade the scarf makes, but at ~12% instead
   of the scarf's 50% thinning, and spread over an arc rather than a junction.
2. **The arm must still sweep slot 0 to slot S.** Dropping the "extra" revolution to
   avoid retracing the surface rings looked right on the scale (16.40 g, near-perfect)
   and was actually a VOID: 5 crossings per ray instead of 6, worst gap 0.91 mm. Weight
   alone cannot distinguish correct fill from underfill - always check the ray crossings.

Ramp azimuths advance monotonically across the wall so the transition arcs tile the
circle; non-monotone scattering would let two revolutions occupy one slot at an azimuth
while the next slot went empty. The bulk of a ramped wall is measurably MORE regular
than a spiral's: uniform 0.40/0.41 mm gaps everywhere except the one column in transit,
where the spiral carries its boundary wedge at every azimuth (inner gap sweeping
0.33-0.53, outer 0.17-0.40).

## Ramp 20 (width-compensated) PRINTED: peel WORSE - never buy volume with contact (2026-08-18)

The compensated tube peeled worse than any spiral revision, and specifically at the
BORE-ADJACENT interface. Cause is the compensation itself: paying for the ramp's extra
material by scaling arm widths thinned every buried bead ~15% (0.40 -> ~0.34 mm), which
shrinks the contact between neighbouring columns - and the bore-side bead is thinned
TWICE, since the merge clearance already tapers the arm there. Interface contact is the
quantity this entire workstream is trying to buy; trading it for volume conservation is
backwards. **Design rule: never pay for a path-length surplus by thinning buried beads.**

Then the ramp sweep falsified the model the compensation rested on. Uncompensated
weights at ramp 6 / 12 / 20: 19.11 / 19.01 / 18.87 g - all ~+15% against ring mode and
FLAT in the ramp fraction. So the surplus is not the transit arcs (which would scale with
the arc); it is one structural extra revolution. Two facts explain it:

- The arm sweeps slot 0 to slot S, i.e. S revolutions of path serving S-1 buried columns,
  plus the two surface rings. The pure spiral absorbs that surplus by thinning both tips
  over a full revolution each (its merge clearance), landing at +4.8%; a ramped arm
  crosses the clearance zones quickly, so it pays the whole +15%.
- DISPERSING the handoffs is what costs the revolution: a column held from one transition
  to the next spans 1 + (delta azimuth) turns, so dispersion over the full circle buys
  exactly one extra revolution of material. Ring mode stacks its junctions at one azimuth
  and is volumetrically exact; the spiral disperses fully and pays. **Dispersing junctions
  on an annulus is not free - it costs material in proportion to the dispersion.**

Final configuration (aligned handoffs, buried columns only, no width compensation):

| | ring mode | aligned ramp 12 | ramp 20 dispersed | pure spiral |
|---|---|---|---|---|
| predicted weight | 16.55 g | **16.50 g** | 18.87 g | 17.34 g |
| constant-radius run | closed ring | **1.01 turn** | 1.10 turn | 0.45 turn |
| bead widths | planned | **planned** | planned | planned |
| rays at full crossing count | all | 15/16 | 13/16 | 13/16 (wedge) |

That is volumetric parity with ring mode at full planned widths, with each column holding
a constant radius for a full turn, and one narrow azimuth band - the seam, already the
sacrificial line - where the columns hand off. What the ramp buys over ring mode is a
junction spread into a shallow ramp instead of an abrupt radial jump.

`sliced/T5_ramp12.gcode.3mf` is the print candidate. Predicted 16.50 g; a tube that
weighs right AND peels would mean the ramp itself is not the answer either.

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
