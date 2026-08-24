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

## Aligned ramp 12 PRINTED: "significantly better", residual peel IS the handoff band (2026-08-19)

The volume-parity tube is much better than every compensated/dispersed predecessor - no
general delamination - but a ~6 mm peel band remains immediately next to the seam. That
is the design's own handoff band, and the agreement is quantitative: measuring the
azimuths whose crossing set deviates from the modal spacing gives a disturbed band of
**7.0 mm at ramp 12**, against the ~6 mm peel observed. The metric predicts the defect.

Mechanism: with the handoffs ALIGNED, every column transits within the same azimuth
window, so that window is the one place where the entire wall thickness is simultaneously
out of register - a through-thickness cleavage plane with no column left in normal
register to arrest a crack. The band scales linearly with the ramp and costs nothing to
shrink (aligned handoffs are volumetrically exact at any ramp):

| ramp | disturbed band | weight | constant-radius run |
|---|---|---|---|
| 12% | 7.0 mm | 16.50 g | 1.01 turn |
| 8% | 4.6 mm | 16.51 g | 1.01 turn |
| **5%** | **3.1 mm** | **16.51 g** | 1.01 turn |

Why narrowing (rather than dispersing) is the right move: ring mode ALSO stacks its
junctions at essentially one azimuth and does not peel, so a stacked through-thickness
disturbance is tolerable when it is narrow - ring mode's own junction footprint is
~1-2 mm. Dispersing instead would cost material in proportion to the dispersion (the law
measured above) for a defect that ring mode shows is survivable when short.

`sliced/T5_ramp5.gcode.3mf` is the print candidate: 3.1 mm band, 16.51 g, columns still
holding a full turn. If the peel disappears, the ramp is viable and the band width is the
governing parameter. If a ~3 mm band still peels while ring mode's ~1-2 mm one never has,
then any spread handoff is a liability, the abrupt junction is what holds a
continuous-shell wall together, and the interior-spiral line should be retired rather
than tuned further.

## Ramp 5 PRINTED: THE PEEL IS GONE (2026-08-21)

`T5_ramp5` printed with no delamination anywhere - the first continuous-shell tube with
a spiral-derived interior that holds together. The section shows the regular triangular
corner-void lattice (the same T-junction packing ring mode produces), and the handoff
band reads as a clean diagonal line on the outside rather than a cleavage plane. Seven
revisions to isolate it; the chain that got here:

| rev | change | print result |
|---|---|---|
| v1 | double-arm spiral | peeled into two helical sheets |
| v2 | single arm (warm radial welds) | still peeled - falsified the thermal story |
| v3 | fixed traversal | still peeled - falsified reverse-join doubling |
| v4 | sine Z-stagger | NO PEEL - interdigitation confirmed, but void rows |
| v5 | golden staircase | mixed: pinpricks in part, splits elsewhere |
| v6 | ramp-tilt compensation | mostly solid, residual splits |
| ramp 12 | ring-like columns, aligned handoffs | much better; ~6 mm peel in the handoff band |
| **ramp 5** | **3.1 mm handoff band** | **no peel** |

Standing verdict on the mechanism: radial integrity comes from COLUMN TOPOLOGY (a column
that wraps past a half turn cannot be pulled out radially) plus a NARROW handoff band.
Interdigitation depth and thermal adjacency are secondary; both were falsified as primary
causes along the way.

Residual observation from the print: slightly less column-to-column bonding than wanted.
That is the Z-stagger's amplitude - at half a layer, adjacent columns overlap only half
their height. With ring-like columns the topology carries the anti-cleavage load that
half a layer was buying, so the amplitude drops to a QUARTER layer (the brick interlock's
proven envelope, and what the first peel-free spiral used). Measured on the re-slice:
adjacent-column offsets now cluster at 0.17-0.20 of a layer, mean overlap 72% of bead
height (was ~62%), still 16.51 g. One interface per ray still sits near 0.48 - a
systematic feature of the golden level sequence worth watching in the section.

Two real bugs fixed at the same time, both found by the suite once the ramped path became
the default: the arm left the outer ring at FULL stagger (a 0.12 mm vertical step at the
entry hop - the stagger now tapers to zero at the arm's own ends), and the ramp-tilt
compensation used the spiral's mid-slot radius per revolution, which mis-estimates the
rise for a ramped arm that sits ON its column's radius - leaving the interdigitation
depths off target, the very quantity that sets overlap.

`continuous_shell_junction_ramp` now DEFAULTS TO 5 rather than 100: defaulting to the
pure spiral meant defaulting to the one configuration print-proven to peel.
`sliced/T5_ramp5.gcode.3mf` is retained as the validated peel-free reference.

## Quarter-stagger PRINTED: outer columns better, one dedicated mid-wall void (2026-08-21)

Halving the stagger did improve the outer columns' bonding, and it exposed the defect the
re-slice measurement had already flagged and I had not chased: ONE interface per ray sat
at ~0.48 of a layer while the rest sat at 0.19. At a fixed radius that prints as a
continuous void line, which is what the section shows mid-wall.

Cause was the end-taper added in the same commit. It tapered the stagger envelope in the
SLOT coordinate v - but during a hold v is CONSTANT at that column's slot, so the taper
never eases in: it zeroed the stagger for the whole of the first and last buried column
while their neighbours kept theirs, leaving those interfaces at a mismatched offset.
**A quantity that is constant during a hold cannot be used to taper anything along a
ramped arm.** The taper now runs along the PATH (tau) over a quarter revolution, which
still removes the vertical step at the entry and exit hops without flattening a column.

Measured per-interface offsets after the fix (fraction of a layer, folded mod h):

| | interface 1 | interface 2 | interface 3 | mean overlap | weight |
|---|---|---|---|---|---|
| half stagger | 0.38-ish | 0.38-ish | 0.38-ish | ~62% | 16.51 g |
| quarter, v-taper | 0.19 | 0.19 | **0.48** | 72% | 16.51 g |
| **quarter, path-taper** | **0.19** | **0.19** | **0.31** | **77%** | 16.51 g |

0.19 / 0.19 / 0.31 are exactly the golden gaps the level sequence is designed to produce,
so the interdigitation is now on target at every interface with no outlier.

## Path-taper PRINTED, and the staircase turns out to be redundant (2026-08-21)

The void line is gone: the section is a uniform corner-void lattice with no continuous
channel at any radius. The voids are distinctly TRIANGULAR though - the T-junction
signature of deliberate staggering, and the same signature the interlock A/B produced on
the tube where a tear could be started (the plain one needed a scalpel). That prompted
measuring how much of the stagger is actually ours:

**The arc-length ramp already offsets radially adjacent columns by its own rise per
revolution - measured h x (one revolution / total layer path) = 0.25 x 93/558 = 0.17 of a
layer on the T5 tube.** Against designed offsets of 0.19 / 0.19 / 0.31, the imposed
staircase was therefore contributing +0.02 / +0.02 / +0.14: nothing on two interfaces,
and on the third an outlier that is the LEAST overlapping interface in the wall.

The staircase exists because a pure spiral is one helical thread with nothing else
resisting radial cleavage. A ramped arm's ring-like columns are topologically blocked
already, so for a ramped arm the staircase now rides the existing
`continuous_shell_interlock` knob - default 0, i.e. keep the natural nesting, the same
verdict that knob already carries for ring mode. A pure spiral still staggers
unconditionally. Measured per-interface offsets (fraction of a layer, folded mod h):

| | interface 1 | 2 | 3 | mean overlap | weight |
|---|---|---|---|---|---|
| half stagger | 0.38 | 0.38 | 0.38 | ~62% | 16.51 g |
| quarter, slot-taper | 0.19 | 0.19 | 0.48 | 72% | 16.51 g |
| quarter, path-taper | 0.19 | 0.19 | 0.31 | 77% | 16.51 g |
| **natural nesting (interlock 0)** | **0.16** | **0.17** | **0.17** | **83%** | 16.51 g |

Uniform across every interface, no outlier, the most bead-to-bead contact of any
configuration tried, and it deletes a knob's worth of machinery from the default path.
`sliced/T5_ramp5c.gcode.3mf` is the candidate. Prediction: smaller and more uniform
voids, less pronounced T-junction rows, no return of the peel (topology carries that).

## Natural nesting PRINTED: TEAR TEST PASSED - strength gate cleared (2026-08-22)

The natural-nesting tube (ramp 5, interlock 0) **could not be torn without a scalpel** -
exactly the result the plain adaptive RING tube gave, and the result the interlocked ring
tube FAILED (a tear could be started and propagated on that one). No peel, no
delamination, uniform section. On the one test that has ever discriminated between
configurations in this campaign, the interior spiral now matches ring mode.

State of the two modes on the same tube, same profile, same printer:

| | ring mode | ramp 5, natural nesting |
|---|---|---|
| predicted weight | 16.55 g | 16.51 g |
| print time | 75 min | 75 min |
| junctions per layer | 6 | **3** |
| tear initiation | scalpel required | **scalpel required** |
| radial peel | none | none |

Caveat worth keeping: "cannot start a tear" is a CEILING result. Both configurations sit
at the ceiling, so this shows no difference rather than proving equality - ranking them
would need a starter notch and a propagation measurement, or a pull test. n=1 as always.

That leaves exactly one question unanswered, and it is the one the whole workstream
exists for: **is the seam actually better than ring mode's?** Never A/B'd on a finished
tube. `sliced/T5_ring.gcode.3mf` is now the control - same geometry, same profile, only
`continuous_shell_interior_spiral` off - so the comparison is a single print away.
Protocol: backlit photograph of both seams, fingernail across the seam, and flex each
tube with the seam on the outside of the bend (the test that exposed the scarf as a
crease initiator). If the seam is not visibly better, the feature has no justification
left and should stay off regardless of the strength parity.

## Non-tube shapes: slicer validation (2026-08-22)

Everything in the campaign had been validated on ONE circular tube, which is the
friendliest case for a ramped arm: constant radius, constant wall, one bore, no corners,
no solid sections. Four shapes each break one of those (`tools/gen_shape_suite.py`,
all 30 mm tall, same 2.485 mm wall as T5):

- **S1_ellipse** radius varies 20 -> 13 around the section
- **S2_rsquare** superellipse, four corners (the case that historically under-counted
  max beads on square annuli)
- **S3_lobed** three lobes, concave-ish flanks, still star-shaped
- **S4_capped** circular tube closing with a solid cap at 22 mm, so the plan must
  decline mid-print and hand back to ring mode

Sliced both ways and compared against ring mode on the same geometry, with the same
72-ray metric applied to the validated tube as a baseline:

| shape | weight ring / spiral | max radial gap ring / spiral | disturbed rays (spiral) |
|---|---|---|---|
| T5 tube (validated) | 16.55 / 16.51 g | 0.41 / 0.78 mm | 3 / 72 |
| S1 ellipse | 8.67 / 8.66 g | 0.52 / 0.79 mm | 2 / 72 |
| S2 rounded square | 9.33 / 9.31 g | 0.50 / 0.80 mm | 2 / 72 |
| S3 three-lobed | 8.06 / 8.04 g | 0.50 / 0.80 mm | 3 / 72 |
| S4 capped | 14.08 / 14.07 g | 0.41 / 0.76 mm | 3 / 72 |

Every shape behaves like the tube: volume parity within 0.02 g, and a handoff band of
2-3 rays in 72 (~4% of the circumference) with the same ~0.8 mm worst gap the printed,
peel-free, untearable tube already carries. No shape hard-failed or declined
unexpectedly. **S4's 30 solid layers above the cap are move-for-move IDENTICAL between
ring and spiral mode** - the plan declines on solid sections and hands back cleanly,
which is the mid-print handover working exactly as designed.

`sliced/PlateS_shapes.gcode.3mf` (built by `tools/make_shape_plate.py`) prints all four
shapes in spiral mode PLUS the legacy 60 mm T5 tube in RING mode as the seam control -
56.4 g, 243 min, five objects printed by object. Running the control on the same plate
puts the seam A/B on one filament, one session and one machine state, and its 60 mm
height matches the spiral tube already in hand.

Two things the script asserts rather than assumes, because both are silent failure modes:

- **The 60 mm control must print LAST.** It exceeds the H2D's 47.4 mm rod clearance
  (`extruder_clearance_height_to_rod`), so it is only safe alongside the 30 mm shapes if
  nothing is printed after it. The arranger does put it last - the script fails the build
  if it ever stops doing so, since getting this wrong drives the gantry through a printed
  part. (The old suite note that 60 mm objects "exceed the by-object height limit" is too
  strong: they are fine as the final object, and OrcaSlicer's own sequential validation
  passes.)
- **The per-object override must actually apply.** Measured on the sliced plate: the four
  shapes run 3.0 junctions per layer (ramped spiral) and the control runs 6.0 (ring mode),
  which is the override working rather than silently inheriting the plate's spiral setting.

`sliced/T5_ring.gcode.3mf` remains as a standalone control (16.55 g, 75 min) if only the
seam comparison is wanted without the 4-hour plate.

## Non-tube shapes PRINTED: cross-sections are shape-INDEPENDENT (2026-08-23)

(A power outage killed the first plate after three objects; the two survivors were
re-plated with `make_shape_plate.py --out PlateS_resume S4_capped.stl`, identifying the
lost pair by matching G-code segments to the project's build-item transforms.)

**The three non-circular shapes came out with practically identical cross-sections** -
the same regular corner-void lattice, indistinguishable between an ellipse whose radius
swings 20 -> 13, a superellipse with four corners, and a three-lobed section with concave
flanks. That is a stronger result than the slicer metric predicted. The metric said each
shape carries the same 2-3 disturbed rays in 72 as the tube; the prints say the packing
is shape-INDEPENDENT, including through the corners where the per-angle bead plan varies
width fastest. Whatever the adaptive plan and the ramped arm are doing, section geometry
does not perturb it.

Taken with the tear and peel results on the tube, the interior spiral is now validated on
geometry as well as strength: constant-radius circular, radius-varying, cornered and
concave. Annulus-to-solid is NOT yet validated - see below.

## S4_capped v1 was a BAD MODEL, not a slicer result (2026-08-23)

The capped shape printed with a clean outside and a mess inside: material drooping into
the bore, blobs around the rim. That is a modelling error of mine, not a continuous-shell
finding. v1 closed the bore with a FLAT HORIZONTAL DISC at z=22, over a bore 27.0 mm
across - so the first solid layer was a 27 mm unsupported bridge printed into open air,
with no fillet to ease into it, and in continuous-shell mode no bridging flow, no anchor
infill and no support to help. It could only droop.

It also means the test measured the wrong thing. The intent was "does the plan decline
cleanly at annulus-to-solid and hand back to ring mode" - which the slicer answered
(the solid layers are move-for-move identical between ring and spiral mode) - but the
print was dominated by a bridging failure that says nothing about the transition.

v2 closes the bore with a 45-degree CONE instead: cylindrical bore to z=13, then the
radius sheds exactly one layer height per layer up to an apex at z=26.5, solid above.
Every layer is self-supporting, and the section now walks the wall from 2.06 mm to
15.4 mm before going solid - so it exercises the plan across nearly its whole range
rather than stepping over it, which makes it a considerably better test than v1 would
have been even if v1 had printed.

`sliced/S4_capped_v2.gcode.3mf` (17.2 g, 78 min) is the reprint. Bore radius traced in
the sliced G-code: 13.72 mm held to z=12.5, then 11.91 / 9.40 / 6.90 / 4.40 / 1.90 at
2.5 mm intervals, closing just past z=26 - a clean 45-degree taper with no bridge.

## Sparse core + continuous shell: S5, and a mode-coupling finding (2026-08-23)

S4 v2's coned closure printed well (outside unchanged, inside clean and tapering to a
smooth point), so S5 generalises it into a part that stacks three regimes in one object
to test the sparse "woven" infill against the newer modes:

    z 0 .. 10     annular wall, 5 beads (2.035 mm)   -> stays dense
    z 10 .. 24    45-degree cone closing             -> wall thickens
    z 24 .. 36    solid column, 27.9 mm core         -> sparse woven core
    top layers    lid over the sparse core

Wall thickness is the whole trick: the sparsifier only thins a chain with buried rings
past its shells (annular needs n > 2*wall_loops+1, i.e. six beads at wall_loops 2), so
2.035 mm is the thickest wall that stays fully dense. No modifiers and no height ranges -
continuous shell rejects objects whose extra regions could split a layer, and that guard
is right for modifier volumes even if it also blocks height ranges.

Sliced at density 25% against a density-0 reference, the core does exactly what it should:
**18.3 -> 7.2 beads per ray** in the solid column. But the annular wall also changed,
5.0 -> 4.0 beads per ray, and that is NOT sparsification - decimation never fires there
(`thick` is false at five beads). Isolating it:

| config | wall beads/ray at z=5 |
|---|---|
| adaptive ON, density 0% | 5.0 |
| adaptive OFF, density 0% | 4.0 |
| adaptive ON, density 25% | **4.0** (identical to adaptive OFF) |

**Enabling `continuous_shell_infill_density` silently disables ADAPTIVE BEADING (C-full)
for the whole object**, because the plan is gated on `!woven.enabled` - and the interior
spiral goes with it, for the same reason. So switching on sparse infill quietly regresses
every wall in the part to the pre-C-full marching path, the one whose printed failure
modes (periodic medial void line, rough transition band) this campaign spent its first
half fixing.

Whether the spiral should compose with sparse infill is arguable - sparse fill decimates
the very buried rings the spiral reorganises, so they are plausibly alternatives for the
same material. Adaptive beading is NOT arguable: it decides bead counts and widths, which
is orthogonal to whether interior rings are later decimated, and losing it is a straight
regression. Worth decoupling, but that is its own change with its own validation, not a
side effect of this test.

`sliced/S5_sparse.gcode.3mf` (18.9 g, 89 min) is the print. Read its wall as a C-lite
wall, not as one of the validated C-full/spiral walls - it is not comparable to the tubes.

### S5 printed (2026-08-23): the regression is now photographed, and it owns every defect

Printed 18.9 g / 89 min. Outside decent; the three zones are visibly distinct (they run at
very different path counts per layer, so layer times differ). Two defects, both reported by
eye and both confirmed in the G-code:

**1. A void up the middle of the 2.035 mm base wall.** C-lite lays FOUR beads there and
stretches the interior pair to 0.555 mm - 1.40x nominal, far past the ~1.1x stretch cap the
T4 microscopy set - and even then they do not close:

| | beads | radii (mm) | widths | flank overlap | medial gap |
|---|---|---|---|---|---|
| density 0 (C-full) | 5 | 14.187 / 14.583 / 14.979 / 15.376 / 15.773 | all 0.396-0.398 | normal | none |
| density 25 (C-lite) | 4 | 14.187 / 14.663 / 15.296 / 15.773 | 0.396 / **0.555** / **0.555** / 0.396 | **zero at every interface** | **0.078 mm** |

Beads laid exactly tangent, nowhere overlapping, with a residual gap at the medial axis:
the missing fifth bead, printed.

**2. Dangling filament inside the cone.** Not stringing - travels never enter the bore
(measured: zero bore-crossing travel moves over the whole cone). It is a DOUBLED BEAD on the
overhanging bore lip: C-lite puts two full-width (0.396) beads 0.203 mm apart, ~33% local
over-extrusion, on every cone layer, exactly where 63% of the lip hangs over air. C-full
spaces the same lip at a uniform 0.398 with no doubling.

Both defects reproduce EXACTLY at density 0 with adaptive beading off (same 4-bead stretched
wall, same 0.203 mm doubled lip). So sparse infill causes neither directly: **every wall
defect in this part is the adaptive-beading regression above.** The sparsifier's own job was
done correctly - a plain thick tube (OD40/ID8) at the same density decimates cleanly across
its entire buried span.

### Why the cone wedge printed solid, and the lever that sparsifies it

The wedge of material forming the cone stays fully dense, growing from 2.4 mm (z=13) to
7.3 mm (z=22) outward from the bore. It is the bottom-shell-over-overhang region, which takes
CS's `force_solid` branch: stride 1 (nothing decimated) plus the brick stagger - confirmed by
consecutive cone layers alternating +-0.20 mm on top of the 0.25 mm cone step, which only
that branch applies.

The lever is `bottom_shell_layers`, and it is **binary, not graded**:

| bottom_shell_layers | dense wedge at z=13 / z=18 / z=22 |
|---|---|
| 3 (default) | 2.37 / 4.87 / 7.33 mm |
| 2 | 2.37 / 4.87 / 7.33 mm |
| 1 | 2.37 / 4.87 / 7.33 mm |
| **0** | **0.21 / 0.22 / 0.21 mm** |

`detect_floating_vertical_shell` has no effect on it.

`sliced/S5_conesparse.gcode.3mf` (**11.97 g, 57 min**) is the follow-up print: identical
geometry and density, `bottom_shell_layers = 0`, single variable. Cone beads/ray fall
11.5 -> 5.4 (z=13), 21.8 -> 8.8 (z=18), 26.0 -> 8.0 (z=22); the base wall and first layer are
untouched. Open question it asks: the cone's overhanging lip now has sparse rings rather than
solid behind it, so the underside may get better (less over-extrusion feeding the droop) or
worse (less backing). The base-wall void and the doubled lip WILL recur - they need the
decoupling, not a config change.

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

### PlateE printed (2026-08-23): the pinch fix is validated; the MILD case is now the worst

E1_ecc (bore offset 2.0 mm, 4/8 mm wall - the case that pinches) and E2_mild (offset
0.5 mm, 5.5/6.5 mm - never pinches), 30 mm tall, sliced density 0, spiral off.

| | predicted medial gap | material | printed |
|---|---|---|---|
| E1 (pinches, fix applies) | 0.82 mm | 101% | near-perfect bricklaying, extremely few internal voids on either the thick or the thin side |
| E2 (never pinches) | 1.29 mm | 94% | the void down the middle of the thick side; otherwise essentially fully solid, barely any voids under the scope |

Both predictions held, including the counter-intuitive ordering.

**The risk the print existed to test did not materialise.** E1 is where the marching fix
(c17d5e36fd) recovered 31% of the part's material by handing the pinched remainder on as a
solid region - so a SEPARATE CHAIN is now joined at a place that used to be a void. That
interface could have been a weak plane. It is invisible in the cut: no delamination, no
void row, bricklaying continuous across it.

**The fix left the extreme case cleaner than the mild one.** E2 never pinches, so it never
reaches the new code path: its band vanishes normally, the scalar deficit correction fires,
and ~2 beads of gap survive at the medial. A 0.5 mm bore offset is a very ordinary
geometry, so this is now the worst case in the family and the target for W1b - which is a
statement about the DEFICIT CORRECTION on a vanishing band, not only about widening
C-full's scope.

### PlateF printed (2026-08-24): W1b/W1 validated, and the strength result of the campaign

E2_mild re-sliced from the SAME STL as the PlateE part, plus a 5-lobe tube, both density 0.

**Strength: the parts were extremely difficult to cut even with a scalpel** - David's
comparison, across the many prints of this campaign, is that this is an enormous improvement
by itself. Against the north star (CS is a strength optimisation) that outranks everything
else on this plate. Two changes are candidates and this print cannot separate them: filling
walls that were previously 6-31% short (W1b), and putting the material where it tiles rather
than leaving overlap-plus-gap pairs (W1). A tear/cut comparison against a concentric tube,
which neither change touches, would separate them.

**The void ring is gone.** The eccentric tube reads "almost identical internally to the star
shape" - solid, no medial line - against the printed PlateE part which carries one. Predicted
from the G-code: rays with a gap over 0.55 mm went 48/72 -> 1/72, median 0.811 -> 0.399 mm
(one bead spacing).

**The lobed wall prints like the round one.** "Incredibly solid, almost no voids at all."
Per-angle bead counts and widths on a wall that varies 5.2-6.8 mm produce the same internal
structure as a wall that does not vary - which is what C-full is for, now shown in plastic.

**One defect found: a column of voids near the outer edge of the lobed part**, and it was
already visible in the G-code as beads below the 0.30 mm transitional floor. The floor
protects the transitional bead, but the REBATE that pays for it could push a recipient under.
Measured on this geometry: a 0.05 mm width ladder running down to 0.246 mm, 164 mm of
extrusion per layer (15.9%) below the floor, concentrated in a 100-200 degree sector - the
eccentric tube had none anywhere, which is why only the star showed the column.

Fixed by paying the floor only out of beads that can afford it (see the commit). Sub-floor
extrusion falls to 70 mm (6.75%), which is almost exactly one ring's worth of a 15-ring wall:
every STRUCTURAL bead is now at or above the floor, and what remains is the transitional ring
itself, in windows where the interior has no headroom to lift it. That residue wants the count
law, not the rebate, and is recorded as follow-up rather than forced here.
