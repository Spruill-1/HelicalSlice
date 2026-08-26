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

### PlateG printed (2026-08-24): the strength claim, corrected

**The broad claim from PlateF is WITHDRAWN.** G1_conc is a plain concentric tube - geometry
none of the wall-filling work touches, byte-identical in the density-0 fence across every one
of those commits - and it is ALSO difficult to cut. So most of that scalpel resistance is not
from W1b/W1: it belongs to machinery that predates them (the continuous bead, the arc-length
Z ramp, C-full on uniform walls) or to something this campaign never varied. Stating it as a
result of this work was wrong, and the control is what caught it.

**What survives is a differential, and it has a measured mechanism.** The control cut more
easily than the varying-wall parts and showed noticeably LESS bricklaying. Two candidate
mechanisms were measured and BOTH ruled out:

- Z offset between radial neighbours (the arc-length ramp's nesting): identical. Median 0.244
  of a layer on the control vs 0.248 on the lobed part; well-nested pairs 47.5% vs 49.3%.
- Radial stagger of a ring between passes: near zero on both. Median scatter 0.1 um and 0.8 um;
  ~97% of ring-columns tighter than 20 um on each. Rings stack in columns in both parts.

What does differ is the RING INTERFACE. Per-angle bead widths make the boundary between two
rings wander with azimuth instead of forming a cylinder, and the wander tracks wall variation
rather than external shape (sd, micron, by ring depth below the surface):

| part | outside | wall | r1 | r3 | r5 | r7 |
|---|---|---|---|---|---|---|
| G1_conc | round | uniform | 13 | 56 | 82 | 102 |
| E2_mild | round | varying | 13 | 70 | 121 | 176 |
| F2_lobed | lobed | varying | 24 | 109 | 196 | 285 |

E2_mild is the discriminator: a round outside like the control but a varying wall, and it sits
between them. So this is not the lobed profile resisting a blade. On a uniform wall every
interface is close to a cylindrical parting surface a crack can follow; on a varying wall the
same interface wanders by up to 285 um - three quarters of a bead - and there is no continuous
surface to part along.

**Honest limits**: three parts, one filament, and a subjective cutting assessment with no
instrumented force. This is a mechanism supported by a measurement and a consistent ordering,
not a proven causal result. It is a good candidate for the strength campaign's instrumented
work ([[cs-paper-plan]]), where the wander is a knob that can be dialled independently of
geometry.

**The floor guard is validated, and did better than predicted.** The lobed part was expected to
show a REDUCED void column; it came off "basically perfect as one continuous solid". The column
traced to beads below the 0.30 mm transitional floor is gone at the resolution of a cut.

### PlateD printed (2026-08-24, pre-W3 slice): one fix justified, one fix not working

Material note: this campaign prints in TRANSLUCENT MAGENTA TPU. An earlier memory entry said
matte black PETG and was carried forward stale; TPU is why the 0.30 mm w_term floor exists at
all (a TPU bead die-swells and gaps below ~1.2x layer height).

**D1_ladder, the 5-bead step: the medial void is REAL.** Microscope section of the 21-28 mm
band shows a continuous void running along the middle of the ~2 mm wall - the C-lite stretch
defect, four beads laid where five belong. This is the printed justification for W3, which
takes that step from 4.0 to 5.0 beads per ray. A post-W3 reprint is a clean one-variable A/B.

**D2_weave: W2's weave clamp is NOT holding.** The section shows the sparse lattice merging
into the inner face of the outer shell rather than standing a bead clear of it. Measured on
the printed file, and controlled against the detector's own error rate by re-slicing the same
geometry with the weave off:

| same binary, same geometry | median shell1 -> next gap | below 0.30 mm |
|---|---|---|
| weave 0% (plain concentric rings) | 0.400 mm | 13 / 360 azimuths |
| weave 100% | 0.538 mm | **101 / 360 azimuths** |

The 13/360 is the ray-crossing detector's false-positive floor on a wavy lattice; 101/360 is
four times that, so the population is real. At roughly a quarter of azimuths the lattice sits
within 0.30 mm of the inner shell ring, where the clamp is supposed to hold it at one spacing
(0.396 mm). W2 shipped with the caveat that #5 could not be verified from G-code at all and #7
only indirectly - the indirect evidence (ladder bead counts rising) was real but did not mean
the clamp was sufficient.

Open as **W2b**: find why the bound is not binding. Note the clamp caps the OUTWARD swing at
`outer_prof - spacing`, so a lobe that would exceed it should flatten against the bound at
exactly 0.396 - a gap of 0.12 is not explainable by the clamp working as written, so this is a
mechanism question, not a tuning one.

#### Correction (same day): the W2 clamp violation is ~1% of rays, not 28%

The 101/360 figure above is wrong and the conclusion drawn from it - "W2's clamp is NOT
holding" - overstates the defect by about 25x. Two measurement faults stacked in it, both
mine, and the corrected detector is validated by a control that now reads exactly zero:

1. **Junction hops counted as ring crossings.** The continuous bead connects rings with
   EXTRUDED radial moves, so a ray crossing a hop registers material at an intermediate
   radius no ring occupies. Excluding segments whose endpoints differ in radius by more than
   0.15 mm removes them: at weave 0 the sub-0.30 population drops 13/360 -> 0/360.
2. **Pooling layers.** The 8-12 mm band spans ~16 layers and the weave rotates 9 degrees per
   layer, so rings at different heights sit at different radii; pooling them into one ray
   fabricates gaps that exist in no single layer. Measuring per layer is the only valid form.

Per layer, hops excluded:

| | min shell1 -> next gap | rays below 0.30 mm |
|---|---|---|
| weave 0% | 0.399 mm | **0 / 5758** |
| weave 100% | 0.121 mm | **63 / 5814 (1.1%)** |

So there IS a real residue - roughly one ray in ninety puts the lattice ~0.12 mm from the
inner shell ring where the clamp should hold 0.396 - but it is rare and local, not the
systematic failure claimed above. The clamp as written cannot produce it (it caps the outward
swing at `outer_prof - spacing`), so the mechanism is still unidentified.

The printed observation is also better explained now. The lattice IS physically joined to the
shell, by the junction hops, and that is the mode working as designed - one continuous bead
has to connect the lattice to the shells somehow. A cut section cannot distinguish "bonded by
a hop" from "over-extruded against the shell" by eye, which is why the G-code measurement had
to be the arbiter, and why it had to be a per-layer one.

### PlateD reprint printed (2026-08-25): W3 confirmed, and the bore-side gap found

Three objects: D1_ladder (post-W3), D2_weave at density 25, D3_dense at density 50 via a
per-object override.

**Denser is better, by eye.** David's call on D2 vs D3: the denser small model is the better of
the two. That is the dynamic-density question answered in the direction the measurement
predicted (16.0% of the circumference tied at stride 4 against 22.5% at stride 2), so the
ring-count floor is worth building rather than just recording.

**The lattice joined the outer wall and floated free of the bore.** David: "the thickest two
layers have the infill join into the outer wall correctly, it leaves a gap between the 'foam'
and the inner wall". Measured on the printed file, and the cause is exact: kept rings were
chosen by `(i - sl) % K == 0`, which ANCHORS the pattern at the outer shell and dumps the whole
remainder at the bore.

| wall | outer shell -> infill | infill -> bore shell |
|---|---|---|
| 8 beads | 0.508 mm | **1.265 mm** |
| 7 beads | 0.504 mm | **0.702 mm** |
| 6 beads | 0.516 mm | 0.495 mm |
| 5 beads | 0.378 mm | 0.378 mm (below the gate, solid) |

The asymmetry tracks `n_buried mod K` exactly, which is why it shrinks as the wall thins.

Fixed by distributing the kept rings across the buried span inclusive of BOTH ends, with a
floor of two so each shell always gets an adjacent ring:

| wall | outer -> infill | infill -> bore |
|---|---|---|
| 8 beads | 0.523 mm | **0.397 mm** |
| 7 beads | 0.529 mm | **0.400 mm** |
| 6 beads | 0.397 mm | 0.396 mm |

Both faces now tie at one bead spacing. Cost on this plate: 17.09 g -> 17.75 g (+3.9%), all of
it the extra ring the floor adds on the thick steps. Density-0 fence byte-identical on all 20
rows; full suite green.

**Still unanswered, and it needs a different print:** how well the rings lock to each other
BETWEEN layers. A section shows in-plane structure and the top shows one layer; David's point
stands that a de-lidded print - top shells suppressed so the lattice is open to view from above
through several layers - is what would actually show the inter-layer lock.

### PlateL de-lid printed (2026-08-25): the lattice is a bundle of columns, and why

David: "the infill foam has become a set of 'columns' where multiple adjacent tubes connect
properly to form the foam. However, the foam does not connect these different columns together
- either in the less or more dense print."

The de-lid was the right instrument: a section shows in-plane structure and a top view shows
one layer, and neither could have shown this.

**Measured cause.** Adjacent rings alternate phase by exactly pi, so EVERY pair reaches its
closest approach at the same set of azimuths - pair (0,1) where the sine peaks, pair (1,2)
where it troughs, repeating each lobe. The untied gaps therefore LINE UP across every radius,
and those aligned gaps are the column walls. One layer, 360 azimuth samples:

| | azimuths with a tie | tie bands | band width | untied gaps |
|---|---|---|---|---|
| L1 density 25 | 39% | 9 | six at ~22 deg | mean 24 deg |
| L2 density 50 | 66% | 20 | six at ~27-30 deg | mean 6 deg |

Six wide bands, one per lobe, exactly as the phase rule predicts.

**Fix: a per-ring phase DRIFT of 2*pi/m on top of the anti-phase alternation**, so each
successive ring's tie band walks around the circumference and a gap at one radius is covered by
a tie at the next. Measured over the same layers:

| | azimuths with a tie | tie bands/layer |
|---|---|---|
| L1 density 25 | 36.7% -> **66.7%** | 9.6 -> 14.1 |
| L2 density 50 | 61.3% -> **70.2%** | 21.1 -> 18.5 |

The sparse case nearly doubles and its ties break into more, smaller, radially staggered bands.
It costs nothing: 5.18 g -> 5.19 g, same 28 min, because it only moves phase. Density-0 fence
byte-identical on all 20 rows; full suite green.

Note what this does NOT prove: tie fraction is a G-code measure of where beads come within a
bead width of each other. Whether the staggered ties actually bond in plastic, and whether the
columns disappear, needs the reprinted de-lid part.

## PlateL reprint (2026-08-25): the amplitude law was wrong, not the phase

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

## The peel test (2026-08-25): the anchor ring is the peel plane

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

## Pinning the anchor ring, and honouring the density (2026-08-25)

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

## The peel again (2026-08-25): the pin worked and the peel did not move

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

## Radial-traverse interior (2026-08-25): the generator that has radial members

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

## PlateN (2026-08-25): the serpentine's first physical test, with a control

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

## Correction (2026-08-25): the tie profile was a measurement artifact

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

## The 43-vs-16 hunt (2026-08-26): the geometry was right, the target was wrong

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

## Gap width limit (2026-08-26): chosen against the plastic, not the metric

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
