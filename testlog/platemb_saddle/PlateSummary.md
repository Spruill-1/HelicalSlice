<!-- COMPLETE - printed, peeled, weighed, photographed (David) + interpreted. Follows the report
     template at the bottom (carried over from platemb_arcpeel). -->

# PlateMB_saddle - saddle-fill print gate

## Plate layout (which coupon is which)

Positions as the coupons sat on the bed, viewed from above with the printer's FRONT (the door,
low Y) at the BOTTOM and X increasing to the RIGHT. Auto-arrange scattered the pairs, so each
piece is located by position; identify the SHAPE by its bores, then soff vs son by cell (and the
_son of each pair is visibly DENSER between the bores - the whole point of the test).

```
              LEFT (X~50)          CENTER (X~180)        RIGHT (X~305)
  BACK        MB2e_son             --                    --
              (2 equal bores)
  MIDDLE      F7_son               MB2u_soff             --
              (3 unequal bores)    (2 unequal bores)
  FRONT       F7_soff              MB2e_soff             MB2u_son
              (3 unequal bores)    (2 equal bores)       (2 unequal bores)
```

Anchors to orient your desk arrangement:
- The LEFT edge had the **column of three**; center had two; right had a single piece.
- Both 3-bore pieces (F7) were on that left edge - front (soff) and middle (son); the piece BEHIND
  them, back-left, is a 2-EQUAL-bore (MB2e_son).
- The lone piece off on the RIGHT-FRONT is MB2u (2 UNEQUAL bores), _son.
- Within each shape pair, the coupon with the denser inter-bore lattice is the _son (saddle_fill on).

| coupon | column | row | bores | saddle_fill | bed center (mm) |
|---|---|---|---|---|---|
| F7_soff   | left   | front  | 3 unequal | off | (50, 43) |
| F7_son    | left   | middle | 3 unequal | **ON**  | (50, 160) |
| MB2e_son  | left   | back   | 2 equal   | **ON**  | (53, 277) |
| MB2e_soff | center | front  | 2 equal   | off | (180, 43) |
| MB2u_soff | center | middle | 2 unequal | off | (180, 160) |
| MB2u_son  | right  | front  | 2 unequal | **ON**  | (305, 51) |

---

| Name | Measured Weight | Pictures | Qualitative |
|------|-----------------|----------|-------------|
| MB2e_soff | 7.5 g | [cut](MB2e_soff_cut.jpg) [top](MB2e_soff_top.jpg) [peel](MB2e_soff_peel.jpg) | Yup - matches expectations, the center is under dense. |
| MB2e_son  | 7.9 g | [cut](MB2e_son_cut.jpg) [top](MB2e_son_top.jpg) [peel](MB2e_son_peel.jpg) | Yes this does resist the peel, but at a cost of the overal print's 'squishiness' |
| MB2u_soff | 7.4 g | [cut](MB2u_soff_cut.jpg) [top](MB2u_soff_top.jpg) [peel](MB2u_soff_peel.jpg) | very solid outer rings, appropriately dense center between the bores |
| MB2u_son  | 8.3 g | [cut](MB2u_son_cut.jpg) [top](MB2u_son_top.jpg) [peel](MB2u_son_peel.jpg) | Yeah this one is pretty good, fairly resilient but retains some squishiness. |
| F7_soff   | 4.9 g | [cut](F7_soff_cut.jpg) [top](F7_soff_top.jpg) [peel](F7_soff_peel.jpg) | peels like crazy around the bores |
| F7_son    | 7.3 g | [cut](F7_son_cut.jpg) [top](F7_son_top.jpg) [peel](F7_son_peel.jpg) | doesn't peel - but appears basically fully solid and has virtually no squish. |

---

## Sim pre-screen (platecheck, BEFORE printing)

Screened from the sliced gcode at layer 24. Weld-gap = largest bead-free pocket's inscribed
diameter (mm); the weld distance is ~1.19 mm, so > 1.2 mm is a WELD-GAP that peels. Weight is the
per-object sliced total (the densification pre-screen matched plastic within noise on the last plate).

| coupon | sliced wt (g) | weld-gap (mm) | bridges | verdict |
|---|---|---|---|---|
| MB2e_soff | 8.82 | **1.20** | 6 | WELD-GAP (the peel line) |
| MB2e_son  | 9.32 | **0.80** | 6 | ok - closes, **+6 % weight** |
| MB2u_soff | 8.56 | 0.40 | 5 | ok |
| MB2u_son  | 9.68 | 0.80 | 3 | ok - holds, fewer bridges, +13 % |
| F7_soff   | 5.71 | 0.80 | 5 | ok |
| F7_son    | 8.78 | 0.80 | 5 | ok - void unchanged, **+54 % weight** |

Gate result: PASS - every _son coupon welds (< 1.2 mm), no coupon gains bridges.

---

## Verdict

**Bottom line:** saddle_fill WORKS - it closes the inter-bore peel on all three shapes (MB2e_son
resists where MB2e_soff separated, MB2u_son is resilient, F7_son does not peel at all). This is the
first multibore fix this campaign has proposed in sim that SURVIVED plastic - arc inverted, this held.
But the plate also surfaced the axis the sim never measured: COMPLIANCE. TPU's whole point is squish,
and saddle_fill trades squish for peel resistance in proportion to how much it densifies. MB2u_son
(+12 %) is the sweet spot - resilient AND still squishy. F7_son (+49 %) went too far - it does not
peel because it printed basically SOLID, "virtually no squish," ruined as a flexible part. So the
mechanism is validated and the target is now clear: close the peel with the MINIMUM density, not the
weld-distance blanket. The F7 over-fire I hoped was "harmless denser" is NOT harmless on a flexible
material - which decides the deferred question: build the post-weave (minimum-density) targeting.

**Setup:** PlateMB_saddle, TPU 0.4 mm, cs_all25 (weld 25, interpolated rings, de-lidded top+bottom).
Six coupons = {MB2e 2-equal-bore, MB2u 2-unequal-bore, F7 3-unequal-bore} x {saddle-off control,
saddle-on}, one variable changed per pair (continuous_shell_saddle_fill). 50.93 g, ~3 h 36 m.

### Prediction vs. reality

| Axis | Sim said | Plastic showed |
|---|---|---|
| MB2e weld-gap (peel line) | closes 1.20 -> 0.80 mm at +6 % weight | CONFIRMED - MB2e_son resists the peel where MB2e_soff separates; cut face shows the strip filled |
| MB2u (already good) | stays welded, no peel | CONFIRMED - both hold; son "fairly resilient but retains some squishiness" (the sweet spot) |
| F7 (+54 % weight) | void unchanged; extra density in the strips | INVERTED - soff "peels like crazy around the bores" (sim called it OK), son "doesn't peel but basically fully solid, no squish". The fill was NEEDED, but the blanket target over-shot to solid |
| Compliance / squish | not measured | NEW AXIS - saddle_fill costs squish in proportion to density: MB2u_son keeps it, MB2e_son less, F7_son loses it entirely |
| Bridges | none rise on any coupon | CONFIRMED - no coupon peeled from bridging; the fix added no travel |
| Densification (son/soff weight) | MB2e 1.06, MB2u 1.13, F7 1.54 | MATCHED - measured 7.5->7.9 (1.05), 7.4->8.3 (1.12), 4.9->7.3 (1.49); pre-screen accurate again |

### Mechanisms (what saddle_fill does)

- **saddle_fill (son)** - the interpolated rings are level sets whose count is set by the widest
  place, which holds the weld distance everywhere EXCEPT the inter-bore saddle: there the field has a
  critical point, the level sets go singular and spread apart however high the count goes (why the
  void was density-invariant). `family_saddle` threads extra full level sets through any corridor
  still wider than the weld distance, so the strip densifies without cutting the section or adding a
  bridge - the weave stays fully connected. CONFIRMED in the cut faces: MB2e_son_cut shows the
  centre strip filled with lattice where MB2e_soff_cut shows an open, bubbly foam void; MB2u_son_cut
  shows the same, still a lattice (retains squish). At the extreme, F7_son_cut is nearly glassy-SOLID
  through the section - the level sets were threaded so densely they merged into solid plastic.
- **soff (control)** - the shipped fraction path: one ring count for the whole section, welds the
  bores but leaves the saddle strip as evenly-spread foam void. CONFIRMED: MB2e_soff_cut and
  F7_soff_cut are loose open foam in the inter-bore region - exactly the PlateMB_arcpeel behaviour,
  and where both peeled. MB2u_soff was the exception again ("appropriately dense centre") - its
  UNEQUAL bores sit farther apart, so its saddle corridor is narrower and the base count already
  reaches it, the same reason MB2u_frac was the one acceptable coupon last plate.

### Quantitative

- Weights, measured (g): MB2e 7.5 -> 7.9 (son/soff 1.05) / MB2u 7.4 -> 8.3 (1.12) /
  F7 4.9 -> 7.3 (1.49). Sliced pre-screen said 1.06 / 1.13 / 1.54 - matched within noise AGAIN, so the
  densification pre-screen is trustworthy; the physical parts run ~15 % lighter than the sliced number
  across the board (a fixed offset, not a per-coupon error).
- The compliance/peel trade tracks weight monotonically: F7 (heaviest son, 1.49x) = solid/no-squish;
  MB2u (1.12x) = resilient + squishy; MB2e (1.05x) = resists but a little stiff. The lighter the
  added fill, the better the balance - which is the argument for minimum-density targeting.
- Weld-gap (mm, layer 24): MB2e 1.20 -> 0.80, MB2u 0.40 -> 0.80, F7 0.80 -> 0.80. Note F7's 0.80
  "OK" number did NOT predict its "peels like crazy" - the metric missed the defect (below).
- Bridges: none rise son vs soff; no coupon peeled from bridging.

### Why prediction diverged (sim blind spots to watch)

1. **void_w missed F7_soff's peel entirely.** platecheck read F7_soff "ok, 0.80 mm" yet it "peels
   like crazy around the bores." The single-largest-inscribed-disc metric cannot see a distributed
   under-dense strip wrapping three bores - the void is spread thin, not pooled in one pocket. So the
   gate PASSED a coupon that fails by hand. -> HARNESS FIX: a coverage/solid-FRACTION metric over the
   inter-bore band (not a single disc), and always peel by hand. The good news: the fix that void_w
   under-motivated (densify F7) was still the right call - F7_son does not peel.
2. **The sim never measured COMPLIANCE, and that turned out to be half the answer.** Every metric this
   campaign built - void, bridges, weld-gap, weight - scores a part on strength/connectedness. None
   scores whether it stays SQUISHY, which for TPU is the point. F7_son passes every sim gate and is a
   useless brick. -> HARNESS FIX: a compliance proxy from the slice - solid-fraction of the cross
   section, or lattice-vs-solid area - so "closed the peel but went solid" is visible BEFORE printing.
3. **Pre-weave over-fire IS harmful (my earlier call was wrong).** I predicted F7's +54 % was "not
   harmful, denser only." On a flexible material, denser-to-solid IS the harm. `family_saddle`
   measures ring centrelines before the weave displaces them, so it fires on F7's wide-but-weave-
   weldable corridors as hard as on MB2e's genuine hole. The plate decided the deferred question:
   build POST-weave detection (a 2-pass slice) so only the corridors the weave truly leaves open get
   filled - F7 would get little, MB2e the full fix.

### Decision & next actions

- **KEEP `continuous_shell_saddle_fill`** - the mechanism is validated (closes the peel on all three
  shapes, no bridges, weight as predicted). Commit the config code. But keep it EXPERIMENTAL / off,
  NOT toward default yet: at the weld-distance target it over-densifies weave-weldable sections (F7)
  into solid, which ruins a flexible part.
- **BUILD post-weave (minimum-density) targeting next** - the plate turned "maybe wasteful" into
  "definitely over-shoots." Detect on the assembled/woven spiral, not the pre-weave ring centrelines,
  and only thread extra levels where the weave actually leaves the strip open. MB2u_son proves the
  balance exists (resilient + squishy); the job is to hit it on every shape, not just the one whose
  bores happen to sit far apart.
- **ADD a compliance metric to the harness** - solid-fraction / lattice-vs-solid over the inter-bore
  band - so "closed the peel but went solid" (F7_son) is caught in sim, and re-screen with it before
  the next plate. Pair it with a strip-coverage metric so void_w stops passing coupons that peel
  (F7_soff). Both are the PlateMB_arcpeel standing lesson: the sim needs measures the plate taught it.
- **Re-plate** once minimum-density targeting lands: the same six coupons should then close every
  peel AND keep every coupon squishy - F7_son is the one to watch.

<!-- REPORT TEMPLATE (keep this shape for every result plate):
  1. Observations table: Name | Measured Weight | Pictures | Qualitative
  2. Verdict: bottom line + setup
  3. Prediction vs. reality table (per axis: what was expected vs what the print showed)
  4. Mechanisms per mode, each tied to a photo or a metric
  5. Quantitative (weights vs prediction, plus any defect metric)
  6. Why prediction diverged (name the blind spots -> concrete harness fixes)
  7. Decision & next actions (what to retire/keep/build next)
-->
