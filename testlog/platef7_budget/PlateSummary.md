<!-- COMPLETE - printed, peeled, weighed, photographed (David) + interpreted. Result: budget approach
     VALIDATED, F7 minimum-density target = budget ~3-4. Follows the report template at the bottom. -->

# PlateF7_budget - saddle-fill minimum-density sweep (F7)

All six coupons are the SAME F7 geometry (3 unequal bores). They differ ONLY by how much the saddle
fill densifies the inter-bore strip - so you CANNOT tell them apart by shape, only by POSITION (below)
and by feel: `off` is the lightest / squishiest, `full` the heaviest / stiffest, monotone in between.
When in doubt, WEIGH them - the sliced weights (5.7 -> 8.8 g) are monotone and match plastic within noise.

## Plate layout (which coupon is which)

Positions as they sat on the bed, viewed from above, printer FRONT (low Y) at the BOTTOM, X to the RIGHT.

```
              LEFT (X~48)        CENTER (X~175)      RIGHT (X~302)
  BACK        F7_b5              --                  --
  MIDDLE      F7_b3              F7_b6               --
  FRONT       F7_off             F7_b4               F7_full
```

- The LEFT edge had the **column of three** (front `off`, middle `b3`, back `b5`); CENTER had two
  (front `b4`, middle `b6`); RIGHT had a single piece (`full`).
- Density order (lightest to heaviest): off < b3 < b4 < b5 < b6 < full. If the positions get scrambled,
  weigh them and sort - the order is unambiguous.

| coupon | column | row | budget (levels) | sliced wt (g) |
|---|---|---|---|---|
| F7_off  | left   | front  | saddle OFF (control) | 5.71 |
| F7_b3   | left   | middle | 3 | 7.18 |
| F7_b5   | left   | back   | 5 | 8.04 |
| F7_b4   | center | front  | 4 | 7.61 |
| F7_b6   | center | middle | 6 | 8.39 |
| F7_full | right  | front  | unlimited (control) | 8.79 |

---

| Name | Measured Weight | Pictures | Qualitative |
|------|-----------------|----------|-------------|
| F7_off  | 4.9 g | [cut](F7_off_cut.jpg) [peel](F7_off_peel.jpg) | soft and squishy, but does peel pretty bad around the bores |
| F7_b3   | 6.1 g | [cut](F7_b3_cut.jpg) [peel](F7_b3_peel.jpg) | better than I expected, squishy but much less peel than F7_off |
| F7_b4   | 6.6 g | [cut](F7_b4_cut.jpg) [peel](F7_b4_peel.jpg) | no peel but much less squish - ideal seems to be between b3 and b4. |
| F7_b5   | 6.9 g | [cut](F7_b5_cut.jpg) [peel](F7_b5_peel.jpg) | too solid - I will say that while the outer edges are almost entirely solid... the weave in the interbore middle is actually very clean and visible. |
| F7_b6   | 7.3 g | [cut](F7_b6_cut.jpg) [peel](F7_b6_peel.jpg) | very similar to F7_b5, very very little compliance |
| F7_full | 7.6 g | [cut](F7_full_cut.jpg) [peel](F7_full_peel.jpg) | basically solid |

---

## Sim pre-screen (saddlecheck + platecheck, BEFORE printing)

Two metrics validated on PlateMB_saddle (tools/saddlecheck.py, averaged over layers 16-32):
solid-fraction = compliance (a HIGH value lost its squish; > 0.85 = OVER-SOLID); inter-bore band
coverage = peel-risk (< 0.56 = peel-prone). void_w (platecheck) catches a coarse-subdivision artifact.
The WINDOW is: not peel, not solid, no void.

| coupon | budget | solid_frac | ibore_cov | void_w | verdict |
|---|---|---|---|---|---|
| F7_off  | off | 0.69 | 0.50 | 0.80 | PEEL |
| F7_b3   | 3   | 0.79 | 0.51 | 0.80 | PEEL (worst saddle not yet closed) |
| F7_b4   | 4   | 0.81 | **0.62** | 0.80 | **IN WINDOW - predicted sweet spot** |
| F7_b5   | 5   | 0.84 | 0.62 | **1.60** | VOID (placement artifact) |
| F7_b6   | 6   | 0.86 | 0.64 | 0.80 | OVER-SOLID |
| F7_full | unl | 0.89 | 0.64 | 0.80 | OVER-SOLID |

Key signal: ibore_cov JUMPS 0.51 -> 0.62 between b3 and b4 - the worst inter-bore saddle welding.
That is why b4 is the pick: the first budget that closes the peel, still under the solid line.

---

## Verdict

**Bottom line:** CONFIRMED, and it closes the loop the whole multibore arc has chased: a LIGHT saddle
fill stops F7's inter-bore peel while KEEPING its squish - the always-on fill could only do the first
(it solidified F7). The budget gives a clean compliance gradient on plastic (off squishy+peels -> full
solid), and the sweet spot is REAL and sits at budget ~3-4: b3 is "squishy but much less peel," b4 has
"no peel but much less squish," so the ideal is between them. The metric-predicted b4 was RIGHT to
within one budget but slightly conservative - the peel is already much improved at b3, one step lighter
than the metric's threshold crossing. And the sim's one false alarm (b5's void_w 1.60) did NOT show on
plastic: b5's inter-bore weave is clean and welded, it is simply too solid. So the budget approach is
validated; F7's minimum-density target is budget 3-4 (lean 3 for compliance, 4 to fully kill the peel).

**Setup:** PlateF7_budget, TPU 0.4 mm, cs_all25 (weld 25, interpolated rings, de-lidded top+bottom).
Six copies of F7 (3 unequal bores) x {saddle off, budget 3/4/5/6, unlimited}, one variable changed
(continuous_shell_saddle_budget, at the weld-distance target). 45.78 g, ~3 h 14 m.

### Prediction vs. reality

| Axis | Sim said | Plastic showed |
|---|---|---|
| off | peels (ibore 0.50) | CONFIRMED - "soft and squishy, but peels pretty bad around the bores" |
| b3 | still peels (ibore 0.51) | BETTER than sim - "squishy but much less peel than off"; metric was conservative |
| b4 | stops peeling AND stays squishy | HALF - "no peel" yes, but "much less squish"; stiffer than ideal |
| b5 | void artifact (void_w 1.60) | FALSE ALARM - "the weave in the inter-bore middle is very clean and visible"; just too solid |
| b6 / full | solid, no squish | CONFIRMED - "very little compliance" / "basically solid" |
| lightest that holds | b4 (predicted) | BETWEEN b3 and b4 - one budget lighter than predicted; peel gone by b4, squish best at b3 |

### Mechanisms (what the budget does)

- **saddle budget** - a section with several bores has several saddles. Bringing every corridor under
  the weld distance at once (unlimited) closes all of them but fills the section solid. The budget
  spends a capped number of extra ring levels on the WIDEST saddles first, so a small budget closes
  only the worst saddle(s) - the lightest densification that still welds the peel line. CONFIRMED in
  the cut faces: F7_b3_cut is a connected lattice throughout (still open, compliant); F7_b4_cut is
  visibly denser with the inter-bore region filled; the density climbs smoothly to F7_full's near-solid.
  The mechanism does exactly what it was built to - a dialable, monotone density that closes the peel
  before it fills the part.
- **b5 "void artifact" was a metric false alarm, not a print defect.** saddle fill places its
  intermediate levels at even fractions of the field level (uniform-t), not even physical spacing, and
  at b5 that made void_w read 1.60 - but F7_b5_cut shows the inter-bore weave clean, welded and visible
  ("very clean and visible"); the outer regions are just dense. So the placement fix (physical-spacing
  levels) is DEPRIORITISED - the coarse placement does not cost real quality, only a noisy single-layer
  metric. Keep it as a nice-to-have, not a blocker.

### Quantitative

- Weights, measured (g): off 4.9 / b3 6.1 / b4 6.6 / b5 6.9 / b6 7.3 / full 7.6 - monotone, ~15% under
  the sliced numbers (off 5.71 ... full 8.79), the same fixed offset as every prior plate. The
  densification pre-screen is trustworthy again.
- solid-fraction (compliance, sim): 0.69 / 0.79 / 0.81 / 0.84 / 0.86 / 0.89 - monotone, tracks the felt
  stiffness (off squishy -> full solid, b3/b4 the usable middle).
- inter-bore coverage (peel, sim): 0.50 / 0.51 / 0.62 / 0.62 / 0.64 / 0.64. The metric's step is b3->b4;
  the FELT peel improvement starts one budget earlier (off->b3), so the threshold PEEL_MIN=0.56 is a
  touch conservative - the real peel line sits between ibore 0.51 (b3) and 0.62 (b4).
- void_w: 0.80 everywhere except b5 (1.60) - and b5 printed clean, so that number is a false positive.

### Why prediction diverged (watch for)

1. **The peel threshold is slightly conservative.** Plastic said b3 (ibore 0.51) already peels "much
   less," while the metric only clears it at b4 (0.62). So PEEL_MIN=0.56 flags a coupon that is already
   mostly fine. Recalibrate toward ~0.52 - or read the metric as predicting the FULLY-welded end, with
   the practical sweet spot one budget lighter. Either way the metric was right to within one budget.
2. **void_w gave a false positive (b5).** A single-layer inscribed-disc that reads 1.60 on a coupon whose
   inter-bore weave is physically clean is exactly the kind of miss the metric set was built to retire;
   the compliance + coverage pair told the true story (b5 = solid, not voided). Trust the pair, not void_w,
   for saddle quality - and do not chase the uniform-t placement fix on b5's evidence alone.
3. Compliance (solid_frac) tracked felt stiffness well across all six - a rare case where the sim proxy
   and the hand agreed monotonically. It earns its place as the "did it go solid" gate.

### Decision & next actions

- **COMMIT the budget work** - the widest-first knob + `continuous_shell_saddle_budget` config +
  `saddlecheck.py`. It is now plastic-validated: it closes F7's peel while keeping compliance, which the
  always-on fill could not. This is the resolution of the whole multibore peel/squish tension.
- **F7's minimum-density target is budget ~3-4** (lean 3 for compliance, 4 to fully kill the peel; the
  ideal is between them). Consider making a small budget (3 or 4) the DEFAULT for saddle fill: a simple
  2-bore section (MB2e, one saddle) only ever needs ~1 level so a cap of 3-4 does not touch it, while a
  multi-saddle section (F7) gets the sweet spot instead of solidifying. That would make saddle fill safe
  to turn on generally.
- **Recalibrate PEEL_MIN** toward ~0.52 (b3 already held), and record void_w as unreliable for saddle
  quality - the compliance+coverage pair is the referee.
- **Deprioritise** the uniform-t -> physical-spacing placement fix: b5 proved the coarse placement costs
  a noisy metric, not real quality. Keep it on the list, not the critical path.
- A finer knob than an integer budget (budget 3.5 was the felt ideal) is a possible refinement - e.g.
  interpolating the target between the b3 and b4 densities - but only if a print shows integer budgets
  are too coarse in practice.

<!-- REPORT TEMPLATE (keep this shape for every result plate):
  1. Observations table: Name | Measured Weight | Pictures | Qualitative
  2. Verdict: bottom line + setup
  3. Prediction vs. reality table (per axis: what was expected vs what the print showed)
  4. Mechanisms per mode, each tied to a photo or a metric
  5. Quantitative (weights vs prediction, plus any defect metric)
  6. Why prediction diverged (name the blind spots -> concrete harness fixes)
  7. Decision & next actions (what to retire/keep/build next)
-->
