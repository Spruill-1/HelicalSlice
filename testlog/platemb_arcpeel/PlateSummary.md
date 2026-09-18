| Name | Measured Weight | Pictures | Qualitative |
|------|-----------------|----------|-------------|
| MB2e_frac | 7.4g | [cut](MB2e_frac_cut.jpg) [top](MB2e_frac_top.jpg) [middle peel](MB2e_frac_middlepeel.jpg) | The parts closer to the outside edges came out very well, but the area in between the two bores is much less dense and less uniform and peels too easily. |
| F7_frac | 4.8g | [cut](F7_frac_cut.jpg) [top](F7_frac_top.jpg) | very similar to MB2e_frac, the outer walls are properly latticed (if possibly too dense) but space directly between the bores is not dense enough. |
| F7_arc | 7.6g | [cut](F7_arc_cut.jpg) [top](F7_arc_top.jpg) [bore gap](F7_arc_boregap.jpg) | Overall this print feels much denser, the internals mostly appear solid, however there is a noticeable gap between the arcs around the bores and the surrounding lattice. There also appear to be bridges moving through the bores. |
| MB2u_arc | 7.9g | [cut](MB2u_arc_cut.jpg) [top](MB2u_arc_top.jpg) | This one is interesting... it's basically a bunch of solid 'rings' of various shapes separated by near voids. The bridging in this one is particularly bad. |
| MB2e_arc | 7.2g | [cut](MB2e_arc_cut.jpg) [top](MB2e_arc_top.jpg) | Same with MB2u_arc, this is solid rings and voids - with terrible bridging. |
| MB2u_frac | 7.3g | [cut](MB2u_frac_cut.jpg) [top](MB2u_frac_top.jpg) | This is the only one of the bunch whose quality I would call acceptable. The top/bottom surfaces are clean and the weave is fully connected, although it does seem under-dense in the center and over-dense close to the edges. |
---

## Verdict

**Bottom line:** the plate INVERTS the sim. Arc was predicted to beat fraction on multibore (10.2%
vs 31.7% uncovered); on plastic arc is the WORSE mode. The per-bore-cell partition prints solid
rings separated by voids, an unwelded gap at every cell cut, and heavy bridging. Fraction is the
better foundation; its one real flaw is a narrow under-dense strip between bores.

**Setup:** PlateMB_arcpeel, TPU 0.4mm, cs_all25 (weld 25, interpolated rings, de-lidded top+bottom).
Six coupons = {MB2e 2-equal-bore, MB2u 2-unequal-bore, F7 3-unequal-bore} x {arc, fraction control},
one variable changed per pair (continuous_shell_arc_rings).

### Prediction vs. reality

| Axis | Sim said | Plastic showed |
|---|---|---|
| Densification (arc/frac weight) | F7 1.57, MB2e 0.98, MB2u 1.12 | F7 1.58, MB2e 0.97, MB2u 1.08 — MATCHED |
| Coverage / fill quality | arc fills denser, "beats" fraction | arc = solid columns + void stripes + bridging — WORSE |
| Bore welding | arc welds all bores (0 "not reached") | arc leaves a GAP at every cell cut (F7_arc_boregap) |
| Continuity | 0 teleports, both modes clean | arc bridges ~4x more (extruding, so the check missed it) |

### Mechanisms (what each mode actually does)

- **arc** — distance rings pack into dense columns per Voronoi cell but don't knit to each other
  (MB2e_arc_cut: vertical solid stripes with dark voids between); the cell cut abuts the wall lattice
  rather than fusing, leaving a void seam (F7_arc_boregap); the per-cell chains and inter-cell
  connectors span gaps and bores with extrusion (bridging).
- **fraction** — smooth u-level-sets wrap all bores and join cleanly at the saddle, fully connected
  weave (MB2u_frac_top is genuinely good). The saddle strip between two bores is where the level sets
  space out, so it is under-dense and is the peel line (MB2e_frac_middlepeel).

### Quantitative

- Weights (measured, g): MB2e_frac 7.4 / F7_frac 4.8 / F7_arc 7.6 / MB2u_arc 7.9 / MB2e_arc 7.2 /
  MB2u_frac 7.3. Arc/frac ratios matched the sim within noise (densification pre-screen was accurate).
- Bridge metric (extruding segments >2mm): arc ~1100-1180 per object vs fraction ~290-330 (~4x),
  ~7-8 m of bead in long spans vs ~2-2.5 m.

### Why the prediction diverged (sim blind spots to fix)

1. Coverage % is 2D, one layer — a bridge stretched across a void reads as "covered".
2. The continuity check only counts NON-extruding travels; every bridge is an extruding move.
3. "Welded bores" was a chain-assembly property (a ring encloses the bore skin), not physical adhesion
   across the cut.
   -> Add to the harness before the next print: long-extruding-span (bridge) detection, a true 3D
   across-layers void check, and a cut/interface weld measure.

### Decision & next actions

- RETIRE the per-bore-cell (arc) approach for multibore - geometric cuts that don't weld + bridging
  is inherent, not tunable. continuous_shell_arc_rings stays experimental and OFF.
- Real target: densify the fraction path's inter-bore SADDLE within the fraction path (bias the ring
  count / the percentile-count knob toward that strip), without ever cutting the section - keep
  MB2u_frac's connected weave, close its weak line.
- Print-before-default did its job: a sim-passing / plastic-failing approach was caught before it
  could become the default.

<!-- REPORT TEMPLATE (keep this shape for future plates):
  1. Observations table: Name | Measured Weight | Pictures | Qualitative
  2. Verdict: bottom line + setup
  3. Prediction vs. reality table (per axis: what was expected vs what the print showed)
  4. Mechanisms per mode, each tied to a photo or a metric
  5. Quantitative (weights vs prediction, plus any defect metric)
  6. Why prediction diverged (name the blind spots -> concrete harness fixes)
  7. Decision & next actions (what to retire/keep/build next)
-->
