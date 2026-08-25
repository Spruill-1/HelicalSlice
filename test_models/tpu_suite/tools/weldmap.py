"""Weld count per (azimuth, layer), measured on the printed G-code.

The per-pair tie fraction used earlier is the wrong statistic: it asks whether ring j and ring
j+1 touch SOMEWHERE, which says nothing about whether any given azimuth is braced. The wall
reads as columns when an azimuth has zero welds anywhere across its depth, and reads as a
fracture plane when those zero-weld azimuths stack in Z instead of helixing away.
"""
import math, sys
import numpy as np
from tiechain import segs, crossings

W_TOL = 1.0

def weldmap(path, name_filter, zlo=4.0, zhi=10.0, naz=360):
    lh, objs = segs(path)
    out = {}
    for name, sg in objs:
        if name_filter not in name: continue
        zs = sorted({int(round(s[4] / lh)) for s in sg if zlo <= s[4] <= zhi})
        grid = np.zeros((len(zs), naz), dtype=int)
        for li, z in enumerate(zs):
            mid = [s for s in sg if int(round(s[4] / lh)) == z]
            cx = sum(s[0] for s in mid) / len(mid); cy = sum(s[1] for s in mid) / len(mid)
            for k, c in enumerate(crossings(mid, cx, cy)):
                m = []
                for r, w in c:
                    if m and r - m[-1][0] < 0.12: continue
                    m.append((r, w))
                grid[li, k] = sum(1 for i in range(len(m) - 1)
                                  if m[i+1][0] - m[i][0] <= W_TOL * 0.5 * (m[i][1] + m[i+1][1]))
        out[name] = grid
    return out

def describe(tag, grid):
    dead = grid == 0
    # vertical stacking of dead azimuths: longest run of consecutive layers dead at one azimuth
    tall = 0
    for k in range(grid.shape[1]):
        run = best = 0
        for li in range(grid.shape[0]):
            run = run + 1 if dead[li, k] else 0
            best = max(best, run)
        tall = max(tall, best)
    # azimuthal width of dead bands, in degrees, averaged over layers
    widths = []
    for li in range(grid.shape[0]):
        row = dead[li]
        run = 0
        for k in range(2 * len(row)):
            if row[k % len(row)]: run += 1
            elif run: widths.append(run); run = 0
        if run: widths.append(run)
    deg = 360.0 / grid.shape[1]
    print(f"{tag:>34s} | {grid.mean():9.2f} {100*dead.mean():9.1f}% "
          f"{(np.mean(widths)*deg if widths else 0):11.1f}d {tall:9d} ly")

if __name__ == "__main__":
    print(f"{'file / object':>34s} | {'mean weld':>9s} {'azim w/ 0':>10s} "
          f"{'dead band':>12s} {'tallest':>12s}")
    for f, tag in (("delid_predrift.3mf", "PRE "), ("delid_drift.3mf", "POST")):
        for obj in ("L1_open25", "L2_open50"):
            for n, g in weldmap(f, obj).items():
                describe(f"{tag} {n}", g)
