"""Radial tie CHAINS in the woven lattice, and how they stack in Z.

The earlier metric counted, per adjacent ring pair, the azimuths where the pair comes within a
bead width. That is a per-PAIR measure and it improved. But the printed complaint is about the
wall acting as a bundle of columns, and a column is defined by what happens ACROSS the radius:
whether, at a given azimuth, the ties chain all the way from the outer shell to the bore shell.
A wall where every pair is 67% tied can still have zero full bridges if the tied azimuths never
line up radially. So measure the chain, not the links.
"""
import math, re, sys, zipfile
import numpy as np

NRAY = 360

def segs(path):
    with zipfile.ZipFile(path) as z:
        g = z.read("Metadata/plate_1.gcode").decode("utf-8", "replace")
    lh = float(re.search(r"; layer_height = ([\d.]+)", g).group(1))
    fd = float(re.search(r"; filament_diameter = ([\d.]+)", g).group(1))
    area = math.pi * fd * fd / 4.0
    objs, order, cur = {}, [], None
    x = y = zz = None
    for line in g.splitlines():
        m = re.match(r"; printing object (\S+)", line)
        if m:
            cur = m.group(1)
            if cur not in objs: objs[cur] = []; order.append(cur)
            continue
        if line.startswith("; stop printing object"): cur = None; continue
        if not (line.startswith("G1 ") or line.startswith("G0 ")): continue
        nx = re.search(r" X(-?[\d.]+)", line); ny = re.search(r" Y(-?[\d.]+)", line)
        nz = re.search(r" Z(-?[\d.]+)", line); ne = re.search(r" E(-?[\d.]+)", line)
        px, py = x, y
        if nx: x = float(nx.group(1))
        if ny: y = float(ny.group(1))
        if nz: zz = float(nz.group(1))
        if ne and float(ne.group(1)) > 0 and cur is not None and px is not None and zz is not None:
            L = math.hypot(x - px, y - py)
            if L > 1e-9:
                objs[cur].append((px, py, x, y, zz, float(ne.group(1)) * area / (L * lh), L))
    return lh, [(n, objs[n]) for n in order if objs[n]]


def crossings(mid, cx, cy):
    """(radius, width) of every bead crossing each of NRAY rays, sorted outward-in per ray."""
    out = [[] for _ in range(NRAY)]
    az = -math.pi + (np.arange(NRAY) + 0.5) * (2 * math.pi / NRAY)
    for (x0, y0, x1, y1, _z, w, L) in mid:
        if L < 1e-6: continue
        a0 = math.atan2(y0 - cy, x0 - cx); a1 = math.atan2(y1 - cy, x1 - cx)
        d0 = np.arctan2(np.sin(a0 - az), np.cos(a0 - az))
        d1 = np.arctan2(np.sin(a1 - az), np.cos(a1 - az))
        hit = ((d0 < 0) != (d1 < 0)) & (np.abs(d0 - d1) < math.pi)
        for k in np.nonzero(hit)[0]:
            t = d0[k] / (d0[k] - d1[k])
            out[k].append((math.hypot(x0 + t*(x1-x0) - cx, y0 + t*(y1-y0) - cy), w))
    return [sorted(c) for c in out]


def bridges(cl, weld=1.0):
    """Per ray: (full bridge outer->inner?, longest tied run, n beads). Tied when the centreline
    gap is at most `weld` x the mean of the two bead widths - i.e. the beads actually touch."""
    res = []
    for c in cl:
        # merge duplicates from the same bead (a polyline can cross a ray twice within a bead)
        m = []
        for r, w in c:
            if m and r - m[-1][0] < 0.12: continue
            m.append((r, w))
        if len(m) < 2:
            res.append((False, len(m), len(m))); continue
        tied = [(m[i+1][0] - m[i][0]) <= weld * 0.5 * (m[i][1] + m[i+1][1]) for i in range(len(m)-1)]
        full = all(tied)
        run, best = 1, 1
        for t in tied:
            run = run + 1 if t else 1
            best = max(best, run)
        res.append((full, best, len(m)))
    return res


def report(path, label, zlo, zhi, weld=1.0):
    lh, objs = segs(path)
    print(f"\n=== {label} (weld tol {weld:.2f}) ===")
    print(f"{'object':>10s} {'layers':>6s} {'beads/ray':>9s} {'full bridge':>12s} "
          f"{'longest run':>12s} {'dead azim':>10s} {'tallest dead col':>17s}")
    for name, sg in objs:
        # Bin by NOMINAL layer index. The Z-wave gives every woven point its own Z, so
        # grouping on the raw value produced one "layer" per extrusion and every metric came
        # back empty. lh/2 tolerance recovers the whole layer including its Z excursion.
        zs = sorted({int(round(s[4] / lh)) for s in sg if zlo <= s[4] <= zhi})
        if not zs: continue
        per_layer_bridge, runs, beads = [], [], []
        grid = np.zeros((len(zs), NRAY), dtype=bool)   # True = full radial bridge
        for li, z in enumerate(zs):
            mid = [s for s in sg if int(round(s[4] / lh)) == z]
            if len(mid) < 20: continue
            cx = sum(s[0] for s in mid) / len(mid); cy = sum(s[1] for s in mid) / len(mid)
            b = bridges(crossings(mid, cx, cy), weld)
            grid[li] = [t[0] for t in b]
            per_layer_bridge.append(np.mean(grid[li]))
            runs.append(np.mean([t[1] for t in b])); beads.append(np.mean([t[2] for t in b]))
        # A DEAD azimuth is one never bridged in any layer of the band: a fracture plane that
        # runs the full height. A dead COLUMN is the longest stack of consecutive layers in
        # which one azimuth stays unbridged - the thing that reads as a column boundary.
        dead = np.mean(~grid.any(axis=0))
        tall = 0
        for k in range(NRAY):
            run = best = 0
            for li in range(grid.shape[0]):
                run = 0 if grid[li, k] else run + 1
                best = max(best, run)
            tall = max(tall, best)
        print(f"{name:>10s} {len(zs):6d} {np.mean(beads):9.1f} {100*np.mean(per_layer_bridge):11.1f}% "
              f"{np.mean(runs):12.2f} {100*dead:9.1f}% {tall:14d} ly")


if __name__ == "__main__":
    for p, lab in ((sys.argv[1], sys.argv[2]), (sys.argv[3], sys.argv[4])):
        report(p, lab, 4.0, 10.0, weld=float(sys.argv[5]) if len(sys.argv) > 5 else 1.0)
