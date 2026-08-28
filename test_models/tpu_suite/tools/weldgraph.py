"""Which rings are welded to which - across chains, not only within one.

A peel is a radial fault, and per-azimuth statistics are blind to it: a cylinder of poor bond is
equally unbraced at every azimuth. This asks the question the peel asks. Every ring in the section
against every other, whichever chain each belongs to, so a chain BOUNDARY is visible as an
interface like any other.
"""
import sys
import numpy as np
from scipy.spatial import cKDTree

SP = 0.3963
TIE = 1.15 * SP


def load(p):
    out = []
    for line in open(p, encoding="utf-8"):
        t = line.split()
        if t and t[0] == "CHAIN":
            out.append((int(t[1]), int(t[3]),
                        np.array([float(v) for v in t[10:]]).reshape(-1, 2)))
    return out


def report(path, tag):
    rings = load(path)
    n = len(rings)
    seglen = [np.linalg.norm(np.roll(P, -1, 0) - P, axis=1) for _, _, P in rings]
    print(f"\n=== {tag}: {n} rings ===")
    print("      ring        len   welded to (fraction of ITS length within "
          f"{TIE:.3f} mm)")
    frac = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = cKDTree(rings[j][2]).query(rings[i][2])[0]
            m = d <= TIE
            frac[i, j] = seglen[i][m].sum() / seglen[i].sum()
    for i in range(n):
        ci, ri, P = rings[i]
        parts = [f"c{rings[j][0]}r{rings[j][1]}:{100*frac[i,j]:.0f}%"
                 for j in range(n) if frac[i, j] > 0.02]
        print(f"  c{ci}r{ri:<2d} {seglen[i].sum():9.2f}   " + (" ".join(parts) or "NOTHING"))

    # connectivity of the whole section at the tie threshold
    adj = (frac > 0.02) | (frac.T > 0.02)
    seen, stack, comp = {0}, [0], []
    while stack:
        k = stack.pop()
        comp.append(k)
        for j in range(n):
            if adj[k, j] and j not in seen:
                seen.add(j); stack.append(j)
    print(f"  -> the section falls into "
          f"{'ONE connected piece' if len(comp) == n else str(n - len(comp)) + ' ring(s) DETACHED'}")

    # the chain boundary specifically
    chains = sorted({c for c, _, _ in rings})
    for a in chains:
        for b in chains:
            if a >= b:
                continue
            ia = [i for i in range(n) if rings[i][0] == a]
            ib = [i for i in range(n) if rings[i][0] == b]
            best = max(((frac[i, j], i, j) for i in ia for j in ib), default=(0, 0, 0))
            print(f"  -> chain {a} to chain {b}: best interface "
                  f"c{rings[best[1]][0]}r{rings[best[1]][1]} -> c{rings[best[2]][0]}r{rings[best[2]][1]}"
                  f" at {100*best[0]:.1f}% of its length")
    return frac


for path, tag in ((sys.argv[1], sys.argv[2]), (sys.argv[3], sys.argv[4])):
    report(path, tag)
