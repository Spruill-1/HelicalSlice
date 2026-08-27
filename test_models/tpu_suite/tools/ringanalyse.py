"""Tie quality per interface, from the slicer's own ring polygons.

Every instrument used before this needed rays, a centroid, or a ring correspondence, and none of
those survive a concave outline - four attempts to measure the crescent produced numbers that were
provably wrong. This reads the rings the slicer actually built, so ring identity is exact and the
outline's shape is irrelevant.

Gaps are nearest-point, which is what "across the wall" means locally, and each point of the outer
ring is labelled CONVEX or CONCAVE by the sign of its discrete curvature - so the question the
print raised ("the concave-side outer ring is not tied to the foam") can be answered directly.
"""
import math, sys
import numpy as np
from scipy.spatial import cKDTree

SP = 0.3963
TIE = 1.15 * SP


def load(path):
    rings = []
    for line in open(path, encoding="utf-8"):
        t = line.split()
        if not t or t[0] != "CHAIN":
            continue
        ci, ri = int(t[1]), int(t[3])
        bore, pinch, n = int(t[5]), int(t[7]), int(t[9])
        xy = np.array([float(v) for v in t[10:]]).reshape(-1, 2)
        rings.append(dict(chain=ci, idx=ri, bore=bore, pinch=pinch, pts=xy))
    return rings


def curvature_sign(P):
    """+1 where the outline turns the way the bulk of it does (convex), -1 where it turns back."""
    n = len(P)
    s = np.zeros(n)
    for i in range(n):
        a, b, c = P[i - 1], P[i], P[(i + 1) % n]
        cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
        s[i] = np.sign(cross)
    return s if np.sum(s) >= 0 else -s


if __name__ == "__main__":
    rings = load(sys.argv[1] if len(sys.argv) > 1 else "cs_rings.txt")
    by_chain = {}
    for r in rings:
        by_chain.setdefault(r["chain"], []).append(r)
    for ci, rs in sorted(by_chain.items()):
        rs.sort(key=lambda r: r["idx"])
        print(f"\n=== chain {ci}: {len(rs)} rings, bore={rs[0]['bore']}, "
              f"outer end is {'a PINCH' if rs[0]['pinch'] else 'a real face'} ===")
        outer = rs[0]["pts"]
        sgn = curvature_sign(outer)
        concave_tree = cKDTree(outer[sgn < 0]) if (sgn < 0).any() else None
        print(f"  outline: {(sgn<0).sum()} of {len(outer)} points are CONCAVE")
        print(f"\n  {'iface':>5s} {'pts':>6s} {'med gap':>8s} {'tied':>7s} "
              f"{'tied CONVEX':>12s} {'tied CONCAVE':>13s}")
        for j in range(len(rs) - 1):
            A, B = rs[j]["pts"], rs[j + 1]["pts"]
            if len(A) < 8 or len(B) < 8:
                continue
            d, _ = cKDTree(B).query(A)
            tied = d <= TIE
            # label each point of ring A by whether it lies near the outline's concave stretch
            if concave_tree is not None:
                dc, _ = concave_tree.query(A)
                near_concave = dc < 2.0
            else:
                near_concave = np.zeros(len(A), dtype=bool)
            cvx = ~near_concave
            f_cvx = 100 * tied[cvx].mean() if cvx.any() else float("nan")
            f_ccv = 100 * tied[near_concave].mean() if near_concave.any() else float("nan")
            flag = ""
            if near_concave.any() and cvx.any() and f_ccv < 0.5 * f_cvx:
                flag = "  <-- concave side much worse"
            print(f"  {j:5d} {len(A):6d} {np.median(d):7.3f}mm {100*tied.mean():6.1f}% "
                  f"{f_cvx:11.1f}% {f_ccv:12.1f}%{flag}")
