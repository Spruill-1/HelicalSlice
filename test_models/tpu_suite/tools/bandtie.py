"""Tie fraction per interface, with rings assigned by DISJOINT INTERVAL rather than by rank.

rank_rows matched rings by their order along a ray and kept only rays showing the modal number of
crossings. On a saturated weave that filter discards more than half the circumference - a ray
that grazes a steep transition can cross one ring twice or miss it - and the rays it keeps are
not a fair sample of the ones that matter. It reported an alternating tie profile on geometry
whose adjacent rings correlate at -0.87 to -0.98 and whose every interface reaches a minimum gap
of exactly one spacing, which those two facts together say should be near-uniform.

Each ring's excursion is bounded by its own caps, so its radius lives in a known interval, and on
this geometry the intervals do not overlap. Assigning by interval needs no filter at all: a ring
is whichever crossing falls in its band, and a ray contributes to an interface whenever both of
its rings are found. Nothing is discarded for being inconvenient.
"""
import math, sys
import numpy as np
from tiechain import segs
from interface import crossings_z

SP = 0.3963


def nominal_rings(nw_path, obj, lh):
    _, objs = segs(nw_path)
    for name, sg in objs:
        if obj not in name:
            continue
        zs = sorted({int(round(s[4] / lh)) for s in sg if 4.0 <= s[4] <= 10.0})
        z = zs[len(zs) // 2]
        mid = [s for s in sg if int(round(s[4] / lh)) == z]
        cx = sum(s[0] for s in mid) / len(mid)
        cy = sum(s[1] for s in mid) / len(mid)
        m = []
        for r, w, _ in crossings_z(mid, cx, cy)[0]:
            if m and r - m[-1] < 0.12:
                continue
            m.append(r)
        return sorted(m)
    return []


def bands(nom, shell_n=2, pinned=True):
    """Interval each ring's radius can occupy, from the same cap law the slicer uses."""
    n = len(nom)
    static = set(range(shell_n)) | set(range(n - shell_n, n))
    if pinned:
        static |= {shell_n, n - shell_n - 1}
    out = []
    for j, r in enumerate(nom):
        if j in static:
            out.append((r - 0.12, r + 0.12)); continue
        ci = co = 0.0
        if j > 0:
            g = r - nom[j - 1]
            ci = max(0.0, (1.0 if (j - 1) in static else 0.5) * (g - SP))
        if j + 1 < n:
            g = nom[j + 1] - r
            co = max(0.0, (1.0 if (j + 1) in static else 0.5) * (g - SP))
        out.append((r - ci - 0.12, r + co + 0.12))
    return out


def tie_profile(path, nw_path, obj, zlo=4.0, zhi=10.0, naz=720):
    lh, objs = segs(path)
    nom = nominal_rings(nw_path, obj, lh)
    bd = bands(nom)
    overlap = sum(1 for i in range(len(bd) - 1) if bd[i][1] > bd[i + 1][0])
    for name, sg in objs:
        if obj not in name:
            continue
        zs = sorted({int(round(s[4] / lh)) for s in sg if zlo <= s[4] <= zhi})
        hit = np.zeros(len(nom) - 1)
        seen = np.zeros(len(nom) - 1)
        for z in zs:
            mid = [s for s in sg if int(round(s[4] / lh)) == z]
            cx = sum(s[0] for s in mid) / len(mid)
            cy = sum(s[1] for s in mid) / len(mid)
            for c in crossings_z(mid, cx, cy):
                m = []
                for r, w, zz in c:
                    if m and r - m[-1][0] < 0.12:
                        continue
                    m.append((r, w, zz))
                pick = [None] * len(nom)
                for r, w, zz in m:
                    for j, (lo, hi) in enumerate(bd):
                        if lo <= r <= hi:
                            if pick[j] is None or abs(r - nom[j]) < abs(pick[j] - nom[j]):
                                pick[j] = r
                            break
                for i in range(len(nom) - 1):
                    if pick[i] is not None and pick[i + 1] is not None:
                        seen[i] += 1
                        if pick[i + 1] - pick[i] <= 1.15 * SP:
                            hit[i] += 1
        return nom, overlap, 100.0 * hit / np.maximum(seen, 1), seen
    return [], 0, [], []


if __name__ == "__main__":
    nw = "delid_noweave.3mf"
    print(f"{'config':>24s}  tie% per interface                     worst   mean   coverage")
    for f, tag in ((sys.argv[1], sys.argv[2]),) if len(sys.argv) > 2 else (
            ("sy_base", "pin + drift 2 (printed)"), ("pin_nodrift", "pin + drift 0"),
            ("sy_nopin_nodrift", "no pin + drift 0")):
        nom, ov, tie, seen = tie_profile(f + ".3mf", nw, "L2_open50")
        if ov:
            print(f"{tag:>24s}  {ov} band overlaps - assignment unsafe, skipping")
            continue
        inner = tie[2:-2]
        print(f"{tag:>24s}  " + " ".join(f"{t:4.0f}" for t in tie) +
              f"  {min(inner):5.1f}% {np.mean(inner):5.1f}%  {seen.min():.0f} rays min")
