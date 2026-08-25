"""Where does the wall separate?

A part that peels into concentric sections separates at RADII, not at azimuths, so the weld map
cannot see it: that measures whether a given azimuth is braced, and a cylinder of non-contact is
unbraced at every azimuth equally.

Coverage answers it directly and needs no ring identity. For each radius, over every azimuth and
every layer, what fraction of samples has material there? A woven ring sweeps radially, so it
paints a BAND of radii - low coverage at a radius is normal and is where the sparse voids live.
A radius covered at NO azimuth is different in kind: a continuous cylinder of air, which is
exactly the surface a peel finds.

Z is reported alongside because the woven infill carries a per-angle Z wave that the shells do
not, so infill and shell beads meant to fuse can sit up to half a layer apart vertically - a
weakness no in-plane measure in this campaign has ever looked at.
"""
import math, sys
import numpy as np
from tiechain import segs
from interface import crossings_z

R0, R1, DR = 2.4, 10.0, 0.02


def profile(path, obj_filter, zlo=4.0, zhi=10.0):
    lh, objs = segs(path)
    out = {}
    nb = int((R1 - R0) / DR)
    for name, sg in objs:
        if obj_filter not in name:
            continue
        zs = sorted({int(round(s[4] / lh)) for s in sg if zlo <= s[4] <= zhi})
        cover = np.zeros(nb)          # samples with material at this radius
        zdev  = np.zeros(nb)          # sum |z - layer nominal| of the covering bead
        zn    = np.zeros(nb)
        total = 0
        for z in zs:
            mid = [s for s in sg if int(round(s[4] / lh)) == z]
            cx = sum(s[0] for s in mid) / len(mid)
            cy = sum(s[1] for s in mid) / len(mid)
            znom = np.median([s[4] for s in mid])
            for c in crossings_z(mid, cx, cy):
                total += 1
                for r, w, zz in c:
                    lo = int((r - 0.5 * w - R0) / DR)
                    hi = int((r + 0.5 * w - R0) / DR)
                    for b in range(max(0, lo), min(nb, hi + 1)):
                        cover[b] += 1
                        zdev[b] += abs(zz - znom)
                        zn[b] += 1
        out[name] = (cover / max(total, 1), np.where(zn > 0, zdev / np.maximum(zn, 1), 0.0), lh)
    return out


if __name__ == "__main__":
    path, tag = sys.argv[1], sys.argv[2]
    for obj in ("L2_open50", "L1_open25"):
        for name, (cov, zd, lh) in profile(path, obj).items():
            print(f"\n=== {tag}  {name} ===")
            print("  radius   coverage                                   mean |dZ|")
            for b in range(len(cov)):
                r = R0 + (b + 0.5) * DR
                if b % 2:
                    continue
                bar = "#" * int(round(cov[b] * 40))
                flag = "   <-- CYLINDER OF AIR" if cov[b] < 0.02 else ""
                print(f"  {r:6.2f}  {cov[b]*100:5.1f}% {bar:<41s} {zd[b]:.3f}mm{flag}")
            dead = [R0 + (b + 0.5) * DR for b in range(len(cov)) if cov[b] < 0.02]
            # group contiguous dead radii into bands
            bands, cur = [], []
            for r in dead:
                if cur and r - cur[-1] > 1.5 * DR:
                    bands.append((cur[0], cur[-1])); cur = []
                cur.append(r)
            if cur:
                bands.append((cur[0], cur[-1]))
            print(f"\n  separation surfaces (<2% coverage at every azimuth and height):")
            for a, b2 in bands:
                print(f"     r = {a:.2f} .. {b2:.2f} mm   (width {b2-a+DR:.2f} mm)")
