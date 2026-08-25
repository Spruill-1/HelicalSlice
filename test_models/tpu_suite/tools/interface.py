"""Which radial interface is the weak one?

A part that peels into concentric sections has separation surfaces at particular RADII, not at
particular azimuths. The weld map answers "is this azimuth braced" and cannot see a cylinder of
non-contact; this asks, for every radially adjacent bead pair, how it does around the whole
circumference and up the whole height. Two numbers matter per interface:

  in-plane gap - centreline distance, against a bead width; and
  Z offset     - the woven infill carries a per-angle Z wave and the SHELLS DO NOT, so an
                 infill bead can sit up to half a layer above or below the shell bead it is
                 supposed to fuse with. That is invisible to every in-plane measure taken so far.
"""
import math, sys
import numpy as np
from tiechain import segs, crossings

NAZ = 360


def interfaces(path, obj_filter, zlo=4.0, zhi=10.0):
    lh, objs = segs(path)
    out = {}
    for name, sg in objs:
        if obj_filter not in name:
            continue
        zs = sorted({int(round(s[4] / lh)) for s in sg if zlo <= s[4] <= zhi})
        gap_acc, zoff_acc = None, None
        for z in zs:
            mid = [s for s in sg if int(round(s[4] / lh)) == z]
            cx = sum(s[0] for s in mid) / len(mid)
            cy = sum(s[1] for s in mid) / len(mid)
            # crossings() returns (radius, width); recompute with Z carried along
            per_ray = crossings_z(mid, cx, cy)
            rows_g, rows_z = [], []
            for c in per_ray:
                m = []
                for r, w, zz in c:
                    if m and r - m[-1][0] < 0.12:
                        continue
                    m.append((r, w, zz))
                rows_g.append([m[i + 1][0] - m[i][0] for i in range(len(m) - 1)])
                rows_z.append([abs(m[i + 1][2] - m[i][2]) for i in range(len(m) - 1)])
            n = min(len(r) for r in rows_g)
            g = np.array([r[:n] for r in rows_g])          # [ray][interface] from OUTSIDE in
            zo = np.array([r[:n] for r in rows_z])
            gap_acc = g if gap_acc is None else np.concatenate([gap_acc, g[:, :gap_acc.shape[1]]])
            zoff_acc = zo if zoff_acc is None else np.concatenate([zoff_acc, zo[:, :zoff_acc.shape[1]]])
        out[name] = (gap_acc, zoff_acc, lh)
    return out


def crossings_z(mid, cx, cy):
    out = [[] for _ in range(NAZ)]
    az = -math.pi + (np.arange(NAZ) + 0.5) * (2 * math.pi / NAZ)
    for (x0, y0, x1, y1, zz, w, L) in mid:
        if L < 1e-6:
            continue
        a0 = math.atan2(y0 - cy, x0 - cx)
        a1 = math.atan2(y1 - cy, x1 - cx)
        d0 = np.arctan2(np.sin(a0 - az), np.cos(a0 - az))
        d1 = np.arctan2(np.sin(a1 - az), np.cos(a1 - az))
        hit = ((d0 < 0) != (d1 < 0)) & (np.abs(d0 - d1) < math.pi)
        for k in np.nonzero(hit)[0]:
            t = d0[k] / (d0[k] - d1[k])
            out[k].append((math.hypot(x0 + t * (x1 - x0) - cx, y0 + t * (y1 - y0) - cy), w, zz))
    return [sorted(c) for c in out]


if __name__ == "__main__":
    for f, tag in ((sys.argv[1], sys.argv[2]),):
        for obj in ("L2_open50", "L1_open25"):
            for name, (g, zo, lh) in interfaces(f, obj).items():
                print(f"\n=== {tag} {name}  ({g.shape[0]} ray-layers, layer height {lh}) ===")
                print(f"{'iface':>6s} {'median gap':>11s} {'min gap':>9s} {'welded':>8s} "
                      f"{'median dZ':>10s} {'dZ > .3lh':>10s}")
                for i in range(g.shape[1]):
                    col, zc = g[:, i], zo[:, i]
                    print(f"{i:6d} {np.median(col):10.3f}mm {col.min():8.3f}mm "
                          f"{100*np.mean(col <= 0.44):6.1f}% {np.median(zc):9.3f}mm "
                          f"{100*np.mean(zc > 0.3*lh):9.1f}%")
                print("  interface 0 is the outermost bead pair; the last is at the bore.")
