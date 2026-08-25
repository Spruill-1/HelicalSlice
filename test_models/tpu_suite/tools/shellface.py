"""Is the peel surface the infill-to-SHELL interface?

The de-lid part separated into three concentric sections, which is [outer shell | infill | bore
shell] if the two infill-to-shell interfaces are the weak ones. Coverage already showed there is
no void there, so the suspect is Z: the woven infill carries a per-angle Z wave of up to 0.45 of
a layer and THE SHELLS DO NOT, so an infill bead and the shell bead it is meant to fuse with sit
at different heights and share only part of their walls. Nothing in this campaign has measured
that - every metric so far has been in-plane.

Shell beads are identified from a weave-off slice, where they sit at fixed radii; they do not
move when the weave is switched on, so the same radii identify them in the woven run.
"""
import sys
import numpy as np
from tiechain import segs
from interface import crossings_z


def shell_radii(nominal_path, obj_filter, lh, zlo=4.0, zhi=10.0):
    """Radii of beads that never move: the outermost/innermost runs of tightly spaced rings."""
    _, objs = segs(nominal_path)
    for name, sg in objs:
        if obj_filter not in name:
            continue
        zs = sorted({int(round(s[4] / lh)) for s in sg if zlo <= s[4] <= zhi})
        z = zs[len(zs) // 2]
        mid = [s for s in sg if int(round(s[4] / lh)) == z]
        cx = sum(s[0] for s in mid) / len(mid)
        cy = sum(s[1] for s in mid) / len(mid)
        m = []
        for r, w, _ in crossings_z(mid, cx, cy)[0]:
            if m and r - m[-1] < 0.12:
                continue
            m.append(r)
        # Shells are the contiguous runs at each end whose spacing is ~one bead.
        outer, inner = [m[-1]], [m[0]]
        for i in range(len(m) - 1, 0, -1):
            if m[i] - m[i - 1] < 0.55:
                outer.append(m[i - 1])
            else:
                break
        for i in range(len(m) - 1):
            if m[i + 1] - m[i] < 0.55:
                inner.append(m[i + 1])
            else:
                break
        return sorted(set(inner + outer)), m
    return [], []


if __name__ == "__main__":
    woven, nominal = sys.argv[1], sys.argv[2]
    lh, objs = segs(woven)
    for obj in ("L2_open50", "L1_open25"):
        sh, allr = shell_radii(nominal, obj, lh)
        print(f"\n=== {obj} ===")
        print("  rings (weave off):", " ".join(f"{r:.2f}" for r in allr))
        print("  shell rings:      ", " ".join(f"{r:.2f}" for r in sh))
        for name, sg in objs:
            if obj not in name:
                continue
            zs = sorted({int(round(s[4] / lh)) for s in sg if 4.0 <= s[4] <= 10.0})
            cats = {"shell-shell": [], "infill-SHELL": [], "infill-infill": []}
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
                    for i in range(len(m) - 1):
                        a_sh = any(abs(m[i][0] - s) < 0.12 for s in sh)
                        b_sh = any(abs(m[i + 1][0] - s) < 0.12 for s in sh)
                        key = ("shell-shell" if a_sh and b_sh
                               else "infill-infill" if not a_sh and not b_sh else "infill-SHELL")
                        cats[key].append((m[i + 1][0] - m[i][0], abs(m[i + 1][2] - m[i][2])))
            print(f"\n  {'pair type':>14s} {'n':>7s} {'med gap':>9s} {'med dZ':>9s} "
                  f"{'dZ/layer':>9s} {'dZ > .3 layer':>14s}")
            for k, v in cats.items():
                if not v:
                    continue
                a = np.array(v)
                print(f"  {k:>14s} {len(a):7d} {np.median(a[:,0]):8.3f}mm "
                      f"{np.median(a[:,1]):8.3f}mm {np.median(a[:,1])/lh:9.2f} "
                      f"{100*np.mean(a[:,1] > 0.3*lh):13.1f}%")
