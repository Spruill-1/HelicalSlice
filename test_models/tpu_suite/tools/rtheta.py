"""Actual r(theta) per ring against the model's prediction, ring by ring.

The model of the shipped cap law says a uniform 43% tie at every interface; the slicer delivers
16-44% on the same geometry. Both cannot be right, and the difference has to show up as either a
different WAVE (wrong amplitude, phase, or shape) or a correct wave measured wrongly. So extract
each ring's radius as a function of azimuth from the printed G-code and lay the prediction on top.

Ring order matters and is easy to get backwards: cc.loops runs OUTERMOST first, so a ring's
"cap_out" is its budget toward the LARGER radius, and the per-ring phase counter skips pinned
rings without consuming an index.
"""
import math, sys
import numpy as np
from tiechain import segs
from interface import crossings_z
from bandtie import nominal_rings, bands

SP = 0.3963
OVERSHOOT = 4.0
DRIFT_SLOTS = 2.0
DPHASE = 0.30 * math.pi          # continuous_shell_infill_weave_shift = 30%
LOBES = 6.0
SHELL_N = 2


def actual(path, nw, obj, layer_z, naz=720):
    lh, objs = segs(path)
    nom = nominal_rings(nw, obj, lh)
    bd = bands(nom, SHELL_N, pinned=True)
    li = int(round(layer_z / lh))
    for name, sg in objs:
        if obj not in name:
            continue
        mid = [s for s in sg if int(round(s[4] / lh)) == li]
        cx = sum(s[0] for s in mid) / len(mid)
        cy = sum(s[1] for s in mid) / len(mid)
        out = [dict() for _ in nom]
        for k, c in enumerate(crossings_z(mid, cx, cy)):
            th = -math.pi + (k + 0.5) * (2 * math.pi / len(crossings_z(mid, cx, cy)))
            m = []
            for r, w, zz in c:
                if m and r - m[-1] < 0.12:
                    continue
                m.append(r)
            for r in m:
                for j, (lo, hi) in enumerate(bd):
                    if lo <= r <= hi:
                        out[j][k] = r
                        break
        return nom, out, lh, li
    return [], [], 0.25, 0


def model(nom, li, naz):
    """Predicted radius per ring, from the law as written in ContinuousShell.cpp."""
    n = len(nom)
    # loops run outermost first; nominal_rings is sorted ascending, so index-map it
    outer_first = list(reversed(nom))
    shells = set(range(SHELL_N)) | set(range(n - SHELL_N, n))
    kept = [j for j in range(n) if j not in shells]
    pinned = {kept[0], kept[-1]}
    n_infill = len(kept)
    drift = min(2.0 * math.pi / max(2, n_infill), DRIFT_SLOTS * 2.0 * math.asin(1.0 / OVERSHOOT))
    layer_phase = li * DPHASE
    th = -math.pi + (np.arange(naz) + 0.5) * (2 * math.pi / naz)
    pred = {}
    ring_i = 0
    for j in range(n):
        r0 = outer_first[j]
        if j in shells or j in pinned:
            pred[j] = (np.full(naz, r0), 0.0, 0.0)
            continue
        co = ci = 0.0
        if j > 0:
            g = abs(outer_first[j - 1] - r0)
            co = max(0.0, ((1.0 if (j - 1) in pinned or (j - 1) in shells else 0.5)) * (g - SP))
        if j + 1 < n:
            g = abs(r0 - outer_first[j + 1])
            ci = max(0.0, ((1.0 if (j + 1) in pinned or (j + 1) in shells else 0.5)) * (g - SP))
        ph = layer_phase + ring_i * (math.pi + drift)
        ring_i += 1
        e = np.clip(OVERSHOOT * np.sin(LOBES * th + ph), -1, 1)
        pred[j] = (r0 + np.where(e >= 0, e * co, e * ci), co, ci)
    return pred, outer_first, shells, pinned


if __name__ == "__main__":
    path, nw, obj = "landed.3mf", "delid_noweave.3mf", "L2_open50"
    nom, act, lh, li = actual(path, nw, obj, 7.0)
    naz = 720
    pred, outer_first, shells, pinned = model(nom, li, naz)
    print(f"layer {li} (z {li*lh:.2f}), {len(nom)} rings, ascending: "
          + " ".join(f"{r:.2f}" for r in nom))
    print(f"\n{'ring r':>8s} {'kind':>7s} {'ACTUAL min..max':>18s} {'MODEL min..max':>18s} "
          f"{'amp act':>8s} {'amp pred':>9s}")
    for j in range(len(nom)):
        r0 = outer_first[j]
        ai = len(nom) - 1 - j                     # ascending index of the same ring
        vals = np.array(sorted(act[ai].values())) if act[ai] else np.array([r0])
        p = pred[j][0]
        kind = "shell" if j in shells else ("pinned" if j in pinned else "weave")
        print(f"{r0:7.2f}mm {kind:>7s} {vals.min():8.3f}..{vals.max():<8.3f} "
              f"{p.min():8.3f}..{p.max():<8.3f} {vals.max()-vals.min():7.3f} "
              f"{p.max()-p.min():8.3f}")
