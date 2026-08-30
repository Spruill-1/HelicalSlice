"""Sweep the weave's drive SHAPE, from a pure sinusoid to nearly a square wave.

    sweep_overshoot.py [object ...]

The weave displaces a ring by `clamp(overshoot * sin(phase), -1, 1)` times its local amplitude.
At overshoot 1 that is a pure sinusoid and never clamps: rings meet at POINTS, once per lobe,
which is the shape the foam hypothesis actually described. At 8 it is very nearly a square wave:
rings run alongside their neighbours for long arcs and meet over them, which is where the
honeycomb look comes from. 8 was chosen on one number - the bonded fraction went 22% to 29% -
and nothing else was measured at the time.

Measured here, per layer, from the emitted G-code rather than from ring dumps:

  weld      fraction of extruded length lying within a bead width of ANOTHER part of the same
            layer's path (at least 40 points away along it, so a bead's own neighbours on the
            path do not count as a weld). This is the bonded fraction the campaign tracks.
  void      the largest enclosed void in the LATTICE, in mm2 and at its widest point in mm.
            The part's own bore is excluded and printed separately: it is an enclosed region
            too, it is much the widest, and it has been mistaken for a lattice void three times
            in this campaign.
            Widest untied void is the number that has predicted every peel: 2.251 mm held,
            3.013 parted.
  sub-floor share of length below a printable bead - beads the volumetric cap had to thin.
"""
import math, os, re, subprocess, sys, zipfile
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

SUITE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CFG   = os.path.join(SUITE, "config")
REPO  = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE   = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE  = os.path.dirname(os.path.abspath(__file__))
OUT   = os.path.join(HERE, "_verify")
PIX, FLOOR = 0.04, 0.34


def slice_at(tag, obj, env_extra):
    name = f"ov_{tag}_{obj}.gcode.3mf"
    dst = os.path.join(OUT, name)
    if os.path.exists(dst):
        os.remove(dst)
    env = dict(os.environ); env.update(env_extra)
    cmd = [EXE, "--load-settings", f"{os.path.join(CFG,'machine.json')};{os.path.join(CFG,'cs_interp.json')}",
           "--load-filaments", os.path.join(CFG, "filament.json"),
           "--slice", "0", "--export-3mf", name, "--outputdir", OUT,
           os.path.join(SUITE, obj + ".stl")]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, env=env)
    return dst if r.returncode == 0 and os.path.exists(dst) else None


def segments(path, layer_index):
    with zipfile.ZipFile(path) as z:
        g = z.read("Metadata/plate_1.gcode").decode("utf-8", "replace")
    lh = float(re.search(r"; layer_height = ([\d.]+)", g).group(1))
    fd = float(re.search(r"; filament_diameter = ([\d.]+)", g).group(1))
    area = math.pi * fd * fd / 4.0
    segs, allsegs = [], []
    x = y = None; cur = None; lay = -1
    for line in g.splitlines():
        if line.startswith("; CHANGE_LAYER"):
            lay += 1; continue
        if re.match(r"; printing object", line):
            cur = 1; continue
        if line.startswith("; stop printing object"):
            cur = None; continue
        if not line.startswith("G1 "):
            continue
        nx = re.search(r" X(-?[\d.]+)", line); ny = re.search(r" Y(-?[\d.]+)", line)
        ne = re.search(r" E(-?[\d.]+)", line)
        px, py = x, y
        if nx: x = float(nx.group(1))
        if ny: y = float(ny.group(1))
        if not (cur and px is not None and ne and float(ne.group(1)) > 0):
            continue
        d = math.hypot(x - px, y - py)
        if d <= 1e-9:
            continue
        w = float(ne.group(1)) * area / (d * lh)
        allsegs.append((px, py, x, y, w, d))
        if lay == layer_index:
            segs.append((px, py, x, y, w, d))
    return segs, allsegs


def weld_fraction(segs, bead):
    """Length lying within a bead of another part of the same path, 40+ points away along it."""
    if len(segs) < 100:
        return 0.0
    P = np.array([[s[0], s[1]] for s in segs])
    L = np.array([s[5] for s in segs])
    t = cKDTree(P)
    pairs = t.query_ball_point(P, r=bead)
    hit = np.zeros(len(P), bool)
    for i, js in enumerate(pairs):
        for j in js:
            if abs(i - j) >= 40:
                hit[i] = True
                break
    return 100.0 * L[hit].sum() / L.sum()


def void_cells(segs):
    if not segs:
        return 0.0, 0.0
    xs = [s[0] for s in segs] + [s[2] for s in segs]
    ys = [s[1] for s in segs] + [s[3] for s in segs]
    x0, y0 = min(xs) - 1, min(ys) - 1
    W = int((max(xs) + 1 - x0) / PIX) + 1
    H = int((max(ys) + 1 - y0) / PIX) + 1
    yy, xx = np.mgrid[0:H, 0:W]
    px = x0 + xx * PIX; py = y0 + yy * PIX
    mat = np.zeros((H, W), bool)
    for ax, ay, bx, by, w, _ in segs:
        dx, dy = bx - ax, by - ay
        L2 = dx * dx + dy * dy
        t = np.clip(((px - ax) * dx + (py - ay) * dy) / L2, 0.0, 1.0)
        mat |= (px - (ax + t * dx)) ** 2 + (py - (ay + t * dy)) ** 2 <= (0.5 * w) ** 2
    out = np.zeros_like(mat)
    out[0, :] = out[-1, :] = out[:, 0] = out[:, -1] = True
    out = ndimage.binary_propagation(out & ~mat, mask=~mat)
    void = (~mat) & (~out)
    lab, n = ndimage.label(void)
    if n == 0:
        return 0.0, 0.0, 0.0
    dist = ndimage.distance_transform_edt(void, sampling=PIX)
    sizes = np.bincount(lab.ravel())[1:] * PIX * PIX
    # THE PART'S OWN BORE IS NOT A VOID. This detector has counted it three times in this
    # campaign, and it announces itself the same way every time: the enclosed region with much
    # the largest inscribed circle. Report it separately rather than letting it masquerade as
    # the lattice's worst cell - a 5 mm bore beside 0.8 mm cells is not a finding.
    widest = np.array([dist[lab == (i + 1)].max() * 2.0 for i in range(n)])
    bore = int(np.argmax(widest))
    others = [i for i in range(n) if i != bore]
    if not others:
        return float(widest[bore]), 0.0, 0.0
    j = max(others, key=lambda i: widest[i])
    return float(widest[bore]), float(sizes[j]), float(widest[j])


def main():
    objs = sys.argv[1:] or ["P3_ecc"]
    print(__doc__)
    for obj in objs:
        print(f"\n{obj}")
        print(f"  {'overshoot':>9s} {'weld %':>7s} {'lattice void mm2':>17s} {'widest mm':>10s} "
              f"{'sub-floor':>10s} {'path mm':>9s} {'(bore mm)':>10s}")
        for ov in ("1.0", "1.5", "2.0", "4.0", "6.0", "8.0"):
            p = slice_at(ov.replace(".", "_"), obj, {"CS_OVERSHOOT": ov})
            if not p:
                print(f"  {ov:>9s}   slice failed"); continue
            segs, allsegs = segments(p, 24)
            bead = float(np.median([s[4] for s in segs])) if segs else 0.4
            wf = weld_fraction(segs, bead)
            bore_w, va, vw = void_cells(segs)
            tot = sum(s[5] for s in allsegs)
            sub = sum(s[5] for s in allsegs if s[4] < FLOOR)
            print(f"  {ov:>9s} {wf:>6.1f}% {va:>17.2f} {vw:>10.2f} "
                  f"{100*sub/max(tot,1e-9):>9.2f}% {tot:>9.0f} {bore_w:>10.2f}")


if __name__ == "__main__":
    main()
