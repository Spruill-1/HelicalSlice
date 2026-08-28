"""Void CELLS, measured the way a cut face is read: how big is the biggest hole.

At density 25 the section is supposed to be mostly void - that is what foam is. The question the
eye answers on a cut face is not how much void there is but whether it is broken into cells or
runs as one region. So: rasterise the beads, take the complement inside the section, and report
the connected components by size.
"""
import math, re, sys, zipfile
import numpy as np
from scipy import ndimage

PIX = 0.04  # mm


def layer_segments(path, want, layer_index):
    with zipfile.ZipFile(path) as z:
        g = z.read("Metadata/plate_1.gcode").decode("utf-8", "replace")
    lh = float(re.search(r"; layer_height = ([\d.]+)", g).group(1))
    fd = float(re.search(r"; filament_diameter = ([\d.]+)", g).group(1))
    area = math.pi * fd * fd / 4.0
    segs, layer, cur = [], -1, None
    x = y = None
    for line in g.splitlines():
        if line.startswith("; CHANGE_LAYER"):
            layer += 1
            continue
        m = re.match(r"; printing object (\S+)", line)
        if m:
            cur = m.group(1)
            continue
        if line.startswith("; stop printing object"):
            cur = None
            continue
        if not (line.startswith("G1 ") or line.startswith("G0 ")):
            continue
        nx = re.search(r" X(-?[\d.]+)", line)
        ny = re.search(r" Y(-?[\d.]+)", line)
        ne = re.search(r" E(-?[\d.]+)", line)
        px, py = x, y
        if nx: x = float(nx.group(1))
        if ny: y = float(ny.group(1))
        if ne and float(ne.group(1)) > 0 and cur and want in cur and px is not None \
           and layer == layer_index:
            L = math.hypot(x - px, y - py)
            if L > 1e-9:
                segs.append((px, py, x, y, float(ne.group(1)) * area / (L * lh)))
    return segs


def analyse(path, want, layer_index, tag):
    segs = layer_segments(path, want, layer_index)
    if not segs:
        print(f"{tag}: no extrusion"); return
    xs = [s[0] for s in segs] + [s[2] for s in segs]
    ys = [s[1] for s in segs] + [s[3] for s in segs]
    x0, x1, y0, y1 = min(xs) - 1, max(xs) + 1, min(ys) - 1, max(ys) + 1
    W = int((x1 - x0) / PIX) + 1
    H = int((y1 - y0) / PIX) + 1
    mat = np.zeros((H, W), dtype=bool)
    yy, xx = np.mgrid[0:H, 0:W]
    px = x0 + xx * PIX
    py = y0 + yy * PIX
    for ax, ay, bx, by, w in segs:
        dx, dy = bx - ax, by - ay
        L2 = dx * dx + dy * dy
        t = np.clip(((px - ax) * dx + (py - ay) * dy) / L2, 0.0, 1.0)
        d2 = (px - (ax + t * dx)) ** 2 + (py - (ay + t * dy)) ** 2
        mat |= d2 <= (0.5 * w) ** 2

    # The SECTION is the material plus every void it encloses: fill from outside.
    outside = np.zeros_like(mat)
    outside[0, :] = outside[-1, :] = outside[:, 0] = outside[:, -1] = True
    outside = ndimage.binary_propagation(outside & ~mat, mask=~mat)
    void = (~mat) & (~outside)

    lab, n = ndimage.label(void)
    if n == 0:
        print(f"{tag}: no void at all"); return
    sizes = np.bincount(lab.ravel())[1:] * PIX * PIX
    order = np.argsort(sizes)[::-1]
    sect = (mat | void).sum() * PIX * PIX
    # widest point of the biggest void: twice the largest distance-to-material inside it
    dist = ndimage.distance_transform_edt(void, sampling=PIX)
    big = dist[lab == (order[0] + 1)].max() * 2.0
    print(f"{tag}")
    print(f"   section {sect:6.1f} mm2, material {mat.sum()*PIX*PIX/sect*100:4.1f}%, "
          f"{n} void cells")
    print(f"   biggest cell {sizes[order[0]]:6.2f} mm2 = {100*sizes[order[0]]/sect:4.1f}% of the "
          f"section, {big:.2f} mm across at its widest")
    print(f"   top 6 cells (mm2): " + " ".join(f"{sizes[k]:.2f}" for k in order[:6]))
    print(f"   median cell {np.median(sizes):.3f} mm2")


if __name__ == "__main__":
    path = sys.argv[1]
    layer = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    for want in sys.argv[3:]:
        analyse(path, want, layer, f"  {want} @ layer {layer}")
