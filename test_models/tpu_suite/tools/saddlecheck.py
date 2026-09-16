"""Inter-bore SADDLE-STRIP density, the defect platecheck's single-pocket void_w is blind to.

The printed plate's fraction coupons failed the same way: "the area in between the two bores is
much less dense... and peels too easily" (MB2e_frac), "space directly between the bores is not
dense enough" (F7_frac). That is a distributed under-dense STRIP where the u-level-sets spread
(|grad u|->0 at the saddle), not a single fillable pocket - so void_w (largest inscribed bead-free
disc) barely moves while the strip stays weak.

This measures bead COVERAGE in a band around the inter-bore saddle line and compares it to coverage
in the rest of the material. The saddle line is taken between the two nearest bore centroids; the
band is +/-BAND_MM of it. A well-knit saddle scores near the body coverage (ratio ~1.0); the peel
line shows as ratio < 1. Reports coverage_saddle, coverage_body, and their ratio per object.

    python saddlecheck.py <plate.gcode.3mf> [layer=24] [obj_substr ...]
"""
import sys, zipfile, math
import numpy as np
import platecheck as pc

BAND_MM = 1.6     # half-width of the saddle band (a couple of bead spacings each side of the line)


def bore_centroids(cov, xs, ys):
    """centroids (mm) and cell-sets of the large enclosed holes = bores, reusing platecheck's rule."""
    out = []
    for c in pc.components(~cov):
        if pc.touches_border(c, cov.shape):
            continue
        if len(c) * pc.G * pc.G > 2.0:      # a bore (same threshold as platecheck's domain)
            cx = float(np.mean([xs[a] for a, b in c])); cy = float(np.mean([ys[b] for a, b in c]))
            out.append(((cx, cy), c))
    return out


def analyze_saddle(name, segs):
    if not segs:
        return None
    xs0 = min(min(s[0], s[2]) for s in segs); xs1 = max(max(s[0], s[2]) for s in segs)
    ys0 = min(min(s[1], s[3]) for s in segs); ys1 = max(max(s[1], s[3]) for s in segs)
    xs = np.arange(xs0 - 1, xs1 + 1, pc.G); ys = np.arange(ys0 - 1, ys1 + 1, pc.G)
    cov = pc.raster(segs, only_beads=False, xs=xs, ys=ys)
    # material domain = cov with small enclosed foam gaps filled, bores kept open (platecheck's rule)
    domain = cov.copy(); bore = np.zeros(cov.shape, bool)
    for c in pc.components(~cov):
        if pc.touches_border(c, cov.shape):
            continue
        if len(c) * pc.G * pc.G > 2.0:
            for a, b in c: bore[a, b] = True
        else:
            for a, b in c: domain[a, b] = True
    bores = bore_centroids(cov, xs, ys)
    if len(bores) < 2:
        return {"n_bore": len(bores)}
    # the saddle line = between the two nearest bore centroids
    cens = [c for c, _ in bores]
    best = None
    for i in range(len(cens)):
        for j in range(i + 1, len(cens)):
            d = math.hypot(cens[i][0] - cens[j][0], cens[i][1] - cens[j][1])
            if best is None or d < best[0]:
                best = (d, cens[i], cens[j])
    _, A, B = best
    A = np.array(A); B = np.array(B); AB = B - A; L2 = float(AB @ AB)
    # per material cell: distance to the saddle SEGMENT; band = within BAND_MM
    XX, YY = np.meshgrid(xs, ys, indexing='ij')
    P = np.stack([XX, YY], axis=-1)
    t = np.clip(((P - A) @ AB) / max(L2, 1e-9), 0.0, 1.0)
    proj = A + t[..., None] * AB
    dist = np.hypot(P[..., 0] - proj[..., 0], P[..., 1] - proj[..., 1])
    mat = domain & ~bore                      # solid material cells (exclude bores)
    band = mat & (dist <= BAND_MM)
    body = mat & (dist > BAND_MM)
    cov_saddle = float(cov[band].mean()) if band.any() else float('nan')
    cov_body = float(cov[body].mean()) if body.any() else float('nan')
    ratio = cov_saddle / cov_body if cov_body and cov_body == cov_body else float('nan')
    return {"n_bore": len(bores), "cov_saddle": cov_saddle, "cov_body": cov_body,
            "ratio": ratio, "band_mm2": int(band.sum()) * pc.G * pc.G,
            "saddle_mid": ((A + B) / 2).tolist()}


def main():
    path = sys.argv[1]
    layer = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 24
    wants = [a for a in sys.argv[3:]]
    with zipfile.ZipFile(path) as z:
        nm = next(n for n in z.namelist() if n.endswith('plate_1.gcode'))
        g = z.read(nm).decode('utf-8', 'replace')
    per = pc.layer_segments(g, layer)
    print(f"layer {layer}  |  BAND_MM={BAND_MM}  (saddle coverage vs body coverage)\n")
    hdr = f"{'object':14s} {'bores':>5s} {'cov_saddle':>10s} {'cov_body':>9s} {'ratio':>6s}  verdict"
    print(hdr); print("-" * len(hdr))
    for name in sorted(per):
        if wants and not any(w in name for w in wants):
            continue
        r = analyze_saddle(name, per[name])
        if r is None:
            continue
        if r.get("n_bore", 0) < 2:
            print(f"{name:14s} {r.get('n_bore',0):>5d}   (need >=2 bores)")
            continue
        verdict = "SADDLE-THIN" if r["ratio"] < 0.85 else "ok"
        print(f"{name:14s} {r['n_bore']:>5d} {r['cov_saddle']:>10.2f} {r['cov_body']:>9.2f} "
              f"{r['ratio']:>6.2f}  {verdict}")


if __name__ == '__main__':
    main()
