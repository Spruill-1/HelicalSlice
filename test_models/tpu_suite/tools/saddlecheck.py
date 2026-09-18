"""Two print-quality metrics for the inter-bore region, validated against PlateMB_saddle's six
printed coupons (testlog/platemb_saddle/PlateSummary.md). void_w (platecheck's single-inscribed-disc)
could not tell a coupon that PEELS from one that holds, nor one that went SOLID from one that stayed
squishy - the two axes the plate actually judged on. These two, together, reproduce the plate ranking:

  COMPLIANCE  = solid-fraction: the fraction of the material cross-section covered by beads. TPU's
                value is squish, so a HIGH solid-fraction is a part that lost its compliance. F7_son
                (printed "basically solid, no squish") scores 0.89; every squishy coupon is <= 0.79.
                > SOLID_MAX flags "closed the peel but went solid."

  PEEL-RISK   = inter-bore band coverage: bead coverage in a BAND_MM band around the line joining the
                two nearest bores - the strip the plate peeled along. The two coupons that peeled
                (F7_soff 0.47, MB2e_soff 0.51) are the two lowest; every coupon that held is >= 0.54.
                < PEEL_MIN flags a peel-prone under-dense strip - which void_w missed on F7_soff.

Bores are the enclosed holes that are BIG (> BORE_MM2) AND ROUND (compactness > 0.5); foam gaps are
small or stringy and are counted as material. Metrics are averaged over several mid-height layers to
damp single-layer raster noise. The thresholds are the midpoints of the plate's peel/hold and
solid/squishy gaps - re-calibrate them if the material or nozzle changes.

    python saddlecheck.py <plate.gcode.3mf> [obj_substr ...]
"""
import sys, zipfile, math
import numpy as np
import platecheck as pc

BORE_MM2 = 8.0      # an enclosed hole bigger than this AND round is a bore, not a foam gap
BAND_MM  = 1.6      # half-width of the inter-bore band
LAYERS   = [16, 20, 24, 28, 32]   # averaged, to damp per-layer raster noise
# Thresholds calibrated on PlateMB_saddle (averaged over LAYERS): peel F7_soff 0.50 / MB2e_soff 0.55,
# hold MB2u_soff 0.58, so PEEL_MIN 0.56 flags both peels and passes the holds. The metric is ROBUST at
# the extremes (F7's dramatic peel at 0.50, and the solid F7_son at 0.89 vs squishy <=0.80) and THIN in
# the middle (mild-peel MB2e_soff 0.55 vs hold MB2u_soff 0.58 is within raster noise) - for a borderline
# coupon the physical peel stays the arbiter. Re-calibrate if the material/nozzle/bead width changes.
PEEL_MIN = 0.56     # inter-bore band coverage below this = peel-prone
SOLID_MAX = 0.85    # solid-fraction above this = lost its squish


def bore_split(cov, xs, ys):
    """material mask (foam gaps kept, bores removed) and the bore centroids."""
    domain = cov.copy(); bore = np.zeros(cov.shape, bool); bores = []
    for c in pc.components(~cov):
        if pc.touches_border(c, cov.shape):
            continue
        a = len(c) * pc.G * pc.G
        cells = set(c); perim = 0
        for (i, j) in c:
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                if (i + di, j + dj) not in cells:
                    perim += pc.G
        comp = 4 * math.pi * a / (perim * perim) if perim > 0 else 0.0
        if a > BORE_MM2 and comp > 0.5:                     # big and round = a bore
            for i, j in c: bore[i, j] = True
            bores.append((float(np.mean([xs[i] for i, j in c])),
                          float(np.mean([ys[j] for i, j in c]))))
        else:
            for i, j in c: domain[i, j] = True              # small/stringy = material
    return (domain & ~bore), bores


def layer_metrics(segs):
    if not segs:
        return None
    xs0 = min(min(s[0], s[2]) for s in segs); xs1 = max(max(s[0], s[2]) for s in segs)
    ys0 = min(min(s[1], s[3]) for s in segs); ys1 = max(max(s[1], s[3]) for s in segs)
    xs = np.arange(xs0 - 1, xs1 + 1, pc.G); ys = np.arange(ys0 - 1, ys1 + 1, pc.G)
    cov = pc.raster(segs, only_beads=False, xs=xs, ys=ys)
    mat, bores = bore_split(cov, xs, ys)
    if not mat.any():
        return None
    cov_all = float(cov[mat].mean())                        # solid-fraction (compliance)
    cov_ib = float('nan')
    if len(bores) >= 2:
        XX, YY = np.meshgrid(xs, ys, indexing='ij')
        best = None
        for i in range(len(bores)):
            for j in range(i + 1, len(bores)):
                A = np.array(bores[i]); B = np.array(bores[j]); AB = B - A; L2 = float(AB @ AB)
                if L2 < 1e-9:
                    continue
                t = np.clip(((np.stack([XX, YY], -1) - A) @ AB) / L2, 0, 1)
                proj = A + t[..., None] * AB
                dist = np.hypot(XX - proj[..., 0], YY - proj[..., 1])
                band = mat & (dist <= BAND_MM)
                if band.any():
                    c = float(cov[band].mean())
                    if best is None or c < best:            # weakest inter-bore line
                        best = c
        cov_ib = best if best is not None else float('nan')
    return cov_all, cov_ib, len(bores)


def analyze(name, per_layer):
    accA, accI, nb = [], [], 0
    for segs in per_layer:
        m = layer_metrics(segs)
        if m is None:
            continue
        accA.append(m[0])
        if m[1] == m[1]:
            accI.append(m[1])
        nb = max(nb, m[2])
    if not accA:
        return None
    return {"cov_all": float(np.mean(accA)),
            "cov_ibore": float(np.mean(accI)) if accI else float('nan'),
            "n_bore": nb}


def main():
    path = sys.argv[1]
    wants = [a for a in sys.argv[2:]]
    with zipfile.ZipFile(path) as z:
        nm = next(n for n in z.namelist() if n.endswith('plate_1.gcode'))
        g = z.read(nm).decode('utf-8', 'replace')
    per = {L: pc.layer_segments(g, L) for L in LAYERS}
    names = sorted({n for L in LAYERS for n in per[L]})
    print(f"layers {LAYERS}  |  PEEL_MIN={PEEL_MIN}  SOLID_MAX={SOLID_MAX}\n")
    hdr = f"{'object':14s} {'bores':>5s} {'solid_frac':>10s} {'ibore_cov':>9s}  verdict"
    print(hdr); print("-" * len(hdr))
    for name in names:
        if wants and not any(w in name for w in wants):
            continue
        r = analyze(name, [per[L].get(name, []) for L in LAYERS])
        if r is None:
            continue
        flags = []
        if r['n_bore'] >= 2 and r['cov_ibore'] == r['cov_ibore'] and r['cov_ibore'] < PEEL_MIN:
            flags.append("PEEL-RISK")
        if r['cov_all'] > SOLID_MAX:
            flags.append("OVER-SOLID")
        verdict = " ".join(flags) if flags else "ok"
        ib = f"{r['cov_ibore']:.2f}" if r['cov_ibore'] == r['cov_ibore'] else "  - "
        print(f"{name:14s} {r['n_bore']:>5d} {r['cov_all']:>10.2f} {ib:>9s}  {verdict}")


if __name__ == '__main__':
    main()
