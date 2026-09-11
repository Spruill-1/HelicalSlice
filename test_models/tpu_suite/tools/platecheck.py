"""Per-object print-quality metrics from a sliced plate, tuned to catch the two failure modes the
arc validation plate (PlateMB_arcpeel) exposed on plastic that the earlier sim metrics missed:

  BRIDGES  - extruding spans that jump a gap (unsupported filament). Arc bridged ~4x the fraction
             path; the old continuity check only counted NON-extruding travels, so every bridge
             (E>0) passed it clean and arc scored "0 teleports". Here a bridge is an extruding
             segment whose straight span exceeds BRIDGE_MM. Bore-crossing bridges are called out
             separately - those are the worst (a bead strung across an open bore).

  WELD-GAP - the largest bead-free pocket inside the material (its inscribed diameter is the gap a
             bead there must weld across). Measured against ALL beads, bridges INCLUDED: a bridge is
             a real bead and bridging is its own metric above; the earlier version excluded bridges,
             which turned every bridge SPAN into a phantom void that beat the real inter-bore pocket.
             The domain is the bead footprint WITHOUT expanding past the outermost bead (a `close`
             rim otherwise reads as a spurious edge void). Flags fraction's inter-bore pinch and
             (with the bridge metric) arc's inter-column holes; passes the coupons that welded.

Reads the sliced gcode (what we screen BEFORE printing). Validated to reproduce the plate verdict:
arc = heavy bridges + big voids; fraction = few bridges but a central saddle void.

    python platecheck.py <plate.gcode.3mf> [layer=24] [obj_substr ...]
"""
import sys, re, zipfile, math
import numpy as np

BRIDGE_MM   = 2.0      # extruding span above this is an unsupported bridge, not a bead
G           = 0.20     # raster cell (mm)
BEAD_R      = 0.225    # half a 0.45 bead
CLOSE_MM    = 1.20     # morphological close radius for the material footprint (< bore radius, so bores survive)
FIL_D       = 1.75
FIL_AREA    = math.pi * (FIL_D / 2.0) ** 2
DENSITY     = 1.24e-3  # g/mm^3 TPU


def layer_segments(gcode, want_layer):
    """{obj: [(x0,y0,x1,y1, extruding, fil_mm)]} for one OBJECT-LOCAL layer. By-object sequencing
    makes CHANGE_LAYER a global counter, so each object's layer index is offset from its first
    layer; want_layer is the index within each object. Relative-E (M83) gcode."""
    cur = None; x = y = None; lay = -1; out = {}; first = {}
    for line in gcode.splitlines():
        if line.startswith('; CHANGE_LAYER'): lay += 1; continue
        mo = re.match(r'; printing object (.+?) id:\d+ copy \d+', line)
        if mo:
            cur = mo.group(1).replace('.stl', ''); out.setdefault(cur, [])
            first.setdefault(cur, lay); continue
        if line.startswith('; stop printing object'): cur = None; continue
        if not line.startswith('G1 '): continue
        nx = re.search(r' X(-?[\d.]+)', line); ny = re.search(r' Y(-?[\d.]+)', line); ne = re.search(r' E(-?[\d.]+)', line)
        px, py = x, y
        if nx: x = float(nx.group(1))
        if ny: y = float(ny.group(1))
        if cur is None or px is None: continue
        if lay - first[cur] != want_layer: continue
        moved = bool(nx or ny)
        e = float(ne.group(1)) if ne else 0.0
        if moved and e > 0:
            out[cur].append((px, py, x, y, math.hypot((x or px) - px, (y or py) - py), e))
    return out


def raster(segs, only_beads, xs, ys):
    """boolean coverage grid; only_beads=True rasterizes segments shorter than BRIDGE_MM only."""
    cov = np.zeros((len(xs), len(ys)), bool)
    for (x0, y0, x1, y1, span, _e) in segs:
        if only_beads and span >= BRIDGE_MM: continue
        n = max(2, int(span / (G * 0.7)))
        for t in np.linspace(0, 1, n):
            cx = x0 + t * (x1 - x0); cy = y0 + t * (y1 - y0)
            i0 = max(0, int((cx - BEAD_R - xs[0]) / G)); i1 = min(len(xs), int((cx + BEAD_R - xs[0]) / G) + 1)
            j0 = max(0, int((cy - BEAD_R - ys[0]) / G)); j1 = min(len(ys), int((cy + BEAD_R - ys[0]) / G) + 1)
            for i in range(i0, i1):
                for j in range(j0, j1):
                    if (xs[i] - cx) ** 2 + (ys[j] - cy) ** 2 <= BEAD_R * BEAD_R: cov[i, j] = True
    return cov


def dilate(m, r_cells):
    for _ in range(r_cells):
        g = m.copy(); g[1:] |= m[:-1]; g[:-1] |= m[1:]; g[:, 1:] |= m[:, :-1]; g[:, :-1] |= m[:, 1:]; m = g
    return m


def erode(m, r_cells):
    return ~dilate(~m, r_cells)


def components(mask):
    """connected-component sizes+cells of a boolean mask (4-connectivity)."""
    from collections import deque
    lbl = np.zeros(mask.shape, int); comps = []; nxt = 0
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            if mask[i, j] and lbl[i, j] == 0:
                nxt += 1; q = deque([(i, j)]); lbl[i, j] = nxt; cells = []
                while q:
                    a, b = q.popleft(); cells.append((a, b))
                    for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        na, nb = a + da, b + db
                        if 0 <= na < mask.shape[0] and 0 <= nb < mask.shape[1] and mask[na, nb] and lbl[na, nb] == 0:
                            lbl[na, nb] = nxt; q.append((na, nb))
                comps.append(cells)
    return comps


def touches_border(cells, shape):
    return any(a == 0 or b == 0 or a == shape[0] - 1 or b == shape[1] - 1 for a, b in cells)


def void_width_mm(cells, shape):
    """diameter of the largest inscribed disc in the void = 2x the max distance-to-edge, via erosion peel."""
    m = np.zeros(shape, bool)
    for a, b in cells: m[a, b] = True
    k = 0
    while m.any() and k < 60:
        m = erode(m, 1); k += 1
    return 2 * k * G


def analyze(name, segs):
    if not segs: return None
    xs0 = min(min(s[0], s[2]) for s in segs); xs1 = max(max(s[0], s[2]) for s in segs)
    ys0 = min(min(s[1], s[3]) for s in segs); ys1 = max(max(s[1], s[3]) for s in segs)
    xs = np.arange(xs0 - 1, xs1 + 1, G); ys = np.arange(ys0 - 1, ys1 + 1, G)
    cov = raster(segs, only_beads=False, xs=xs, ys=ys)    # ALL beads (bridges included - see below)
    # DOMAIN = the material cross-section WITHOUT expanding past the outermost bead. Fill the small
    # foam gaps the beads enclose (they are material); keep the bores (large enclosed holes) and the
    # outside open. A footprint from `close` pushes the boundary ~CLOSE_MM past the beads and its rim
    # then reads as a phantom edge void.
    domain = cov.copy(); bore = np.zeros(cov.shape, bool)
    for c in components(~cov):
        if touches_border(c, cov.shape): continue         # the outside
        if len(c) * G * G > 2.0:                           # a bore: keep it as a hole
            for a, b in c: bore[a, b] = True
        else:
            for a, b in c: domain[a, b] = True             # a small foam gap: it is material
    # WELD-GAP = the largest pocket the beads leave inside the material, measured against ALL beads.
    # Bridges are INCLUDED here: a bridge IS a real bead, and bridging is already its own metric;
    # excluding bridges (the old cov_bead) turned every bridge span into a phantom void that beat the
    # real inter-bore pocket. The pocket's inscribed diameter is the gap a bead there must weld across.
    void = domain & ~cov
    vcomps = sorted((c for c in components(void) if not touches_border(c, void.shape)), key=len, reverse=True)
    big = vcomps[0] if vcomps else []
    big_ctr = (float(np.mean([xs[a] for a, b in big])), float(np.mean([ys[b] for a, b in big]))) if big else (0.0, 0.0)
    # bridges
    bridges = [s for s in segs if s[4] >= BRIDGE_MM]
    bore_cross = 0
    for (x0, y0, x1, y1, span, _e) in bridges:
        n = max(2, int(span / 0.3)); hit = False
        for t in np.linspace(0.15, 0.85, n):        # sample the interior of the span
            cx = x0 + t * (x1 - x0); cy = y0 + t * (y1 - y0)
            i = int((cx - xs[0]) / G); j = int((cy - ys[0]) / G)
            if 0 <= i < len(xs) and 0 <= j < len(ys) and bore[i, j]: hit = True; break
        if hit: bore_cross += 1
    fil_mm = sum(s[5] for s in segs)
    return {
        "weight_g":      fil_mm * FIL_AREA * DENSITY,
        "n_bridge":      len(bridges),
        "bridge_mm":     sum(s[4] for s in bridges),
        "bridge_max":    max((s[4] for s in bridges), default=0.0),
        "bore_cross":    bore_cross,
        "worst_void_mm2": len(big) * G * G,
        "worst_void_w":  void_width_mm(big, void.shape) if big else 0.0,
        "worst_void_at": big_ctr,
        "n_voids":       len(vcomps),
        "domain_mm2":    int(domain.sum()) * G * G,
    }


def main():
    path = sys.argv[1]
    layer = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 24
    wants = [a for a in sys.argv[3:]]
    with zipfile.ZipFile(path) as z:
        nm = next(n for n in z.namelist() if n.endswith('plate_1.gcode'))
        g = z.read(nm).decode('utf-8', 'replace')
    perobj = layer_segments(g, layer)
    print(f"layer {layer}  |  BRIDGE_MM={BRIDGE_MM}  bead={2*BEAD_R}mm\n")
    hdr = f"{'object':14s} {'wt_g':>5s} {'bridges':>8s} {'br_mm':>6s} {'br_max':>6s} {'boreX':>5s} {'worst_void':>10s} {'void_w':>6s}  verdict"
    print(hdr); print("-" * len(hdr))
    for name in sorted(perobj):
        if wants and not any(w in name for w in wants): continue
        r = analyze(name, perobj[name])
        if r is None: continue
        flags = []
        if r['n_bridge'] > 10: flags.append("BRIDGY")
        if r['bore_cross'] > 0: flags.append("BORE-BRIDGE")
        if r['worst_void_w'] > 1.2: flags.append("WELD-GAP")     # gap wider than the weld distance
        verdict = " ".join(flags) if flags else "ok"
        print(f"{name:14s} {r['weight_g']:5.2f} {r['n_bridge']:8d} {r['bridge_mm']:6.0f} {r['bridge_max']:6.1f} "
              f"{r['bore_cross']:5d} {r['worst_void_mm2']:8.1f}mm2 {r['worst_void_w']:5.2f}  {verdict}  "
              f"@({r['worst_void_at'][0]:.1f},{r['worst_void_at'][1]:.1f})")


if __name__ == '__main__':
    main()
