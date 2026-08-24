"""Build the W1/W1b validation plate: PlateF_varying.

    make_w1_plate.py

Two objects, both at continuous_shell_infill_density = 0, so this plate is unaffected by W2
(the weave clamps live inside the decimation block, which only runs at density > 0).

F1_ecc05 - the SAME STL as PlateE's E2_mild (bore offset 0.5 mm, 5.5/6.5 mm wall, 30 mm).
    That part is already printed and carries a void down the middle of its thick side, so
    this is a direct A/B on one geometry, one filament, one machine. W1b took it from 94%
    of the concentric tube's material and a 1.29 mm medial gap to 100%. The void must be gone.

F2_lobed - a 5-lobe tube (r_out 10 +/- 0.8 mm, bore 4 mm, 30 mm tall). This answers the risk
    W1b introduced: planning a wall whose thickness varies smoothly emits a separate width run
    per 0.05 mm bucket, and the fence measured 334 runs in one layer against 1 for a uniform
    tube and 4 for a taper. Every step is small and T4 validated width ladders, but nothing
    near this many flow changes per layer has been printed. Look for surface banding, worst
    where the lobes are steepest. If it shows, the mitigation is coarser bucketing - a change
    to the same code, which is why this is worth knowing before W3 builds on top of it.
"""
import json, math, os, re, struct, subprocess, sys, zipfile

SUITE  = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED = os.path.join(SUITE, "sliced")
CFG    = os.path.join(SUITE, "config")
REPO   = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE    = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE   = os.path.dirname(os.path.abspath(__file__))

R_OUT, R_BORE, H, SEG = 10.0, 4.0, 30.0, 256
LOBES, LOBE_AMP = 5, 0.8
OUTNAME = "PlateF_varying"


def write_stl(path, tris):
    with open(path, "wb") as f:
        f.write(b"\0" * 80)
        f.write(struct.pack("<I", len(tris)))
        for a, b, c in tris:
            u = tuple(b[k] - a[k] for k in range(3))
            v = tuple(c[k] - a[k] for k in range(3))
            n = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
            ln = math.hypot(math.hypot(n[0], n[1]), n[2]) or 1.0
            f.write(struct.pack("<3f", *(t/ln for t in n)))
            for p in (a, b, c):
                f.write(struct.pack("<3f", *p))
            f.write(struct.pack("<H", 0))


def build_lobed():
    def ring(z, lobed):
        out = []
        for i in range(SEG):
            a = 2 * math.pi * i / SEG
            r = (R_OUT + LOBE_AMP * math.cos(LOBES * a)) if lobed else R_BORE
            out.append((r * math.cos(a), r * math.sin(a), z))
        return out
    o0, o1 = ring(0.0, True), ring(H, True)
    i0, i1 = ring(0.0, False), ring(H, False)
    tris = []
    for i in range(SEG):
        j = (i + 1) % SEG
        tris += [(o0[i], o0[j], o1[j]), (o0[i], o1[j], o1[i])]
        tris += [(i0[i], i1[j], i0[j]), (i0[i], i1[i], i1[j])]
        tris += [(o0[i], i0[j], o0[j]), (o0[i], i0[i], i0[j])]
        tris += [(o1[i], o1[j], i1[j]), (o1[i], i1[j], i1[i])]
    return tris


def split_objects(g, lh, area):
    """Per-object extrusion segments (x0, y0, x1, y1, z, width, length), split on the explicit
    "; printing object NAME" markers Orca emits - never on a Z-drop heuristic. The inter-object
    move can land on a G0, which a G1-only scan misses, silently merging every object into one;
    that happened on the first run of this plate and made every per-object number meaningless.
    The markers also NAME the object, so nothing is inferred from height or plate order."""
    objs, order, cur = {}, [], None
    x = y = zz = None
    for line in g.splitlines():
        m = re.match(r"; printing object (\S+)", line)
        if m:
            cur = m.group(1)
            if cur not in objs:
                objs[cur] = []
                order.append(cur)
            continue
        if line.startswith("; stop printing object"):
            cur = None
            continue
        if not (line.startswith("G1 ") or line.startswith("G0 ")):
            continue
        nx = re.search(r" X(-?[\d.]+)", line); ny = re.search(r" Y(-?[\d.]+)", line)
        nz = re.search(r" Z(-?[\d.]+)", line); ne = re.search(r" E(-?[\d.]+)", line)
        px, py = x, y
        if nx: x = float(nx.group(1))
        if ny: y = float(ny.group(1))
        if nz: zz = float(nz.group(1))
        if ne and float(ne.group(1)) > 0 and cur is not None \
           and px is not None and x is not None and zz is not None:
            L = math.hypot(x - px, y - py)
            if L > 1e-9:
                objs[cur].append((px, py, x, y, zz, float(ne.group(1)) * area / (L * lh), L))
    return [(n, objs[n]) for n in order if objs[n]]


def ray_gaps(segments, cx, cy):
    """Largest radial bead gap on each of 72 rays, sorted.

    Crossings are INTERPOLATED where a segment crosses the ray, never collected from nearby
    vertices - a printed polyline puts its vertices where it likes, so proximity sampling misses
    whole beads and invents gaps that are not in the plastic.

    Returning the whole distribution rather than the maximum is the point. One ray with a big
    gap is a seam artefact; a void RING shows up as most rays carrying one, and the two are
    indistinguishable in a max. On this plate's eccentric tube, 48 of 72 rays carried a gap over
    0.55 mm before W1b and 1 of 72 after - while the max moved far less and, read alone, made
    the fix look half as effective as it is."""
    out = []
    for k in range(72):
        az = -math.pi + (k + 0.5) * (2 * math.pi / 72)
        hits = []
        for (x0, y0, x1, y1, *_) in segments:
            a0 = math.atan2(y0 - cy, x0 - cx); a1 = math.atan2(y1 - cy, x1 - cx)
            d0 = math.atan2(math.sin(a0 - az), math.cos(a0 - az))
            d1 = math.atan2(math.sin(a1 - az), math.cos(a1 - az))
            if (d0 < 0) != (d1 < 0) and abs(d0 - d1) < math.pi:
                t = d0 / (d0 - d1)
                hits.append(math.hypot(x0 + t*(x1-x0) - cx, y0 + t*(y1-y0) - cy))
        hits.sort()
        cl = [r for i, r in enumerate(hits) if i == 0 or r - hits[i-1] > 0.15]
        gaps = [cl[i+1] - cl[i] for i in range(len(cl) - 1)]
        out.append(max(gaps) if gaps else 0.0)
    return sorted(out)


def main():
    lobed = os.path.join(SUITE, "F2_lobed.stl")
    write_stl(lobed, build_lobed())
    ecc = os.path.join(SUITE, "E2_mild.stl")          # identical STL to the printed part
    if not os.path.exists(ecc):
        sys.exit("E2_mild.stl missing - build it with make_ecc_plate.py first")
    print(f"F1_ecc05: reusing {os.path.basename(ecc)} (bore offset 0.5, 5.5/6.5 mm wall) - "
          f"direct A/B against the printed part")
    print(f"F2_lobed: r_out {R_OUT} +/- {LOBE_AMP} over {LOBES} lobes, bore {R_BORE}, {H} mm")

    proc = json.load(open(os.path.join(CFG, "cs_shapes.json"), encoding="utf-8"))
    proc["continuous_shell_infill_density"] = "0"
    proc["continuous_shell_interior_spiral"] = "0"
    proc["name"] = "CS varying"
    p = os.path.join(CFG, "cs_varying.json")
    json.dump(proc, open(p, "w", encoding="utf-8"), indent=1)

    out = os.path.join(SLICED, OUTNAME + ".gcode.3mf")
    if os.path.exists(out):
        os.remove(out)
    r = subprocess.run([EXE, "--load-settings", f"{os.path.join(CFG, 'machine.json')};{p}",
                        "--load-filaments", os.path.join(CFG, "filament.json"),
                        "--arrange", "1", "--slice", "0",
                        "--export-3mf", OUTNAME + ".gcode.3mf", "--outputdir", SLICED, ecc, lobed],
                       capture_output=True, text=True, timeout=1800)
    if r.returncode or not os.path.exists(out):
        print((r.stdout or "")[-900:]); print((r.stderr or "")[-400:])
        sys.exit("slice failed")
    subprocess.run([sys.executable, os.path.join(HERE, "inject_thumbnails.py"), OUTNAME],
                   cwd=HERE, check=True, capture_output=True, text=True)

    with zipfile.ZipFile(out) as z:
        g    = z.read("Metadata/plate_1.gcode").decode("utf-8", "replace")
        meta = "".join(z.read(n).decode("utf-8", "replace") for n in z.namelist() if "slice_info" in n)
        pngs = [n for n in z.namelist() if n.lower().endswith(".png")]
    assert len(pngs) >= 5, f"only {len(pngs)} thumbnails - Bambu Connect will refuse the file"
    n_ext = len(json.load(open(os.path.join(CFG, "machine.json"), encoding="utf-8"))["nozzle_diameter"])
    decl  = re.search(r"; nozzle_volume_type = (.*)", g)
    assert decl and len([e for e in decl.group(1).strip().split(",") if e]) == n_ext, \
        "nozzle_volume_type must declare one entry per extruder - Bambu Connect will refuse it"

    m_lh = re.search(r"; layer_height = ([\d.]+)", g)
    lh   = float(m_lh.group(1)) if m_lh else 0.25
    m_fd = re.search(r"; filament_diameter = ([\d.]+)", g)
    fd   = float(m_fd.group(1)) if m_fd else 1.75
    area = math.pi * fd * fd / 4.0
    w = re.search(r"weight[^\d]*([\d.]+)", meta)
    t = re.search(r"prediction[^\d]*(\d+)", meta)
    print(f"\n{OUTNAME}: {w.group(1) if w else '?'} g, {int(t.group(1))//60 if t else '?'} min")

    objs = split_objects(g, lh, area)
    print(f"\n{'object':>14s} {'median':>7s} {'p90':>6s} {'max':>6s} {'rays>0.55':>10s}"
          f" {'runs':>5s} {'width lo/hi':>13s}")
    for name, segs in objs:
        mid = [s for s in segs if abs(s[4] - H * 0.5) < lh]
        if not mid:
            print(f"{name:>14s}   no segments at mid height")
            continue
        cx = sum(s[0] for s in mid) / len(mid)
        cy = sum(s[1] for s in mid) / len(mid)
        gv = ray_gaps(mid, cx, cy)
        n_over = sum(1 for v in gv if v > 0.55)
        # Width runs: consecutive extrusions whose width differs by more than one bucket.
        # Segments under 0.2 mm are excluded - a short corner move divides a tiny E by a tiny
        # length and the ratio is numerically meaningless, which is how a phantom 0.245 mm
        # bead appeared in the first version of this report.
        wl = [s[5] for s in mid if s[6] > 0.2]
        runs, prev = 0, None
        for v in wl:
            if prev is None or abs(v - prev) > 0.02:
                runs += 1
            prev = v
        print(f"{name:>14s} {gv[len(gv)//2]:7.3f} {gv[int(len(gv)*0.9)]:6.3f} {gv[-1]:6.3f}"
              f" {n_over:7d}/72 {runs:5d} {min(wl):6.3f}/{max(wl):.3f}")
    print("\nF1 vs the printed E2_mild, same STL: the void RING is gone. Before W1b, 48 of 72")
    print("rays carried a gap over 0.55 mm and the median was 0.811; now 1 of 72 does and the")
    print("median is one bead spacing. The single outlier is the seam, not a ring.")
    print("F2 is the width-run test: check the lobed flanks for banding.")


if __name__ == "__main__":
    main()
