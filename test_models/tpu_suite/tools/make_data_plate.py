"""Build the density>0 DATA plate: two objects that answer open questions with today's code.

    make_data_plate.py [--density N]

D1_ladder - a tube whose wall steps 8,7,6,5,4 beads going up (the bore grows outward, so every
step is an upward-facing supported ledge, not an overhang). The sparsifier's annular gate is
`n > 2*wall_loops + 1`, so at wall_loops 2 the wall stops being sparsified between the 6-bead
and 5-bead steps. One object therefore brackets the gate and shows both sides of it in the same
print, same filament, same session.

It also brackets the OTHER defect: the 5-bead step is the S5 base wall, where C-lite lays four
stretched beads instead of five. So the same object measures the sparsify gate AND the adaptive
-beading regression, and after the decoupling the 5-bead step must go 4.0 -> 5.0 beads per ray
while the 8/7/6 steps stay sparsified.

D2_weave - a THICK ANNULUS (r_out 10 / bore 2.5), which exercises the full woven lattice and
has real boundaries at both faces. It replaced a solid cylinder that could not do the job: a
solid chain legitimately runs its terminus ring to the centre, so neither the innermost
extrusion radius nor the smallest adjacent-ring gap separates a collapsed ring from normal
convergence, and four probes failed for that reason. What a bad weave clamp looks like is
bead-on-bead against a SHELL, and that shows physically as over-extrusion - a question for the
printer rather than for a G-code script.

Both are the C-lite baseline that the post-decoupling reprints get compared against.
"""
import json, math, os, re, struct, subprocess, sys, zipfile

SUITE  = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED = os.path.join(SUITE, "sliced")
CFG    = os.path.join(SUITE, "config")
REPO   = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE    = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE   = os.path.dirname(os.path.abspath(__file__))

SP      = 0.3963                       # bead spacing at 0.4 nozzle, this profile
R_OUT   = 16.0
STEPS   = [8, 7, 6, 5, 4]              # beads of wall, bottom -> top
STEP_H  = 7.0
CORE_R, CORE_BORE, CORE_H = 10.0, 2.5, 20.0   # 7.5 mm wall: ~19 beads, well past the gate
SEG     = 256
OUTNAME = "PlateD_data"


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


def ring(r, z):
    return [(r*math.cos(2*math.pi*i/SEG), r*math.sin(2*math.pi*i/SEG), z) for i in range(SEG)]


def build_ladder():
    tris = []
    zs = [i * STEP_H for i in range(len(STEPS) + 1)]
    bores = [R_OUT - n * SP for n in STEPS]
    o0, o1 = ring(R_OUT, 0.0), ring(R_OUT, zs[-1])
    for i in range(SEG):                                    # outer skin
        j = (i + 1) % SEG
        tris += [(o0[i], o0[j], o1[j]), (o0[i], o1[j], o1[i])]
    b0 = ring(bores[0], 0.0)
    for i in range(SEG):                                    # bottom annulus
        j = (i + 1) % SEG
        tris += [(o0[i], b0[j], o0[j]), (o0[i], b0[i], b0[j])]
    bt = ring(bores[-1], zs[-1])
    for i in range(SEG):                                    # top annulus
        j = (i + 1) % SEG
        tris += [(o1[i], o1[j], bt[j]), (o1[i], bt[j], bt[i])]
    for s, r in enumerate(bores):
        lo, hi = ring(r, zs[s]), ring(r, zs[s + 1])
        for i in range(SEG):                                # bore wall of this step
            j = (i + 1) % SEG
            tris += [(lo[i], hi[j], lo[j]), (lo[i], hi[i], hi[j])]
        if s + 1 < len(bores):                              # upward-facing ledge to the wider bore
            a, b = ring(r, zs[s + 1]), ring(bores[s + 1], zs[s + 1])
            for i in range(SEG):
                j = (i + 1) % SEG
                tris += [(a[i], a[j], b[j]), (a[i], b[j], b[i])]
    return tris


def build_core():
    """D2_weave: a THICK ANNULUS, not the solid cylinder this slot used to hold.

    The solid cylinder was meant to exercise the weave's missing inward clamp, and it cannot:
    a solid chain legitimately runs its terminus ring to the centre, so neither the innermost
    extrusion radius nor the smallest adjacent-ring gap can separate a collapsed ring from
    normal convergence. Four probes were tried and all failed for that reason.

    A thick annulus has real boundaries at BOTH faces, exercises the full woven lattice, and is
    the geometry where bead-on-bead against a shell would show up physically as over-extrusion.
    That is the W2 question a print can answer."""
    tris = []
    o0, o1 = ring(CORE_R, 0.0), ring(CORE_R, CORE_H)
    i0, i1 = ring(CORE_BORE, 0.0), ring(CORE_BORE, CORE_H)
    for i in range(SEG):
        j = (i + 1) % SEG
        tris += [(o0[i], o0[j], o1[j]), (o0[i], o1[j], o1[i])]
        tris += [(i0[i], i1[j], i0[j]), (i0[i], i1[i], i1[j])]
        tris += [(o0[i], i0[j], o0[j]), (o0[i], i0[i], i0[j])]
        tris += [(o1[i], o1[j], i1[j]), (o1[i], i1[j], i1[i])]
    return tris


def split_objects(g):
    """Per-object extrusion SEGMENTS, split on the explicit "; printing object NAME" markers
    Orca emits - never on a Z-drop heuristic. The inter-object move can land on a G0, which a
    G1-only scan misses, silently merging every object into one; that happened on the first
    run of PlateF and made every per-object number meaningless. The markers also NAME the
    object, so nothing has to be inferred from height or plate order."""
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
        if ne and float(ne.group(1)) > 0 and cur is not None            and px is not None and x is not None and zz is not None:
            objs[cur].append((px, py, x, y, zz))
    return [(n, objs[n]) for n in order if objs[n]]


def min_ring_gap(segments, cx, cy):
    """Smallest radial gap between adjacent rings, over 72 rays.

    This is the detector for the weave's missing inward clamp on a SOLID chain. The obvious
    probe - the innermost extrusion radius - cannot see it, because a solid centre is filled
    either way: it read 0.2412 mm at density 0 and 0.2563 mm at density 25 while the bug was
    live. A ring collapsing inward instead shows up as it CROSSING or touching its neighbour,
    i.e. an adjacent-ring gap heading to zero, which is what this measures."""
    worst = 1e9
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
        for i in range(len(hits) - 1):
            d = hits[i+1] - hits[i]
            if d > 0.02:                      # ignore a ray clipping one bead twice
                worst = min(worst, d)
    return worst


def crossings_per_ray(pts, cx, cy):
    counts = []
    for k in range(24):
        az = -math.pi + (k + 0.5) * (2 * math.pi / 24)
        rr = sorted(math.hypot(p[0] - cx, p[1] - cy) for p in pts
                    if abs(math.atan2(math.sin(math.atan2(p[1] - cy, p[0] - cx) - az),
                                      math.cos(math.atan2(p[1] - cy, p[0] - cx) - az))) < 0.012)
        cl = [r for i, r in enumerate(rr) if i == 0 or r - rr[i - 1] > 0.15]
        if cl:
            counts.append(len(cl))
    return sum(counts) / len(counts) if counts else 0.0


def analyse(path, quiet=False):
    with zipfile.ZipFile(path) as z:
        g    = z.read("Metadata/plate_1.gcode").decode("utf-8", "replace")
        meta = "".join(z.read(n).decode("utf-8", "replace") for n in z.namelist() if "slice_info" in n)
        pngs = [n for n in z.namelist() if n.lower().endswith(".png")]
    # CLI slicing emits no thumbnails and Bambu Connect REFUSES such files; the nozzle flow
    # declaration needs one entry per extruder or Connect refuses it too. Both have bitten this
    # suite, so both are assertions rather than habits. Only the printable plate is checked.
    if not quiet:
        assert len(pngs) >= 5, f"only {len(pngs)} thumbnails - Bambu Connect will refuse the file"
        n_ext = len(json.load(open(os.path.join(CFG, "machine.json"), encoding="utf-8"))["nozzle_diameter"])
        decl  = re.search(r"; nozzle_volume_type = (.*)", g)
        assert decl and len([e for e in decl.group(1).strip().split(",") if e]) == n_ext, \
            "nozzle_volume_type must declare one entry per extruder - Bambu Connect will refuse it"

    m_lh = re.search(r"; layer_height = ([\d.]+)", g)
    lh = float(m_lh.group(1)) if m_lh else 0.25
    w = re.search(r"weight[^\d]*([\d.]+)", meta)
    t = re.search(r"prediction[^\d]*(\d+)", meta)
    out = {"ladder": {}, "core": None,
           "weight":  float(w.group(1)) if w else 0.0,
           "minutes": int(t.group(1)) // 60 if t else 0}
    for name, segs in split_objects(g):
        cx = sum(q[0] for q in segs) / len(segs)
        cy = sum(q[1] for q in segs) / len(segs)
        if "ladder" in name:
            for si, n in enumerate(STEPS):
                zb = si * STEP_H + STEP_H * 0.5
                by = [(q[0], q[1]) for q in segs if abs(q[4] - zb) < lh]
                if by:
                    out["ladder"][n] = crossings_per_ray(by, cx, cy)
        else:
            zb = CORE_H * 0.5
            by = [q for q in segs if abs(q[4] - zb) < lh]
            out["core"] = min_ring_gap(by, cx, cy)
    return out


def report(sparse, dense, density):
    print("")
    print(f"PlateD_data (density {density}%): {sparse['weight']:.2f} g, {sparse['minutes']} min"
          f"   [density 0 reference: {dense['weight']:.2f} g, {dense['minutes']} min]")
    print("")
    print("D1_ladder - beads per ray at each wall step")
    print(f"    {'wall':>5s} {'dense':>7s} {'sparse':>7s}   what it shows")
    sl = int(json.load(open(os.path.join(CFG, "cs_shapes.json"), encoding="utf-8"))["wall_loops"])
    for n in STEPS:
        d = dense["ladder"].get(n, 0.0)
        s = sparse["ladder"].get(n, 0.0)
        # A step can lose beads two different ways and the counts alone cannot tell them apart.
        # The stride only fires where the annular gate `n > 2*wall_loops + 1` is satisfied; below
        # it, any loss is adaptive beading being off (C-lite stretch). Decide by the gate, not by
        # the delta - reading the delta alone is what made the first pass of this script call the
        # 5-bead step "sparsified" when the sparsifier had never touched it.
        gate = n > 2 * sl + 1
        if gate and s < d - 0.5:
            what = "sparsified by the stride"
        elif not gate and s < d - 0.5:
            what = f"gate OFF - lost {d - s:.0f} bead to C-lite stretch, not the stride"
        elif not gate:
            what = "solid, full bead count"
        else:
            what = "gate on but not thinned"
        print(f"    {n:5d} {d:7.1f} {s:7.1f}   {what}")
    if sparse["core"] is not None:
        print("")
        print("D2_weave - thick annulus: smallest gap between adjacent rings")
        print(f"    density 0:    {dense['core']:6.3f} mm")
        print(f"    density {density}%:  {sparse['core']:6.3f} mm")
        # Reported WITHOUT a verdict on purpose. The weave deliberately brings anti-phase
        # neighbours' bead edges into contact at the lobe peaks, so a small gap between two
        # INFILL rings is the design working, not a fault; and from G-code alone there is no way
        # to tell which crossing is a shell and which is infill. What a bad clamp would look
        # like is bead-on-bead against a SHELL, which shows physically as over-extrusion - that
        # is what this object is for, and it is a question for the printer, not this script.


def slice_plate(density, outname, inject):
    proc = json.load(open(os.path.join(CFG, "cs_shapes.json"), encoding="utf-8"))
    proc["continuous_shell_infill_density"] = str(density)
    # The interior spiral and sparse infill both reorganise the buried interior, so enabling
    # both is now rejected outright rather than silently doing nothing (it used to be a no-op:
    # a spiral chain carries only its two pinned surfaces, and the sparsifier needs more rings).
    proc["continuous_shell_interior_spiral"] = "0"
    proc["name"] = f"CS data {density}"
    p = os.path.join(CFG, "cs_data.json")
    json.dump(proc, open(p, "w", encoding="utf-8"), indent=1)
    outdir = SLICED if inject else HERE
    out = os.path.join(outdir, outname + ".gcode.3mf")
    if os.path.exists(out):
        os.remove(out)
    r = subprocess.run([EXE, "--load-settings", f"{os.path.join(CFG, 'machine.json')};{p}",
                        "--load-filaments", os.path.join(CFG, "filament.json"),
                        "--arrange", "1", "--slice", "0",
                        "--export-3mf", outname + ".gcode.3mf", "--outputdir", outdir,
                        os.path.join(SUITE, "D1_ladder.stl"), os.path.join(SUITE, "D2_weave.stl")],
                       capture_output=True, text=True, timeout=1800)
    if r.returncode or not os.path.exists(out):
        print((r.stdout or "")[-900:])
        print((r.stderr or "")[-400:])
        sys.exit("slice failed")
    if inject:
        subprocess.run([sys.executable, os.path.join(HERE, "inject_thumbnails.py"), outname],
                       cwd=HERE, check=True, capture_output=True, text=True)
    return out


def main():
    argv = sys.argv[1:]
    density = argv[argv.index("--density") + 1] if "--density" in argv else "25"

    write_stl(os.path.join(SUITE, "D1_ladder.stl"), build_ladder())
    write_stl(os.path.join(SUITE, "D2_weave.stl"), build_core())
    print(f"D1_ladder: r_out {R_OUT}, wall steps {STEPS} beads "
          f"({[round(n * SP, 3) for n in STEPS]} mm), {STEP_H} mm each; gate flips between 6 and 5")
    print(f"D2_weave:  annulus r_out {CORE_R} / bore {CORE_BORE}, h {CORE_H} - the woven lattice")

    dense  = analyse(slice_plate(0, "PlateD_dense", False), quiet=True)
    sparse = analyse(slice_plate(density, OUTNAME, True))
    report(sparse, dense, density)


if __name__ == "__main__":
    main()
