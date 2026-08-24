"""Build the eccentric-bore validation plate (E1/E2).

    make_ecc_plate.py

Validates the marching-pinch fix in plastic. A wall of varying thickness used to stop
marching the moment its THINNEST azimuth closed, abandoning the thick side: a 4/8 mm wall
kept 69% of its material and carried a 0.9 mm void ring up the middle of the thick side.
Marching now hands the pinched remainder on as a solid region, and the same tube lays 101%.

The G-code question is already settled. What a print has to answer is different: the
recovered material arrives as a SEPARATE CHAIN, joined to the annular chain at the pinch
boundary, so there is a new interface at a place that used to be a void. Does it bond, or
has a void been traded for a weak plane?

    E1_ecc   bore offset 2.0 mm (4 mm thin wall / 8 mm thick) - the fix's main case; the
             pinch boundary runs up the thick side. Section it, and tear it across.
    E2_mild  bore offset 0.5 mm (5.5 / 6.5) - does NOT pinch, and still lays only 94%.
             That residual is a different, much smaller mechanism and is NOT fixed here;
             this object is what it looks like, for comparison against E1's filled wall.

Both are 30 mm tall so they can be cut and torn, and both are sliced at
continuous_shell_infill_density = 0 with the interior spiral off, so the wall structure is
the only variable between them.
"""
import json, math, os, re, struct, subprocess, sys, zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_common import resolve_outdir, inject_env

SUITE  = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED, PUBLISHING = resolve_outdir(SUITE)
CFG    = os.path.join(SUITE, "config")
REPO   = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE    = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE   = os.path.dirname(os.path.abspath(__file__))

R_OUT, R_BORE, H, SEG = 10.0, 4.0, 30.0, 256
OBJECTS = [("E1_ecc", 2.0), ("E2_mild", 0.5)]
OUTNAME = "PlateE_eccentric"


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


def build_tube(dx):
    ring = lambda r, z, off=0.0: [(off + r*math.cos(2*math.pi*i/SEG), r*math.sin(2*math.pi*i/SEG), z)
                                  for i in range(SEG)]
    o0, o1 = ring(R_OUT, 0.0), ring(R_OUT, H)
    i0, i1 = ring(R_BORE, 0.0, dx), ring(R_BORE, H, dx)
    tris = []
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


def worst_gap(segments, cx, cy):
    """Largest radial gap between adjacent bead crossings, over many rays. A filled wall keeps
    every gap near one bead spacing; the pre-fix tube carried a 0.907 mm void ring.

    Crossings are INTERPOLATED where a segment crosses the ray, not collected from endpoints
    near it: a printed polyline can put its vertices anywhere, so sampling by proximity silently
    misses whole beads and invents gaps that are not in the plastic. A first cut did exactly
    that and reported 2.1 mm gaps in a wall that is in fact filled."""
    worst, worst_az = 0.0, 0.0
    for k in range(72):
        az = -math.pi + (k + 0.5) * (2 * math.pi / 72)
        hits = []
        for (x0, y0, x1, y1, _z) in segments:
            a0 = math.atan2(y0 - cy, x0 - cx); a1 = math.atan2(y1 - cy, x1 - cx)
            d0 = math.atan2(math.sin(a0 - az), math.cos(a0 - az))
            d1 = math.atan2(math.sin(a1 - az), math.cos(a1 - az))
            if (d0 < 0) != (d1 < 0) and abs(d0 - d1) < math.pi:
                t = d0 / (d0 - d1)
                hits.append(math.hypot(x0 + t*(x1-x0) - cx, y0 + t*(y1-y0) - cy))
        hits.sort()
        cl = [r for i, r in enumerate(hits) if i == 0 or r - hits[i-1] > 0.15]
        for i in range(len(cl) - 1):
            if cl[i+1] - cl[i] > worst:
                worst, worst_az = cl[i+1] - cl[i], az
    return worst, worst_az


def main():
    stls = []
    for name, dx in OBJECTS:
        p = os.path.join(SUITE, name + ".stl")
        write_stl(p, build_tube(dx))
        stls.append(p)
        print(f"{name}: r_out {R_OUT}, bore {R_BORE} offset {dx} mm "
              f"(wall {R_OUT-R_BORE-dx:.1f} thin / {R_OUT-R_BORE+dx:.1f} thick), {H} mm tall")

    proc = json.load(open(os.path.join(CFG, "cs_shapes.json"), encoding="utf-8"))
    proc["continuous_shell_infill_density"] = "0"
    proc["continuous_shell_interior_spiral"] = "0"   # wall structure is the only variable
    proc["name"] = "CS eccentric"
    p = os.path.join(CFG, "cs_ecc.json")
    json.dump(proc, open(p, "w", encoding="utf-8"), indent=1)

    out = os.path.join(SLICED, OUTNAME + ".gcode.3mf")
    if os.path.exists(out):
        os.remove(out)
    r = subprocess.run([EXE, "--load-settings", f"{os.path.join(CFG, 'machine.json')};{p}",
                        "--load-filaments", os.path.join(CFG, "filament.json"),
                        "--arrange", "1", "--slice", "0",
                        "--export-3mf", OUTNAME + ".gcode.3mf", "--outputdir", SLICED] + stls,
                       capture_output=True, text=True, timeout=1800)
    if r.returncode or not os.path.exists(out):
        print((r.stdout or "")[-900:]); print((r.stderr or "")[-400:])
        sys.exit("slice failed")
    subprocess.run([sys.executable, os.path.join(HERE, "inject_thumbnails.py"), OUTNAME],
                   cwd=HERE, check=True, capture_output=True, text=True, env=inject_env(SLICED))

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
    lh = float(m_lh.group(1)) if m_lh else 0.25
    w = re.search(r"weight[^\d]*([\d.]+)", meta)
    t = re.search(r"prediction[^\d]*(\d+)", meta)
    print(f"\n{OUTNAME}: {w.group(1) if w else '?'} g, {int(t.group(1))//60 if t else '?'} min")
    segs = split_objects(g)
    print(f"\n{'object':>9s} {'wall thin/thick':>16s} {'crossings':>10s} {'worst gap':>10s}")
    for s, (name, dx) in zip(segs, OBJECTS):
        sg = s["seg"]
        cx = sum(q[0] for q in sg) / len(sg)
        cy = sum(q[1] for q in sg) / len(sg)
        mid = [q for q in sg if abs(q[4] - H * 0.5) < lh]
        gap, az = worst_gap(mid, cx, cy)
        n = 0
        for k in range(24):
            az2 = -math.pi + (k + 0.5) * (2 * math.pi / 24)
            hits = []
            for (x0, y0, x1, y1, _z) in mid:
                a0 = math.atan2(y0 - cy, x0 - cx); a1 = math.atan2(y1 - cy, x1 - cx)
                d0 = math.atan2(math.sin(a0 - az2), math.cos(a0 - az2))
                d1 = math.atan2(math.sin(a1 - az2), math.cos(a1 - az2))
                if (d0 < 0) != (d1 < 0) and abs(d0 - d1) < math.pi:
                    t = d0 / (d0 - d1)
                    hits.append(math.hypot(x0 + t*(x1-x0) - cx, y0 + t*(y1-y0) - cy))
            hits.sort()
            n += len([r for i, r in enumerate(hits) if i == 0 or r - hits[i-1] > 0.15])
        print(f"{name:>9s} {f'{R_OUT-R_BORE-dx:.1f} / {R_OUT-R_BORE+dx:.1f} mm':>16s} "
              f"{n/24.0:10.1f} {gap:9.3f} mm  at azimuth {math.degrees(az):6.1f} deg")
    print("\nA filled wall keeps every radial gap near one bead spacing (0.396 mm).")
    print("Pre-fix, the 2.0 mm-offset tube carried a 0.907 mm void ring up its thick side.")
    print("\nWhat to look for when it comes off: section E1 across the THICK side and check the")
    print("pinch boundary - the recovered material is a separate chain joined there, so that")
    print("interface is the thing the fix could plausibly have gotten wrong.")


if __name__ == "__main__":
    main()
