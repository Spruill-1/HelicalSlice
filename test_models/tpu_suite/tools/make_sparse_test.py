"""Build the sparse-core interaction test: S5_sparse.

Generalises the earlier "bored solid with a sparse core" to a part that stacks all three
regimes in one object, using GEOMETRY plus one global setting - no modifiers, no height
ranges (continuous shell rejects objects whose extra regions could split a layer).

    z 0 .. 10     annular wall, 5 beads      -> stays DENSE
    z 10 .. 24    45-degree cone closing     -> wall thickens, core sparsifies
    z 24 .. 36    solid column               -> SPARSE core inside solid shells
    top layers    lid                        -> forced solid over the sparse core

The wall thickness is the whole trick. The sparsifier only thins a chain that has buried
rings past its shells: annular needs n > 2*wall_loops + 1 (six beads at wall_loops 2),
solid needs n > wall_loops + 1. A 2.035 mm wall is exactly five beads - the thickest that
still stays fully dense - so the tube section keeps a solid wall while the core above the
cone, 27.9 mm across, is thick enough to merit sparse fill.

Note the mode interaction this exposes: `continuous_shell_infill_density` above zero
disables the interior spiral AND adaptive beading for the WHOLE object (the plan is gated
on !woven.enabled), so this part prints in woven ring mode throughout. Excluding the spiral
is defensible - it and sparse infill are alternatives for the same buried rings - but losing
adaptive beading is a straight regression, and the first print showed it: a 0.078 mm void up
the middle of the 2.035 mm wall (C-lite lays 4 beads stretched to 0.555 mm with zero overlap
where C-full lays 5 at 0.396), plus a doubled bead 0.203 mm from its neighbour on the cone's
overhanging bore lip. Both reproduce exactly at density 0 with adaptive beading off.

--bottom-shell 0 sparsifies the wedge of material that forms the cone. That wedge is solid
because bottom_shell_layers makes it force_solid, which takes the stride-1 branch (fully
solid, brick-staggered) and is never decimated. The lever is BINARY, not graded: 1, 2 and 3
all give the same 7.33 mm dense wedge at z=22; only 0 collapses it. detect_floating_vertical_
shell has no effect on it.
"""
import json, math, os, re, struct, subprocess, sys, zipfile

SUITE   = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED  = os.path.join(SUITE, "sliced")
CFG     = os.path.join(SUITE, "config")
REPO    = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE     = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")

R_OUT, WALL, SEG = 16.0, 2.035, 256          # 2.035 = five beads: the thickest dense wall
CYL_TOP = 10.0
APEX    = CYL_TOP + (R_OUT - WALL)           # 45 degrees, self-supporting
TALL    = 36.0
SPARSE_PCT = 25
NAME    = "S5_sparse"
SP      = 0.3963

def write_stl(path, tris):
    with open(path, "wb") as f:
        f.write(b"\0" * 80); f.write(struct.pack("<I", len(tris)))
        for a, b, c in tris:
            u = tuple(b[k] - a[k] for k in range(3)); v = tuple(c[k] - a[k] for k in range(3))
            n = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
            ln = math.hypot(math.hypot(n[0], n[1]), n[2]) or 1.0
            f.write(struct.pack("<3f", *(x/ln for x in n)))
            for p in (a, b, c): f.write(struct.pack("<3f", *p))
            f.write(struct.pack("<H", 0))

def build_mesh():
    r_in = R_OUT - WALL
    ring = lambda r, z: [(r*math.cos(2*math.pi*i/SEG), r*math.sin(2*math.pi*i/SEG), z) for i in range(SEG)]
    o0, o1 = ring(R_OUT, 0.0), ring(R_OUT, TALL)
    i0, i1 = ring(r_in, 0.0), ring(r_in, CYL_TOP)
    tris = []
    for i in range(SEG):
        j = (i + 1) % SEG
        tris += [(o0[i], o0[j], o1[j]), (o0[i], o1[j], o1[i])]
        tris += [(i0[i], i1[j], i0[j]), (i0[i], i1[i], i1[j])]
        tris += [(o0[i], i0[j], o0[j]), (o0[i], i0[i], i0[j])]
        tris += [(o1[i], o1[j], (0.0, 0.0, TALL))]
        tris += [(i1[i], (0.0, 0.0, APEX), i1[j])]
    return tris

def slice_it(density, outname, outdir, overrides=None):
    proc = json.load(open(os.path.join(CFG, "cs_shapes.json"), encoding="utf-8"))
    proc["continuous_shell_infill_density"] = str(density)
    proc.update(overrides or {})
    proc["name"] = f"CS sparse {density}"
    p = os.path.join(CFG, "cs_sparse.json")
    json.dump(proc, open(p, "w", encoding="utf-8"), indent=1)
    out = os.path.join(outdir, outname + ".gcode.3mf")
    if os.path.exists(out):
        os.remove(out)
    r = subprocess.run([EXE, "--load-settings", f"{os.path.join(CFG,'machine.json')};{p}",
                        "--load-filaments", os.path.join(CFG, "filament.json"),
                        "--arrange", "1", "--slice", "0",
                        "--export-3mf", outname + ".gcode.3mf", "--outputdir", outdir,
                        os.path.join(SUITE, NAME + ".stl")],
                       capture_output=True, text=True, timeout=1800)
    if r.returncode or not os.path.exists(out):
        print((r.stdout or "")[-900:]); print((r.stderr or "")[-400:]); sys.exit("slice failed")
    # Bambu Connect refuses a file whose nozzle flow declaration does not carry one entry per
    # extruder ("The nozzle flow setting does not match ... Standardx1"). Pinned in machine.json;
    # asserted here so a future edit cannot quietly reintroduce it.
    with zipfile.ZipFile(out) as z:
        g = z.read("Metadata/plate_1.gcode").decode("utf-8", "replace")
    n_ext = len(json.load(open(os.path.join(CFG, "machine.json"), encoding="utf-8"))["nozzle_diameter"])
    decl = re.search(r"; nozzle_volume_type = (.*)", g)
    msg = f"nozzle_volume_type must declare {n_ext} entries - Bambu Connect will refuse the file"
    assert decl and len([e for e in decl.group(1).strip().split(",") if e]) == n_ext, msg
    if os.path.abspath(outdir) == os.path.abspath(SLICED):
        inject_thumbnails(out, outname)
    return out


def inject_thumbnails(out, outname):
    """CLI slicing emits no thumbnail PNGs and Bambu Connect REFUSES such files.
    Every builder here must inject them; forgetting produces a file that slices, measures
    and reports perfectly and then will not load on the printer."""
    here = os.path.dirname(os.path.abspath(__file__))
    subprocess.run([sys.executable, os.path.join(here, "inject_thumbnails.py"), outname],
                   cwd=here, check=True, capture_output=True, text=True)
    with zipfile.ZipFile(out) as z:
        pngs = [n for n in z.namelist() if n.lower().endswith(".png")]
    assert len(pngs) >= 5, f"{outname}: only {len(pngs)} thumbnails - Bambu Connect will refuse the file"

def profile(path, label):
    with zipfile.ZipFile(path) as z:
        g = z.read("Metadata/plate_1.gcode").decode("utf-8", "replace")
        meta = "".join(z.read(n).decode("utf-8", "replace") for n in z.namelist() if "slice_info" in n)
    w = re.search(r'weight[^\d]*([\d.]+)', meta); t = re.search(r'prediction[^\d]*(\d+)', meta)
    x = y = zz = None; per = {}
    for line in g.split("\n"):
        if line.startswith("G1 "):
            nx = re.search(r" X(-?[\d.]+)", line); ny = re.search(r" Y(-?[\d.]+)", line)
            nz = re.search(r" Z(-?[\d.]+)", line); ne = re.search(r" E(-?[\d.]+)", line)
            if nx: x = float(nx.group(1))
            if ny: y = float(ny.group(1))
            if nz: zz = float(nz.group(1))
            if ne and float(ne.group(1)) > 0 and nx and ny and zz is not None:
                per.setdefault(round(zz / 0.25) * 0.25, []).append((x, y))
    allp = [q for v in per.values() for q in v]
    cx = sum(q[0] for q in allp) / len(allp); cy = sum(q[1] for q in allp) / len(allp)
    print(f"\n{label}: {w.group(1) if w else '?'} g, {int(t.group(1))//60 if t else '?'} min")
    print(f"{'z':>6s} {'beads/ray':>10s} {'junct/layer':>11s}   regime")
    for b in (2.5, 5.0, 9.0, 13.0, 18.0, 22.0, 26.0, 30.0, 34.0, 35.5):
        if b not in per: continue
        lp = per[b]
        counts = []
        for k in range(24):
            az = -math.pi + (k + 0.5) * (2*math.pi/24)
            rr = sorted(math.hypot(q[0]-cx, q[1]-cy) for q in lp
                        if abs(math.atan2(math.sin(math.atan2(q[1]-cy, q[0]-cx)-az),
                                          math.cos(math.atan2(q[1]-cy, q[0]-cx)-az))) < 0.012)
            cl = [r for i, r in enumerate(rr) if i == 0 or r - rr[i-1] > 0.15]
            if cl: counts.append(len(cl))
        pr = None; j = 0
        for (px, py) in lp:
            rr = math.hypot(px-cx, py-cy)
            if pr is not None and abs(rr-pr) > 0.6*SP: j += 1
            pr = rr
        reg = ("annular wall (want DENSE)" if b < CYL_TOP else
               ("cone" if b < APEX else "solid core (want SPARSE)"))
        print(f"{b:6.1f} {sum(counts)/len(counts):10.1f} {j:11d}   {reg}")

def main():
    # Usage: make_sparse_test.py [--bottom-shell N] [--out NAME]
    argv = sys.argv[1:]
    over, outname = {}, NAME
    if "--bottom-shell" in argv:
        i = argv.index("--bottom-shell"); over["bottom_shell_layers"] = argv[i + 1]; del argv[i:i + 2]
    if "--out" in argv:
        i = argv.index("--out"); outname = argv[i + 1]; del argv[i:i + 2]

    write_stl(os.path.join(SUITE, NAME + ".stl"), build_mesh())
    print(f"{NAME}: r_out {R_OUT}, wall {WALL} (5 beads), bore closes 45deg {CYL_TOP} -> {APEX:.1f}, "
          f"solid to {TALL}; core above the cone is {2*(R_OUT-WALL):.1f} mm across")
    dense  = slice_it(0, NAME + "_dense", os.path.dirname(os.path.abspath(__file__)))
    if over:
        print("overrides: " + ", ".join(f"{k}={v}" for k, v in over.items()))
    sparse = slice_it(SPARSE_PCT, outname, SLICED, over)
    profile(dense,  "density 0 (reference)")
    profile(sparse, f"density {SPARSE_PCT}% (the test)")
    print("\nWant: the annular wall unchanged between the two runs, the solid core thinned.")

if __name__ == "__main__":
    main()
