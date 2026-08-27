"""Build PlateR_w4w6: the geometry that exercises W4 and W6.

    make_wplate.py [--publish]

Neither change is print-validated, and neither is exercised by any existing plate - every one of
them has its lid off and no internal surface, which is exactly the geometry that avoids
force_solid.

    R1_tube       plain tube, LID ON. Its top layers are forced solid AND annular, which is the
                  case W4 changes: they used to fall back to C-lite for their bead counts while
                  the wall below was planned, so the count stepped at the boundary. On the S5 cone
                  that printed as a doubled lip at the bore. Look down the bore at the top face -
                  there should be no step where the cap meets the wall.

    R2_blindbore  a cylinder with a bore that starts part way up, so the bore's FLOOR is an
                  internal surface. That layer is what W6 is for: it needs solid over the bore's
                  footprint, and the wall around it is buried and should stay foam. Before W6 the
                  whole layer printed solid. Cut it at the floor - solid disc in the middle, foam
                  in the ring around it.

                  The bore opens UPWARD deliberately. A bore closed at the top would put a ceiling
                  over air, and asking TPU to bridge an 8 mm span would test the bridge rather
                  than the mask.

    R3_control    the validated annulus at density 25, unchanged by either change. Same filament,
                  same session - if this one has gone wrong, something unrelated has moved.
"""
import json, math, os, re, subprocess, sys, zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_common import resolve_outdir, inject_env
from make_w1_plate import write_stl
from make_serp_plate import tube, validate, circle

SUITE  = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED, PUBLISHING = resolve_outdir(SUITE)
CFG    = os.path.join(SUITE, "config")
REPO   = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE    = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE   = os.path.dirname(os.path.abspath(__file__))

SEG = 256
OUTNAME = "PlateR_w4w6"


def capped_tube(r_out, r_bore, h):
    """A tube with both ends open - its top layers are annular AND forced solid."""
    return tube(circle(r_out), circle(r_bore), h)


def blind_bore(r_out, r_bore, h, floor_z):
    """Cylinder with a bore from `floor_z` to the top: the bore's floor is an internal surface."""
    o, b = circle(r_out), circle(r_bore)
    tris = []
    n = len(o)
    for k in range(n):
        j = (k + 1) % n
        O0 = (o[k][0], o[k][1], 0.0); O0n = (o[j][0], o[j][1], 0.0)
        O1 = (o[k][0], o[k][1], h);   O1n = (o[j][0], o[j][1], h)
        tris += [(O0, O0n, O1n), (O0, O1n, O1)]                       # outer wall
        # solid bottom, fanned from the axis
        tris += [((0.0, 0.0, 0.0), O0n, O0)]
        # bore wall, floor_z -> h
        B0 = (b[k][0], b[k][1], floor_z); B0n = (b[j][0], b[j][1], floor_z)
        B1 = (b[k][0], b[k][1], h);       B1n = (b[j][0], b[j][1], h)
        tris += [(B0, B1n, B0n), (B0, B1, B1n)]                       # inward-facing
        # the bore floor, fanned from the axis
        tris += [((0.0, 0.0, floor_z), B0, B0n)]
        # top face: annulus between bore and outer
        tris += [(O1, O1n, B1n), (O1, B1n, B1)]
    return tris


def main():
    objs = []

    def add(name, tris):
        validate(tris, name)
        p = os.path.join(SUITE, name + ".stl")
        write_stl(p, tris)
        objs.append(p)

    add("R1_tube", capped_tube(10.0, 3.0, 12.0))
    add("R2_blindbore", blind_bore(10.0, 4.0, 14.0, 6.0))
    add("R3_control", capped_tube(10.0, 2.5, 12.0))

    proc = json.load(open(os.path.join(CFG, "cs_shapes.json"), encoding="utf-8"))
    proc["continuous_shell_infill_density"]      = "25"
    proc["continuous_shell_interior_spiral"]     = "0"
    proc["continuous_shell_interior_serpentine"] = "0"
    # LIDS ON, unlike every previous plate. Without a top shell there are no forced-solid layers
    # at all, and W4 and W6 are exactly about those.
    proc["name"] = "CS w4w6"
    cfgp = os.path.join(CFG, "cs_w4w6.json")
    json.dump(proc, open(cfgp, "w", encoding="utf-8"), indent=1)

    out = os.path.join(SLICED, OUTNAME + ".gcode.3mf")
    if os.path.exists(out):
        os.remove(out)
    r = subprocess.run([EXE, "--load-settings", f"{os.path.join(CFG, 'machine.json')};{cfgp}",
                        "--load-filaments", os.path.join(CFG, "filament.json"),
                        "--arrange", "1", "--slice", "0",
                        "--export-3mf", OUTNAME + ".gcode.3mf", "--outputdir", SLICED] + objs,
                       capture_output=True, text=True, timeout=3600)
    if r.returncode or not os.path.exists(out):
        print((r.stdout or "")[-1500:]); print((r.stderr or "")[-400:]); sys.exit("slice failed")
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

    w = re.search(r"weight[^\d]*([\d.]+)", meta)
    t = re.search(r"prediction[^\d]*(\d+)", meta)
    print(f"\n{OUTNAME}: {w.group(1) if w else '?'} g, "
          f"{int(t.group(1)) // 60 if t else '?'} min")

    # Per-layer path length: W6 should show the bore-floor layer costing much less than a solid
    # layer, and W4 should show no step at the cap boundary.
    from make_w1_plate import split_objects
    lh = float(re.search(r"; layer_height = ([\d.]+)", g).group(1))
    fd = float(re.search(r"; filament_diameter = ([\d.]+)", g).group(1))
    print(f"\n{'object':>14s}  path per layer by height (mm)")
    for name, segs in split_objects(g, lh, math.pi * fd * fd / 4.0):
        by = {}
        for q in segs:
            by.setdefault(int(round(q[4] / lh)), 0.0)
            by[int(round(q[4] / lh))] += q[6]
        ks = sorted(by)
        line = " ".join(f"{by[k]:.0f}" for k in ks[::4])
        print(f"{name.split('.')[0]:>14s}  {line}")
    print("\nR2's spike is the bore floor. Before W6 that layer was solid edge to edge; now it is")
    print("solid only over the bore, so the spike should be well short of a full solid layer.")


if __name__ == "__main__":
    main()
