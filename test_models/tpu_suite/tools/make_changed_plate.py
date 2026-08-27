"""Build PlateQ_changed: only the objects the last two fixes actually moved.

    make_changed_plate.py [--publish]

PlateP's full battery would be eight objects and 94 minutes, and five of them would be printing
geometry that did not change. Measured against the plate that came off the printer, centred so
plate arrangement cannot confuse it, the median and 99th-percentile displacement per object is:

    P3_ecc       0.060 / 0.384 mm     <- the pinch-boundary shell fix
    P4_lobed     0.002 / 0.477 mm     <- the phase field
    P5_crescent  0.048 / 0.505 mm     <- the phase field
    P1_sparse    0.0007 / 0.032 mm    unchanged (coordinate rounding)
    P2_dense     0.0003 / 0.026 mm    unchanged
    P6_solid     0.001 / 0.085 mm     unchanged
    P7_taper     0.0002 / 0.036 mm    unchanged
    P8_thin      0.000 / 0.000 mm     bit-identical

Only three moved by anything a nozzle can resolve, and they are exactly the three the mechanisms
predict: the pinch fix touches only a section that produces a second chain, and the phase field
touches only sections whose rings differ in length. The peel-validated annuli reduce to their old
formula on a circle, so the parts already in hand ARE the control for this plate - it does not
have to carry its own.

PlateL proved the foam holds on ONE geometry - a concentric annulus, at two densities. Everything
that made it hold is geometry-dependent: the gap cap decides ring counts from the buried span,
the per-gap amplitude sizes every excursion from a measured gap, and the anchor pin depends on a
kept ring landing next to a shell. None of that is exercised by a round tube.

So this is the same fix on walls that vary, curve the wrong way, have no bore at all, or change
their ring count partway up:

    P1_sparse     annulus r10 / bore 2.5 at density 25   <- the validated baseline, repeated
    P2_dense      the same STL at density 50             <- ditto
    P3_ecc        bore offset 2 mm: wall thickness varies all the way round
    P4_lobed      5-lobe outer wall against a round bore
    P5_crescent   concave, and its centroid lies OUTSIDE the material
    P6_solid      disc r8, NO bore: no second face, so no anchor pair
    P7_taper      bore widens with height: the ring count CHANGES partway up, so the gap cap
                  re-decides mid-part and the layers either side must still tie to each other
    P8_thin       a 1 mm wall: nothing to sparsify, must decline and print sanely

Lids off, so the foam is visible from above without cutting. Peel P3 through P7 the way P1 and P2
were peeled - the bottom shell comes off first.

Meshes are validated, not assumed: every edge shared by exactly two triangles, no degenerate
triangles, and a fan cap only where its anchor provably sees the whole boundary. That check has
already caught one silently broken mesh in this suite.
"""
import json, math, os, re, subprocess, sys, zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_common import resolve_outdir, inject_env
from make_w1_plate import write_stl
from make_serp_plate import tube, solid, validate, star_ok, circle, lobed, ecc_pair, crescent

SUITE  = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED, PUBLISHING = resolve_outdir(SUITE)
CFG    = os.path.join(SUITE, "config")
REPO   = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE    = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE   = os.path.dirname(os.path.abspath(__file__))

H, SEG = 12.0, 256
OUTNAME = "PlateQ_changed"


def tube_z(outer_fn, inner_fn, h, levels=24):
    """Tube whose loops vary with height - so the wall thickness, and with it the ring count the
    sparsifier chooses, changes partway up the part."""
    zs = [h * k / levels for k in range(levels + 1)]
    rings = [(outer_fn(z), inner_fn(z)) for z in zs]
    n = len(rings[0][0])
    tris = []
    for L in range(levels):
        (o0, i0), (o1, i1) = rings[L], rings[L + 1]
        z0, z1 = zs[L], zs[L + 1]
        for k in range(n):
            j = (k + 1) % n
            O0 = lambda t: (t[0], t[1], z0)
            O1 = lambda t: (t[0], t[1], z1)
            tris += [(O0(o0[k]), O0(o0[j]), O1(o1[j])), (O0(o0[k]), O1(o1[j]), O1(o1[k]))]
            tris += [(O0(i0[k]), O1(i1[j]), O0(i0[j])), (O0(i0[k]), O1(i1[k]), O1(i1[j]))]
    o0, i0 = rings[0]
    o1, i1 = rings[-1]
    for k in range(n):
        j = (k + 1) % n
        b = lambda t: (t[0], t[1], 0.0)
        t = lambda q: (q[0], q[1], h)
        tris += [(b(o0[k]), b(i0[j]), b(o0[j])), (b(o0[k]), b(i0[k]), b(i0[j]))]
        tris += [(t(o1[k]), t(o1[j]), t(i1[j])), (t(o1[k]), t(i1[j]), t(i1[k]))]
    return tris


def main():
    objs = []

    def add(name, tris, dens=None):
        validate(tris, name)
        p = os.path.join(SUITE, name + ".stl")
        write_stl(p, tris)
        objs.append((p, dens))

    eo, ei = ecc_pair(10.0, 3.0, 2.0)
    add("P3_ecc", tube(eo, ei, H))
    add("P4_lobed", tube(lobed(10.0, 1.2, 5), circle(3.0), H))
    cre = crescent(12.0, 9.0, 8.0)
    if not star_ok(cre, (-6.0, 0.0)):
        sys.exit("P5: the fan anchor does not see the whole crescent boundary")
    add("P5_crescent", solid(cre, (-6.0, 0.0), H))

    proc = json.load(open(os.path.join(CFG, "cs_shapes.json"), encoding="utf-8"))
    proc["continuous_shell_infill_density"]      = "25"
    proc["continuous_shell_interior_spiral"]     = "0"
    proc["continuous_shell_interior_serpentine"] = "0"      # the ring lattice: the validated line
    proc["top_shell_layers"]    = "0"
    proc["top_shell_thickness"] = "0"
    proc["name"] = "CS battery"
    cfgp = os.path.join(CFG, "cs_battery.json")
    json.dump(proc, open(cfgp, "w", encoding="utf-8"), indent=1)

    base = [EXE, "--load-settings", f"{os.path.join(CFG, 'machine.json')};{cfgp}",
            "--load-filaments", os.path.join(CFG, "filament.json")]
    probe = os.path.join(SLICED, OUTNAME + "_probe.gcode.3mf")
    if os.path.exists(probe):
        os.remove(probe)
    r = subprocess.run(base + ["--arrange", "1", "--slice", "0", "--export-3mf",
                               OUTNAME + "_probe.gcode.3mf", "--outputdir", SLICED]
                       + [p for p, _ in objs], capture_output=True, text=True, timeout=3600)
    if r.returncode or not os.path.exists(probe):
        print((r.stdout or "")[-1500:]); sys.exit("probe slice failed")

    over = {os.path.basename(p)[:-4]: d for p, d in objs if d}
    zin = zipfile.ZipFile(probe)
    ms  = zin.read("Metadata/model_settings.config").decode("utf-8")
    hits = [0]
    def tag(m):
        d = over.get(m.group(2).replace(".stl", ""))
        if d is None:
            return m.group(0)
        hits[0] += 1
        return m.group(0) + ('\n    <metadata key="continuous_shell_infill_density" value="%d"/>' % d)
    ms2, _ = re.subn(r'<object id="(\d+)">\s*<metadata key="name" value="([^"]+)"/>', tag, ms)
    if hits[0] != len(over):
        sys.exit(f"expected to tag {len(over)} objects, tagged {hits[0]}")
    project = os.path.join(SUITE, OUTNAME + ".3mf")
    with zipfile.ZipFile(project, "w", zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():
            nm = info.filename
            if nm.startswith("Metadata/plate_") or "thumbnail" in nm.lower() \
               or nm.endswith(".gcode") or nm.endswith(".gcode.md5"):
                continue
            zout.writestr(nm, ms2 if nm == "Metadata/model_settings.config" else zin.read(nm))
    zin.close()
    os.remove(probe)

    out = os.path.join(SLICED, OUTNAME + ".gcode.3mf")
    if os.path.exists(out):
        os.remove(out)
    r = subprocess.run(base + ["--slice", "0", "--export-3mf", OUTNAME + ".gcode.3mf",
                               "--outputdir", SLICED, project],
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
    t2 = re.search(r"prediction[^\d]*(\d+)", meta)
    print(f"\n{OUTNAME}: {w.group(1) if w else '?'} g, "
          f"{int(t2.group(1)) // 60 if t2 else '?'} min")
    print("\nPeel every object the way PlateL was peeled: bottom shell off first.")
    print("P7 is the one to look at along its HEIGHT - the ring count changes partway up, so the")
    print("gap cap re-decides mid-part, and the layers either side of that change must still tie.")


if __name__ == "__main__":
    main()
