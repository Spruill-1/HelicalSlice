"""Build the DE-LID plate: the inter-layer lock, open to view from above.

    make_delid_plate.py [--publish]

Everything the woven infill has been judged on so far shows the wrong thing. A section shows
in-plane structure; a top view shows one layer. Neither shows whether consecutive layers LOCK
to each other, which is what the Z-wave and the per-layer phase advance exist to do - and that
is the one claim in this subsystem with no physical evidence behind it at all.

Suppressing the top shells leaves the lattice open from above, so several layers are visible
through each other and the inter-layer diagonals can be seen (or not) directly.

    L1_open25   annulus, density 25 - stride 4, ~4 kept rings across the wall
    L2_open50   the SAME annulus at density 50 via a per-object override - stride 2, ~8 rings

Density is the comparison because density is what actually moved the needle: on PlateD the
denser object was preferred by eye, and a ring-count floor turns out NOT to be what made the
difference (it would have left the density-25 object completely unchanged - its natural ring
count already sits at the floor). So this plate asks the question that is still open, on the
axis that is known to matter.

Short at 12 mm: enough layers to see the lock, little enough to print quickly.
"""
import json, math, os, re, struct, subprocess, sys, zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_common import resolve_outdir, inject_env
from make_w1_plate import split_objects, write_stl

SUITE  = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED, PUBLISHING = resolve_outdir(SUITE)
CFG    = os.path.join(SUITE, "config")
REPO   = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE    = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE   = os.path.dirname(os.path.abspath(__file__))

R_OUT, R_BORE, H, SEG = 10.0, 2.5, 12.0, 256
DENSE_PCT = 50
OUTNAME = "PlateL_delid"


def build_annulus():
    ring = lambda r, z: [(r*math.cos(2*math.pi*i/SEG), r*math.sin(2*math.pi*i/SEG), z)
                         for i in range(SEG)]
    o0, o1 = ring(R_OUT, 0.0), ring(R_OUT, H)
    i0, i1 = ring(R_BORE, 0.0), ring(R_BORE, H)
    tris = []
    for i in range(SEG):
        j = (i + 1) % SEG
        tris += [(o0[i], o0[j], o1[j]), (o0[i], o1[j], o1[i])]
        tris += [(i0[i], i1[j], i0[j]), (i0[i], i1[i], i1[j])]
        tris += [(o0[i], i0[j], o0[j]), (o0[i], i0[i], i0[j])]
        # NO top annulus ring here: the mesh is still closed by the slicer's own top surface
        # handling, and top_shell_layers = 0 is what actually leaves the lattice open.
        tris += [(o1[i], o1[j], i1[j]), (o1[i], i1[j], i1[i])]
    return tris


def main():
    a = os.path.join(SUITE, "L1_open25.stl")
    b = os.path.join(SUITE, "L2_open50.stl")
    write_stl(a, build_annulus())
    write_stl(b, build_annulus())
    print(f"L1_open25 / L2_open50: annulus r_out {R_OUT} / bore {R_BORE}, {H} mm, top shells OFF")
    print(f"   L1 at the process density, L2 overridden to {DENSE_PCT}%")

    proc = json.load(open(os.path.join(CFG, "cs_shapes.json"), encoding="utf-8"))
    proc["continuous_shell_infill_density"]  = "25"
    proc["continuous_shell_interior_spiral"] = "0"
    # The whole point: no top skin, so the lattice is open from above and several layers are
    # visible through each other. Bottom shells stay on so the part still sticks to the plate.
    proc["top_shell_layers"]    = "0"
    proc["top_shell_thickness"] = "0"
    proc["name"] = "CS delid"
    p = os.path.join(CFG, "cs_delid.json")
    json.dump(proc, open(p, "w", encoding="utf-8"), indent=1)

    base = [EXE, "--load-settings", f"{os.path.join(CFG, 'machine.json')};{p}",
            "--load-filaments", os.path.join(CFG, "filament.json")]
    probe = os.path.join(SLICED, OUTNAME + "_probe.gcode.3mf")
    if os.path.exists(probe):
        os.remove(probe)
    r = subprocess.run(base + ["--arrange", "1", "--slice", "0", "--export-3mf",
                               OUTNAME + "_probe.gcode.3mf", "--outputdir", SLICED, a, b],
                       capture_output=True, text=True, timeout=1800)
    if r.returncode or not os.path.exists(probe):
        print((r.stdout or "")[-900:]); sys.exit("probe slice failed")

    zin = zipfile.ZipFile(probe)
    ms  = zin.read("Metadata/model_settings.config").decode("utf-8")
    hits = [0]
    def tag(m):
        if "L2_open50" not in m.group(2):
            return m.group(0)
        hits[0] += 1
        return m.group(0) + ('\n    <metadata key="continuous_shell_infill_density" '
                             'value="%d"/>' % DENSE_PCT)
    ms2, _ = re.subn(r'<object id="(\d+)">\s*<metadata key="name" value="([^"]+)"/>', tag, ms)
    if hits[0] != 1:
        sys.exit(f"expected to tag exactly one L2_open50 object, tagged {hits[0]}")
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

    lh = float(re.search(r"; layer_height = ([\d.]+)", g).group(1))
    fd = float(re.search(r"; filament_diameter = ([\d.]+)", g).group(1))
    w = re.search(r"weight[^\d]*([\d.]+)", meta)
    t = re.search(r"prediction[^\d]*(\d+)", meta)
    print(f"\n{OUTNAME}: {w.group(1) if w else '?'} g, {int(t.group(1))//60 if t else '?'} min")

    # The top layer must NOT be a solid skin, or the plate answers nothing.
    print(f"\n{'object':>14s} {'top-layer paths':>16s} {'rings at mid height':>20s}")
    for name, segs in split_objects(g, lh, math.pi*fd*fd/4.0):
        zmax = max(s[4] for s in segs)
        top  = [s for s in segs if s[4] > zmax - 1.5 * lh]
        mid  = [s for s in segs if abs(s[4] - H * 0.5) < lh]
        cx = sum(s[0] for s in mid) / len(mid)
        cy = sum(s[1] for s in mid) / len(mid)
        rings = 0
        for k in range(24):
            az = -math.pi + (k + 0.5) * (2 * math.pi / 24)
            hits = []
            for (x0, y0, x1, y1, *_) in mid:
                a0 = math.atan2(y0-cy, x0-cx); a1 = math.atan2(y1-cy, x1-cx)
                d0 = math.atan2(math.sin(a0-az), math.cos(a0-az))
                d1 = math.atan2(math.sin(a1-az), math.cos(a1-az))
                if (d0 < 0) != (d1 < 0) and abs(d0-d1) < math.pi:
                    t2 = d0/(d0-d1)
                    hits.append(math.hypot(x0+t2*(x1-x0)-cx, y0+t2*(y1-y0)-cy))
            hits.sort()
            rings += len([r for i, r in enumerate(hits) if i == 0 or r - hits[i-1] > 0.12])
        top_len = sum(s[6] for s in top)
        mid_len = sum(s[6] for s in mid)
        print(f"{name:>14s} {top_len/max(mid_len,1e-9):15.2f}x {rings/24.0:19.1f}")
    print("\nA top layer close to 1x the mid layer means the lid really is off; a solid skin")
    print("would run several times longer. Look DOWN into each: do consecutive layers lock,")
    print("or is it stacked independent rings? L2 has twice the rings of L1.")


if __name__ == "__main__":
    main()
