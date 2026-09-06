"""Build PlateMB_arcpeel: the arc-capable multibore validation plate (6 objects).

    make_mbpeel_plate.py [--publish]

Validates the per-bore-cell partition (commit b9f2959880) on the two axes it introduces, with
an ON-PLATE control for each - the campaign's attribution rule: same STL, one variable changed.

Three de-lidded multibore coupons (gen_mbpeel_shapes.py), each in an ARC copy and a FRACTION
copy that differ by ONE variable - continuous_shell_arc_rings:

    *_arc    continuous_shell_arc_rings = 1: distance-anchored rings, and a section with several
             bores is partitioned by nearest bore into single-bore cells (each threads and welds).
    *_frac   the base (arc_rings = 0): the fraction level-set path, one ring count for the whole
             multibore section - the current default, which welds the bores but leaves the
             inter-bore regions as evenly-spread foam void.

Both are de-lidded (cs_all25 is top_shell_layers = 0 AND bottom_shell_layers = 0), so the woven
lattice is open top and bottom - hand-peel it and see whether the arc coupon separates along a
cut plane that the fraction coupon has no cut to follow (PEEL), and weigh arc vs fraction at
matched geometry (DENSIFICATION - the slicer's per-object mm3 pre-screens it below).

Only the three _arc objects are patched, by name, in the probe's model_settings.config; the
fraction copies take the base untouched. Then the plate is re-sliced from that project.
"""
import json, os, re, shutil, subprocess, sys, zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_common import resolve_outdir, inject_env

SUITE  = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED, PUBLISHING = resolve_outdir(SUITE)
CFG    = os.path.join(SUITE, "config")
REPO   = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE    = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE   = os.path.dirname(os.path.abspath(__file__))

# (source stl, short tag). Each becomes {tag}_arc.stl and {tag}_frac.stl - identical geometry,
# distinct names so the arc override targets one of each pair without touching the control.
#
# The 3-bore case is F7_multibore (the object the fix was validated on: three UNEQUAL bores). A
# purpose-built symmetric-equal 3-bore (MB3_tri3, gen_mbpeel_shapes.py) is NOT used: its three
# equal bores put a degenerate medial vertex at the centroid (equidistant to all three), and the
# foam path blows up there - two copies never finish export. That is a pre-existing limit of the
# shipped fraction path too (fraction MB3 blows up identically), not the arc partition, so it is
# left out of the plate and F7's asymmetric bores stand in for the realistic multi-bore.
SHAPES  = [("MB2e_equal2", "MB2e"), ("MB2u_uneq2", "MB2u"), ("F7_multibore", "F7")]
OUTNAME = "PlateMB_arcpeel"


def main():
    stls = []
    n_arc = 0
    for src, tag in SHAPES:
        srcpath = os.path.join(SUITE, src + ".stl")
        if not os.path.exists(srcpath):
            sys.exit(f"missing coupon {srcpath} - run gen_mbpeel_shapes.py first")
        for kind in ("arc", "frac"):
            dst = os.path.join(SUITE, f"{tag}_{kind}.stl")
            shutil.copyfile(srcpath, dst)
            stls.append(dst)
            if kind == "arc":
                n_arc += 1
    print(f"{len(stls)} objects: {[os.path.splitext(os.path.basename(s))[0] for s in stls]}")
    print(f"   _arc -> continuous_shell_arc_rings=1; _frac -> base (fraction)")

    # Base = the de-lidded weld-25 coupon preset (top+bottom shell 0, interpolated_rings 1).
    proc = os.path.join(CFG, "cs_all25.json")
    base = [EXE, "--load-settings", f"{os.path.join(CFG, 'machine.json')};{proc}",
            "--load-filaments", os.path.join(CFG, "filament.json")]

    probe = os.path.join(SLICED, OUTNAME + "_probe.gcode.3mf")
    if os.path.exists(probe):
        os.remove(probe)
    r = subprocess.run(base + ["--arrange", "1", "--slice", "0", "--export-3mf",
                               OUTNAME + "_probe.gcode.3mf", "--outputdir", SLICED] + stls,
                       capture_output=True, text=True, timeout=3600)
    if r.returncode or not os.path.exists(probe):
        print((r.stdout or "")[-1600:]); print((r.stderr or "")[-400:]); sys.exit("probe slice failed")

    zin = zipfile.ZipFile(probe)
    ms  = zin.read("Metadata/model_settings.config").decode("utf-8")
    hits = [0]
    def arctag(m):
        # Names carry the .stl extension (MB2e_arc.stl). endswith("_arc.stl") matches the arc
        # copies and NOT the fraction ones (suffix "frac.stl").
        if not m.group(2).endswith("_arc.stl"):
            return m.group(0)
        hits[0] += 1
        return m.group(0) + '\n    <metadata key="continuous_shell_arc_rings" value="1"/>'
    ms2, _ = re.subn(r'(<object id="\d+">\s*<metadata key="name" value="([^"]+)"/>)', arctag, ms)
    if hits[0] != n_arc:
        sys.exit(f"expected to tag exactly {n_arc} arc objects, tagged {hits[0]}")

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
        print((r.stdout or "")[-1600:]); print((r.stderr or "")[-400:]); sys.exit("slice failed")

    subprocess.run([sys.executable, os.path.join(HERE, "inject_thumbnails.py"), OUTNAME],
                   cwd=HERE, check=True, capture_output=True, text=True, env=inject_env(SLICED))

    with zipfile.ZipFile(out) as z:
        g    = z.read("Metadata/plate_1.gcode").decode("utf-8", "replace")
        pngs = [n for n in z.namelist() if n.lower().endswith(".png")]
    assert len(pngs) >= 5, f"only {len(pngs)} thumbnails - Bambu Connect will refuse the file"
    n_ext = len(json.load(open(os.path.join(CFG, "machine.json"), encoding="utf-8"))["nozzle_diameter"])
    decl  = re.search(r"; nozzle_volume_type = (.*)", g)
    assert decl and len([e for e in decl.group(1).strip().split(",") if e]) == n_ext, \
        "nozzle_volume_type must declare one entry per extruder - Bambu Connect will refuse it"

    # Densification is pre-screened post-slice by reading each object's deposited length/mm3 from
    # the plate gcode (the analysis step); the whole-plate filament total is printed here.
    print(f"\nOK: {out}")
    print(f"   objects: {len(stls)} | thumbnails: {len(pngs)}")
    tm  = re.search(r"; model printing time: (.+?)(?:\n|$)", g) or re.search(r"; total estimated time: (.+?)(?:\n|$)", g)
    fil = re.search(r"; filament used \[g\] = (.+)", g)
    print(f"   print time: {tm.group(1) if tm else '?'} | filament: {fil.group(1) if fil else '?'} g")
    print(f"   {'PUBLISHED to sliced/' if PUBLISHING else 'verify only - pass --publish to write sliced/'}")


if __name__ == "__main__":
    main()
