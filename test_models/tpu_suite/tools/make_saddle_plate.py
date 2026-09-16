"""Build PlateMB_saddle: the saddle-fill print gate (6 objects).

    make_saddle_plate.py [--publish]

Validates continuous_shell_saddle_fill on the axis it changes, with an on-plate control for each
- the campaign's attribution rule: same STL, one variable changed. The fraction path leaves the
inter-bore SADDLE under-dense (the strip PlateMB_arcpeel peeled along); saddle_fill threads extra
level sets through any corridor still wider than the weld distance to densify that strip without
cutting the section.

Three de-lidded multibore coupons (the same sources as PlateMB_arcpeel), each in a saddle-OFF copy
and a saddle-ON copy that differ by ONE variable - continuous_shell_saddle_fill:

    *_soff   the base (saddle_fill = 0): the fraction level-set path, one ring count for the whole
             section - welds the bores but leaves the inter-bore strip as evenly-spread foam void.
    *_son    continuous_shell_saddle_fill = 1: extra level sets threaded through the saddle strip.

Both are de-lidded (cs_all25: top+bottom shell 0), so the woven lattice is open top and bottom -
hand-peel the strip BETWEEN THE BORES and see whether the _son coupon resists where the _soff coupon
separates (PEEL), and weigh son vs soff at matched geometry (the added density - the slicer's
per-object mm3 pre-screens it, and platecheck screens the weld-gap, below).

Only the three _son objects are patched, by name, in the probe's model_settings.config; the _soff
copies take the base untouched. Then the plate is re-sliced from that project.
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

# (source stl, short tag). Each becomes {tag}_soff.stl and {tag}_son.stl - identical geometry,
# distinct names so the saddle_fill override targets the ON copy without touching the control.
# Same sources as PlateMB_arcpeel: MB2e (2 equal bores), MB2u (2 unequal), F7 (3 unequal).
SHAPES  = [("MB2e_equal2", "MB2e"), ("MB2u_uneq2", "MB2u"), ("F7_multibore", "F7")]
OUTNAME = "PlateMB_saddle"


def main():
    stls = []
    n_on = 0
    for src, tag in SHAPES:
        srcpath = os.path.join(SUITE, src + ".stl")
        if not os.path.exists(srcpath):
            sys.exit(f"missing coupon {srcpath} - run gen_mbpeel_shapes.py first")
        for kind in ("soff", "son"):
            dst = os.path.join(SUITE, f"{tag}_{kind}.stl")
            shutil.copyfile(srcpath, dst)
            stls.append(dst)
            if kind == "son":
                n_on += 1
    print(f"{len(stls)} objects: {[os.path.splitext(os.path.basename(s))[0] for s in stls]}")
    print(f"   _son -> continuous_shell_saddle_fill=1; _soff -> base (fraction)")

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
    def ontag(m):
        # Names carry the .stl extension (MB2e_son.stl). endswith("_son.stl") matches the ON
        # copies and NOT the OFF ones (suffix "soff.stl").
        if not m.group(2).endswith("_son.stl"):
            return m.group(0)
        hits[0] += 1
        return m.group(0) + '\n    <metadata key="continuous_shell_saddle_fill" value="1"/>'
    ms2, _ = re.subn(r'(<object id="\d+">\s*<metadata key="name" value="([^"]+)"/>)', ontag, ms)
    if hits[0] != n_on:
        sys.exit(f"expected to tag exactly {n_on} saddle-on objects, tagged {hits[0]}")

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

    print(f"\nOK: {out}")
    print(f"   objects: {len(stls)} | thumbnails: {len(pngs)}")
    tm  = re.search(r"; model printing time: (.+?)(?:\n|$)", g) or re.search(r"; total estimated time: (.+?)(?:\n|$)", g)
    fil = re.search(r"; filament used \[g\] = (.+)", g)
    print(f"   print time: {tm.group(1) if tm else '?'} | filament: {fil.group(1) if fil else '?'} g")
    print(f"   {'PUBLISHED to sliced/' if PUBLISHING else 'verify only - pass --publish to write sliced/'}")


if __name__ == "__main__":
    main()
