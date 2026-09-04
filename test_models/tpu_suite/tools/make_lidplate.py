"""Build PlateM_lidset: F1, F7, F8 at weld 25, each lidded AND delidded (6 objects).

    make_lidplate.py [--publish]

The morning-of-print validation plate for the 2026-09-02 lone-ring work. Three suite shapes,
each in two copies that differ by ONE variable - the lid:

    *_delid   the campaign's coupon convention: top_shell_layers = 0, so the woven lattice is
              open from above and several layers read through each other. This is exactly what
              the bead-width inspector shows.
    *_lid     the same object with top_shell_layers = 4 (a ~1 mm solid cap over the lattice).
              The body stays continuous-shell (verified: body/bottom path length is identical
              to the delidded copy layer for layer; only the top four layers carry the cap's
              solid fill, ~1.6x the lattice length).

Bottom stays at the base's 0 for both, matching every coupon this campaign has printed - the
part prints lattice-to-bed and the only difference across each pair is the top. The three
shapes are the highest-signal set for this fix: F1 to see the kept tip rings in plastic, F7
for the multibore closing tour, F8 for the cone taper.

Delidded copies take the base config untouched; only the three _lid objects are patched, by
name, in the probe's model_settings.config, then the plate is re-sliced from that project.
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

# (source stl, short tag, can_lid). All three lid now: the multibore cap that once failed at its
# first cap layer (march_bands gave up on a band with several bores) is built by the marcher's
# Voronoi-cell partition as of commit 083551e34b, with floor-width connectors and the ordered
# geo tour keeping it travel-free and out of the bores. 6 objects.
SHAPES   = [("F1_crescent", "F1", True), ("F7_multibore", "F7", True), ("F8_cone_trunc", "F8", True)]
# A lidded object is FULLY ENCASED: a top AND a bottom solid cap over the lattice. A delidded
# object keeps the base's 0/0 - open both ends, so the foam reads through it. The multibore's
# bottom cap uses the same marched Voronoi-cell partition as the top (validated: F7 encased
# builds clean, both caps solid, zero faux-travel).
LID_TOP_LAYERS    = 4
LID_BOTTOM_LAYERS = 4
# Every object gets one extra wall loop (thicker shell): base wall_loops 2 -> 3.
WALL_LOOPS = 3
OUTNAME  = "PlateM_lidset"


def main():
    # Distinct-named copies so each object carries its own name in model_settings.config,
    # and the _lid / _delid override can target one without touching the other. Same geometry.
    stls = []
    n_lids = 0
    for src, tag, can_lid in SHAPES:
        srcpath = os.path.join(SUITE, src + ".stl")
        if not os.path.exists(srcpath):
            sys.exit(f"missing suite shape {srcpath}")
        kinds = ("lid", "delid") if can_lid else ("delid",)
        for kind in kinds:
            dst = os.path.join(SUITE, f"{tag}_{kind}.stl")
            shutil.copyfile(srcpath, dst)
            stls.append(dst)
            if kind == "lid":
                n_lids += 1
    print(f"{len(stls)} objects: {[os.path.splitext(os.path.basename(s))[0] for s in stls]}")
    print(f"   _lid -> top_shell_layers={LID_TOP_LAYERS}; _delid -> base (open top)")

    # Base process config = weld-25 suite preset (already top_shell_layers=0 = delidded).
    proc = json.load(open(os.path.join(CFG, "cs_all25.json"), encoding="utf-8"))
    proc["name"] = "CS lidset w25"
    proc["wall_loops"] = str(WALL_LOOPS)         # +1 thicker shell on every object
    p = os.path.join(CFG, "cs_lidset.json")
    json.dump(proc, open(p, "w", encoding="utf-8"), indent=1)

    base = [EXE, "--load-settings", f"{os.path.join(CFG, 'machine.json')};{p}",
            "--load-filaments", os.path.join(CFG, "filament.json")]

    probe = os.path.join(SLICED, OUTNAME + "_probe.gcode.3mf")
    if os.path.exists(probe):
        os.remove(probe)
    r = subprocess.run(base + ["--arrange", "1", "--slice", "0", "--export-3mf",
                               OUTNAME + "_probe.gcode.3mf", "--outputdir", SLICED] + stls,
                       capture_output=True, text=True, timeout=3600)
    if r.returncode or not os.path.exists(probe):
        print((r.stdout or "")[-1200:]); print((r.stderr or "")[-400:]); sys.exit("probe slice failed")

    zin = zipfile.ZipFile(probe)
    ms  = zin.read("Metadata/model_settings.config").decode("utf-8")
    hits = [0]
    def lidtag(m):
        # Names carry the .stl extension (F1_lid.stl). endswith("_lid.stl") matches the lids
        # and NOT the delids, whose suffix is "elid.stl".
        if not m.group(2).endswith("_lid.stl"):
            return m.group(0)
        hits[0] += 1
        return m.group(0) + (
            '\n    <metadata key="top_shell_layers" value="%d"/>'
            '\n    <metadata key="top_shell_thickness" value="0"/>'
            '\n    <metadata key="bottom_shell_layers" value="%d"/>'
            '\n    <metadata key="bottom_shell_thickness" value="0"/>'
            % (LID_TOP_LAYERS, LID_BOTTOM_LAYERS))
    ms2, _ = re.subn(r'(<object id="\d+">\s*<metadata key="name" value="([^"]+)"/>)', lidtag, ms)
    if hits[0] != n_lids:
        sys.exit(f"expected to lid exactly {n_lids} objects, tagged {hits[0]}")

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
        print((r.stdout or "")[-1200:]); print((r.stderr or "")[-400:]); sys.exit("slice failed")

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

    tm = re.search(r"; model printing time: (.+?)(?:\n|$)", g) or re.search(r"; total estimated time: (.+?)(?:\n|$)", g)
    fil = re.search(r"; filament used \[g\] = (.+)", g)
    print(f"\nOK: {out}")
    print(f"   objects: {len(stls)} | thumbnails: {len(pngs)}")
    print(f"   print time: {tm.group(1) if tm else '?'} | filament: {fil.group(1) if fil else '?'} g")
    print(f"   {'PUBLISHED to sliced/' if PUBLISHING else 'verify only - pass --publish to write sliced/'}")


if __name__ == "__main__":
    main()
