"""Build PlateF7_budget: the F7 minimum-density sweep (6 objects, one print).

    make_f7budget_plate.py [--publish]

PlateMB_saddle showed saddle fill at the weld distance closes F7's peel but prints it SOLID (no
squish). The widest-first budget (continuous_shell_saddle_budget) spends a capped number of extra
ring levels on the worst saddles first, so a small budget is a LIGHT densification. saddlecheck
predicts cap4 lands F7 in the window (peel closed, still squishy); this plate brackets it on plastic.

Six F7_multibore copies, differing only by the saddle knobs (per-object model_settings override):

    F7_off    saddle_fill = 0                    - the peel control (fraction path, under-dense strip)
    F7_b3..b6 saddle_fill = 1, saddle_budget = K - K widest saddles closed (the density gradient)
    F7_full   saddle_fill = 1, saddle_budget = 0 - unlimited (the solid control, no squish)

De-lidded (cs_all25), so peel the inter-bore strip on each and find the LIGHTEST that stops the peel
while still squishing - that budget is the minimum-density target. saddlecheck screens it below first.
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

SRC     = "F7_multibore"
OUTNAME = "PlateF7_budget"
# (tag, saddle_fill, saddle_budget or None). None budget = key omitted (off control).
VARIANTS = [("off", 0, None), ("b3", 1, 3), ("b4", 1, 4), ("b5", 1, 5), ("b6", 1, 6), ("full", 1, 0)]


def main():
    srcpath = os.path.join(SUITE, SRC + ".stl")
    if not os.path.exists(srcpath):
        sys.exit(f"missing {srcpath} - run gen_mbpeel_shapes.py first")
    stls, over = [], {}
    for tag, fill, budget in VARIANTS:
        dst = os.path.join(SUITE, f"F7_{tag}.stl")
        shutil.copyfile(srcpath, dst)
        stls.append(dst)
        over[f"F7_{tag}.stl"] = (fill, budget)
    print(f"{len(stls)} objects: {[os.path.splitext(os.path.basename(s))[0] for s in stls]}")

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
    def tag(m):
        v = over.get(m.group(2))
        if v is None:
            return m.group(0)
        hits[0] += 1
        fill, budget = v
        add = '\n    <metadata key="continuous_shell_saddle_fill" value="%d"/>' % fill
        if budget is not None:
            add += '\n    <metadata key="continuous_shell_saddle_budget" value="%d"/>' % budget
        return m.group(0) + add
    ms2, _ = re.subn(r'(<object id="\d+">\s*<metadata key="name" value="([^"]+)"/>)', tag, ms)
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
        "nozzle_volume_type must declare one entry per extruder"

    print(f"\nOK: {out}")
    print(f"   objects: {len(stls)} | thumbnails: {len(pngs)}")
    tm  = re.search(r"; model printing time: (.+?)(?:\n|$)", g) or re.search(r"; total estimated time: (.+?)(?:\n|$)", g)
    fil = re.search(r"; filament used \[g\] = (.+)", g)
    print(f"   print time: {tm.group(1) if tm else '?'} | filament: {fil.group(1) if fil else '?'} g")
    print(f"   {'PUBLISHED to sliced/' if PUBLISHING else 'verify only - pass --publish to write sliced/'}")


if __name__ == "__main__":
    main()
