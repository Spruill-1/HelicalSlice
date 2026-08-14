"""Build a single-plate T5 trio project: three tubes front/middle/back with per-object
continuous_shell_scarf_length overrides (0/6/12), then slice it."""
import zipfile, re, os, subprocess, sys

SUITE = r"C:\Users\david\source\HelicalSlice\test_models\tpu_suite"
SLICED = os.path.join(SUITE, "sliced")
PROBE = os.path.join(SLICED, "T5_trio_probe.gcode.3mf")
TRIO  = os.path.join(SUITE, "T5_trio.3mf")
EXE   = r"C:\Users\david\source\HelicalSlice\build-arm64\src\Release\orca-slicer.exe"

zin = zipfile.ZipFile(PROBE)
model = zin.read("3D/3dmodel.model").decode("utf-8")
ms    = zin.read("Metadata/model_settings.config").decode("utf-8")

# The probe's arranged layout is already the front/middle/back column we want
# (x=175, y=32/160/288 - the sequential arranger's own clearance-safe pitch for
# 60mm-tall objects). Leave the model geometry COMPLETELY untouched (items carry
# required p:UUID attrs); only tag per-object scarf overrides in model_settings.
# Build-item order maps y=32 (front), 160 (middle), 288 (back) to ids 2, 4, 6.
items = re.findall(r'<item objectid="(\d+)"', model)
print("build items (object ids):", items)
model2 = model

# --- per-object scarf overrides + names, front->back = 0/6/12.
SCARF = {items[0]: ("T5_scarf0_front", "0"), items[1]: ("T5_scarf6_middle", "6"), items[2]: ("T5_scarf12_back", "12")}
def fix_object(m):
    oid = m.group(1)
    name, scarf = SCARF[oid]
    return (m.group(0)
            .replace('<metadata key="name" value="T5_seam.stl"/>',
                     f'<metadata key="name" value="{name}"/>\n'
                     f'    <metadata key="continuous_shell_scarf_length" value="{scarf}"/>', 1))
ms2, n = re.subn(r'<object id="(\d+)">\s*<metadata key="name" value="T5_seam\.stl"/>',
                 lambda m: fix_object(m), ms)
assert n == 3, f"expected 3 objects, patched {n}"

with zipfile.ZipFile(TRIO, "w", zipfile.ZIP_DEFLATED) as zout:
    for info in zin.infolist():
        n = info.filename
        if n.startswith("Metadata/plate_") or "thumbnail" in n.lower() or n.endswith(".gcode") or n.endswith(".gcode.md5"):
            continue                       # strip sliced artifacts; keep pure project
        if n == "3D/3dmodel.model":
            zout.writestr(n, model2)
        elif n == "Metadata/model_settings.config":
            zout.writestr(n, ms2)
        else:
            zout.writestr(n, zin.read(n))
print("wrote", TRIO)

cmd = [EXE,
       "--load-settings", os.path.join(SUITE, "config", "machine.json") + ";" + os.path.join(SUITE, "config", "cs_t5_scarf0.json"),
       "--load-filaments", os.path.join(SUITE, "config", "filament.json"),
       "--slice", "0",
       "--export-3mf", "T5_trio_fmb.gcode.3mf",
       "--outputdir", SLICED,
       TRIO]
r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
out = os.path.join(SLICED, "T5_trio_fmb.gcode.3mf")
print("rc:", r.returncode, "exists:", os.path.exists(out), "size:", os.path.getsize(out)//1024 if os.path.exists(out) else 0, "KB")
if r.returncode != 0:
    print("stdout:", (r.stdout or "")[-500:])
    print("stderr:", (r.stderr or "")[-500:])
    sys.exit(1)
