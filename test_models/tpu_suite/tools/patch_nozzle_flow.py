"""Fix the nozzle flow-type declaration in already-sliced .gcode.3mf files.

Bambu Connect refused the suite's files with "The nozzle flow setting does not
match". Root cause (established empirically, 2026-08-14): the CLI slice left
the machine config's nozzle_volume_type at its 1-element default, so the
G-code CONFIG_BLOCK declared `; nozzle_volume_type = Standard` — ONE entry for
a TWO-extruder printer. Connect could not resolve a flow type for the right
nozzle and failed the check even though the value itself matched the hardware
(David's H2D reports BOTH nozzles as 0.4 hardened steel, STANDARD flow — the
High Flow hotends are an optional upgrade he does not have). Notes:
  - Connect reads the flow type from the G-code CONFIG_BLOCK comments; the 3MF
    metadata (slice_info.config / project_settings.config) alone is not enough.
  - An intermediate attempt declaring the right nozzle High Flow was refused
    as a genuine value mismatch, confirming Connect does receive and compare
    the printer's nozzle state.

This normalizes every sliced/*.gcode.3mf — from the original 1-element state
or from the intermediate High Flow state — to exactly what a slice with
nozzle_volume_type = ["Standard", "Standard"] would have emitted:
  - Metadata/plate_N.gcode: `; nozzle_volume_type = Standard,Standard` and
    `; printer_extruder_variant = "Direct Drive Standard";"Direct Drive
    Standard"` (comment lines only — no motion/command bytes change), then the
    stored Metadata/plate_N.gcode.md5 is recomputed.
  - Metadata/slice_info.config: plate metadata vector -> "0 0";
    <filament>/<nozzle> tags keep/regain volume_type="Standard".
  - Metadata/project_settings.config: "nozzle_volume_type" ->
    ["Standard", "Standard"].
Idempotent; fails loudly on anything unexpected.
"""
import hashlib, os, re, shutil, tempfile, zipfile

SUITE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED = os.path.join(SUITE, "sliced")

SLICE_INFO = "Metadata/slice_info.config"
PROJ_SETTINGS = "Metadata/project_settings.config"

PLATE_FINAL = '<metadata key="nozzle_volume_type" value="0 0"/>'
PLATE_OLDS = ['<metadata key="nozzle_volume_type" value="0"/>',
              '<metadata key="nozzle_volume_type" value="0 1"/>']
TAG_HF = 'volume_type="High Flow"'
TAG_STD = 'volume_type="Standard"'
PROJ_FINAL = '"nozzle_volume_type": [\r\n\t\t"Standard",\r\n\t\t"Standard"\r\n\t]'
PROJ_OLDS = ['"nozzle_volume_type": [\r\n\t\t"Standard"\r\n\t]',
             '"nozzle_volume_type": [\r\n\t\t"Standard",\r\n\t\t"High Flow"\r\n\t]']
GCODE_NVT_FINAL = b"\n; nozzle_volume_type = Standard,Standard\n"
GCODE_NVT_OLDS = [b"\n; nozzle_volume_type = Standard\n",
                  b"\n; nozzle_volume_type = Standard,High Flow\n"]
GCODE_PEV_FINAL = b'\n; printer_extruder_variant = "Direct Drive Standard";"Direct Drive Standard"\n'
GCODE_PEV_OLDS = [b'\n; printer_extruder_variant = "Direct Drive Standard";"Direct Drive High Flow"\n']


def normalize(content, final, olds, label, name):
    if final in content:
        return content
    for old in olds:
        if content.count(old) == 1:
            return content.replace(old, final, 1)
    raise SystemExit(f"{name}: '{label}' not in any expected form")


def patch_3mf(path):
    name = os.path.basename(path)
    with zipfile.ZipFile(path) as zin:
        entries = zin.infolist()
        data = {e.filename: zin.read(e.filename) for e in entries}

    si = data[SLICE_INFO].decode("utf-8")
    si = normalize(si, PLATE_FINAL, PLATE_OLDS, "plate nozzle_volume_type", name)
    si = si.replace(TAG_HF, TAG_STD)
    if si.count(TAG_STD) < 2:  # one <filament> + one <nozzle> tag minimum
        raise SystemExit(f"{name}: expected >=2 volume_type=\"Standard\" tags")
    data[SLICE_INFO] = si.encode("utf-8")

    ps = data[PROJ_SETTINGS].decode("utf-8")
    data[PROJ_SETTINGS] = normalize(ps, PROJ_FINAL, PROJ_OLDS, "project nozzle_volume_type", name).encode("utf-8")

    gcode_names = [n for n in data if re.fullmatch(r"Metadata/plate_\d+\.gcode", n)]
    if not gcode_names:
        raise SystemExit(f"{name}: no plate gcode entry")
    for gn in gcode_names:
        g = data[gn]
        g = normalize(g, GCODE_NVT_FINAL, GCODE_NVT_OLDS, "gcode nozzle_volume_type", name)
        g = normalize(g, GCODE_PEV_FINAL, GCODE_PEV_OLDS, "gcode printer_extruder_variant", name)
        data[gn] = g
        md5_name = gn + ".md5"
        if md5_name not in data:
            raise SystemExit(f"{name}: missing {md5_name}")
        old_md5 = data[md5_name].decode("ascii")
        digest = hashlib.md5(g).hexdigest()
        if old_md5.strip().isupper():
            digest = digest.upper()
        data[md5_name] = (digest + old_md5[len(old_md5.rstrip()):]).encode("ascii")

    fd, tmp = tempfile.mkstemp(suffix=".3mf", dir=os.path.dirname(path))
    os.close(fd)
    with zipfile.ZipFile(tmp, "w") as zout:
        for e in entries:  # keep original entry order and per-entry compression
            zout.writestr(e, data[e.filename], compress_type=e.compress_type)
    shutil.move(tmp, path)

    with zipfile.ZipFile(path) as zchk:  # verify the rewrite took
        si = zchk.read(SLICE_INFO).decode("utf-8")
        ok = PLATE_FINAL in si and TAG_HF not in si
        ok = ok and PROJ_FINAL in zchk.read(PROJ_SETTINGS).decode("utf-8")
        for gn in gcode_names:
            g = zchk.read(gn)
            ok = ok and GCODE_NVT_FINAL in g and GCODE_PEV_FINAL in g and b"High Flow\n" not in g
            ok = ok and hashlib.md5(g).hexdigest().lower() == zchk.read(gn + ".md5").decode("ascii").strip().lower()
    if not ok:
        raise SystemExit(f"{name}: verification after rewrite failed")
    print(f"{name:28s} normalized to Standard,Standard (metadata + gcode + md5)")


files = sorted(f for f in os.listdir(SLICED) if f.endswith(".gcode.3mf"))
if not files:
    raise SystemExit(f"no .gcode.3mf files under {SLICED}")
for f in files:
    patch_3mf(os.path.join(SLICED, f))
print(f"done: {len(files)} files")
