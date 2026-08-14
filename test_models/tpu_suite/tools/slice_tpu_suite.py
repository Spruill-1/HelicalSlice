"""Flatten BBL profile inherit chains and slice the TPU tube suite via orca-slicer CLI.

Setup: Bambu H2D, right-side 0.4 hardened nozzle (TPU maps to master extruder = right),
Bambu TPU 95A HF (translucent), Engineering Plate, 0.45 line width / 0.25 layer height.
"""
import json, os, subprocess, sys

REPO   = r"C:\Users\david\source\HelicalSlice"
BASE   = os.path.join(REPO, "resources", "profiles", "BBL")
SUITE  = os.path.join(REPO, "test_models", "tpu_suite")
CFG    = os.path.join(SUITE, "config")
OUT    = os.path.join(SUITE, "sliced")
EXE    = os.path.join(REPO, "build-arm64", "src", "Release", "orca-slicer.exe")

os.makedirs(CFG, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

def flatten(kind, name):
    merged = {}
    stack = []
    n = name
    while n:
        d = json.load(open(os.path.join(BASE, kind, n + ".json"), encoding="utf-8"))
        stack.append(d)
        n = d.get("inherits")
    for d in reversed(stack):          # base first, child overrides
        merged.update(d)
    merged.pop("inherits", None)
    merged["name"] = name
    return merged

machine  = flatten("machine",  "Bambu Lab H2D 0.4 nozzle")
proc     = flatten("process",  "0.20mm Standard @BBL H2D")

# Generic third-party 95A TPU spool: rated nozzle 190-235C, bed 0-80C. The stock
# Generic TPU profile runs 240C - above the spool's max - so pin the temps in-spec.
filament = flatten("filament", "Generic TPU @BBL H2D")
filament.update({
    "nozzle_temperature":               ["230", "230"],
    "nozzle_temperature_initial_layer": ["230", "230"],
    "nozzle_temperature_range_low":     ["190"],
    "nozzle_temperature_range_high":    ["235"],
})

# Nozzle: right-side hardened steel. TPU maps to the master extruder (right) by itself,
# but pin the nozzle type explicitly if the key exists in the machine profile.
for key in ("nozzle_type",):
    if key in machine:
        v = machine[key]
        machine[key] = ["hardened_steel"] * len(v) if isinstance(v, list) else "hardened_steel"

common = {
    "layer_height": "0.25",
    "line_width": "0.45",
    "outer_wall_line_width": "0.45",
    "inner_wall_line_width": "0.45",
    "wall_loops": "2",
    "sparse_infill_density": "0%",
    "continuous_shell_mode": "1",
    "curr_bed_type": "Engineering Plate",
    # David's hardened 0.4 nozzle is on the RIGHT side of the H2D. filament_map is
    # 1-based (1 = left, 2 = right); the CLI does not run the auto grouping that the
    # GUI would use, so pin it manually.
    "filament_map_mode": "Manual",
    "filament_map": ["2"],
    # By-object sequencing: required by continuous shell whenever a plate has more
    # than one object; harmless on single-object plates.
    "print_sequence": "by object",
}

variants = {
    "cs_t1_t4":       {**common, "continuous_shell_interlock": "0",   "continuous_shell_scarf_length": "0"},
    "cs_t5_scarf0":   {**common, "continuous_shell_interlock": "100", "continuous_shell_scarf_length": "0"},
    "cs_t5_scarf6":   {**common, "continuous_shell_interlock": "100", "continuous_shell_scarf_length": "6"},
    "cs_t5_scarf12":  {**common, "continuous_shell_interlock": "100", "continuous_shell_scarf_length": "12"},
}

def write_json(path, data):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=1)

write_json(os.path.join(CFG, "machine.json"),  machine)
write_json(os.path.join(CFG, "filament.json"), filament)
for vname, overrides in variants.items():
    p = dict(proc)
    p.update(overrides)
    p["name"] = f"CS {vname}"
    write_json(os.path.join(CFG, vname + ".json"), p)

jobs = [
    # One combined plate for T1-T4 (same config; by-object sequencing; all under the
    # H2D's 47.4mm gantry clearance). The T5 scarf variants need different process
    # configs AND exceed the by-object height limit at 60mm -> separate prints.
    (["T1_tile_exact", "T2_deficit", "T3_thin", "T4_taper"], "cs_t1_t4", "PlateA_T1-T4"),
    (["T5_seam"],      "cs_t5_scarf0",  "T5_seam_scarf0"),
    (["T5_seam"],      "cs_t5_scarf6",  "T5_seam_scarf6"),
    (["T5_seam"],      "cs_t5_scarf12", "T5_seam_scarf12"),
]

failures = []
for stls, variant, outname in jobs:
    cmd = [
        EXE,
        "--load-settings", f"{os.path.join(CFG, 'machine.json')};{os.path.join(CFG, variant + '.json')}",
        "--load-filaments", os.path.join(CFG, "filament.json"),
        "--arrange", "1",
        "--slice", "0",
        "--export-3mf", outname + ".gcode.3mf",
        "--outputdir", OUT,
    ] + [os.path.join(SUITE, s + ".stl") for s in stls]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    out_path = os.path.join(OUT, outname + ".gcode.3mf")
    ok = r.returncode == 0 and os.path.exists(out_path)
    size = os.path.getsize(out_path) if os.path.exists(out_path) else 0
    print(f"{outname:22s} rc={r.returncode} exists={os.path.exists(out_path)} size={size//1024}KB")
    if not ok:
        failures.append((outname, r.returncode, (r.stdout or "")[-400:], (r.stderr or "")[-400:]))

for name, rc, so, se in failures:
    print(f"\n--- FAILED {name} rc={rc}\nstdout: {so}\nstderr: {se}")
sys.exit(1 if failures else 0)
