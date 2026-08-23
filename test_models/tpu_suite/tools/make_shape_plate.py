"""Build a shape + seam-control validation plate.

    make_shape_plate.py                                  # the full four-shape suite
    make_shape_plate.py --out PlateS_resume S4_capped.stl  # just these, e.g. to resume

Objects print by object: the named shapes with the interior spiral ON, plus the legacy
60 mm T5 tube with it OFF (per-object override) as the RING-MODE SEAM CONTROL. Running
the control on the same plate puts the seam A/B on one filament, one session and one
machine state - better than comparing tubes printed days apart - and its 60 mm height
matches the spiral tube already in hand.

Three things are asserted rather than assumed, because each fails silently:

  * The 60 mm control prints LAST. It exceeds the H2D's 47.4 mm rod clearance, so it is
    only safe alongside 30 mm shapes if nothing follows it. Getting this wrong drives
    the gantry through a printed part.
  * The per-object override applies - the control must slice in ring mode (~k junctions
    per layer) rather than inheriting the plate's spiral setting (~3).
  * The nozzle flow declaration carries one entry per extruder, or Bambu Connect
    refuses the file.
"""
import json, os, re, subprocess, sys, zipfile

SUITE   = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED  = os.path.join(SUITE, "sliced")
CFG     = os.path.join(SUITE, "config")
REPO    = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE     = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
CONTROL = "T5_seam.stl"                       # the legacy 60 mm tube, sliced in ring mode
ALL_SHAPES = ["S1_ellipse.stl", "S2_rsquare.stl", "S3_lobed.stl", "S4_capped.stl"]
ROD_CLEARANCE = 47.4

# Usage: make_shape_plate.py [--out NAME] [shape.stl ...]
# With no shapes listed the full suite is built. Naming a subset re-plates just those
# (e.g. to resume after a failed print) - the ring-mode control is always appended.
argv   = sys.argv[1:]
OUTNAME = "PlateS_shapes"
if "--out" in argv:
    i = argv.index("--out"); OUTNAME = argv[i + 1]; del argv[i:i + 2]
WITH_CONTROL = "--no-control" not in argv        # drop the control once it has been printed
if not WITH_CONTROL:
    argv.remove("--no-control")
SHAPES  = argv or list(ALL_SHAPES)
PROJECT = os.path.join(SUITE, OUTNAME + ".3mf")

def write_process():
    proc = json.load(open(os.path.join(CFG, "cs_t1_t4.json"), encoding="utf-8"))
    proc["continuous_shell_adaptive_beading"] = "1"
    proc["continuous_shell_interior_spiral"]  = "1"
    proc["name"] = "CS shapes"
    path = os.path.join(CFG, "cs_shapes.json")
    json.dump(proc, open(path, "w", encoding="utf-8"), indent=1)
    return path

def segments(gcode):
    """Split a sequential plate into per-object segments on big Z resets."""
    segs = []
    cur  = {"zmax": 0.0, "pts": []}
    prevz = None; x = y = None
    for line in gcode.split("\n"):
        if not line.startswith("G1 "):
            continue
        nz = re.search(r" Z(-?[\d.]+)", line); nx = re.search(r" X(-?[\d.]+)", line)
        ny = re.search(r" Y(-?[\d.]+)", line); ne = re.search(r" E(-?[\d.]+)", line)
        if nx: x = float(nx.group(1))
        if ny: y = float(ny.group(1))
        if nz:
            z = float(nz.group(1))
            if prevz is not None and z < prevz - 5.0:
                segs.append(cur); cur = {"zmax": 0.0, "pts": []}
            prevz = z; cur["zmax"] = max(cur["zmax"], z)
        if ne and float(ne.group(1)) > 0 and x is not None and y is not None:
            cur["pts"].append((x, y, prevz))
    segs.append(cur)
    return [s for s in segs if s["pts"]]

def junctions_per_layer(seg, sp=0.3963):
    """Radial jumps per layer - ring mode runs ~k, a ramped spiral ~3."""
    import math
    by = {}
    for (x, y, z) in seg["pts"]:
        by.setdefault(round(z / 0.25) * 0.25, []).append((x, y))
    mid = [b for b in sorted(by) if 0.3 * seg["zmax"] < b < 0.7 * seg["zmax"]]
    if not mid:
        return 0.0
    cx = sum(p[0] for b in mid for p in by[b]) / sum(len(by[b]) for b in mid)
    cy = sum(p[1] for b in mid for p in by[b]) / sum(len(by[b]) for b in mid)
    tot = 0
    for b in mid:
        pr = None
        for (x, y) in by[b]:
            r = math.hypot(x - cx, y - cy)
            if pr is not None and abs(r - pr) > 0.6 * sp:
                tot += 1
            pr = r
    return tot / len(mid)

def main():
    proc    = write_process()
    machine = os.path.join(CFG, "machine.json")
    fil     = os.path.join(CFG, "filament.json")
    stls    = [os.path.join(SUITE, f) for f in SHAPES + ([CONTROL] if WITH_CONTROL else [])]

    probe = os.path.join(SLICED, OUTNAME + "_probe.gcode.3mf")
    r = subprocess.run([EXE, "--load-settings", f"{machine};{proc}", "--load-filaments", fil,
                        "--arrange", "1", "--slice", "0",
                        "--export-3mf", OUTNAME + "_probe.gcode.3mf", "--outputdir", SLICED] + stls,
                       capture_output=True, text=True, timeout=1800)
    if r.returncode or not os.path.exists(probe):
        print((r.stdout or "")[-900:]); sys.exit("probe slice failed")

    zin   = zipfile.ZipFile(probe)
    ms    = zin.read("Metadata/model_settings.config").decode("utf-8")
    hits  = [0]
    def fix(m):
        if m.group(2) != CONTROL:
            return m.group(0)
        hits[0] += 1
        return (m.group(0).replace(f'value="{CONTROL}"', 'value="T5_ring_control"', 1) +
                '\n    <metadata key="continuous_shell_interior_spiral" value="0"/>')
    ms2, _ = re.subn(r'<object id="(\d+)">\s*<metadata key="name" value="([^"]+)"/>', fix, ms)
    if hits[0] != (1 if WITH_CONTROL else 0):
        sys.exit(f"expected {1 if WITH_CONTROL else 0} {CONTROL} object(s) to tag, tagged {hits[0]}")

    with zipfile.ZipFile(PROJECT, "w", zipfile.ZIP_DEFLATED) as zout:
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
    r = subprocess.run([EXE, "--load-settings", f"{machine};{proc}", "--load-filaments", fil,
                        "--slice", "0", "--export-3mf", OUTNAME + ".gcode.3mf",
                        "--outputdir", SLICED, PROJECT],
                       capture_output=True, text=True, timeout=1800)
    if r.returncode or not os.path.exists(out):
        print((r.stdout or "")[-900:]); sys.exit("plate slice failed")

    # CLI slicing emits no thumbnail PNGs and Bambu Connect REFUSES such files.
    subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "inject_thumbnails.py"), OUTNAME],
                   cwd=os.path.dirname(os.path.abspath(__file__)),
                   check=True, capture_output=True, text=True)

    with zipfile.ZipFile(out) as z:
        g    = z.read("Metadata/plate_1.gcode").decode("utf-8", "replace")
        pngs = [n for n in z.namelist() if n.lower().endswith(".png")]
        meta = "".join(z.read(n).decode("utf-8", "replace") for n in z.namelist() if "slice_info" in n)
    assert len(pngs) >= 5, f"only {len(pngs)} thumbnails - Bambu Connect will refuse the file"
    w = re.search(r'weight[^\d]*([\d.]+)', meta); t = re.search(r'prediction[^\d]*(\d+)', meta)
    segs = segments(g)
    print(f"plate: {w.group(1) if w else '?'} g, {int(t.group(1))//60 if t else '?'} min, "
          f"{len(segs)} objects by object")
    for i, s in enumerate(segs):
        print(f"  {i+1}. peak Z {s['zmax']:5.1f} mm   junctions/layer {junctions_per_layer(s):.1f}")

    tall = [i for i, s in enumerate(segs) if s["zmax"] > ROD_CLEARANCE]
    assert tall in ([], [len(segs) - 1]), \
        f"objects over {ROD_CLEARANCE} mm rod clearance must print LAST; got {tall} of {len(segs)}"
    if WITH_CONTROL:
        assert junctions_per_layer(segs[-1]) > 4.5, \
            "the control object is not slicing in ring mode - per-object override did not apply"

    # Bambu Connect validates the nozzle flow declared in the G-code CONFIG_BLOCK against
    # the printer and refuses the file on mismatch ("The nozzle flow setting does not
    # match ... Standardx1"). It needs ONE ENTRY PER EXTRUDER; a single entry on a
    # two-extruder machine leaves the right nozzle unresolvable. This has bitten the suite
    # twice, because the fix used to live only as a runtime patch inside a slicing script -
    # it is now pinned in config/machine.json, and asserted here so it cannot regress again.
    n_extruders = len(json.load(open(machine, encoding="utf-8"))["nozzle_diameter"])
    decl = re.search(r"; nozzle_volume_type = (.*)", g)
    assert decl, "no nozzle_volume_type declared - Bambu Connect will refuse the file"
    entries = [s for s in decl.group(1).strip().split(",") if s]
    assert len(entries) == n_extruders, \
        (f"nozzle_volume_type declares {len(entries)} entry/entries ({decl.group(1).strip()}) "
         f"for {n_extruders} extruders - Bambu Connect will refuse the file")
    print(f"nozzle flow declaration: {decl.group(1).strip()} ({n_extruders} extruders) OK")
    if WITH_CONTROL:
        print(f"\nOK: the {segs[-1]['zmax']:.0f} mm control prints last, and its junction count "
              f"confirms ring mode against the shapes' ramped spiral.")
    else:
        print("\nOK: shapes only, no ring-mode control on this plate.")

if __name__ == "__main__":
    main()
