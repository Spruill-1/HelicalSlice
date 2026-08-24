"""Build PlateG_control: the strength control, plus the floor-guard reprint.

    make_control_plate.py

G1_conc - a plain CONCENTRIC tube, r_out 10 / bore 4, 30 mm. Geometry that NONE of the
    2026-08 wall-filling work touches: W1b only widens the plan's scope on walls of varying
    thickness, and W1 (plus the floor-affordability guard) only fires where a transitional
    bead is being floored, which a uniform wall of this thickness never does. Confirmed by the
    density-0 fence, where the concentric rows are byte-identical across every one of those
    commits.

    It exists because PlateF came off "extremely difficult to cut even with a scalpel" and
    that result is currently unattributed. If this control cuts like the older campaign parts,
    the gain belongs to the wall-filling work. If it is equally scalpel-resistant, the cause is
    something this campaign did not change - filament, spool, machine state - and the strength
    claim has to be withdrawn. Cut it the same way as F2 and compare honestly.

G2_lobed - the same 5-lobe STL as PlateF, re-sliced with the floor-affordability guard. That
    part printed with a column of voids near its outer edge, traced to beads laid below the
    0.30 mm transitional floor: the rebate paying for the floor was pushing recipients under
    it. Sub-floor extrusion is now 15.9% -> 6.75% of the layer, which is one transitional
    ring's worth, so every structural bead is at or above the floor. The column should be
    reduced, not gone - the residue is the transitional ring itself and belongs to the count
    law rather than the rebate.
"""
import json, math, os, re, struct, subprocess, sys, zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_common import resolve_outdir, inject_env

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_w1_plate import split_objects, ray_gaps, write_stl   # one implementation, one set of bugs

SUITE  = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED, PUBLISHING = resolve_outdir(SUITE)
CFG    = os.path.join(SUITE, "config")
REPO   = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE    = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE   = os.path.dirname(os.path.abspath(__file__))

R_OUT, R_BORE, H, SEG = 10.0, 4.0, 30.0, 256
OUTNAME = "PlateG_control"


def build_concentric():
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
        tris += [(o1[i], o1[j], i1[j]), (o1[i], i1[j], i1[i])]
    return tris


def main():
    conc = os.path.join(SUITE, "G1_conc.stl")
    write_stl(conc, build_concentric())
    lobed = os.path.join(SUITE, "F2_lobed.stl")
    if not os.path.exists(lobed):
        sys.exit("F2_lobed.stl missing - build it with make_w1_plate.py first")
    print(f"G1_conc:  concentric r_out {R_OUT} / bore {R_BORE}, {H} mm - untouched by W1b/W1/guard")
    print(f"G2_lobed: reusing {os.path.basename(lobed)} - the floor-guard reprint")

    proc = json.load(open(os.path.join(CFG, "cs_shapes.json"), encoding="utf-8"))
    proc["continuous_shell_infill_density"] = "0"
    proc["continuous_shell_interior_spiral"] = "0"
    proc["name"] = "CS control"
    p = os.path.join(CFG, "cs_control.json")
    json.dump(proc, open(p, "w", encoding="utf-8"), indent=1)

    out = os.path.join(SLICED, OUTNAME + ".gcode.3mf")
    if os.path.exists(out):
        os.remove(out)
    r = subprocess.run([EXE, "--load-settings", f"{os.path.join(CFG, 'machine.json')};{p}",
                        "--load-filaments", os.path.join(CFG, "filament.json"),
                        "--arrange", "1", "--slice", "0",
                        "--export-3mf", OUTNAME + ".gcode.3mf", "--outputdir", SLICED, conc, lobed],
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

    lh   = float(re.search(r"; layer_height = ([\d.]+)", g).group(1))
    fd   = float(re.search(r"; filament_diameter = ([\d.]+)", g).group(1))
    area = math.pi * fd * fd / 4.0
    w = re.search(r"weight[^\d]*([\d.]+)", meta)
    t = re.search(r"prediction[^\d]*(\d+)", meta)
    print(f"\n{OUTNAME}: {w.group(1) if w else '?'} g, {int(t.group(1))//60 if t else '?'} min")
    print(f"\n{'object':>14s} {'median':>7s} {'p90':>6s} {'rays>0.55':>10s} {'sub-floor':>10s} {'width lo/hi':>13s}")
    for name, segs in split_objects(g, lh, area):
        mid = [s for s in segs if abs(s[4] - H * 0.5) < lh]
        cx = sum(s[0] for s in mid) / len(mid)
        cy = sum(s[1] for s in mid) / len(mid)
        gv = ray_gaps(mid, cx, cy)
        band = [s for s in segs if 5.0 < s[4] < 25.0 and s[6] > 0.2]
        tot  = sum(s[6] for s in band)
        sub  = sum(s[6] for s in band if s[5] < 0.29)
        wl   = [s[5] for s in band]
        print(f"{name:>14s} {gv[len(gv)//2]:7.3f} {gv[int(len(gv)*0.9)]:6.3f} "
              f"{sum(1 for v in gv if v > 0.55):7d}/72 {100*sub/tot:9.2f}% {min(wl):6.3f}/{max(wl):.3f}")
    print("\nG1 is the attribution control: cut it exactly as F2 was cut. If it resists the")
    print("scalpel too, the strength gain is NOT from this campaign's wall-filling work.")
    print("G2 should show a REDUCED void column, not an absent one.")


if __name__ == "__main__":
    main()
