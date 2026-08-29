"""Build PlateQ_changed: only the objects a fix actually moves.

    make_changed_plate.py [--publish]

PlateP's full battery would be eight objects and 94 minutes, and most of them would be printing
geometry that did not change. The same three objects have carried two rounds of validation now.

ROUND 2 (2026-08-27) fixed the crescent and broke the eccentric. The crescent stopped peeling -
its concave side went from 5.7% of that interface welded to 38.1%, against 41.5% convex - and the
lobed object came out as the foam is meant to look. The eccentric came apart into three pieces.

ROUND 5 (2026-08-29) - weave the medial terminus. Round 4 passed the peel test on both objects,
which closed the eccentric's peel problem after four rounds. What was left on both was a solid in
the middle of the foam, and rendering the material present on ALL of twelve consecutive layers
named it: the MEDIAL CORE, the terminus plus the ring the distribution puts beside it, one slot
apart and both immovable, so identical on every layer. On the crescent that is 10.7 mm2 of a
155 mm2 section, a millimetre wide, running its whole length.

The terminus is the one held ring here with no face behind it - a shell backs onto the printed
surface, a pocket face onto another chain's material, a terminus is just where the marching ran
out - and the pin rule freezes its neighbour along with it. So weave it: kept whatever the density
asks, but tagged infill, which releases the neighbour too.

    object         material  interior solid  worst void  moves/layer
    P5_crescent      155.1      10.7 -> 4.6    1.33 (=)   28.0 -> 29.7%
    P3_ecc           205.9         49.7 (=)    0.59 (=)      20.6% (=)
    P6_solid         115.0          2.8 (=)    0.07 (=)      38.7% (=)

ONE object changes. P3_ecc is unchanged by it - its pocket chain is too short to qualify - so it
rides along as a same-session control and as the peel re-check.

The obvious alternative, dropping the ring beside the terminus, was built and measured and is
worse: it removes the ridge and puts a VOID in the same place an order of magnitude bigger
(crescent worst void 1.33 -> 14.72 mm2). A gap next to something immovable sits at a fixed radius
and stacks. The medial wants a ring there; it does not want a frozen one.

ROUND 4 (2026-08-28) - pin ONE face of the pocket, not both. Round 3 held together but had
voids, and the reason turned out to be Z, not the layer: the eccentric's lattice barely moves
between layers, so its voids sit at the same radius on every one and stack into channels. Of the
six battery objects it is the least mobile, and the only one whose voids stack:

    fraction of ring length the weave moves    P6 44.0  P1 39.1  P4 38.6  P5 30.5  P7 24.7
                                               P3_ecc 17.5   <- the one with voids
    void area present on ALL of ten layers     P3_ecc 47%    P5_crescent 23%

A pin buys a weld that holds all the way round and costs the ring all of its motion. The chain
filling a pocket welds to whichever face is held, so connectivity needs ONE pin, not two:

    eccentric      material  void/layer  through-void  moves/layer  worst channel
    both pinned      204.8      107.5      40.0 (37%)     17.3%       1.79 mm2
    outer only       205.9      106.4      34.2 (32%)     20.6%       0.59 mm2

Same material; the worst persistent channel falls to a third. Still one connected piece, the two
chains still meeting over 58.5% of that interface (68.7% with both pinned; 15.4% in the version
that printed acceptably), and the pocket chain's ring now welds to four rings rather than two.

Only chains with a pocket change, and the eccentric is the only object that has one.

ROUND 3 (2026-08-28) - the pocket a pinching wall leaves behind. Band marching walks a chain down
one face to the medial and back up the other; on a wall of varying thickness those two middle
rings do not meet, and the remainder is handed on as a separate chain that fills the space between
them. Those middle rings are that pocket's faces and the chain filling it has nothing else to weld
to - but the gap cap looks only within a chain, and so does the weave, so the stride stepped over
both. Measured on this object: middle pair 7.41 slots apart, the next chain's outermost ring one
slot inside EACH of them, the two chains four slots apart with a void between them.

                                    detached rings   chain-to-chain welded
    before the excursions bounded                0                   15.4%
    what you printed (round 2)                   7                    0.0%
    this plate                                   0                   68.7%

ONE object changes. P4_lobed and P5_crescent are byte-identical to the parts from round 2 - same
ring count, same path to 0.1 mm - so those parts ARE this plate's control and are not reprinted.
P5 is kept only as a peel reference in the same filament and session: it held last time and must
hold again.

PlateL proved the foam holds on ONE geometry - a concentric annulus, at two densities. Everything
that made it hold is geometry-dependent: the gap cap decides ring counts from the buried span,
the per-gap amplitude sizes every excursion from a measured gap, and the anchor pin depends on a
kept ring landing next to a shell. None of that is exercised by a round tube.

So this is the same fix on walls that vary, curve the wrong way, have no bore at all, or change
their ring count partway up:

    P1_sparse     annulus r10 / bore 2.5 at density 25   <- the validated baseline, repeated
    P2_dense      the same STL at density 50             <- ditto
    P3_ecc        bore offset 2 mm: wall thickness varies all the way round
    P4_lobed      5-lobe outer wall against a round bore
    P5_crescent   concave, and its centroid lies OUTSIDE the material
    P6_solid      disc r8, NO bore: no second face, so no anchor pair
    P7_taper      bore widens with height: the ring count CHANGES partway up, so the gap cap
                  re-decides mid-part and the layers either side must still tie to each other
    P8_thin       a 1 mm wall: nothing to sparsify, must decline and print sanely

Lids off, so the foam is visible from above without cutting. Peel P3 through P7 the way P1 and P2
were peeled - the bottom shell comes off first.

Meshes are validated, not assumed: every edge shared by exactly two triangles, no degenerate
triangles, and a fan cap only where its anchor provably sees the whole boundary. That check has
already caught one silently broken mesh in this suite.
"""
import json, math, os, re, subprocess, sys, zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_common import resolve_outdir, inject_env
from make_w1_plate import write_stl
from make_serp_plate import tube, solid, validate, star_ok, circle, lobed, ecc_pair, crescent

SUITE  = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED, PUBLISHING = resolve_outdir(SUITE)
CFG    = os.path.join(SUITE, "config")
REPO   = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE    = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE   = os.path.dirname(os.path.abspath(__file__))

H, SEG = 12.0, 256
OUTNAME = "PlateQ_changed"


def tube_z(outer_fn, inner_fn, h, levels=24):
    """Tube whose loops vary with height - so the wall thickness, and with it the ring count the
    sparsifier chooses, changes partway up the part."""
    zs = [h * k / levels for k in range(levels + 1)]
    rings = [(outer_fn(z), inner_fn(z)) for z in zs]
    n = len(rings[0][0])
    tris = []
    for L in range(levels):
        (o0, i0), (o1, i1) = rings[L], rings[L + 1]
        z0, z1 = zs[L], zs[L + 1]
        for k in range(n):
            j = (k + 1) % n
            O0 = lambda t: (t[0], t[1], z0)
            O1 = lambda t: (t[0], t[1], z1)
            tris += [(O0(o0[k]), O0(o0[j]), O1(o1[j])), (O0(o0[k]), O1(o1[j]), O1(o1[k]))]
            tris += [(O0(i0[k]), O1(i1[j]), O0(i0[j])), (O0(i0[k]), O1(i1[k]), O1(i1[j]))]
    o0, i0 = rings[0]
    o1, i1 = rings[-1]
    for k in range(n):
        j = (k + 1) % n
        b = lambda t: (t[0], t[1], 0.0)
        t = lambda q: (q[0], q[1], h)
        tris += [(b(o0[k]), b(i0[j]), b(o0[j])), (b(o0[k]), b(i0[k]), b(i0[j]))]
        tris += [(t(o1[k]), t(o1[j]), t(i1[j])), (t(o1[k]), t(i1[j]), t(i1[k]))]
    return tris


def main():
    objs = []

    def add(name, tris, dens=None):
        validate(tris, name)
        p = os.path.join(SUITE, name + ".stl")
        write_stl(p, tris)
        objs.append((p, dens))

    eo, ei = ecc_pair(10.0, 3.0, 2.0)
    add("P3_ecc", tube(eo, ei, H))
    # P4_lobed is byte-identical to the part in hand and is not reprinted. P5 is, as the peel
    # reference: it is also byte-identical, so if it peels this time the machine moved, not the code.
    cre = crescent(12.0, 9.0, 8.0)
    if not star_ok(cre, (-6.0, 0.0)):
        sys.exit("P5: the fan anchor does not see the whole crescent boundary")
    add("P5_crescent", solid(cre, (-6.0, 0.0), H))

    proc = json.load(open(os.path.join(CFG, "cs_shapes.json"), encoding="utf-8"))
    proc["continuous_shell_infill_density"]      = "25"
    proc["continuous_shell_interior_spiral"]     = "0"
    proc["continuous_shell_interior_serpentine"] = "0"      # the ring lattice: the validated line
    proc["top_shell_layers"]    = "0"
    proc["top_shell_thickness"] = "0"
    proc["name"] = "CS battery"
    cfgp = os.path.join(CFG, "cs_battery.json")
    json.dump(proc, open(cfgp, "w", encoding="utf-8"), indent=1)

    base = [EXE, "--load-settings", f"{os.path.join(CFG, 'machine.json')};{cfgp}",
            "--load-filaments", os.path.join(CFG, "filament.json")]
    probe = os.path.join(SLICED, OUTNAME + "_probe.gcode.3mf")
    if os.path.exists(probe):
        os.remove(probe)
    r = subprocess.run(base + ["--arrange", "1", "--slice", "0", "--export-3mf",
                               OUTNAME + "_probe.gcode.3mf", "--outputdir", SLICED]
                       + [p for p, _ in objs], capture_output=True, text=True, timeout=3600)
    if r.returncode or not os.path.exists(probe):
        print((r.stdout or "")[-1500:]); sys.exit("probe slice failed")

    over = {os.path.basename(p)[:-4]: d for p, d in objs if d}
    zin = zipfile.ZipFile(probe)
    ms  = zin.read("Metadata/model_settings.config").decode("utf-8")
    hits = [0]
    def tag(m):
        d = over.get(m.group(2).replace(".stl", ""))
        if d is None:
            return m.group(0)
        hits[0] += 1
        return m.group(0) + ('\n    <metadata key="continuous_shell_infill_density" value="%d"/>' % d)
    ms2, _ = re.subn(r'<object id="(\d+)">\s*<metadata key="name" value="([^"]+)"/>', tag, ms)
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
        print((r.stdout or "")[-1500:]); print((r.stderr or "")[-400:]); sys.exit("slice failed")
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

    w = re.search(r"weight[^\d]*([\d.]+)", meta)
    t2 = re.search(r"prediction[^\d]*(\d+)", meta)
    print(f"\n{OUTNAME}: {w.group(1) if w else '?'} g, "
          f"{int(t2.group(1)) // 60 if t2 else '?'} min")
    print("\nP5 is the question this time, and it is about the SOLID down the middle.")
    print("  P5: cut it across and look at the centre of the foam. The solid lens running along")
    print("      the medial should be roughly half the width it was - 10.7 mm2 of section down to")
    print("      4.6 - at the same material, with the cells either side of it unchanged. Peel it")
    print("      as well: the whole point is that halving it costs nothing, so it must still hold.")
    print("  P3: the control, and a real one this time - BIT-FOR-BIT identical to the part in")
    print("      hand (0.0000 mm displacement at every point), because its pocket chain is too")
    print("      short to qualify for the change. It is here to say whether the machine or the")
    print("      filament moved, and to re-check the peel that round 4 won. Its own middle solid")
    print("      is the pocket band and is NOT expected to change.")
    print("\nWhat would falsify this: P5's middle looking the same, which would mean what the eye")
    print("reads as solid is not the persistent-material measure; or P5 peeling at the medial,")
    print("which would mean the terminus was load-bearing frozen and the ring beside it was doing")
    print("work the weave cannot replace.")


if __name__ == "__main__":
    main()
