"""Build PlateR_interp: the first plastic for interpolated rings.

    make_interp_plate.py [--publish]

ROUND 8 (2026-08-30) - stop marching the rings and interpolate them.

Every round up to now has tuned how a MARCHED ring stack is decimated and woven. This one
changes where the rings are. Measured on real slices, consecutive undecimated offsets sit
1.00 bead slots apart at the median and 11.58 at the maximum, so a slot is not a distance -
and the medial pocket, the tip voids, the solid core and the gap cap are all that one fact.
Rings are now level sets of

    u(p) = d_out(p) / (d_out(p) + d_in(p))

between the inner edge of the peeled skin and the bore wall, so the gap is W(p)/(n+1):
even in the FRACTION, by construction. Nothing is decimated - the count is chosen so that
the widest place in the region stays inside the distance two rings can still weld across.

Three objects, and each is on the plate to answer a different question.

  P1_sparse   The tube, and the object rounds 6 and 7 were about. Its failure is a ~1 mm
              void in the middle of the wall that "more or less forms a tube of its own",
              crossed only at the helical seam. Interpolation does NOT close that void -
              the anti-weld gap is inherent to a woven lattice and belongs to the amplitude
              law, which has not changed. What it does is make every gap in the wall the
              same, so the void should read as an even lattice rather than one fat ring
              flanked by tight ones. THE TEST IS STILL THE PEEL, and the number that has
              predicted every peel in this campaign is the widest untied void: 2.251 mm
              held, 3.013 parted.

  P5_crescent The object whose tips coalesce. Its wall runs from 9.6 mm at the back to a
              fraction of a millimetre at the tips, so a fixed ring count crowds at the
              tips by construction. The bench measurement said "still a slight solid bar in
              the middle". Interpolation removes the pocket that bar sat in; whether the
              tips now over-pack instead is exactly what the degradation ladder is for, and
              this plate is what says whether the ladder is needed yet.

  P3_ecc      The object that took three rounds to stop coming apart, and the one this
              change should transform rather than improve. Marched, this section fragments
              into four chains; interpolated it is one nested stack with no chain boundary
              anywhere, and chain boundaries are what the pockets were made of. (Both print
              as ONE continuous run per layer - the connector stitching already saw to that,
              and it is checked below - so the difference is in where the rings are, not in
              whether the extrusion is continuous.) If the solid crescent of material
              between the bore and the wall is gone, that is this change working.

What to look for, in order:
  1. Peel each object. P1 first.
  2. Cut P3 and look for the solid crescent. It should not be there.
  3. Cut P5 and look at the tips: coalesced (over-solid) or open (under-filled)?
  4. Outer surface finish on all three. The skin is peeled and never subsectioned, so it
     should be at least as good as round 7 - a visible seam change would mean the invariant
     is not holding in practice.

The requested density is 25%, and it is NOT what gets printed. How far apart two rings may
sit and still weld is the bound a peel test set at 3 bead slots - 3 held, 4 parted the tube -
and 25% asks for 4. Marched, that bound was a cap on a decimator working over ring positions
marching had already fixed, so it landed at 2.7 slots and filled at 38% while reporting
neither number honestly. Placing rings instead, the bound is honoured exactly: 3 rings, a
widest gap of 1.15 mm against the 1.19 mm bound, filling at 39%.

So this plate is NOT a density change dressed up as a geometry change. It prints about the
same amount of material as round 7 did, in the same places measured radially, and what
differs is that every gap in the wall is now the same fraction of it.
"""
import json, os, re, subprocess, sys, zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_common import resolve_outdir, inject_env
from make_w1_plate import write_stl
from make_serp_plate import tube, solid, validate, star_ok, circle, ecc_pair, crescent

SUITE  = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED, PUBLISHING = resolve_outdir(SUITE)
CFG    = os.path.join(SUITE, "config")
REPO   = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE    = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE   = os.path.dirname(os.path.abspath(__file__))

H = 12.0
OUTNAME = "PlateR_interp"


def main():
    objs = []

    def add(name, tris):
        validate(tris, name)
        p = os.path.join(SUITE, name + ".stl")
        write_stl(p, tris)
        objs.append(p)

    add("P1_sparse", tube(circle(10.0), circle(2.5), H))

    cre = crescent(12.0, 9.0, 8.0)
    if not star_ok(cre, (-6.0, 0.0)):
        sys.exit("P5: the fan anchor does not see the whole crescent boundary")
    add("P5_crescent", solid(cre, (-6.0, 0.0), H))

    outer, bore = ecc_pair(10.0, 2.5, 4.0)
    add("P3_ecc", tube(outer, bore, H))

    proc = json.load(open(os.path.join(CFG, "cs_shapes.json"), encoding="utf-8"))
    proc["continuous_shell_infill_density"]      = "25"
    proc["continuous_shell_interpolated_rings"]  = "1"      # the change under test
    proc["continuous_shell_interior_spiral"]     = "0"
    proc["continuous_shell_interior_serpentine"] = "0"
    proc["top_shell_layers"]    = "0"                        # open top: the lattice is the point
    proc["top_shell_thickness"] = "0"
    proc["name"] = "CS interpolated"
    cfgp = os.path.join(CFG, "cs_interp.json")
    json.dump(proc, open(cfgp, "w", encoding="utf-8"), indent=1)

    # --debug 3 raises the log level far enough for the per-layer degrade notes. They go to the
    # CLI's rotating log file, not the console - see the note counters below.
    base = [EXE, "--debug", "3",
            "--load-settings", f"{os.path.join(CFG, 'machine.json')};{cfgp}",
            "--load-filaments", os.path.join(CFG, "filament.json")]

    out = os.path.join(SLICED, OUTNAME + ".gcode.3mf")
    if os.path.exists(out):
        os.remove(out)
    for f in os.listdir(HERE):                       # this run's notes must be this run's
        if re.fullmatch(r"\d+\.log", f):
            os.remove(os.path.join(HERE, f))
    r = subprocess.run(base + ["--arrange", "1", "--slice", "0",
                               "--export-3mf", OUTNAME + ".gcode.3mf",
                               "--outputdir", SLICED] + objs,
                       capture_output=True, text=True, timeout=3600)
    if r.returncode or not os.path.exists(out):
        print((r.stdout or "")[-2500:]); print((r.stderr or "")[-600:]); sys.exit("slice failed")
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

    # The degrade notes are the whole point of this round: a wall that fell back to marching is
    # not testing anything, and a wall still reporting the gap cap has not taken the new path.
    #
    # They do NOT come back on stdout. The CLI logs to a rotating NNNNN.log in its working
    # directory and writes nothing to the console at all, so counting notes in r.stdout scored
    # zero on every run - which reads exactly like "the gap cap never fired" and was really
    # "nothing was listening". Stale logs are cleared above so this run's file is this run's.
    log = "".join(open(os.path.join(HERE, f), encoding="utf-8", errors="replace").read()
                  for f in sorted(os.listdir(HERE)) if re.fullmatch(r"\d+\.log", f))
    if not log:
        sys.exit("no slicer log found - the note counters below would be meaningless")
    layers    = len(re.findall(r"Continuous shell, layer", log))
    fell_back = len(re.findall(r"interpolated rings declined", log))
    capped    = len(re.findall(r"lattice gap cap", log))
    weld      = len(re.findall(r"limited by the weld distance", log))
    rounded   = len(re.findall(r"infill density rounded", log))
    collapse  = len(re.findall(r"interpolated rings collapse", log))
    print(f"\n  layers reporting anything:           {layers}")
    print(f"  layers that fell back to marching:   {fell_back}   (solid layers only)")
    print(f"  layers still reporting the gap cap:  {capped}   (must be 0)")
    print(f"  layers clamped by the weld distance: {weld}")
    print(f"  layers with the ring count rounded:  {rounded}")
    print(f"  layers whose rings collapse:         {collapse}   (the ladder's work, not built)")
    # A solid layer declining is CORRECT: a ring count fixed by the widest place is a lattice
    # rule, and solid fill needs the count to follow the local thickness - which Arachne's
    # per-angle beading already does properly. Any OTHER decline means this plate is not
    # testing the change on that layer.
    reasons    = re.findall(r"interpolated rings declined: ([^;\n]*)", log)
    unexpected = [w.strip() for w in reasons if "this layer is solid" not in w]
    if capped or unexpected:
        for why in sorted(set(unexpected)):
            print("   ", why)
        sys.exit("some sparse layers did not take the interpolated path - this plate would not test it")
    print(f"  (all {len(reasons)} of those are solid bottom layers, which should decline)")

    # The premise of the whole mode: one unbroken extrusion per layer that never stops climbing.
    # Ring generation sits upstream of the Z ramp and the connector stitching, so it can only
    # break this, never improve it - which is why it is checked rather than assumed.
    sys.path.insert(0, HERE)
    from continuity import measure
    print()
    for o, m in measure(out).items():
        print(f"  {o:>13s}: {m['runs_med']} run/layer (max {m['runs_max']}), {m['travels']} travels, "
              f"{m['z_heights']} of {m['z_possible']} Z heights per layer")
        if m["runs_max"] > 1 or m["travels"] > len(objs):
            sys.exit(f"{o} does not print as one continuous extrusion per layer")

    w  = re.search(r"weight[^\d]*([\d.]+)", meta)
    t2 = re.search(r"prediction[^\d]*(\d+)", meta)
    print(f"\n{OUTNAME}: {w.group(1) if w else '?'} g, "
          f"{int(t2.group(1)) // 60 if t2 else '?'} min")
    print("\nWhat to look for, in order:")
    print("  1. Peel each object, P1 first. The number that has predicted every peel in this")
    print("     campaign is the widest untied void: 2.251 mm held, 3.013 parted.")
    print("  2. Cut P3_ecc and look for the solid crescent between bore and wall. Marched, this")
    print("     section fragments into four chains; interpolated it is one, and chain boundaries")
    print("     are what the pockets were made of. If the crescent is gone, that is the change.")
    print("  3. Cut P5_crescent and look at the tips: coalesced (over-solid) or open? A fixed")
    print("     ring count crowds where the wall thins, and that is what the ladder is for.")
    print("  4. Outer surface finish on all three. The skin is peeled and never subsectioned, so")
    print("     a visible seam change would mean that invariant is not holding in practice.")
    print(f"\n  file: {out}")
    if not PUBLISHING:
        print("  (scratch build - re-run with --publish to write the file you intend to print)")


if __name__ == "__main__":
    main()
