"""Build PlateN_serp: the radial serpentine across the shapes the slicer actually meets.

    make_serp_plate.py [--publish]

Everything known about the serpentine so far is geometry. It has never been printed. This plate
is the first physical evidence, and it is built to answer two questions at once.

FIRST, does it beat the ring lattice on the failure that motivated it? N0 and N1 are the SAME
STL, the de-lid annulus, differing only in whether the serpentine is on. Peel both. Without that
control the answer is unattributable - PlateG already taught this campaign that lesson once, when
an uncontrolled "enormous strength improvement" turned out to be present in the untouched part
too.

SECOND, does it survive the shapes it claims to? The generator was rebuilt in ring-index space
precisely so it would not be confined to round tubes, and that claim is untested in plastic:

    N0_rings      annulus r10 / bore 2.5, ring lattice          <- the control
    N1_serp       the same STL with the serpentine on           <- the comparison
    N2_solid      disc r8, NO bore: sweeps to the medial terminus and back
    N3_crescent   concave, and its centroid lies OUTSIDE the material
    N4_ecc        bore offset 2 mm: wall thickness varies all the way round
    N5_lobed      5-lobe outer wall against a round bore
    N6_thin       a wall too thin to traverse: must DECLINE and print sanely, not crash

Lids off, so the interior is open from above and the traverse is visible without cutting.

Mesh note, because this campaign has lost time to it before: nothing here needs a general
polygon-with-holes tessellator. A tube is capped by a STRIP between its two loops sampled on a
common parameter, and a solid section by a FAN from an anchor that is checked to see the whole
boundary. Both are validated below rather than assumed - an unnoticed non-manifold mesh is a
whole wasted print.
"""
import json, math, os, re, subprocess, sys, zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_common import resolve_outdir, inject_env
from make_w1_plate import write_stl

SUITE  = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SLICED, PUBLISHING = resolve_outdir(SUITE)
CFG    = os.path.join(SUITE, "config")
REPO   = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE    = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE   = os.path.dirname(os.path.abspath(__file__))

H, SEG = 12.0, 256
OUTNAME = "PlateN_serp"


# ---------------------------------------------------------------- mesh helpers
def tube(outer, inner, h):
    """Closed mesh between two loops sampled on a COMMON parameter, swept to height h.

    The caps are a strip between outer[i] and inner[i], which is only valid while both loops are
    traversed in the same direction and neither doubles back relative to the other - true for
    every tube here, and checked by validate() rather than trusted."""
    n = len(outer)
    assert len(inner) == n
    o0 = [(p[0], p[1], 0.0) for p in outer]; o1 = [(p[0], p[1], h) for p in outer]
    i0 = [(p[0], p[1], 0.0) for p in inner]; i1 = [(p[0], p[1], h) for p in inner]
    tris = []
    for k in range(n):
        j = (k + 1) % n
        tris += [(o0[k], o0[j], o1[j]), (o0[k], o1[j], o1[k])]      # outer wall
        tris += [(i0[k], i1[j], i0[j]), (i0[k], i1[k], i1[j])]      # bore wall
        tris += [(o0[k], i0[j], o0[j]), (o0[k], i0[k], i0[j])]      # bottom
        tris += [(o1[k], o1[j], i1[j]), (o1[k], i1[j], i1[k])]      # top
    return tris


def solid(loop, anchor, h):
    """Closed mesh of a simple (hole-free) section, capped by a fan from `anchor`."""
    n = len(loop)
    a = anchor
    l0 = [(p[0], p[1], 0.0) for p in loop]; l1 = [(p[0], p[1], h) for p in loop]
    tris = []
    for k in range(n):
        j = (k + 1) % n
        tris += [(l0[k], l0[j], l1[j]), (l0[k], l1[j], l1[k])]      # wall
        tris += [((a[0], a[1], 0.0), l0[j], l0[k])]                 # bottom fan
        tris += [((a[0], a[1], h), l1[k], l1[j])]                   # top fan
    return tris


def validate(tris, name):
    """Every edge must be shared by exactly two triangles. A mesh that fails this slices into
    nonsense and the print is wasted before it starts."""
    from collections import Counter
    edges = Counter()
    for t in tris:
        for k in range(3):
            a, b = t[k], t[(k + 1) % 3]
            edges[(a, b) if a < b else (b, a)] += 1
    bad = [e for e, c in edges.items() if c != 2]
    if bad:
        sys.exit(f"{name}: mesh is not closed - {len(bad)} edges shared by != 2 triangles")
    degenerate = sum(1 for t in tris if len(set(t)) < 3)
    if degenerate:
        sys.exit(f"{name}: {degenerate} degenerate triangles")
    print(f"   {name}: {len(tris)} triangles, closed and manifold")


def star_ok(loop, anchor):
    """Does `anchor` see the whole boundary? A fan cap is only valid if it does."""
    prev = None
    total = 0.0
    for p in loop:
        a = math.atan2(p[1] - anchor[1], p[0] - anchor[0])
        if prev is not None:
            d = a - prev
            while d > math.pi:  d -= 2 * math.pi
            while d < -math.pi: d += 2 * math.pi
            total += d
        prev = a
    return abs(abs(total) - 2 * math.pi) < 0.2


# ---------------------------------------------------------------- the shapes
def circle(r, cx=0.0, cy=0.0, n=SEG, phase=0.0):
    return [(cx + r * math.cos(2 * math.pi * i / n + phase),
             cy + r * math.sin(2 * math.pi * i / n + phase)) for i in range(n)]


def lobed(r, amp, lobes, n=SEG):
    return [(( r + amp * math.cos(lobes * 2 * math.pi * i / n)) * math.cos(2 * math.pi * i / n),
             ( r + amp * math.cos(lobes * 2 * math.pi * i / n)) * math.sin(2 * math.pi * i / n))
            for i in range(n)]


def ecc_pair(r_out, r_bore, dx, n=SEG):
    """Outer and bore loops parameterised by angle about the BORE centre, so the cap strip
    between them never doubles back."""
    o, i = [], []
    for k in range(n):
        a = 2 * math.pi * k / n
        # ray from the bore centre at angle a, out to the outer circle (centred at the origin)
        ux, uy = math.cos(a), math.sin(a)
        b = 2.0 * (dx * ux)
        c = dx * dx - r_out * r_out
        t = (-b + math.sqrt(b * b - 4 * c)) / 2.0
        o.append((dx + t * ux, t * uy))
        i.append((dx + r_bore * ux, r_bore * uy))
    return o, i


def crescent(r_out, r_bite, bite_x, n=SEG):
    """Outer disc minus an off-centre disc that pokes out through its edge, so the result is a
    simple hole-free section - concave, with its centroid outside the material."""
    # circle intersections
    d = bite_x
    x = (r_out * r_out - r_bite * r_bite + d * d) / (2 * d)
    y = math.sqrt(max(r_out * r_out - x * x, 0.0))
    a1 = math.atan2(-y, x)                       # on the outer circle
    a2 = math.atan2(y, x)
    b1 = math.atan2(-y, x - d)                   # the same points on the bite circle
    b2 = math.atan2(y, x - d)
    pts = []
    # the LONG way round the outer circle, from a2 to a1
    span = (a1 - a2) % (2 * math.pi)
    for k in range(n):
        a = a2 + span * k / n
        pts.append((r_out * math.cos(a), r_out * math.sin(a)))
    # Back along the bite arc, from b1 to b2 the way that stays INSIDE the outer disc - through
    # (d - r_bite, 0). The other way round leaves the disc entirely, and since the two endpoints
    # sit near +/-90 deg the two spans are nearly equal, so picking by "shorter" would have been
    # a coin flip.
    span = (b1 - b2) % (2 * math.pi)
    for k in range(n // 2):
        a = b1 - span * k / (n // 2)
        pts.append((d + r_bite * math.cos(a), r_bite * math.sin(a)))
    return pts


def main():
    objs, notes = [], []

    ann_o, ann_i = circle(10.0), circle(2.5)
    t = tube(ann_o, ann_i, H); validate(t, "N0/N1 annulus")
    a0 = os.path.join(SUITE, "N0_rings.stl"); write_stl(a0, t); objs.append(a0)
    a1 = os.path.join(SUITE, "N1_serp.stl");  write_stl(a1, t); objs.append(a1)

    disc = circle(8.0)
    t = solid(disc, (0.0, 0.0), H); validate(t, "N2 solid")
    p = os.path.join(SUITE, "N2_solid.stl"); write_stl(p, t); objs.append(p)

    cre = crescent(12.0, 9.0, 8.0)
    anchor = (-6.0, 0.0)
    if not star_ok(cre, anchor):
        sys.exit("N3: the fan anchor does not see the whole crescent boundary")
    t = solid(cre, anchor, H); validate(t, "N3 crescent")
    p = os.path.join(SUITE, "N3_crescent.stl"); write_stl(p, t); objs.append(p)

    eo, ei = ecc_pair(10.0, 3.0, 2.0)
    t = tube(eo, ei, H); validate(t, "N4 eccentric")
    p = os.path.join(SUITE, "N4_ecc.stl"); write_stl(p, t); objs.append(p)

    t = tube(lobed(10.0, 1.2, 5), circle(3.0), H); validate(t, "N5 lobed")
    p = os.path.join(SUITE, "N5_lobed.stl"); write_stl(p, t); objs.append(p)

    t = tube(circle(6.0), circle(4.6), H); validate(t, "N6 thin")
    p = os.path.join(SUITE, "N6_thin.stl"); write_stl(p, t); objs.append(p)

    proc = json.load(open(os.path.join(CFG, "cs_shapes.json"), encoding="utf-8"))
    proc["continuous_shell_infill_density"]      = "25"
    proc["continuous_shell_interior_spiral"]     = "0"
    proc["continuous_shell_interior_serpentine"] = "1"    # N0 is overridden back off, per object
    proc["top_shell_layers"]    = "0"
    proc["top_shell_thickness"] = "0"
    proc["name"] = "CS serpentine"
    cfgp = os.path.join(CFG, "cs_serp.json")
    json.dump(proc, open(cfgp, "w", encoding="utf-8"), indent=1)

    base = [EXE, "--load-settings", f"{os.path.join(CFG, 'machine.json')};{cfgp}",
            "--load-filaments", os.path.join(CFG, "filament.json")]
    probe = os.path.join(SLICED, OUTNAME + "_probe.gcode.3mf")
    if os.path.exists(probe):
        os.remove(probe)
    r = subprocess.run(base + ["--arrange", "1", "--slice", "0", "--export-3mf",
                               OUTNAME + "_probe.gcode.3mf", "--outputdir", SLICED] + objs,
                       capture_output=True, text=True, timeout=3600)
    if r.returncode or not os.path.exists(probe):
        print((r.stdout or "")[-1500:]); sys.exit("probe slice failed")

    # N0 is the CONTROL: the same STL with the serpentine switched back off.
    zin = zipfile.ZipFile(probe)
    ms  = zin.read("Metadata/model_settings.config").decode("utf-8")
    hits = [0]
    def tag(m):
        if "N0_rings" not in m.group(2):
            return m.group(0)
        hits[0] += 1
        return m.group(0) + ('\n    <metadata key="continuous_shell_interior_serpentine" value="0"/>')
    ms2, _ = re.subn(r'<object id="(\d+)">\s*<metadata key="name" value="([^"]+)"/>', tag, ms)
    if hits[0] != 1:
        sys.exit(f"expected to tag exactly one N0_rings object, tagged {hits[0]}")
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
    for ln in (r.stdout or "").splitlines():
        if "serpentine" in ln or "declined" in ln:
            notes.append(ln.strip()[-170:])
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

    if notes:
        print("\nslicer notes (declines are expected on N6, and are the point of it):")
        for ln in dict.fromkeys(notes):
            print("   ", ln)

    lh = float(re.search(r"; layer_height = ([\d.]+)", g).group(1))
    fd = float(re.search(r"; filament_diameter = ([\d.]+)", g).group(1))
    w = re.search(r"weight[^\d]*([\d.]+)", meta)
    t2 = re.search(r"prediction[^\d]*(\d+)", meta)
    print(f"\n{OUTNAME}: {w.group(1) if w else '?'} g, {int(t2.group(1))//60 if t2 else '?'} min")

    # Per object: how far a single uninterrupted run sweeps through the wall's DEPTH.
    #
    # Counting direction reversals does not separate these two generators - the ring lattice
    # undulates six times a revolution on every one of its rings, so it racks up MORE
    # reversals (121 against 61 on the control pair) while never actually crossing the wall.
    # What separates them is the SPAN of each monotonic run: a woven ring turns back within
    # about twice its weave amplitude, a traverse runs from one face to the other.
    from make_w1_plate import split_objects
    print(f"\n{'object':>16s} {'mid-layer':>11s} {'wall':>8s} {'p90 run':>9s} {'deep runs':>10s}")
    for name, segs in split_objects(g, lh, math.pi * fd * fd / 4.0):
        mid = [q for q in segs if abs(q[4] - H * 0.5) < lh]
        if not mid:
            print(f"{name:>16s}   no segments at mid height")
            continue
        cx = sum(q[0] for q in mid) / len(mid)
        cy = sum(q[1] for q in mid) / len(mid)
        d = [math.hypot(q[0] - cx, q[1] - cy) for q in mid]
        wall = max(d) - min(d)
        spans, run_start, direction = [], 0, 0
        for i in range(1, len(d)):
            step = d[i] - d[i - 1]
            sgn = (step > 1e-6) - (step < -1e-6)
            if sgn == 0:
                continue
            if direction and sgn != direction:
                spans.append(abs(d[i - 1] - d[run_start]))
                run_start = i - 1
            direction = sgn
        spans.append(abs(d[-1] - d[run_start]))
        spans.sort()
        p90 = spans[int(len(spans) * 0.9)] if spans else 0.0
        deep = sum(1 for v in spans if v > 0.5 * wall)
        print(f"{name:>16s} {sum(q[6] for q in mid):9.0f}mm {wall:7.2f}mm "
              f"{p90:8.2f}mm {deep:10d}")
    print("\n'p90 run' is the depth one monotonic run sweeps; 'deep runs' counts those")
    print("crossing more than half the wall. A ring lattice cannot cross it at all.")
    print("\nN0 vs N1 is the control pair - same STL, serpentine off then on. Peel BOTH.")
    print("N6 should have declined; the note above says so if it did.")


if __name__ == "__main__":
    main()
