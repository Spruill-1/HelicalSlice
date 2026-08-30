"""A/B the two answers to a ring count that does not fit everywhere.

    sweep_interp.py

ONE COUNT, ANSWERED TO A PERCENTILE (CS_INTERP_PCT). The count serves the whole region. At 1.0
it answers to the widest place, so no gap exceeds the weld distance and the thin parts crowd.
Lower it and the thin parts stop crowding while the wide parts open past that distance - which
is the peel risk, measured twice in this campaign (2.251 mm held, 3.013 parted).

CELLS, EACH CHOOSING ITS OWN COUNT (CS_INTERP_CELLS). Each cell spans a narrower range of
extents, so its count fits it. The cost is an interface at every cut: two cells that choose
different counts put their level sets in different places, so their rings weld across the
boundary instead of continuing through it.

The columns that decide it: `collapse` is what each is trying to cure, `widest gap` is the peel
risk each is spending to do so, and `sub-floor` is how much of the print comes out below a
printable bead whatever the notes claim.
"""
import json, math, os, re, subprocess, sys, zipfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_common import inject_env

SUITE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CFG   = os.path.join(SUITE, "config")
REPO  = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE   = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
HERE  = os.path.dirname(os.path.abspath(__file__))
OUT   = os.path.join(HERE, "_verify")
OBJS  = ["P1_sparse", "P5_crescent", "P3_ecc"]
FLOOR = 0.34


def measure(path):
    with zipfile.ZipFile(path) as z:
        g = z.read("Metadata/plate_1.gcode").decode("utf-8", "replace")
    lh = float(re.search(r"; layer_height = ([\d.]+)", g).group(1))
    fd = float(re.search(r"; filament_diameter = ([\d.]+)", g).group(1))
    area = math.pi * fd * fd / 4.0
    cur = None; x = y = None; runs = {}; prev = False; L = {}; thin = {}; vol = {}; lay = -1
    for line in g.splitlines():
        if line.startswith("; CHANGE_LAYER"):
            lay += 1; prev = False; continue
        m = re.match(r"; printing object (\S+)", line)
        if m: cur = m.group(1).split(".stl")[0]; prev = False; continue
        if line.startswith("; stop printing object"): cur = None; prev = False; continue
        if not (line.startswith("G1 ") or line.startswith("G0 ")): continue
        nx = re.search(r" X(-?[\d.]+)", line); ny = re.search(r" Y(-?[\d.]+)", line)
        ne = re.search(r" E(-?[\d.]+)", line)
        px, py = x, y
        if nx: x = float(nx.group(1))
        if ny: y = float(ny.group(1))
        if cur is None or px is None: prev = False; continue
        moved = bool(nx or ny) and math.hypot(x - px, y - py) > 1e-9
        ext = bool(ne) and float(ne.group(1)) > 0
        if moved and ext:
            d = math.hypot(x - px, y - py)
            w = float(ne.group(1)) * area / (d * lh)
            L[cur] = L.get(cur, 0.) + d
            vol[cur] = vol.get(cur, 0.) + float(ne.group(1)) * area
            if w < FLOOR: thin[cur] = thin.get(cur, 0.) + d
            if not prev: runs[(cur, lay)] = runs.get((cur, lay), 0) + 1
            prev = True
        elif moved:
            prev = False
    out = {}
    for o in OBJS:
        ks = [k for k in runs if k[0] == o]
        out[o] = {"len": L.get(o, 0.), "vol": vol.get(o, 0.),
                  "thin": 100 * thin.get(o, 0.) / max(L.get(o, 1e-9), 1e-9),
                  "runs": max((runs[k] for k in ks), default=0)}
    return out


def run(tag, env_extra):
    for f in os.listdir(HERE):
        if re.fullmatch(r"\d+\.log", f):
            os.remove(os.path.join(HERE, f))
    name = "sweep_" + tag + ".gcode.3mf"
    dst = os.path.join(OUT, name)
    if os.path.exists(dst): os.remove(dst)
    env = dict(os.environ); env.update(env_extra)
    cmd = [EXE, "--debug", "3",
           "--load-settings", f"{os.path.join(CFG,'machine.json')};{os.path.join(CFG,'cs_interp.json')}",
           "--load-filaments", os.path.join(CFG, "filament.json"),
           "--arrange", "1", "--slice", "0", "--export-3mf", name, "--outputdir", OUT] \
          + [os.path.join(SUITE, o + ".stl") for o in OBJS]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, env=env)
    if r.returncode or not os.path.exists(dst):
        print(f"  {tag}: SLICE FAILED"); return None
    log = "".join(open(os.path.join(HERE, f), encoding="utf-8", errors="replace").read()
                  for f in sorted(os.listdir(HERE)) if re.fullmatch(r"\d+\.log", f))
    col = [float(v) for v in re.findall(r"still collapse over ([\d.]+)%", log)]
    gaps = [float(v) for v in re.findall(r"widest gap at ([\d.]+) mm", log)]
    cells = len(re.findall(r"subsectioned into", log))
    m = measure(dst)
    return {"collapse": max(col) if col else 0.0, "gap": max(gaps) if gaps else 0.0,
            "cells": cells, "m": m}


def main():
    print(__doc__)
    rows = [("baseline  pct=1.0", {}),
            ("pct=0.90", {"CS_INTERP_PCT": "0.90"}),
            ("pct=0.75", {"CS_INTERP_PCT": "0.75"}),
            ("pct=0.50", {"CS_INTERP_PCT": "0.50"}),
            ("cells=2",  {"CS_INTERP_CELLS": "2"}),
            ("cells=4",  {"CS_INTERP_CELLS": "4"}),
            ("cells=6",  {"CS_INTERP_CELLS": "6"})]
    print(f"{'setting':>18s} {'collapse':>9s} {'widest gap':>11s} {'cell layers':>12s} "
          f"{'sub-floor %':>28s} {'runs/layer':>12s} {'mm3':>8s}")
    for tag, env in rows:
        r = run(tag, env)
        if r is None: continue
        sub = " ".join(f"{r['m'][o]['thin']:5.1f}" for o in OBJS)
        runs = " ".join(f"{r['m'][o]['runs']:2d}" for o in OBJS)
        vol = sum(r['m'][o]['vol'] for o in OBJS)
        print(f"{tag:>18s} {r['collapse']:>8.1f}% {r['gap']:>11.2f} {r['cells']:>12d} "
              f"{sub:>28s} {runs:>12s} {vol:>8.0f}")
    print(f"\n  sub-floor % and runs/layer are per object, in the order {', '.join(OBJS)}")
    print("  runs/layer MUST stay 1: more than one means the layer stopped extruding mid-way.")


if __name__ == "__main__":
    main()
