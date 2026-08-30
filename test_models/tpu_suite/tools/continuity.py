"""Is the layer actually ONE continuous extrusion that never stops climbing?

That is the whole premise of continuous-shell mode, and it is a property of the emitted
G-code, not of the ring generator's intentions. Two things can quietly break it:

  * A layer that comes out as several paths. Anything that fragments a section - band
    marching splitting at a pinch, a partition cutting a cell - puts a travel move in the
    middle of the layer, and every travel is a stop, a retract and a seam.

  * A Z ramp that is not there. `apply_continuous_z` spreads the layer height across the
    layer by cumulative arc length, so the nozzle climbs the whole way round.

One caution about reading the numbers. G-code emits Z with three decimals, so on a 0.25 mm
layer the ramp can only ever show 250 distinct heights: most individual segments repeat the
previous Z simply because the true step is below a micron. Counting "segments whose Z
changed" therefore reports about 10% and means nothing. What is meaningful is the number of
DISTINCT heights per layer (250 of a possible 250 is a perfect ramp), that no extruding
segment descends, and that there are no travels at all.

    continuity.py <file.gcode.3mf> [<file.gcode.3mf> ...]
"""
import math, re, statistics, sys, zipfile


def measure(path):
    with zipfile.ZipFile(path) as z:
        g = z.read("Metadata/plate_1.gcode").decode("utf-8", "replace")
    lh = float(re.search(r"; layer_height = ([\d.]+)", g).group(1))

    cur = None; x = y = zc = None; layer = -1; prev_ext = False
    runs, travels, zvals, descents, seglen = {}, {}, {}, {}, {}
    for line in g.splitlines():
        if line.startswith("; CHANGE_LAYER"):
            layer += 1; prev_ext = False; continue
        m = re.match(r"; printing object (\S+)", line)
        if m:
            cur = m.group(1).split(".stl")[0]; prev_ext = False; continue
        if line.startswith("; stop printing object"):
            cur = None; prev_ext = False; continue
        if not (line.startswith("G1 ") or line.startswith("G0 ")):
            continue
        nx = re.search(r" X(-?[\d.]+)", line); ny = re.search(r" Y(-?[\d.]+)", line)
        nz = re.search(r" Z(-?[\d.]+)", line); ne = re.search(r" E(-?[\d.]+)", line)
        px, py, pz = x, y, zc
        if nx: x = float(nx.group(1))
        if ny: y = float(ny.group(1))
        if nz: zc = float(nz.group(1))
        if cur is None or px is None:
            prev_ext = False; continue
        moved = bool(nx or ny) and math.hypot(x - px, y - py) > 1e-9
        ext   = bool(ne) and float(ne.group(1)) > 0
        k = (cur, layer)
        if moved and ext:
            if not prev_ext:
                runs[k] = runs.get(k, 0) + 1
            seglen[k] = seglen.get(k, 0.0) + math.hypot(x - px, y - py)
            if nz:
                zvals.setdefault(k, set()).add(zc)
                if pz is not None and zc < pz - 1e-9:
                    descents[k] = descents.get(k, 0) + 1
            prev_ext = True
        elif moved:
            travels[k] = travels.get(k, 0) + 1
            prev_ext = False

    out = {}
    for o in sorted({k[0] for k in runs}):
        ks = [k for k in runs if k[0] == o]
        r  = sorted(runs[k] for k in ks)
        out[o] = {
            "layers":       len(ks),
            "runs_med":     r[len(r) // 2],
            "runs_max":     r[-1],
            "travels":      sum(travels.get(k, 0) for k in ks),
            "z_heights":    int(statistics.median(len(zvals.get(k, ())) for k in ks)),
            "z_possible":   int(round(lh * 1000)),      # 1 um is all G-code can express
            "descents":     sum(descents.get(k, 0) for k in ks),
            "mm_per_layer": statistics.median(seglen[k] for k in ks),
        }
    return out


def main(argv):
    if not argv:
        print(__doc__); return 1
    for path in argv:
        print(f"\n{path}")
        print(f"  {'object':>13s} {'layers':>6s} {'runs/layer':>11s} {'travels':>8s} "
              f"{'Z heights/layer':>16s} {'descents':>9s} {'mm/layer':>9s}")
        for o, m in measure(path).items():
            runs = "med %d max %d" % (m["runs_med"], m["runs_max"])
            zed  = "%d of %d" % (m["z_heights"], m["z_possible"])
            print(f"  {o:>13s} {m['layers']:>6d} {runs:>11s} {m['travels']:>8d} "
                  f"{zed:>16s} {m['descents']:>9d} {m['mm_per_layer']:>9.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
