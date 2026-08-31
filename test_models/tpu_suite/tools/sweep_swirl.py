"""Does the swirl knob turn the interior, and by how much?

    sweep_swirl.py

`continuous_shell_swirl` is degrees per layer, and it has to move two things that used to be
driven separately and disagreed by about 26x: the ring-to-ring JOINS (a seam-jitter scale) and
the LATTICE (a phase advance whose angle depends on the lobe count).

Three earlier versions of this script reported the knob doing nothing, and all three were
measurement artefacts rather than findings. They are worth naming, because the same traps are
sitting in every other instrument here:

  the first sharp joint's azimuth      not a stable feature - it is whichever joint sorts first,
                                       and which one that is changes between layers
  correlating rasterised layers        P3_ecc's bore is off-centre, so 63% of the radial signal
                                       is k=1 and IDENTICAL every layer because it is the part.
                                       The correlation peak pins to zero no matter what moves.
  the same, with k<3 filtered out      the skin and bore loops are long, high-curvature and also
                                       layer-invariant, and they still dominated the residual.

So the lattice is measured on the rings THEMSELVES, dumped by CS_DUMPRINGS, one ring at a time -
no centroid shared between rings, nothing from the part's outline in the signal. The joins are
measured by correlating the whole SET of joint azimuths between layers rather than one of them.

Measured on P3_ecc, layers 24-25, weave 50%:

    knob      joins    lattice, outer ring -> inner ring
     0.0        1.0    -8.8  -7.8  -6.8  -5.2 -13.2 -11.5  -3.8      (the old, sheared behaviour)
     7.5        7.0    +7.2  +6.0  +5.0  +4.2  +3.2  +2.8  +2.2
    15.0       17.0   +14.5 +13.5 +12.8 +12.0  +7.5  +6.2  +3.2

Surfaces read 0.0 at every setting, which is the invariant that matters most: the outer loop and
the bore loop must not move. The lattice decays inward because the weave is spaced along the
path, not by angle, and equal steps of arc length are smaller angles near a bore. Equalising the
degrees would pull the octave welds apart, so it is documented, not corrected.
"""
import json, math, os, subprocess, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sweep_overshoot import segments

HERE  = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.abspath(os.path.join(HERE, ".."))
CFG   = os.path.join(SUITE, "config")
REPO  = os.path.abspath(os.path.join(SUITE, "..", ".."))
EXE   = os.environ.get("ORCA_EXE") or os.path.join(REPO, "build", "src", "RelWithDebInfo", "orca-slicer.exe")
OUT   = os.path.join(HERE, "_verify")
NA, KMIN, LAYERS = 1440, 3, (24, 25)


def slice_at(deg, layer=None, obj="P3_ecc"):
    """Slice once. With `layer` set, also dump that layer's rings and return them."""
    base = json.load(open(os.path.join(CFG, "cs_w50.json"), encoding="utf-8"))
    base["continuous_shell_swirl"] = str(deg)
    cfgp = os.path.join(CFG, f"cs_swirl{deg}.json")
    json.dump(base, open(cfgp, "w", encoding="utf-8"), indent=1)
    tag = f"swirl{str(deg).replace('.', '_')}" + (f"_L{layer}" if layer is not None else "")
    dst = os.path.join(OUT, tag + ".gcode.3mf")
    cache = os.path.join(OUT, f"rings_{deg}_{layer}.txt")
    env = dict(os.environ); env["CS_OVERSHOOT"] = "1.0"
    if layer is not None:
        # Both halves come from the SAME slice: the rings for the lattice, the G-code for the
        # joins. Reusing a cached dump whose 3mf has been cleaned away would measure two
        # different runs against each other.
        if os.path.exists(cache) and os.path.exists(dst):
            return dst, read_rings(cache)
        env["CS_DUMPRINGS"] = str(layer)
    live = os.path.join(HERE, "cs_rings.txt")
    if os.path.exists(live):
        os.remove(live)
    subprocess.run([EXE, "--load-settings", f"{os.path.join(CFG,'machine.json')};{cfgp}",
                    "--load-filaments", os.path.join(CFG, "filament.json"),
                    "--slice", "0", "--export-3mf", tag + ".gcode.3mf",
                    "--outputdir", OUT, os.path.join(SUITE, obj + ".stl")],
                   capture_output=True, text=True, timeout=3600, env=env, cwd=HERE)
    if layer is None:
        return (dst if os.path.exists(dst) else None), None
    if not os.path.exists(live):
        return dst, []
    os.replace(live, cache)
    return dst, read_rings(cache)


def read_rings(path):
    out = []
    for line in open(path):
        t = line.split()
        if not t or t[0] != "CHAIN":
            continue
        n = int(t[9])
        out.append(np.array([float(x) for x in t[10:10 + 2 * n]]).reshape(-1, 2))
    return out


def shift(a, b, na):
    """angular shift, in degrees, that best maps profile a onto profile b"""
    x = np.fft.irfft(np.conj(np.fft.rfft(a)) * np.fft.rfft(b), n=na)
    k = int(np.argmax(x))
    return (k if k < na / 2 else k - na) * 360.0 / na


def ring_profile(v):
    """One ring's radius against its own azimuth, with the part's own shape removed."""
    c = v.mean(0); d = v - c
    r = np.linalg.norm(d, axis=1)
    a = np.arctan2(d[:, 1], d[:, 0]) % (2 * math.pi)
    b = (a / (2 * math.pi) * NA).astype(int) % NA
    acc = np.zeros(NA); cnt = np.zeros(NA)
    np.add.at(acc, b, r); np.add.at(cnt, b, 1.0)
    ok = cnt > 0
    p = np.interp(np.arange(NA), np.arange(NA)[ok], acc[ok] / cnt[ok], period=NA)
    F = np.fft.rfft(p - p.mean()); F[:KMIN] = 0
    return np.fft.irfft(F, n=NA)


def join_profile(path, layer):
    """Histogram of the azimuths at which the path turns sharply - the ring-to-ring hops."""
    segs, _ = segments(path, layer)
    if not segs:
        return None
    P = np.array([[s[0], s[1]] for s in segs]); C = P.mean(0)
    v = np.diff(np.vstack([P, P[:1]]), axis=0)
    n = np.linalg.norm(v, axis=1); keep = n > 1e-9
    vv = v[keep] / n[keep, None]; Pk = P[keep]
    t = np.degrees(np.arccos(np.clip((vv[:-1] * vv[1:]).sum(axis=1), -1, 1)))
    idx = np.where(t > 45)[0]
    if len(idx) < 4:
        return None
    a = np.arctan2(Pk[idx, 1] - C[1], Pk[idx, 0] - C[0]) % (2 * math.pi)
    h = np.zeros(720)
    np.add.at(h, (a / (2 * math.pi) * 720).astype(int) % 720, 1.0)
    return h - h.mean()


def main():
    print(__doc__)
    print(f"  {'knob':>6s} {'joins':>7s}   lattice, ring by ring (surfaces should read 0.0)")
    for deg in (0.0, 7.5, 15.0):
        pa, A = slice_at(deg, LAYERS[0])
        pb, B = slice_at(deg, LAYERS[1])
        if not A or not B:
            print(f"  {deg:>6.1f}   no ring dump - is CS_DUMPRINGS still in ContinuousShell.cpp?")
            continue
        ja = join_profile(pa, LAYERS[0]); jb = join_profile(pb, LAYERS[1])
        j = shift(ja, jb, 720) if ja is not None and jb is not None else float("nan")
        cells = []
        for i in range(min(len(A), len(B))):
            if len(A[i]) < 64 or len(B[i]) < 64:
                continue
            cells.append(f"{shift(ring_profile(A[i]), ring_profile(B[i]), NA):+6.1f}")
        print(f"  {deg:>6.1f} {j:>7.1f}   " + " ".join(cells))
    print("\n  joins should track the knob (plus ~1.1 deg/layer of pre-existing seam step);")
    print("  the lattice should track it at the outer rings and decay toward a bore.")


if __name__ == "__main__":
    main()
