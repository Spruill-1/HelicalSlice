"""Prototype of the radial-traverse interior, to fix its parameterisation before any C++.

The ring lattice has no radial members: its only radial connections are where two rings happen
to touch, so a continuous cylindrical void always exists between them and a peel needs exactly
one. A serpentine that sweeps ACROSS the wall as it goes round is made of radial members - every
traverse crosses every radius - so no such void can exist by construction.

    r(theta) = R_mid(theta) + A(theta) * tri(N * theta + phase)

with tri a triangle wave in [-1, 1] and A half the buried thickness less a margin, so the
turning points land one spacing inside each shell face and weld to it. Cells are triangles
bounded by two radial members and a shell: a truss, which is what resists peel.

Three things have to be right before implementing, and none of them need the slicer:
  1. the density law - what N delivers a requested fill fraction;
  2. the cost - path length and cusp count against the ring lattice it replaces;
  3. that the tie profile really is flat, measured the same way as the printed parts.
"""
import math
import numpy as np

SP = 0.3963          # bead spacing / width
R_OUT_FACE = 8.95    # innermost outer-shell ring, from the L2 de-lid part
R_IN_FACE = 3.50     # outermost bore-shell ring
NSAMP = 20000


def zigzag(N, phase=0.0, margin=SP):
    """(theta, r) of one revolution of the serpentine."""
    th = np.linspace(0, 2 * math.pi, NSAMP, endpoint=False)
    lo, hi = R_IN_FACE + margin, R_OUT_FACE - margin
    mid, amp = 0.5 * (lo + hi), 0.5 * (hi - lo)
    x = (N * th + phase) / (2 * math.pi)
    tri = 4.0 * np.abs(x - np.floor(x + 0.5)) - 1.0      # triangle wave in [-1, 1]
    return th, mid + amp * tri


def path_length(th, r):
    x, y = r * np.cos(th), r * np.sin(th)
    return float(np.sum(np.hypot(np.diff(x, append=x[0]), np.diff(y, append=y[0]))))


def ring_lattice_length(k):
    """k concentric rings evenly spread across the same buried band."""
    if k < 1:
        return 0.0
    rs = np.linspace(R_IN_FACE + SP, R_OUT_FACE - SP, k)
    return float(np.sum(2 * math.pi * rs))


def tie_profile(th, r, nbins=40, naz=360):
    """Fraction of azimuths at which each radius carries material within a bead width.

    For the ring lattice this is what collapsed to 0.8% at one radius. For a serpentine every
    traverse crosses every radius, so it should be near-uniform - and if it is not, the design
    is wrong before a line of C++ is written."""
    edges = np.linspace(R_IN_FACE + SP, R_OUT_FACE - SP, nbins + 1)
    az = np.linspace(-math.pi, math.pi, naz, endpoint=False)
    thw = np.mod(th + math.pi, 2 * math.pi) - math.pi
    out = []
    for i in range(nbins):
        rl, rh = edges[i], edges[i + 1]
        hit = (r >= rl - SP / 2) & (r <= rh + SP / 2)
        if not hit.any():
            out.append(0.0); continue
        a = thw[hit]
        # an azimuth is covered if some sample sits within half a bead of it, in arc length
        rmid = 0.5 * (rl + rh)
        tol = (SP / 2) / max(rmid, 1e-6)
        covered = np.zeros(naz, dtype=bool)
        idx = np.searchsorted(np.sort(a), az)
        sa = np.sort(a)
        for j, v in enumerate(az):
            k = min(max(idx[j], 0), len(sa) - 1)
            d = min(abs(sa[k] - v), abs(sa[k - 1] - v) if len(sa) > 1 else 9e9)
            covered[j] = d <= tol
        out.append(covered.mean())
    return np.array(out), edges


T = R_OUT_FACE - R_IN_FACE - 2 * SP
R_MID = 0.5 * (R_OUT_FACE + R_IN_FACE)
print(f"buried band: {R_IN_FACE + SP:.2f} .. {R_OUT_FACE - SP:.2f} mm  "
      f"(thickness {T:.2f} mm, mid radius {R_MID:.2f} mm)")
print(f"annulus area per layer: {2 * math.pi * R_MID * T:.1f} mm^2\n")

print(f"{'N':>4s} {'fill':>7s} {'path':>9s} {'cusps':>6s} {'vs rings':>9s} {'ring k':>7s} "
      f"{'tie min':>8s} {'tie mean':>9s}")
area = 2 * math.pi * R_MID * T
for N in (4, 6, 8, 10, 12, 16, 20, 24):
    th, r = zigzag(N)
    L = path_length(th, r)
    fill = L * SP / area
    # the ring lattice that delivers the SAME fill, for a like-for-like cost comparison
    k = max(1, int(round(fill * area / (SP * 2 * math.pi * R_MID))))
    tp, _ = tie_profile(th, r)
    print(f"{N:4d} {100 * fill:6.1f}% {L:8.1f}mm {2 * N:6d} "
          f"{L / max(ring_lattice_length(k), 1e-9):8.2f}x {k:7d} "
          f"{100 * tp.min():7.1f}% {100 * tp.mean():8.1f}%")

print("\nfill is path area over annulus area; 'vs rings' is serpentine length divided by the")
print("length of the concentric lattice delivering the same fill. tie min is the WORST radius -")
print("the statistic the ring lattice failed at 0.8%.")
