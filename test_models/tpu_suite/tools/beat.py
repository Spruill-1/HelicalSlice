"""Per-ring lobe COUNT, chosen from each ring's own circumference - the beat idea.

Every attempt so far has kept one angular lobe count for the whole stack and varied phase. That
locks the whole wall into one pattern: with anti-phase, rings j and j+2 are in phase, so at any
azimuth alternate interfaces close and the ones between them open, and which set closes is fixed
by the amplitudes rather than by the azimuth. Measured on the printed part, the result is
permanent doublets - tie fractions of 74/1/15/63/14/44/21 across the interfaces, and the near-zero
ones are free cylinders a peel walks straight down.

If instead each ring gets its OWN lobe count, adjacent rings BEAT against each other: their
relative phase sweeps through a full turn |dL| times per revolution, so they pass through
opposition |dL| times whatever their amplitudes are. Contact stops depending on amplitude
symmetry, which is what kept failing, and becomes a property of the frequencies. The count also
stops being an angular constant and starts being what it should have been all along - a
WAVELENGTH along the ring, so a bore ring and an outer ring get the same lobe pitch in mm rather
than the same number of lobes across very different circumferences.

    tie zones per interface = |L_j - L_{j+1}|   per revolution, by construction

This scores the three schemes on the de-lid L2 geometry against the statistic that matches the
hand: the WORST interface, since a peel finds the weakest surface.
"""
import math
import numpy as np

SP = 0.3963
RINGS = [2.70, 3.10, 3.50, 4.30, 5.10, 5.90, 6.55, 7.35, 8.15, 8.95, 9.36, 9.75]
SHELL = {0, 1, 2, 9, 10, 11}            # rings that do not weave
NAZ = 3600
OVERSHOOT = 4.0


def caps(rings):
    """Per-ring inward/outward budget: half the excess over one spacing against a mobile
    neighbour, all of it against a static one - the law already in the slicer."""
    out = []
    for j, r in enumerate(rings):
        co = ci = 0.0
        if j + 1 < len(rings):
            g = rings[j + 1] - r
            co = max(0.0, (0.5 if (j + 1) not in SHELL else 1.0) * (g - SP))
        if j > 0:
            g = r - rings[j - 1]
            ci = max(0.0, (0.5 if (j - 1) not in SHELL else 1.0) * (g - SP))
        out.append((co, ci))
    return out


def profile(lobes, phases):
    th = np.linspace(0, 2 * math.pi, NAZ, endpoint=False)
    cp = caps(RINGS)
    r = []
    for j, r0 in enumerate(RINGS):
        if j in SHELL:
            r.append(np.full(NAZ, r0)); continue
        e = np.clip(OVERSHOOT * np.sin(lobes[j] * th + phases[j]), -1, 1)
        co, ci = cp[j]
        r.append(r0 + np.where(e >= 0, e * co, e * ci))
    r = np.array(r)
    ties = []
    for i in range(len(RINGS) - 1):
        gap = r[i + 1] - r[i]
        ties.append(100.0 * np.mean(gap <= 1.15 * SP))
    return ties


def show(name, ties):
    inner = ties[2:-2]
    print(f"{name:>28s} " + " ".join(f"{t:4.0f}" for t in ties) +
          f"   worst {min(inner):5.1f}%  mean {np.mean(inner):5.1f}%")


n = len(RINGS)
print(f"{'scheme':>28s} {'tie % per interface (bore -> out)':>50s}")

# (a) what ships today: one angular lobe count, anti-phase plus a drift
L = [6.0] * n
drift = 2.0 * 2.0 * math.asin(1.0 / OVERSHOOT)
show("today: L=6, anti-phase+drift", profile(L, [j * (math.pi + drift) for j in range(n)]))
show("today: L=6, anti-phase only", profile(L, [j * math.pi for j in range(n)]))

# (b) the beat: a lobe count from each ring's own circumference, at a target wavelength
for wl in (4.0, 6.0, 8.0, 10.0):
    Lb = [max(2.0, round(2 * math.pi * r / wl)) for r in RINGS]
    show(f"beat: wavelength {wl:.0f}mm {Lb[3:9]}", profile(Lb, [0.0] * n))

# (c) the beat with a deliberate step, so neighbours always differ by dL
for dl in (1, 2, 3, 4):
    Lb = [6.0 + dl * j for j in range(n)]
    show(f"beat: forced step dL={dl}", profile(Lb, [0.0] * n))

print("\nA pair with lobe counts differing by dL passes through opposition dL times per")
print("revolution whatever its amplitudes are - which is the property every amplitude-based")
print("attempt failed to deliver.")
