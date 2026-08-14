"""Generate the TPU continuous-shell test tube suite (binary STLs).

Parametrized on the print flow: nozzle 0.4, line width w=0.45, layer height h=0.25.
spacing s = w - h*(1 - pi/4). Wall thicknesses are tuned so each tube exercises one
specific continuous-shell behavior. Regenerate with different W/H if the profile changes.
"""
import math
import struct
import os

W = 0.45           # line width (mm)
H = 0.25           # layer height (mm)
S = W - H * (1 - math.pi / 4)          # bead spacing
DENSITY = 1.21     # TPU g/cm^3
SEG = 192
OUT = r"C:\Users\david\source\HelicalSlice\test_models\tpu_suite"

def t_exact(M):    # wall thickness that tiles exactly with M band pairs (2M rings)
    return 2 * W + 2 * (M - 1) * S

def t_deficit(M, v):  # wall leaving residual void v between innermost bead edges
    return 2 * (W / 2 + (M - 1) * S + (v + W) / 2)

def tube_tris(r_out, r_in_bot, r_in_top, h):
    """Watertight annular tube; inner radius may differ bottom->top (taper)."""
    tris = []
    for i in range(SEG):
        a0 = 2 * math.pi * i / SEG
        a1 = 2 * math.pi * (i + 1) / SEG
        c0, s0, c1, s1 = math.cos(a0), math.sin(a0), math.cos(a1), math.sin(a1)
        ob0 = (r_out * c0, r_out * s0, 0); ob1 = (r_out * c1, r_out * s1, 0)
        ot0 = (r_out * c0, r_out * s0, h); ot1 = (r_out * c1, r_out * s1, h)
        ib0 = (r_in_bot * c0, r_in_bot * s0, 0); ib1 = (r_in_bot * c1, r_in_bot * s1, 0)
        it0 = (r_in_top * c0, r_in_top * s0, h); it1 = (r_in_top * c1, r_in_top * s1, h)
        tris += [(ob0, ob1, ot1), (ob0, ot1, ot0)]      # outer wall (out-facing CCW)
        tris += [(ib0, it1, ib1), (ib0, it0, it1)]      # inner wall (in-facing)
        tris += [(ob0, ib1, ob1), (ob0, ib0, ib1)]      # bottom (down-facing)
        tris += [(ot0, ot1, it1), (ot0, it1, it0)]      # top (up-facing)
    return tris

def write_stl(path, tris):
    with open(path, "wb") as f:
        f.write(b"\0" * 80)
        f.write(struct.pack("<I", len(tris)))
        for a, b, c in tris:
            u = tuple(b[k] - a[k] for k in range(3))
            v = tuple(c[k] - a[k] for k in range(3))
            n = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
            ln = math.sqrt(sum(x*x for x in n)) or 1.0
            f.write(struct.pack("<3f", *(x/ln for x in n)))
            for p in (a, b, c):
                f.write(struct.pack("<3f", *p))
            f.write(struct.pack("<H", 0))

def mesh_volume(tris):  # divergence theorem, mm^3
    vol = 0.0
    for a, b, c in tris:
        vol += (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0])) / 6.0
    return abs(vol)

os.makedirs(OUT, exist_ok=True)
print(f"flow: w={W} h={H} -> spacing s={S:.4f}")
print(f"t_exact(2)={t_exact(2):.4f}  t_exact(3)={t_exact(3):.4f}  "
      f"t_deficit(2, 0.30)={t_deficit(2, 0.30):.4f}  v_max={2*S-W:.4f}")

suite = [
    # name, OD, wall_bot, wall_top, height, purpose
    ("T1_tile_exact",  32.0, t_exact(2),          t_exact(2),          30.0, "control: 4 rings tile exactly, no medial correction fires"),
    ("T2_deficit",     32.0, t_deficit(2, 0.30),  t_deficit(2, 0.30),  30.0, "medial deficit v=0.30: innermost pair widened to ~0.60"),
    ("T3_thin",        32.0, 0.35,                0.35,                30.0, "sub-bead wall: single fitted 0.35 bead"),
    ("T4_taper",       32.0, 2.60,                1.50,                40.0, "wall tapers across ring-count transitions"),
    ("T5_seam",        32.0, t_exact(3),          t_exact(3),          60.0, "6-ring clean wall for scarf A/B: print at scarf 0/6/12"),
]

for name, od, wb, wt, h, purpose in suite:
    r_out = od / 2
    tris = tube_tris(r_out, r_out - wb, r_out - wt, h)
    path = os.path.join(OUT, name + ".stl")
    write_stl(path, tris)
    vol = mesh_volume(tris)
    print(f"{name:15s} OD{od:.0f} wall {wb:.3f}->{wt:.3f} H{h:.0f}  "
          f"vol {vol/1000:.3f} cm^3  weight {vol/1000*DENSITY:.2f} g   ({purpose})")
