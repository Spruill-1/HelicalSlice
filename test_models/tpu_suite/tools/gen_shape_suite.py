"""Generate non-tube shapes for continuous-shell validation (binary STLs).

Everything in the campaign so far was validated on ONE circular tube, which is the
friendliest possible case for the ramped interior arm: constant radius, constant wall,
one bore, no corners, no solid sections. These shapes each break one of those.

  S1_ellipse   radius varies 20 -> 13 around the section (per-angle plan tables, and
               a "revolution" that is no longer a constant radius)
  S2_rsquare   superellipse: four corners, strong curvature variation (the case that
               historically under-counted max beads on square annuli)
  S3_lobed     three-lobed outer wall - still star-shaped, but concave-ish flanks
  S4_capped    circular tube that closes with a 45-degree CONE, so the bore shrinks
               continuously to nothing and the section walks from a 2.5 mm annulus all
               the way to solid - the interior spiral must track that and then decline
               and hand back to ring mode. (v1 closed the bore with a flat disc, which
               made the first solid layer a 27 mm unsupported bridge printed into open
               air; it drooped into the bore and told us nothing about the transition.
               A cone keeps every layer within a 45-degree overhang and tests the
               annulus-to-solid handover instead of the printer's bridging.)

Wall thickness matches T5 (2.485 mm, ~6 beads) so results compare directly, and every
shape is 30 mm tall to stay under the H2D's by-object gantry clearance.
"""
import math, os, struct

W, H = 0.45, 0.25
S    = W - H * (1 - math.pi / 4)
WALL = 2.485                       # same as T5
TALL = 30.0
SEG  = 256
OUT  = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def offset_inward(pts, t):
    """Offset a sampled closed curve inward by t along per-vertex normals."""
    n = len(pts)
    out = []
    for i in range(n):
        x0, y0 = pts[(i - 1) % n]
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        nx, ny = 0.0, 0.0
        for (ax, ay), (bx, by) in (((x0, y0), (x1, y1)), ((x1, y1), (x2, y2))):
            ex, ey = bx - ax, by - ay
            L = math.hypot(ex, ey) or 1.0
            nx += ey / L; ny += -ex / L        # outward normal for CCW input
        L = math.hypot(nx, ny) or 1.0
        out.append((x1 - t * nx / L, y1 - t * ny / L))
    return out

def prism(outer, inner, h, cap_z=None):
    """Watertight prism between two closed loops. cap_z: solidify above that height."""
    tris = []
    def wall(loop, z0, z1, outward):
        for i in range(len(loop)):
            a = loop[i]; b = loop[(i + 1) % len(loop)]
            p = [(a[0], a[1], z0), (b[0], b[1], z0), (b[0], b[1], z1), (a[0], a[1], z1)]
            if outward: tris.extend([(p[0], p[1], p[2]), (p[0], p[2], p[3])])
            else:       tris.extend([(p[0], p[2], p[1]), (p[0], p[3], p[2])])
    def ring_face(lo, li, z, up):
        for i in range(len(lo)):
            a = lo[i]; b = lo[(i + 1) % len(lo)]
            c = li[i]; d = li[(i + 1) % len(li)]
            q = [(a[0],a[1],z),(b[0],b[1],z),(d[0],d[1],z),(c[0],c[1],z)]
            if up: tris.extend([(q[0], q[1], q[2]), (q[0], q[2], q[3])])
            else:  tris.extend([(q[0], q[2], q[1]), (q[0], q[3], q[2])])
    def disc(loop, z, up):
        cx = sum(p[0] for p in loop) / len(loop)
        cy = sum(p[1] for p in loop) / len(loop)
        for i in range(len(loop)):
            a = loop[i]; b = loop[(i + 1) % len(loop)]
            t = ((cx, cy, z), (a[0], a[1], z), (b[0], b[1], z))
            tris.append(t if up else (t[0], t[2], t[1]))
    if cap_z is None:
        wall(outer, 0, h, True); wall(inner, 0, h, False)
        ring_face(outer, inner, 0, False); ring_face(outer, inner, h, True)
    else:
        wall(outer, 0, h, True)                    # full-height outer wall
        wall(inner, 0, cap_z, False)               # bore only up to the cap
        ring_face(outer, inner, 0, False)          # bottom annulus
        disc(inner, cap_z, True)                   # cap floor closes the bore
        disc(outer, h, True)                       # solid top
    return tris

def write_stl(path, tris):
    with open(path, "wb") as f:
        f.write(b"\0" * 80); f.write(struct.pack("<I", len(tris)))
        for a, b, c in tris:
            u = tuple(b[k] - a[k] for k in range(3)); v = tuple(c[k] - a[k] for k in range(3))
            n = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
            ln = math.hypot(math.hypot(n[0], n[1]), n[2]) or 1.0
            f.write(struct.pack("<3f", *(x/ln for x in n)))
            for p in (a, b, c): f.write(struct.pack("<3f", *p))
            f.write(struct.pack("<H", 0))

def revolve_coned(r_out, wall, h, cyl_top):
    """Cylinder with a bore that closes as a 45-degree cone (self-supporting).

    Solid boundary = outer cylinder + bottom annulus + inner surface (cylinder then
    cone to an apex) + top disc. The apex sits inside the solid, so everything above it
    is a full disc - the annulus-to-solid transition happens gradually, one layer's
    worth of bore at a time, instead of as an unsupported bridge.
    """
    r_in  = r_out - wall
    apex  = cyl_top + r_in                      # 45 degrees: radius sheds 1 mm per mm
    assert apex < h, f"cone apex {apex} must sit below the top face {h}"
    tris  = []
    def ring(r, z):
        return [(r * math.cos(2*math.pi*i/SEG), r * math.sin(2*math.pi*i/SEG), z) for i in range(SEG)]
    o0, o1 = ring(r_out, 0.0), ring(r_out, h)
    i0, i1 = ring(r_in, 0.0), ring(r_in, cyl_top)
    for i in range(SEG):
        j = (i + 1) % SEG
        tris += [(o0[i], o0[j], o1[j]), (o0[i], o1[j], o1[i])]      # outer wall, outward
        tris += [(i0[i], i1[j], i0[j]), (i0[i], i1[i], i1[j])]      # bore wall, inward
        tris += [(o0[i], i0[j], o0[j]), (o0[i], i0[i], i0[j])]      # bottom annulus, down
        tris += [(o1[i], o1[j], (0.0, 0.0, h))]                     # top disc, up
        tris += [(i1[i], (0.0, 0.0, apex), i1[j])]                  # cone to apex, inward
    return tris

def curve(fn):
    return [fn(2 * math.pi * i / SEG) for i in range(SEG)]

shapes = {
    # ellipse: 20 x 13 semi-axes
    "S1_ellipse": (curve(lambda a: (20.0 * math.cos(a), 13.0 * math.sin(a))), None),
    # superellipse n=4 -> rounded square, half-width 16
    "S2_rsquare": (curve(lambda a: (16.0 * math.copysign(abs(math.cos(a)) ** 0.5, math.cos(a)),
                                    16.0 * math.copysign(abs(math.sin(a)) ** 0.5, math.sin(a)))), None),
    # three lobes, amplitude 2 on radius 15 (min curvature stays well above the wall)
    "S3_lobed":   (curve(lambda a: ((15.0 + 2.0 * math.cos(3 * a)) * math.cos(a),
                                    (15.0 + 2.0 * math.cos(3 * a)) * math.sin(a))), None),

}

for name, (outer, cap_z) in shapes.items():
    inner = offset_inward(outer, WALL)
    tris  = prism(outer, inner, TALL, cap_z)
    path  = os.path.join(OUT, name + ".stl")
    write_stl(path, tris)
    cx = sum(p[0] for p in outer) / len(outer); cy = sum(p[1] for p in outer) / len(outer)
    rr = [math.hypot(x - cx, y - cy) for (x, y) in outer]
    print(f"{name:12s} {len(tris):5d} tris  outer r {min(rr):.1f}-{max(rr):.1f} mm  wall {WALL} mm  open tube")

# S4: circular tube whose bore closes as a 45-degree cone (see the module docstring)
CYL_TOP = 13.0
tris = revolve_coned(16.0, WALL, TALL, CYL_TOP)
write_stl(os.path.join(OUT, "S4_capped.stl"), tris)
print(f"{'S4_capped':12s} {len(tris):5d} tris  outer r 16.0 mm  wall {WALL} mm  "
      f"bore closes 45deg from z={CYL_TOP} to apex z={CYL_TOP + 16.0 - WALL:.1f}, solid above")
