"""Generate the shapes the interpolated-ring foam has not been looked at on (binary STLs).

    gen_foam_shapes.py

Everything the foam has been tuned on is a round annulus with one bore. Each of these breaks
one assumption the ring generator leans on:

  F1_crescent      no bore at all, so the level sets close onto a MEDIAL AXIS rather than a
                   bore, and the section is deeply non-convex with two thin tips
  F2_square        four sharp corners: the distance field's gradient turns through 90 degrees
                   at each, and the level sets have to stay parallel through that
  F3_crescent_bore a bore inside the thick limb of a crescent, so the field closes on the bore
                   near it and on the medial axis out at the tips, in ONE region
  F4_square_bore   corners outside, a circle inside, so the rings have to morph from square to
                   round across the wall
  F5_hex_twist     a hexagon rotating 90 degrees over its height: the same section at every
                   layer but at a different angle, so anything quietly tied to a fixed azimuth
                   shows up as a drift
  F6_hex_twist_bore  the same, with a bore, so the wall thickness pattern rotates with it
  F7_multibore     three bores of different sizes, which is what the nearest-bore partition was
                   built for: the cut between two bores is their equidistant set, so u is
                   continuous across it and the rings should JOIN there rather than terminate
  F8_cone_trunc    a truncated cone: the section AREA shrinks with height, so the ring count
                   has to fall as the part narrows and the lattice has to re-space itself
                   rather than keep a count chosen at the bottom
  F9_pyramid_trunc the same taper with corners, so the count changes and the corners move at
                   once - the two things F2 and F8 each test alone

12 mm tall - these are for inspection rather than the plate, so they slice quickly and layers
20-27 sit in the middle of the part instead of in its top or bottom shell.
"""
import math, os, struct
import numpy as np
from shapely.geometry import Polygon, Point, MultiPoint
from shapely import delaunay_triangles

TALL, NZ = 12.0, 24            # 24 lofted steps is plenty: only the twist needs them
SEG      = 192
OUT      = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def circle(cx, cy, r, seg=SEG):
    a = np.linspace(0, 2 * math.pi, seg, endpoint=False)
    return Polygon(np.c_[cx + r * np.cos(a), cy + r * np.sin(a)])


def ngon(n, r, rot=0.0, fillet=0.0):
    a = np.linspace(0, 2 * math.pi, n, endpoint=False) + rot
    p = Polygon(np.c_[r * np.cos(a), r * np.sin(a)])
    # A fillet keeps the corner from being a single point the slicer has to resolve exactly;
    # 0.4 mm is well under a bead, so the corner is still a corner.
    return p.buffer(fillet, join_style=1).buffer(-fillet, join_style=1) if fillet else p


def rings_of(poly):
    """exterior first, then holes - each as a closed list of (x, y), no repeated last point"""
    out = [list(poly.exterior.coords)[:-1]]
    out += [list(h.coords)[:-1] for h in poly.interiors]
    return out


def cap_tris(poly, z, up):
    """Triangulate a polygon WITH HOLES for the end caps.

    Delaunay over the boundary vertices plus an interior grid, keeping only triangles whose
    centre lies in the polygon. The grid is what makes it work on a crescent: a Delaunay of
    boundary points alone spans the concave side, and the centre test then throws those away
    and leaves a hole in the cap.
    """
    pts = [p for ring in rings_of(poly) for p in ring]
    x0, y0, x1, y1 = poly.bounds
    step = 0.8
    gx, gy = np.meshgrid(np.arange(x0, x1, step), np.arange(y0, y1, step))
    for px, py in zip(gx.ravel(), gy.ravel()):
        if poly.contains(Point(px, py)) and poly.exterior.distance(Point(px, py)) > step * 0.5:
            pts.append((px, py))
    tris = []
    for t in delaunay_triangles(MultiPoint(pts)).geoms:
        if not poly.contains(t.representative_point()):
            continue
        c = list(t.exterior.coords)[:3]
        q = [(p[0], p[1], z) for p in c]
        tris.append(tuple(q) if up else (q[0], q[2], q[1]))
    return tris


def loft(sections, path):
    """Watertight mesh through a stack of (z, polygon). Every section needs the same rings
    with the same vertex counts, which holds here because a twist only rotates them."""
    tris = []
    for (z0, p0), (z1, p1) in zip(sections, sections[1:]):
        for ri, (r0, r1) in enumerate(zip(rings_of(p0), rings_of(p1))):
            outward = (ri == 0)
            for i in range(len(r0)):
                a, b = r0[i], r0[(i + 1) % len(r0)]
                c, d = r1[i], r1[(i + 1) % len(r1)]
                q = [(a[0], a[1], z0), (b[0], b[1], z0), (d[0], d[1], z1), (c[0], c[1], z1)]
                if outward: tris += [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
                else:       tris += [(q[0], q[2], q[1]), (q[0], q[3], q[2])]
    tris += cap_tris(sections[0][1], sections[0][0], up=False)
    tris += cap_tris(sections[-1][1], sections[-1][0], up=True)
    with open(path, "wb") as f:
        f.write(b"\0" * 80)
        f.write(struct.pack("<I", len(tris)))
        for t in tris:
            f.write(struct.pack("<3f", 0.0, 0.0, 0.0))
            for p in t:
                f.write(struct.pack("<3f", *p))
            f.write(struct.pack("<H", 0))
    return len(tris)


def straight(poly, name):
    return loft([(0.0, poly), (TALL, poly)], os.path.join(OUT, name + ".stl")), name


def twisted(fn, name, turn=math.pi / 2):
    zs = np.linspace(0, TALL, NZ + 1)
    return loft([(float(z), fn(turn * z / TALL)) for z in zs],
                os.path.join(OUT, name + ".stl")), name


def tapered(fn, name):
    """fn(t) for t in 0..1 up the height - the section may change SIZE, not just angle."""
    zs = np.linspace(0, TALL, NZ + 1)
    return loft([(float(z), fn(z / TALL)) for z in zs],
                os.path.join(OUT, name + ".stl")), name


def main():
    print(__doc__)
    made = []
    # 1. a solid crescent: no bore, so the rings close onto the medial axis
    cres = circle(0, 0, 12).difference(circle(5, 0, 10))
    made.append(straight(cres, "F1_crescent"))
    # 2. a solid square
    made.append(straight(ngon(4, 14.1, math.pi / 4, fillet=0.4), "F2_square"))
    # 3. the crescent, bored through its thick limb
    made.append(straight(cres.difference(circle(-8.5, 0, 1.6)), "F3_crescent_bore"))
    # 4. square outside, circle inside
    made.append(straight(ngon(4, 14.1, math.pi / 4, fillet=0.4).difference(circle(0, 0, 5)),
                         "F4_square_bore"))
    # 5. a hexagon that turns 90 degrees on the way up
    made.append(twisted(lambda t: ngon(6, 12, t, fillet=0.4), "F5_hex_twist"))
    # 6. the same, bored, so the wall thickness pattern turns with it
    made.append(twisted(lambda t: ngon(6, 12, t, fillet=0.4).difference(circle(0, 0, 5)),
                        "F6_hex_twist_bore"))
    # 7. three unequal bores for the nearest-bore partition
    slab = ngon(4, 1, 0).buffer(0)          # placeholder, replaced below
    slab = Polygon([(-14, -9), (14, -9), (14, 9), (-14, 9)]).buffer(1.5, join_style=1)
    multi = slab.difference(circle(-8, 0, 3)).difference(circle(2, 3, 2)) \
                .difference(circle(8, -3, 2.5))
    made.append(straight(multi, "F7_multibore"))
    # 8, 9. the section AREA changes with height. Tapering INWARD, so no layer overhangs the
    # one below it and what is being tested is the ring count following the section rather
    # than the printer's bridging.
    made.append(tapered(lambda t: circle(0, 0, 13.0 - 6.5 * t), "F8_cone_trunc"))
    made.append(tapered(lambda t: ngon(4, (15.6 - 7.8 * t), math.pi / 4, fillet=0.4),
                        "F9_pyramid_trunc"))
    for n, name in made:
        p = os.path.join(OUT, name + ".stl")
        print(f"  {name:<20s} {n:>6d} triangles  {os.path.getsize(p)/1024:>7.1f} kB")


if __name__ == "__main__":
    main()
