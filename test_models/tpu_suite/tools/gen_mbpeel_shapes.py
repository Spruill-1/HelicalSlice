"""Generate the multi-bore coupons that validate the arc-capable per-bore-cell partition
(commit b9f2959880) on the two axes it introduces:

  PEEL         - the cells are geometric cuts along the nearest-bore equidistant set, so
                 adjacent cells' rings ABUT the cut rather than continuing through it. Print
                 de-lidded (top+bottom shell 0), hand-peel the exposed lattice, and see whether
                 separation follows a cut plane in the arc coupon that the fraction control has
                 no cut to follow.
  DENSIFICATION- the cells each choose their own ring count, so arc lays down more bead than the
                 single-count fraction path (F7 sim: +56% path). Weigh arc vs fraction at matched
                 geometry; the slicer's per-object mm3 pre-screens it.

Sized to the F7-scale that the suite already slices six-up (make_lidplate): ~34 mm slabs, 14 mm
tall. The first cut at 20 mm / 44x28 mm made each coupon heavy enough that two on a bed thrashed
the per-layer foam generation (2xMB3 never finished); at this scale six-up slices in minutes.

Three coupons, each de-lidded, 14 mm tall (a good peel handle; layers 20-27 sit mid-part):

  MB2e_equal2  two EQUAL bores -> ONE clean central cut at x=0. The clean case: any peel that
               follows x=0 is the cut, not geometry.
  MB2u_uneq2   two UNEQUAL bores -> probes the bore-CENTROID Voronoi mispartition. The centroid
               bisector sits at x=2 (midpoint of centres); the true nearest-BOUNDARY equidistant
               plane is at x=3.25. If the 1.25 mm offset matters, the smaller bore's cell is
               starved and this coupon shows it.
  MB3_tri3     three EQUAL bores in a symmetric triangle. Generated but NOT plated: the three
               equal bores put a degenerate medial vertex at the centroid (equidistant to all
               three), where the foam path blows up - two copies never finish export, in the
               shipped fraction path as much as under arc, so it is a pre-existing CS limit, not
               the arc partition. The plate uses F7_multibore (three UNEQUAL bores) for the
               realistic 3-bore case instead. Kept here as the reproducer for that limit.

Run from anywhere; writes binary STLs next to the other suite shapes.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shapely.geometry import box
from shapely.geometry.polygon import Polygon
from gen_foam_shapes import circle, loft

OUT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
H = 14.0


def rrect(w, h, fillet=1.5):
    # A slab, corners filleted a hair so the slicer never resolves a bare 90-degree point.
    return box(-w / 2.0, -h / 2.0, w / 2.0, h / 2.0).buffer(fillet, join_style=1)


def make(poly, name, bores):
    for (cx, cy, r) in bores:
        poly = poly.difference(circle(cx, cy, r))
    if not isinstance(poly, Polygon):
        raise SystemExit(f"{name}: bores split the slab into a MultiPolygon - move them apart")
    if len(poly.interiors) != len(bores):
        raise SystemExit(f"{name}: {len(poly.interiors)} holes for {len(bores)} bores - a bore "
                         f"reaches the wall or another bore")
    path = os.path.join(OUT, name + ".stl")
    n = loft([(0.0, poly), (H, poly)], path)
    print(f"{name}: {n} tris, {len(poly.interiors)} bores, bbox {tuple(round(v,1) for v in poly.bounds)} -> {path}")


if __name__ == "__main__":
    make(rrect(34, 18), "MB2e_equal2", [(-8, 0, 3.5), (8, 0, 3.5)])
    make(rrect(34, 18), "MB2u_uneq2",  [(-6, 0, 4.5), (9, 0, 2.0)])
    make(rrect(34, 22), "MB3_tri3",    [(-8, -5, 3.0), (8, -5, 3.0), (0, 6, 3.0)])
