"""Prismatic void: void at the same (x,y) on EVERY layer of the part.

The property gyroid has and a concentric lattice may not - no cross-section that, extruded
vertically, stays entirely in empty space. The part's own bores are excluded: a through-hole is
not a fill defect, and counting it is a mistake this instrument has made twice.
"""
import re, sys, zipfile
import numpy as np
from scipy import ndimage
sys.path.insert(0, ".")
import voids
from voids import PIX

def run(path, want, skip_bottom=4, label=None):
    g = zipfile.ZipFile(path).read("Metadata/plate_1.gcode").decode("utf-8", "replace")
    lay = -1; seen = set()
    for line in g.splitlines():
        if line.startswith("; CHANGE_LAYER"): lay += 1
        elif line.startswith("; printing object") and want in line: seen.add(lay)
    lo, hi = min(seen), max(seen)
    s0 = voids.layer_segments(path, want, lo + skip_bottom)
    xs=[s[0] for s in s0]+[s[2] for s in s0]; ys=[s[1] for s in s0]+[s[3] for s in s0]
    x0,x1,y0,y1=min(xs)-1,max(xs)+1,min(ys)-1,max(ys)+1
    W=int((x1-x0)/PIX)+1; H=int((y1-y0)/PIX)+1
    yy,xx=np.mgrid[0:H,0:W]; px=x0+xx*PIX; py=y0+yy*PIX
    stack=None; anymat=None; n=0
    for L in range(lo+skip_bottom, hi+1):
        sg=voids.layer_segments(path,want,L)
        if not sg: continue
        mat=np.zeros((H,W),dtype=bool)
        for ax,ay,bx,by,w in sg:
            dx,dy=bx-ax,by-ay; L2=dx*dx+dy*dy
            t=np.clip(((px-ax)*dx+(py-ay)*dy)/L2,0,1)
            mat |= (px-(ax+t*dx))**2+(py-(ay+t*dy))**2 <= (0.5*w)**2
        o=np.zeros_like(mat); o[0,:]=o[-1,:]=o[:,0]=o[:,-1]=True
        o=ndimage.binary_propagation(o & ~mat, mask=~mat)
        v=(~mat)&(~o)
        stack = v if stack is None else (stack & v)
        anymat = mat if anymat is None else (anymat | mat)
        n += 1
    # a BORE is void on every layer AND never has material in it on any layer
    lab,k=ndimage.label(stack)
    sizes=np.bincount(lab.ravel())[1:]*PIX*PIX if k else np.array([])
    d=ndimage.distance_transform_edt(stack,sampling=PIX)
    keep=[]
    for i in range(len(sizes)):
        comp = lab == i+1
        if (anymat & comp).sum() == 0 and 2*d[comp].max() > 1.5:
            continue                       # the part's own hole: never any material there
        keep.append(i)
    tot = sum(sizes[i] for i in keep)
    keep.sort(key=lambda i: -sizes[i])
    print(f"  {label or want}: {n} layers, prismatic void {tot:.2f} mm2 in {len(keep)} pieces")
    for i in keep[:4]:
        if sizes[i] < 0.01: break
        print(f"      {sizes[i]:7.2f} mm2, {2*d[lab==i+1].max():.2f} mm across")
    return stack, lab, keep, sizes, d

if __name__ == "__main__":
    for path, want, lbl in [tuple(a.split(":")) for a in sys.argv[1:]]:
        run(path, want, label=lbl)
