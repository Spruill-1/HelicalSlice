"""Tie fraction as a function of RADIUS - the profile a peel actually follows.

A peel does not care what the average tie fraction is. It follows the weakest continuous
surface, and a surface is weak when its tie LAPSES, however often it also holds. So the thing
to look at is the tie fraction of each interface in radial order, and in particular where it
STEPS: a solid block tied 100% sitting against a lattice tied 40% has a fracture plane between
them wherever that step is, regardless of which ring is pinned.

Rings are matched by rank on rays showing the full set, which the weld test validated at 715/720.
"""
import sys
import numpy as np
from tiechain import segs
from interface import crossings_z

SP = 0.3963

def rank_rows(path, obj, lh, zlo=4.0, zhi=10.0, naz=360):
    _, objs = segs(path)
    for name, sg in objs:
        if obj not in name: continue
        zs = sorted({int(round(s[4]/lh)) for s in sg if zlo <= s[4] <= zhi})
        rows = []
        for z in zs:
            mid = [s for s in sg if int(round(s[4]/lh)) == z]
            cx = sum(s[0] for s in mid)/len(mid); cy = sum(s[1] for s in mid)/len(mid)
            per = []
            for c in crossings_z(mid, cx, cy):
                m = []
                for r, w, zz in c:
                    if m and r - m[-1][0] < 0.12: continue
                    m.append((r, w, zz))
                per.append(m)
            rows.extend(per)
        # One ring count for the whole band: the modal one. Per-layer maxima gave rows of
        # different widths and made the array ragged.
        from collections import Counter
        n = Counter(len(p) for p in rows).most_common(1)[0][0]
        keep = [p for p in rows if len(p) == n]
        return (np.array([[q[0] for q in p] for p in keep]),
                np.array([[q[2] for q in p] for p in keep]), lh)
    return None, None, lh

if __name__ == "__main__":
  for f, tag in ((sys.argv[1], sys.argv[2]),):
      lh, _ = segs(f)[0], None
      lh = segs(f)[0]
      for obj in ("L2_open50", "L1_open25"):
          R, Z, lh = rank_rows(f, obj, lh)
          if R is None: continue
          print(f"\n=== {tag}  {obj}   {R.shape[0]} rays x {R.shape[1]} rings ===")
          print(f"  {'iface':>5s} {'mean radius':>11s} {'tied':>7s} {'tied+flat':>10s} {'med dZ':>9s}")
          for i in range(R.shape[1]-1):
              gap = R[:, i+1] - R[:, i]
              dz  = np.abs(Z[:, i+1] - Z[:, i])
              tied = gap <= 1.15*SP
              print(f"  {i:5d} {R[:, i:i+2].mean():10.2f}mm {100*tied.mean():6.1f}% "
                    f"{100*np.mean(tied & (dz <= 0.35*lh)):9.1f}% {np.median(dz)/lh:8.2f}L")
  