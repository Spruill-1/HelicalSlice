"""Render a Z-range of a sliced plate as an inspectable toolpath viewer.

    viewgcode.py <file.gcode.3mf> [--z LO HI | --layers A B] [--object NAME]... [--out FILE]

The campaign reads bead width through printf tables, which is fine for a total and useless for
"where". This emits a self-contained HTML page - no network, no libraries - drawing the actual
moves for a slice of Z, coloured by the quantity everything here turns on: the extrusion width
against the narrowest bead the nozzle can lay. Anything under that floor is what "collapse"
means, and on this page it is simply the red part of the picture.

Defaults to six layers from the middle of the part, because the whole plate is far more data
than a browser wants and six layers is enough to see the weave precess.

    viewgcode.py sliced/PlateR_interp.gcode.3mf
    viewgcode.py _verify/sweep_cells=6.gcode.3mf --z 4.0 5.5 --object P3_ecc
"""
import argparse, json, math, os, re, sys, zipfile


def parse(path, want_objects=None):
    with zipfile.ZipFile(path) as z:
        name = next(n for n in z.namelist() if n.endswith("plate_1.gcode"))
        g = z.read(name).decode("utf-8", "replace")

    def head(key, default=None):
        m = re.search(r"; " + key + r" = ([^\n]+)", g)
        return m.group(1).strip() if m else default

    lh = float(head("layer_height", "0.2"))
    fd = float(head("filament_diameter", "1.75").split(",")[0])
    nozzle = float(head("nozzle_diameter", "0.4").split(",")[0])
    nominal = float(head("line_width", "0.45").split(",")[0])
    area = math.pi * fd * fd / 4.0

    runs, travels, objs = [], [], []
    cur = None; x = y = zc = None; layer = -1
    run = None
    for line in g.splitlines():
        if line.startswith("; CHANGE_LAYER"):
            layer += 1; run = None; continue
        m = re.match(r"; printing object (\S+)", line)
        if m:
            cur = m.group(1).split(".stl")[0]
            if cur not in objs: objs.append(cur)
            run = None; continue
        if line.startswith("; stop printing object"):
            cur = None; run = None; continue
        if not (line.startswith("G1 ") or line.startswith("G0 ")):
            continue
        nx = re.search(r" X(-?[\d.]+)", line); ny = re.search(r" Y(-?[\d.]+)", line)
        nz = re.search(r" Z(-?[\d.]+)", line); ne = re.search(r" E(-?[\d.]+)", line)
        px, py = x, y
        if nx: x = float(nx.group(1))
        if ny: y = float(ny.group(1))
        if nz: zc = float(nz.group(1))
        if cur is None or px is None or zc is None:
            run = None; continue
        d = math.hypot(x - px, y - py)
        if d <= 1e-9:
            continue
        if ne and float(ne.group(1)) > 0:
            w = float(ne.group(1)) * area / (d * lh)
            if run is None:
                run = {"o": cur, "l": layer, "z": zc, "p": [round(px, 3), round(py, 3)],
                       "w": [round(w, 4)]}
                runs.append(run)
            run["p"] += [round(x, 3), round(y, 3)]
            run["w"].append(round(w, 4))
            run["z"] = max(run["z"], zc)
        else:
            run = None
            travels.append({"o": cur, "l": layer, "_z": zc,
                            "p": [round(px, 3), round(py, 3), round(x, 3), round(y, 3)]})
    # The plate emits `; CHANGE_LAYER` once per object per layer, so counting them gives an
    # OBJECT-layer index - on a three-object plate that made "layers 20..25" mean two layers of
    # three objects rather than six layers of one. Each object's own layer number is simply the
    # order its runs appear in.
    seen = {}
    for r in runs:
        r["l"] = seen.get(r["o"], 0)
        seen[r["o"]] = r["l"] + 1
    # A travel belongs to the layer of the run nearest it in Z, which for a plate that ramps Z
    # across every layer is just the closest run of the same object.
    zof = {}
    for r in runs:
        zof.setdefault(r["o"], []).append((r["z"], r["l"]))
    for t in travels:
        cand = zof.get(t["o"])
        tz = t.pop("_z", None)
        t["l"] = min(cand, key=lambda zl: abs(zl[0] - tz))[1] if (cand and tz is not None) else 0

    if want_objects:
        runs = [r for r in runs if r["o"] in want_objects]
        travels = [t for t in travels if t["o"] in want_objects]
        objs = [o for o in objs if o in want_objects]
    return {"layer_height": lh, "nozzle": nozzle, "nominal": nominal,
            "floor": round(max(0.30, 0.85 * nozzle), 3),
            "objects": objs, "runs": runs, "travels": travels}


def select(data, zlo, zhi, lo_layer, hi_layer):
    if lo_layer is not None:
        keep = lambda r: lo_layer <= r["l"] <= hi_layer
    elif zlo is not None:
        keep = lambda r: zlo <= r["z"] <= zhi
    else:
        ls = sorted({r["l"] for r in data["runs"]})
        if not ls:
            return data
        mid = ls[len(ls) // 2]
        a, b = max(ls[0], mid - 3), min(ls[-1], mid + 2)
        keep = lambda r: a <= r["l"] <= b
    data["runs"] = [r for r in data["runs"] if keep(r)]
    lset = {r["l"] for r in data["runs"]}
    data["travels"] = [t for t in data["travels"] if t["l"] in lset]
    return data


HTML = r"""<title>Bead Width Inspector</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600;700&display=swap">
<style>
/* A single committed look: a dark instrument ground. The page's whole job is a colour ramp
   over toolpaths, and a ramp reads against one consistent ground or it does not read. Every
   colour is painted explicitly so the page holds whatever the host is set to. */
:root{
  --ground:#0D1317; --panel:#141C22; --panel-2:#1A242B; --line:#26333C;
  --ink:#E6EDF2; --ink-2:#9DAFBC; --ink-3:#6B7E8C;
  --accent:#48C7E8; --accent-dim:#1E5C6E;
  --hot:#FF6A4D; --warn:#F4B942; --cool:#48C7E8; --cold:#A9E9FB;
  --pad:16px; --r:3px;
}
*{box-sizing:border-box}
html,body{height:100%}
body{
  margin:0; background:var(--ground); color:var(--ink);
  font-family:"IBM Plex Sans",system-ui,-apple-system,"Segoe UI",sans-serif;
  font-size:14px; line-height:1.45; overflow:hidden;
}
.wrap{display:grid; grid-template-columns:272px 1fr; height:100vh}
@media (max-width:820px){ .wrap{grid-template-columns:1fr; grid-template-rows:auto 1fr; overflow:auto} }

aside{
  background:var(--panel); border-right:1px solid var(--line);
  padding:var(--pad); display:flex; flex-direction:column; gap:18px; overflow-y:auto;
}
h1{
  font-size:15px; font-weight:700; margin:0; letter-spacing:-.01em;
}
h1 span{display:block; font-size:11px; font-weight:500; color:var(--ink-3);
  letter-spacing:.08em; text-transform:uppercase; margin-top:3px}
.grp{display:flex; flex-direction:column; gap:8px}
.grp > .lbl{
  font-size:10px; font-weight:600; letter-spacing:.11em; text-transform:uppercase;
  color:var(--ink-3);
}
label.row{display:flex; align-items:center; gap:8px; font-size:13px; color:var(--ink-2); cursor:pointer}
label.row:hover{color:var(--ink)}
input[type=checkbox]{accent-color:var(--accent); width:14px; height:14px; cursor:pointer}
input[type=range]{width:100%; accent-color:var(--accent); cursor:pointer}
.pair{display:flex; align-items:center; justify-content:space-between; gap:10px;
  font-family:"IBM Plex Mono",ui-monospace,monospace; font-size:12px; font-variant-numeric:tabular-nums}
.pair b{font-weight:500; color:var(--ink); }
.pair span{color:var(--ink-3)}
button{
  font:inherit; font-size:12px; padding:6px 10px; border-radius:var(--r); cursor:pointer;
  background:var(--panel-2); color:var(--ink-2); border:1px solid var(--line);
}
button:hover{color:var(--ink); border-color:var(--accent-dim)}
button:focus-visible, input:focus-visible{outline:2px solid var(--accent); outline-offset:2px}

/* the ramp legend: the floor is a hard mark, because it is a hard threshold */
.ramp{height:10px; border-radius:2px; border:1px solid var(--line); position:relative}
.ramp .mark{position:absolute; top:-4px; bottom:-4px; width:1px; background:var(--ink)}
.ramp .mark::after{content:"floor"; position:absolute; top:-14px; left:-13px; font-size:9px;
  letter-spacing:.06em; color:var(--ink)}
.scale{display:flex; justify-content:space-between; font-family:"IBM Plex Mono",monospace;
  font-size:10px; color:var(--ink-3); font-variant-numeric:tabular-nums}

main{position:relative; min-height:0}
canvas{display:block; width:100%; height:100%; cursor:crosshair}
.stats{
  position:absolute; left:var(--pad); top:var(--pad); display:flex; gap:1px;
  background:var(--line); border:1px solid var(--line); border-radius:var(--r); overflow:hidden;
}
.stat{background:var(--panel); padding:7px 13px; min-width:82px}
.stat .k{font-size:9px; letter-spacing:.1em; text-transform:uppercase; color:var(--ink-3)}
.stat .v{font-family:"IBM Plex Mono",monospace; font-size:15px; font-variant-numeric:tabular-nums;
  color:var(--ink); margin-top:1px}
.stat.alarm .v{color:var(--hot)}
.hint{position:absolute; right:var(--pad); bottom:var(--pad); font-size:11px; color:var(--ink-3)}
.probe{
  position:absolute; pointer-events:none; background:var(--panel-2); border:1px solid var(--line);
  border-radius:var(--r); padding:6px 9px; font-family:"IBM Plex Mono",monospace; font-size:11px;
  font-variant-numeric:tabular-nums; color:var(--ink); white-space:pre; display:none;
  box-shadow:0 6px 20px rgba(0,0,0,.45);
}
</style>

<div class="wrap">
<aside>
  <h1>Bead Width Inspector<span id="src"></span></h1>

  <div class="grp">
    <div class="lbl">Extrusion width</div>
    <div class="ramp" id="ramp"><div class="mark" id="rampmark"></div></div>
    <div class="scale"><span id="s0"></span><span id="s1"></span></div>
  </div>

  <div class="grp">
    <div class="lbl">Layers</div>
    <div class="pair"><span>from</span><b id="lofrom"></b></div>
    <input type="range" id="lo" min="0" max="0" value="0">
    <div class="pair"><span>to</span><b id="loto"></b></div>
    <input type="range" id="hi" min="0" max="0" value="0">
  </div>

  <div class="grp" id="objs"><div class="lbl">Objects</div></div>

  <div class="grp">
    <div class="lbl">Show</div>
    <label class="row" id="trrow"><input type="checkbox" id="tr"> <span id="trlbl">Travel moves</span></label>
    <label class="row"><input type="checkbox" id="only"> Sub-floor beads only</label>
    <label class="row"><input type="checkbox" id="fade" checked> Fade older layers</label>
  </div>

  <div class="grp"><button id="fit">Reset view</button></div>
</aside>

<main>
  <canvas id="c"></canvas>
  <div class="stats" id="stats"></div>
  <div class="probe" id="probe"></div>
  <div class="hint">scroll to zoom &middot; drag to pan &middot; hover a bead to read it</div>
</main>
</div>

<script>
const D = __DATA__;
const cv = document.getElementById('c'), ctx = cv.getContext('2d');
const probe = document.getElementById('probe');
const FLOOR = D.floor, NOM = D.nominal;
const LO = Math.min(FLOOR * 0.55, NOM * 0.5), HI = NOM * 1.15;

// Below the floor the ramp goes hot: that is not decoration, it is the threshold the whole
// campaign is measured against, so it gets the one alarming colour on the page.
function ramp(w){
  if (w < FLOOR){ const t = Math.max(0, Math.min(1, (w - LO) / Math.max(1e-6, FLOOR - LO)));
                  return mix([255,106,77], [244,185,66], t); }
  const t = Math.max(0, Math.min(1, (w - FLOOR) / Math.max(1e-6, HI - FLOOR)));
  return mix([244,185,66], [169,233,251], t);
}
function mix(a,b,t){ return `rgb(${Math.round(a[0]+(b[0]-a[0])*t)},${Math.round(a[1]+(b[1]-a[1])*t)},${Math.round(a[2]+(b[2]-a[2])*t)})`; }

const layers = [...new Set(D.runs.map(r=>r.l))].sort((a,b)=>a-b);
const objOn = new Map(D.objects.map(o=>[o,true]));
let loI = 0, hiI = Math.max(0, layers.length-1);
let view = {s:1, x:0, y:0}, drag = null;

const el = id => document.getElementById(id);
el('src').textContent = D.title || '';
el('s0').textContent = LO.toFixed(2)+' mm';
el('s1').textContent = HI.toFixed(2)+' mm';
{ const g = []; for(let i=0;i<=24;i++){ const w = LO + (HI-LO)*i/24; g.push(ramp(w)+' '+(i/24*100)+'%'); }
  el('ramp').style.background = 'linear-gradient(90deg,'+g.join(',')+')';
  el('rampmark').style.left = ((FLOOR-LO)/(HI-LO)*100)+'%'; }

const objBox = el('objs');
for (const o of D.objects){
  const l = document.createElement('label'); l.className='row';
  l.innerHTML = `<input type="checkbox" checked> ${o}`;
  l.querySelector('input').addEventListener('change', e=>{ objOn.set(o, e.target.checked); draw(); });
  objBox.appendChild(l);
}
// A dead switch is bad, but this one's emptiness is the finding: a continuous-shell layer that
// never stops extruding has no travels inside an object at all. Say the count rather than
// offering a toggle that does nothing.
{
  const n = D.travels.length;
  el('trlbl').textContent = n ? `Travel moves (${n})` : 'Travel moves — none, layer is unbroken';
  if (!n){ el('tr').disabled = true; el('trrow').style.color = 'var(--ink-3)';
           el('trrow').style.cursor = 'default'; }
}
const loS = el('lo'), hiS = el('hi');
loS.max = hiS.max = Math.max(0, layers.length-1); hiS.value = hiI;
loS.addEventListener('input', ()=>{ loI = Math.min(+loS.value, hiI); loS.value = loI; draw(); });
hiS.addEventListener('input', ()=>{ hiI = Math.max(+hiS.value, loI); hiS.value = hiI; draw(); });
for (const id of ['tr','only','fade']) el(id).addEventListener('change', draw);
el('fit').addEventListener('click', ()=>{ fit(); draw(); });

function visible(r){
  if (!objOn.get(r.o)) return false;
  const i = layers.indexOf(r.l);
  return i >= loI && i <= hiI;
}
function bounds(){
  let x0=1e9,y0=1e9,x1=-1e9,y1=-1e9;
  for (const r of D.runs){ if(!objOn.get(r.o)) continue;
    for (let i=0;i<r.p.length;i+=2){ x0=Math.min(x0,r.p[i]); x1=Math.max(x1,r.p[i]);
                                     y0=Math.min(y0,r.p[i+1]); y1=Math.max(y1,r.p[i+1]); } }
  if (x0>x1) return {x0:0,y0:0,x1:1,y1:1};
  return {x0,y0,x1,y1};
}
function fit(){
  const b = bounds(), w = cv.width/devicePixelRatio, h = cv.height/devicePixelRatio;
  const s = 0.9 * Math.min(w/Math.max(b.x1-b.x0,1e-6), h/Math.max(b.y1-b.y0,1e-6));
  view.s = s;
  view.x = w/2 - s*(b.x0+b.x1)/2;
  view.y = h/2 + s*(b.y0+b.y1)/2;
}
const T = (x,y) => [view.x + x*view.s, view.y - y*view.s];

function resize(){
  const r = cv.getBoundingClientRect(), d = devicePixelRatio||1;
  cv.width = Math.max(1, r.width*d); cv.height = Math.max(1, r.height*d);
  ctx.setTransform(d,0,0,d,0,0);
}
function draw(){
  const w = cv.width/devicePixelRatio, h = cv.height/devicePixelRatio;
  ctx.clearRect(0,0,w,h);
  ctx.fillStyle = getComputedStyle(document.body).backgroundColor;
  ctx.fillRect(0,0,w,h);
  const onlySub = el('only').checked, fade = el('fade').checked;

  if (el('tr').checked){
    ctx.strokeStyle = 'rgba(107,126,140,.5)'; ctx.lineWidth = 1; ctx.setLineDash([3,3]);
    ctx.beginPath();
    for (const t of D.travels){
      if (!objOn.get(t.o)) continue;
      const i = layers.indexOf(t.l); if (i<loI||i>hiI) continue;
      const a = T(t.p[0],t.p[1]), b = T(t.p[2],t.p[3]);
      ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]);
    }
    ctx.stroke(); ctx.setLineDash([]);
  }

  ctx.lineCap = 'round'; ctx.lineJoin = 'round';
  let len=0, sub=0, vol=0, wmin=9, wmax=0;
  for (const r of D.runs){
    if (!visible(r)) continue;
    const li = layers.indexOf(r.l);
    const age = hiI>loI ? (li-loI)/(hiI-loI) : 1;
    const alpha = fade ? 0.25 + 0.75*age : 1;
    for (let i=0;i+3<r.p.length;i+=2){
      const wd = r.w[i/2+1] ?? r.w[r.w.length-1];
      const a = T(r.p[i],r.p[i+1]), b = T(r.p[i+2],r.p[i+3]);
      const d = Math.hypot(r.p[i+2]-r.p[i], r.p[i+3]-r.p[i+1]);
      len += d; vol += d*wd*D.layer_height; if (wd<FLOOR) sub += d;
      wmin = Math.min(wmin,wd); wmax = Math.max(wmax,wd);
      if (onlySub && wd >= FLOOR) continue;
      ctx.globalAlpha = alpha;
      ctx.strokeStyle = ramp(wd);
      ctx.lineWidth = Math.max(0.7, wd*view.s);
      ctx.beginPath(); ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]); ctx.stroke();
    }
  }
  ctx.globalAlpha = 1;
  const pct = len>0 ? 100*sub/len : 0;
  el('stats').innerHTML =
    stat('layers', (hiI-loI+1)) + stat('path', len.toFixed(0)+' mm') +
    stat('volume', vol.toFixed(0)+' mm³') +
    stat('width', wmax>0 ? wmin.toFixed(3)+'–'+wmax.toFixed(3) : '–') +
    stat('below floor', pct.toFixed(2)+'%', pct > 1);
}
function stat(k,v,alarm){
  return `<div class="stat${alarm?' alarm':''}"><div class="k">${k}</div><div class="v">${v}</div></div>`;
}

cv.addEventListener('wheel', e=>{
  e.preventDefault();
  const r = cv.getBoundingClientRect(), mx = e.clientX-r.left, my = e.clientY-r.top;
  const k = Math.exp(-e.deltaY*0.0015);
  view.x = mx - (mx-view.x)*k; view.y = my - (my-view.y)*k; view.s *= k;
  draw();
}, {passive:false});
cv.addEventListener('pointerdown', e=>{ drag={x:e.clientX,y:e.clientY}; cv.setPointerCapture(e.pointerId); });
cv.addEventListener('pointerup', ()=> drag=null);
cv.addEventListener('pointermove', e=>{
  if (drag){ view.x += e.clientX-drag.x; view.y += e.clientY-drag.y; drag={x:e.clientX,y:e.clientY}; draw(); return; }
  const r = cv.getBoundingClientRect(), mx = e.clientX-r.left, my = e.clientY-r.top;
  const wx = (mx-view.x)/view.s, wy = (view.y-my)/view.s;
  let best=null, bd=Infinity;
  for (const run of D.runs){
    if (!visible(run)) continue;
    for (let i=0;i+1<run.p.length;i+=2){
      const d = (run.p[i]-wx)**2 + (run.p[i+1]-wy)**2;
      if (d<bd){ bd=d; best={run, k:i/2}; }
    }
  }
  const tol = (12/view.s)**2;
  if (best && bd<tol){
    probe.style.display='block';
    probe.style.left = (mx+14)+'px'; probe.style.top = (my+14)+'px';
    probe.textContent =
      `${best.run.o}\nwidth  ${best.run.w[best.k].toFixed(3)} mm\nZ      ${best.run.z.toFixed(3)} mm\nlayer  ${best.run.l}`;
  } else probe.style.display='none';
});
cv.addEventListener('pointerleave', ()=> probe.style.display='none');

function boot(){ resize(); fit(); draw(); }
addEventListener('resize', ()=>{ resize(); draw(); });
function sync(){ el('lofrom').textContent = layers[loI] ?? '–'; el('loto').textContent = layers[hiI] ?? '–'; }
const _draw = draw; draw = function(){ sync(); _draw(); };
boot();
</script>
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("--z", nargs=2, type=float, metavar=("LO", "HI"))
    ap.add_argument("--layers", nargs=2, type=int, metavar=("A", "B"))
    ap.add_argument("--object", action="append")
    ap.add_argument("--out")
    a = ap.parse_args()

    data = parse(a.file, a.object)
    data = select(data, a.z[0] if a.z else None, a.z[1] if a.z else None,
                  a.layers[0] if a.layers else None, a.layers[1] if a.layers else None)
    if not data["runs"]:
        sys.exit("no extrusion in that range - widen --z or check --object")
    data["title"] = os.path.basename(a.file)

    out = a.out or os.path.splitext(os.path.splitext(a.file)[0])[0] + ".view.html"
    html = HTML.replace("__DATA__", json.dumps(data, separators=(",", ":")))
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    ls = sorted({r["l"] for r in data["runs"]})
    print(f"{out}  ({len(data['runs'])} runs, layers {ls[0]}-{ls[-1]}, "
          f"{os.path.getsize(out)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
