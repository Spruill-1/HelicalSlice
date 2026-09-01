"""Render a Z-range of a sliced plate as an inspectable toolpath viewer.

    viewgcode.py <file.gcode.3mf>... [--z LO HI | --layers A B] [--object NAME]... [--out FILE]

Pass more than one file to compare them side by side on a shared colour scale and a shared
view: pan and zoom move both, so a difference in the weave is a difference in the picture
rather than a difference in two pictures.

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
    # The width this page derives from E is the Flow SPACING, not the plan width: area/h =
    # width - h*(1-pi/4), i.e. 0.0537 low at these settings. The floor must live in the same
    # units or every bead PLANNED exactly at the 0.34 floor reads as sub-floor - which painted
    # a third of the bored crescent red for material that is actually 2-4x OVERFILLED.
    return {"layer_height": lh, "nozzle": nozzle, "nominal": nominal,
            "floor": round(max(0.30, 0.85 * nozzle) - (1 - math.pi / 4) * lh, 4),
            "objects": objs, "runs": runs, "travels": travels}


def shrink(sets):
    """Round the payload down to what the page can actually show.

    Coordinates carry six significant figures out of the G-code and the page draws them a few
    hundred pixels wide, so the last three are bytes nobody sees. 0.01 mm is a fortieth of a
    bead. A nine-panel, four-variant page came to 16.6 MB and the publish limit is 16."""
    for st in sets:
        for r in st["runs"]:
            r["p"] = [round(v, 2) for v in r["p"]]
            r["w"] = [round(v, 3) for v in r["w"]]
        for t in st["travels"]:
            t["p"] = [round(v, 2) for v in t["p"]]
    return sets


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
.seg{display:flex; gap:1px; background:var(--line); border:1px solid var(--line);
  border-radius:var(--r); overflow:hidden}
.seg button{flex:1; border:0; border-radius:0; background:var(--panel-2)}
.seg button[aria-pressed=true]{background:var(--accent-dim); color:var(--ink)}
.note{font-family:"IBM Plex Mono",ui-monospace,monospace; font-size:11px; line-height:1.5;
  color:var(--ink-2); font-variant-numeric:tabular-nums; white-space:pre; overflow-x:auto;
  margin:0}
.note b{color:var(--ink); font-weight:600}
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
    <div class="lbl">Colour by</div>
    <div class="seg">
      <button id="mw" aria-pressed="true">width</button>
      <button id="ml" aria-pressed="false">layer</button>
    </div>
    <div class="ramp" id="ramp"><div class="mark" id="rampmark"></div></div>
    <div class="scale"><span id="s0"></span><span id="s1"></span></div>
  </div>

  <div class="grp" id="vargrp" hidden>
    <div class="lbl" id="varlbl">Variant</div>
    <div class="pair"><span id="varname"></span><b id="varval"></b></div>
    <input type="range" id="var" min="0" max="0" value="0" step="1">
  </div>

  <div class="grp" id="notegrp" hidden>
    <div class="lbl" id="notelbl">Measured</div>
    <pre class="note" id="note"></pre>
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
// Colouring by LAYER is what makes a between-layer effect visible at all: a rotation from
// one layer to the next is invisible when every layer is the same colour, and it is the
// picture the swirl knob has to be judged on. Deliberately not the width family - violet to
// pale is ordered, and nobody will read it as a bead width.
function layerCol(t){
  t = Math.max(0, Math.min(1, t));
  return t < 0.5 ? mix([108,75,224], [72,199,232], t*2)
                 : mix([72,199,232], [234,246,251], (t-0.5)*2);
}
let cmode = 'width';
function colFor(wd, t){ return cmode === 'width' ? ramp(wd) : layerCol(t); }
function mix(a,b,t){ return `rgb(${Math.round(a[0]+(b[0]-a[0])*t)},${Math.round(a[1]+(b[1]-a[1])*t)},${Math.round(a[2]+(b[2]-a[2])*t)})`; }

// One variant is live at a time; the layer range, the view and the colour mode all survive
// switching, so stepping the slider is a true before/after on the same geometry.
const VARIANTS = D.variants || null;
let vi = 0;
if (VARIANTS){
  let best = 0;
  VARIANTS.forEach((v, i) => {
    if (Math.abs(v.value - D.variant_default) < Math.abs(VARIANTS[best].value - D.variant_default))
      best = i;
  });
  vi = best;
}
let SETS = VARIANTS ? VARIANTS[vi].sets : D.sets;
// Layers come from EVERY variant, so the layer slider does not change under the variant
// slider - two controls that move each other are two controls nobody can reason about.
const ALLRUNS = (VARIANTS ? VARIANTS.flatMap(v=>v.sets) : SETS).flatMap(s=>s.runs);
const layers = [...new Set(ALLRUNS.map(r=>r.l))].sort((a,b)=>a-b);
const objOn = new Map(D.objects.map(o=>[o,true]));
let loI = 0, hiI = Math.max(0, layers.length-1);
let view = {s:1, x:0, y:0}, drag = null;

const el = id => document.getElementById(id);
el('src').textContent = D.title || '';
function legend(){
  const g = [];
  for (let i=0;i<=24;i++){
    const t = i/24;
    g.push((cmode==='width' ? ramp(LO+(HI-LO)*t) : layerCol(t)) + ' ' + (t*100) + '%');
  }
  el('ramp').style.background = 'linear-gradient(90deg,'+g.join(',')+')';
  el('rampmark').style.display = cmode==='width' ? 'block' : 'none';
  el('rampmark').style.left = ((FLOOR-LO)/(HI-LO)*100)+'%';
  el('s0').textContent = cmode==='width' ? LO.toFixed(2)+' mm' : 'layer '+layers[loI];
  el('s1').textContent = cmode==='width' ? HI.toFixed(2)+' mm' : 'layer '+layers[hiI];
}
function setMode(m){
  cmode = m;
  el('mw').setAttribute('aria-pressed', String(m==='width'));
  el('ml').setAttribute('aria-pressed', String(m==='layer'));
  // Fade is an age cue, and in layer mode the colour already is one. Two of them at once just
  // makes the early layers hard to see.
  el('fade').disabled = (m==='layer');
  legend(); draw();
}
el('mw').addEventListener('click', ()=>setMode('width'));
el('ml').addEventListener('click', ()=>setMode('layer'));
// Joined on the python side: a newline escape has to survive a raw string, a JSON dump and
// an HTML template, and it did not.
if (D.notes){ el('notegrp').hidden = false; el('note').textContent = D.notes; }
if (VARIANTS){
  el('vargrp').hidden = false;
  el('varlbl').textContent = D.variant_label || 'Variant';
  el('varname').textContent = D.variant_label || '';
  const vs = el('var');
  vs.max = VARIANTS.length - 1; vs.value = vi;
  const showVar = ()=>{ el('varval').textContent = VARIANTS[vi].value + '%'; };
  vs.addEventListener('input', ()=>{
    vi = +vs.value; SETS = VARIANTS[vi].sets; showVar(); draw();
  });
  showVar();
}

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
  const n = SETS.reduce((a,s)=>a+s.travels.length,0);
  el('trlbl').textContent = n ? `Travel moves (${n})` : 'Travel moves — none, layer is unbroken';
  if (!n){ el('tr').disabled = true; el('trrow').style.color = 'var(--ink-3)';
           el('trrow').style.cursor = 'default'; }
}
const loS = el('lo'), hiS = el('hi');
loS.max = hiS.max = Math.max(0, layers.length-1); hiS.value = hiI;
// Once the two ends meet, the pair reads as ONE slice, and pushing either end then drags the
// other with it. That turns the range control into a single-layer stepper without adding a
// third widget: close the range, then walk it up or down a layer at a time.
loS.addEventListener('input', ()=>{
  const v = +loS.value;
  if (loI === hiI && v > hiI){ hiI = Math.min(v, +hiS.max); hiS.value = hiI; }
  loI = Math.min(v, hiI); loS.value = loI; legend(); draw();
});
hiS.addEventListener('input', ()=>{
  const v = +hiS.value;
  if (hiI === loI && v < loI){ loI = Math.max(v, 0); loS.value = loI; }
  hiI = Math.max(v, loI); hiS.value = hiI; legend(); draw();
});
for (const id of ['tr','only','fade']) el(id).addEventListener('change', draw);
el('fit').addEventListener('click', ()=>{ fit(); draw(); });

function visible(r){
  if (!objOn.get(r.o)) return false;
  const i = layers.indexOf(r.l);
  return i >= loI && i <= hiI;
}
// Across every panel AND every variant: the cell has to hold the largest shape, and the
// layout must not shift when the variant slider moves.
function bounds(){
  let x0=1e9,y0=1e9,x1=-1e9,y1=-1e9;
  const all = (VARIANTS ? VARIANTS.flatMap(v=>v.sets) : SETS).flatMap(s=>s.runs);
  for (const r of all){ if(!objOn.get(r.o)) continue;
    for (let i=0;i<r.p.length;i+=2){ x0=Math.min(x0,r.p[i]); x1=Math.max(x1,r.p[i]);
                                     y0=Math.min(y0,r.p[i+1]); y1=Math.max(y1,r.p[i+1]); } }
  if (x0>x1) return {x0:0,y0:0,x1:1,y1:1};
  return {x0,y0,x1,y1};
}
function grid(n){ const c = Math.ceil(Math.sqrt(n)); return [c, Math.ceil(n / c)]; }
function fit(){
  const b = bounds(), w = cv.width/devicePixelRatio, h = cv.height/devicePixelRatio;
  const [cols, rows] = grid(SETS.length);
  COLS = cols;
  const s = 0.92 * Math.min(w / Math.max((b.x1-b.x0) * 1.12 * cols, 1e-6),
                            h / Math.max((b.y1-b.y0) * 1.12 * rows, 1e-6));
  view.s = s; bx = (b.x0+b.x1)/2;
  PW = (b.x1-b.x0)*s*1.12; PH = (b.y1-b.y0)*s*1.12;
  view.x = w/2 - PW*(cols-1)/2;
  view.y = h/2 + s*(b.y0+b.y1)/2 - PH*(rows-1)/2;
}
// Every set is drawn on the SAME transform, offset into a grid cell, so the panels are directly
// comparable rather than each fitted to itself - a row of nine is a strip nobody can see at
// once. `bx`, `PW` and `PH` are cached per draw: deriving the cell size from bounds() inside
// the transform made it O(points x points).
let bx = 0, PW = 0, PH = 0, COLS = 1;
const T = (x,y,k) => [view.x + (x - bx)*view.s + (k % COLS)*PW,
                      view.y - y*view.s + Math.floor(k / COLS)*PH];

function resize(){
  const r = cv.getBoundingClientRect(), d = devicePixelRatio||1;
  cv.width = Math.max(1, r.width*d); cv.height = Math.max(1, r.height*d);
  ctx.setTransform(d,0,0,d,0,0);
}
function draw(){
  const w = cv.width/devicePixelRatio, h = cv.height/devicePixelRatio;
  { const b = bounds();
    PW = (b.x1-b.x0)*view.s*1.12; PH = (b.y1-b.y0)*view.s*1.12;
    COLS = grid(SETS.length)[0]; }
  ctx.clearRect(0,0,w,h);
  ctx.fillStyle = getComputedStyle(document.body).backgroundColor;
  ctx.fillRect(0,0,w,h);
  const onlySub = el('only').checked, fade = el('fade').checked;

  ctx.lineCap = 'round'; ctx.lineJoin = 'round';
  let len=0, sub=0, vol=0, wmin=9, wmax=0;
  SETS.forEach((S,k)=>{
    if (el('tr').checked){
      ctx.strokeStyle = 'rgba(107,126,140,.5)'; ctx.lineWidth = 1; ctx.setLineDash([3,3]);
      ctx.beginPath();
      for (const t of S.travels){
        if (!objOn.get(t.o)) continue;
        const i = layers.indexOf(t.l); if (i<loI||i>hiI) continue;
        const a = T(t.p[0],t.p[1],k), b = T(t.p[2],t.p[3],k);
        ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]);
      }
      ctx.stroke(); ctx.setLineDash([]);
    }
    for (const r of S.runs){
      if (!visible(r)) continue;
      const li = layers.indexOf(r.l);
      const age = hiI>loI ? (li-loI)/(hiI-loI) : 1;
      const alpha = (fade && cmode==='width') ? 0.25 + 0.75*age : 1;
      for (let i=0;i+3<r.p.length;i+=2){
        const wd = r.w[i/2+1] ?? r.w[r.w.length-1];
        const a = T(r.p[i],r.p[i+1],k), b = T(r.p[i+2],r.p[i+3],k);
        const d = Math.hypot(r.p[i+2]-r.p[i], r.p[i+3]-r.p[i+1]);
        if (k===0){ len += d; vol += d*wd*D.layer_height; if (wd<FLOOR) sub += d;
                    wmin = Math.min(wmin,wd); wmax = Math.max(wmax,wd); }
        if (onlySub && wd >= FLOOR) continue;
        ctx.globalAlpha = alpha;
        ctx.strokeStyle = colFor(wd, age);
        ctx.lineWidth = Math.max(0.7, wd*view.s);
        ctx.beginPath(); ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]); ctx.stroke();
      }
    }
    if (SETS.length>1){
      const b = bounds(), p0 = T(b.x0, b.y0, k);
      ctx.globalAlpha = 1; ctx.fillStyle = '#9DAFBC';
      ctx.font = '600 13px "IBM Plex Mono", monospace'; ctx.textAlign = 'left';
      ctx.fillText(S.label, p0[0], p0[1] + 22);
    }
  });
  ctx.globalAlpha = 1;
  const pct = len>0 ? 100*sub/len : 0;
  el('stats').innerHTML =
    stat('layers', (hiI-loI+1)) + stat('path', len.toFixed(0)+' mm') +
    stat('volume', vol.toFixed(0)+' mm³') +
    stat('width', wmax>0 ? wmin.toFixed(3)+'–'+wmax.toFixed(3) : '–') +
    stat('below floor', pct.toFixed(2)+'%', pct > 1) +
    (SETS.length>1 ? stat('shown', SETS.length+' sets') : '');
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
  // Screen x = view.x + (X - bx)*s + col*PW and screen y = view.y - Y*s + row*PH, so the world
  // point under the cursor in panel k inverts both.
  SETS.forEach((S,k)=>{
    const cx = wx - (k % COLS)*PW/view.s + bx;
    const cy = wy + Math.floor(k / COLS)*PH/view.s;
    for (const run of S.runs){
      if (!visible(run)) continue;
      for (let i=0;i+1<run.p.length;i+=2){
        const d = (run.p[i]-cx)**2 + (run.p[i+1]-cy)**2;
        if (d<bd){ bd=d; best={run, k:i/2, set:S}; }
      }
    }
  });
  const tol = (12/view.s)**2;
  if (best && bd<tol){
    probe.style.display='block';
    probe.style.left = (mx+14)+'px'; probe.style.top = (my+14)+'px';
    probe.textContent =
      `${best.run.o}\nwidth  ${best.run.w[best.k].toFixed(3)} mm\nZ      ${best.run.z.toFixed(3)} mm\nlayer  ${best.run.l}`;
  } else probe.style.display='none';
});
cv.addEventListener('pointerleave', ()=> probe.style.display='none');

function boot(){ resize(); fit(); legend(); draw(); }
addEventListener('resize', ()=>{ resize(); draw(); });
function sync(){ el('lofrom').textContent = layers[loI] ?? '–'; el('loto').textContent = layers[hiI] ?? '–'; }
const _draw = draw; draw = function(){ sync(); _draw(); };
boot();
</script>
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", nargs="*")
    # A VARIANT is the same panels sliced at a different value of one knob. Comparing four
    # weld settings by putting them side by side costs four panels for one shape; as a slider
    # it costs none, and the eye compares a shape against ITSELF between settings, which is
    # the comparison that matters.
    ap.add_argument("--variant", action="append", nargs="+", metavar=("VALUE", "FILE"),
                    default=[], help="a knob value followed by one file per panel; repeatable")
    ap.add_argument("--variant-label", default="", help="what the variant slider is measuring")
    ap.add_argument("--variant-default", type=float, default=None,
                    help="which variant value to open on")
    ap.add_argument("--z", nargs=2, type=float, metavar=("LO", "HI"))
    ap.add_argument("--layers", nargs=2, type=int, metavar=("A", "B"))
    ap.add_argument("--object", action="append")
    ap.add_argument("--out")
    ap.add_argument("--note", action="append", default=[],
                    help="a line of measured context to print in the sidebar; repeatable")
    a = ap.parse_args()

    def load(files):
        out, first = [], None
        for f in files:
            d = parse(f, a.object)
            d = select(d, a.z[0] if a.z else None, a.z[1] if a.z else None,
                       a.layers[0] if a.layers else None, a.layers[1] if a.layers else None)
            if not d["runs"]:
                sys.exit(f"no extrusion in that range for {f} - widen --z or check --object")
            first = first or d
            out.append({"label": os.path.basename(f).replace(".gcode.3mf", "").split("@")[0],
                        "runs": d["runs"], "travels": d["travels"]})
        return out, first

    if a.variant:
        groups, base = [], None
        for grp in a.variant:
            val, files = float(grp[0]), grp[1:]
            sets, first = load(files)
            base = base or first
            groups.append({"value": val, "sets": shrink(sets)})
        n = {len(g["sets"]) for g in groups}
        if len(n) != 1:
            sys.exit(f"every --variant needs the same number of files, got {sorted(n)}")
        # EVERY file's objects, not the first one's. Each shape is sliced on its own, so each
        # file names a different object; taking the list from the first file left `objOn`
        # without an entry for the other eight, `visible()` read undefined, and eight of nine
        # panels drew nothing at all while the sidebar offered one checkbox.
        data = {k: base[k] for k in ("layer_height", "nozzle", "nominal", "floor")}
        data["objects"] = sorted({o for g in groups for st in g["sets"]
                                    for o in {r["o"] for r in st["runs"]}})
        data["variants"] = groups
        data["sets"] = groups[0]["sets"]
        data["variant_label"] = a.variant_label
        data["variant_default"] = (a.variant_default if a.variant_default is not None
                                   else groups[0]["value"])
        data["notes"] = chr(10).join(a.note)
        data["title"] = f"{len(groups[0]['sets'])} shapes"
        out = a.out or "variants.view.html"
        with open(out, "w", encoding="utf-8") as f:
            f.write(HTML.replace("__DATA__", json.dumps(data, separators=(",", ":"))))
        print(f"{out}  ({len(groups)} variants x {len(groups[0]['sets'])} panels, "
              f"{os.path.getsize(out)/1e6:.1f} MB)")
        return

    if not a.file:
        sys.exit("give at least one file, or use --variant")
    sets = []
    base = None
    for f in a.file:
        d = parse(f, a.object)
        d = select(d, a.z[0] if a.z else None, a.z[1] if a.z else None,
                   a.layers[0] if a.layers else None, a.layers[1] if a.layers else None)
        if not d["runs"]:
            sys.exit(f"no extrusion in that range for {f} - widen --z or check --object")
        base = base or d
        sets.append({"label": os.path.basename(f).replace(".gcode.3mf", ""),
                     "runs": d["runs"], "travels": d["travels"]})
    data = {k: base[k] for k in ("layer_height", "nozzle", "nominal", "floor")}
    # Same union as the variant path: files that name different objects each need an entry.
    data["objects"] = sorted({r["o"] for st in sets for r in st["runs"]})
    data["sets"] = shrink(sets)
    data["notes"] = chr(10).join(a.note)
    data["title"] = " vs ".join(s["label"] for s in sets) if len(sets) > 1                     else os.path.basename(a.file[0])

    out = a.out or os.path.splitext(os.path.splitext(a.file[0])[0])[0] + ".view.html"
    html = HTML.replace("__DATA__", json.dumps(data, separators=(",", ":")))
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    ls = sorted({r["l"] for st in sets for r in st["runs"]})
    print(f"{out}  ({len(sets)} set(s), {sum(len(st['runs']) for st in sets)} runs, "
          f"layers {ls[0]}-{ls[-1]}, {os.path.getsize(out)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
