"""
Build a self-contained, interactive 3D structure viewer (one HTML file).
======================================================================
Reads all 29 AlphaFold structures from data/raw/ + the clean dataset, and writes
figures/structure_viewer.html — a rotatable, dark-background viewer colored by
AlphaFold confidence (pLDDT). No internet or libraries needed to open it.

Run:  python3 src/build_viewer.py
"""

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")
DATASET = os.path.join(PROJECT_ROOT, "data", "processed", "autoantigens.json")
OUT = os.path.join(PROJECT_ROOT, "figures", "structure_viewer.html")


def read_ca(acc):
    """Centered alpha-carbon coords + pLDDT for one protein: [[x,y,z,p], ...]."""
    path = os.path.join(RAW_DIR, "{}_alphafold.pdb".format(acc))
    pts = []
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                pts.append([float(line[30:38]), float(line[38:46]),
                            float(line[46:54]), float(line[60:66])])
    cx = sum(p[0] for p in pts) / len(pts)
    cy = sum(p[1] for p in pts) / len(pts)
    cz = sum(p[2] for p in pts) / len(pts)
    return [[round(p[0] - cx, 1), round(p[1] - cy, 1), round(p[2] - cz, 1),
             round(p[3])] for p in pts]


def main():
    records = json.load(open(DATASET))["autoantigens"]
    proteins = []
    for r in records:
        proteins.append({
            "short": r["short_name"], "name": r["protein_name"], "gene": r["gene"],
            "acc": r["accession"], "len": r["length_aa"], "plddt": r["structure_mean_plddt"],
            "role": r["role"], "ca": read_ca(r["accession"]),
        })
    data_json = json.dumps(proteins, separators=(",", ":"))

    html = TEMPLATE.replace("/*__DATA__*/", "const PROTEINS = " + data_json + ";")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        f.write(html)
    kb = os.path.getsize(OUT) // 1024
    print("Wrote {}  ({} proteins, {} KB)".format(OUT, len(proteins), kb))


TEMPLATE = r"""<title>Autoantigen Structure Viewer</title>
<style>
  :root{
    --bg:#0A0B10; --bg2:#12141C; --panel:#141824; --line:#232838;
    --ink:#E6E9F2; --muted:#8A91A8; --accent:#65CBF3;
    --pos:#65CBF3; --ctrl:#FFB454;
    --vh:#2E74FF; --cf:#65CBF3; --lo:#FFDB13; --vl:#FF7D45;
    --mono:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
    --sans:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
  }
  *{box-sizing:border-box}
  html,body{margin:0;height:100%}
  body{
    background:radial-gradient(130% 100% at 70% 0%, #151826 0%, var(--bg) 60%);
    color:var(--ink); font-family:var(--sans);
    height:100vh; overflow:hidden;
  }
  .app{display:grid; grid-template-columns:272px 1fr; height:100vh}

  /* ---- sidebar ---- */
  .side{border-right:1px solid var(--line); background:var(--bg2); display:flex;
        flex-direction:column; min-height:0}
  .brand{padding:18px 18px 14px; border-bottom:1px solid var(--line)}
  .brand h1{margin:0; font-size:13px; letter-spacing:.14em; text-transform:uppercase;
            font-weight:650}
  .brand p{margin:6px 0 0; font-size:12px; color:var(--muted); line-height:1.5}
  .list{overflow-y:auto; padding:10px; min-height:0; flex:1}
  .group{font-size:10px; letter-spacing:.16em; text-transform:uppercase;
         color:var(--muted); padding:12px 8px 6px}
  .item{display:flex; align-items:center; gap:10px; width:100%; text-align:left;
        background:none; border:1px solid transparent; border-radius:9px;
        color:var(--ink); padding:8px 10px; cursor:pointer; font-size:13px;
        font-family:inherit}
  .item:hover{background:#1a1f2e}
  .item[aria-current="true"]{background:#1b2740; border-color:#2b3f66}
  .item .dot{width:8px; height:8px; border-radius:50%; flex:none}
  .item .nm{flex:1; white-space:nowrap; overflow:hidden; text-overflow:ellipsis}
  .item .pl{font-family:var(--mono); font-size:11px; color:var(--muted);
            font-variant-numeric:tabular-nums}

  /* ---- stage ---- */
  .stage{position:relative; min-width:0}
  canvas{display:block; width:100%; height:100%; touch-action:none; cursor:grab}
  canvas:active{cursor:grabbing}

  .hud{position:absolute; pointer-events:none; font-family:var(--mono)}
  .info{top:20px; left:22px}
  .info .sh{font-family:var(--sans); font-size:26px; font-weight:680;
            letter-spacing:-.01em; margin:0}
  .info .nm{font-size:12px; color:var(--muted); margin:4px 0 10px; max-width:42ch;
            line-height:1.5; font-family:var(--sans)}
  .meta{display:flex; gap:16px; font-size:12px; color:var(--ink);
        font-variant-numeric:tabular-nums}
  .meta b{color:var(--muted); font-weight:400}
  .chip{display:inline-block; margin-top:12px; font-family:var(--sans); font-size:11px;
        letter-spacing:.04em; padding:4px 10px; border-radius:999px; border:1px solid}
  .chip.pos{color:var(--pos); border-color:#2b4f66; background:#0e2330}
  .chip.ctrl{color:var(--ctrl); border-color:#5c4420; background:#2a1f0c}

  .legend{bottom:20px; left:22px; display:flex; flex-direction:column; gap:7px}
  .legend .lt{font-size:10px; letter-spacing:.14em; text-transform:uppercase;
              color:var(--muted); margin-bottom:2px; font-family:var(--sans)}
  .legend .row{display:flex; align-items:center; gap:9px; font-size:11.5px; color:var(--ink)}
  .legend .sw{width:12px; height:12px; border-radius:3px}

  .hint{bottom:20px; right:22px; font-size:11px; color:var(--muted); text-align:right;
        line-height:1.7}
  .spin{position:absolute; top:18px; right:22px; pointer-events:auto;
        font-family:var(--sans); font-size:11px; letter-spacing:.04em; color:var(--muted);
        background:#141a28; border:1px solid var(--line); border-radius:8px;
        padding:6px 11px; cursor:pointer}
  .spin:hover{color:var(--ink); border-color:#33405e}
  .spin:focus-visible,.item:focus-visible{outline:2px solid var(--accent);
        outline-offset:2px}

  @media (max-width:760px){
    .app{grid-template-columns:1fr; grid-template-rows:auto 1fr}
    .side{border-right:none; border-bottom:1px solid var(--line)}
    .list{display:flex; flex-wrap:nowrap; overflow-x:auto; overflow-y:hidden;
          gap:6px; padding:10px}
    .group{display:none}
    .item{width:auto; flex:none}
    .info .nm{display:none}
  }
</style>

<div class="app">
  <aside class="side">
    <div class="brand">
      <h1>Autoantigen Atlas</h1>
      <p>AlphaFold-predicted structures of lupus autoantigens, colored by prediction
         confidence (pLDDT).</p>
    </div>
    <nav class="list" id="list" aria-label="Proteins"></nav>
  </aside>

  <main class="stage">
    <canvas id="c"></canvas>
    <button class="spin" id="spin" aria-pressed="true">Auto-rotate: on</button>
    <div class="hud info" id="info"></div>
    <div class="hud legend">
      <div class="lt">AlphaFold confidence</div>
      <div class="row"><span class="sw" style="background:var(--vh)"></span>Very high &gt; 90</div>
      <div class="row"><span class="sw" style="background:var(--cf)"></span>Confident 70&ndash;90</div>
      <div class="row"><span class="sw" style="background:var(--lo)"></span>Low 50&ndash;70</div>
      <div class="row"><span class="sw" style="background:var(--vl)"></span>Very low &lt; 50</div>
    </div>
    <div class="hud hint">drag to rotate<br>scroll to zoom</div>
  </main>
</div>

<script>
/*__DATA__*/

const cv = document.getElementById("c"), ctx = cv.getContext("2d");
const listEl = document.getElementById("list"), infoEl = document.getElementById("info");
const spinBtn = document.getElementById("spin");
const reduce = matchMedia("(prefers-reduced-motion: reduce)").matches;

let cur = 0, auto = !reduce, dragging = false, lastX = 0, lastY = 0;
let rotX = 0.4, rotY = 0.6, zoom = 1;
spinBtn.setAttribute("aria-pressed", String(auto));
spinBtn.textContent = "Auto-rotate: " + (auto ? "on" : "off");

function band(p){ return p>90?[46,116,255] : p>=70?[101,203,243] : p>=50?[255,219,19] : [255,125,69]; }

function buildList(){
  const groups = [["Autoantigens","positive"],["Controls","negative_control"]];
  for(const [title,role] of groups){
    const items = PROTEINS.map((p,i)=>[p,i]).filter(([p])=>p.role===role);
    const h = document.createElement("div"); h.className="group";
    h.textContent = title + " (" + items.length + ")"; listEl.appendChild(h);
    for(const [p,i] of items){
      const b = document.createElement("button");
      b.className="item"; b.dataset.i=i;
      const [r,g,bl]=band(p.plddt);
      b.innerHTML = '<span class="dot" style="background:rgb('+r+','+g+','+bl+')"></span>'
        + '<span class="nm">'+p.short+'</span><span class="pl">'+p.plddt.toFixed(1)+'</span>';
      b.onclick = ()=>select(i);
      listEl.appendChild(b);
    }
  }
}

function select(i){
  cur=i; rotX=0.4; rotY=0.6; zoom=1;
  for(const el of listEl.querySelectorAll(".item"))
    el.setAttribute("aria-current", String(+el.dataset.i===i));
  const p = PROTEINS[i], isPos = p.role==="positive";
  infoEl.innerHTML =
    '<p class="sh">'+p.short+'</p>'
    +'<p class="nm">'+p.name+'</p>'
    +'<div class="meta"><span><b>gene</b> '+p.gene+'</span>'
    +'<span><b>'+p.acc+'</b></span><span>'+p.len+' aa</span>'
    +'<span><b>pLDDT</b> '+p.plddt.toFixed(1)+'</span></div>'
    +'<span class="chip '+(isPos?'pos':'ctrl')+'">'
    +(isPos?'Lupus autoantigen':'Negative control')+'</span>';
}

function resize(){
  const r = cv.getBoundingClientRect(), dpr = Math.min(devicePixelRatio||1, 2);
  cv.width = r.width*dpr; cv.height = r.height*dpr; ctx.setTransform(dpr,0,0,dpr,0,0);
}

function draw(){
  const p = PROTEINS[cur], ca = p.ca, n = ca.length;
  const W = cv.clientWidth, H = cv.clientHeight;
  ctx.clearRect(0,0,W,H);
  const cosX=Math.cos(rotX),sinX=Math.sin(rotX),cosY=Math.cos(rotY),sinY=Math.sin(rotY);

  let maxR=1;
  for(const a of ca){ const d=Math.hypot(a[0],a[1],a[2]); if(d>maxR)maxR=d; }
  const scale = Math.min(W,H)*0.42/maxR*zoom, cxp=W/2, cyp=H/2;

  const pr = new Array(n);
  let zmin=1e9, zmax=-1e9;
  for(let i=0;i<n;i++){
    let x=ca[i][0],y=ca[i][1],z=ca[i][2];
    let x1=x*cosY - z*sinY, z1=x*sinY + z*cosY;       // yaw
    let y2=y*cosX - z1*sinX, z2=y*sinX + z1*cosX;      // pitch
    pr[i]=[cxp+x1*scale, cyp - y2*scale, z2, ca[i][3]];
    if(z2<zmin)zmin=z2; if(z2>zmax)zmax=z2;
  }
  const fog = v => { const t=(v-zmin)/(zmax-zmin+1e-6); return 0.35+0.65*t; }; // far=dim

  // backbone
  ctx.lineCap="round";
  for(let i=0;i<n-1;i++){
    const a=pr[i], b=pr[i+1], al=fog((a[2]+b[2])/2);
    const [r,g,bl]=band(a[3]);
    ctx.strokeStyle="rgba("+r+","+g+","+bl+","+(al*0.55).toFixed(3)+")";
    ctx.lineWidth=2.1*al;
    ctx.beginPath(); ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]); ctx.stroke();
  }
  // residues, depth-sorted (far first)
  const order = pr.map((_,i)=>i).sort((i,j)=>pr[i][2]-pr[j][2]);
  for(const i of order){
    const a=pr[i], al=fog(a[2]), [r,g,bl]=band(a[3]);
    ctx.beginPath(); ctx.arc(a[0],a[1], 2.3*al+1.1, 0, 6.283);
    ctx.fillStyle="rgba("+r+","+g+","+bl+","+al.toFixed(3)+")"; ctx.fill();
  }
}

function frame(){ if(auto && !dragging) rotY += 0.004; draw(); requestAnimationFrame(frame); }

// interaction
cv.addEventListener("pointerdown", e=>{dragging=true; lastX=e.clientX; lastY=e.clientY; cv.setPointerCapture(e.pointerId);});
cv.addEventListener("pointermove", e=>{ if(!dragging)return;
  rotY += (e.clientX-lastX)*0.01; rotX += (e.clientY-lastY)*0.01;
  rotX = Math.max(-1.5, Math.min(1.5, rotX)); lastX=e.clientX; lastY=e.clientY; });
cv.addEventListener("pointerup", e=>{dragging=false;});
cv.addEventListener("wheel", e=>{ e.preventDefault();
  zoom *= e.deltaY<0 ? 1.08 : 0.925; zoom=Math.max(0.4,Math.min(4,zoom)); }, {passive:false});
spinBtn.onclick = ()=>{ auto=!auto; spinBtn.setAttribute("aria-pressed",String(auto));
  spinBtn.textContent="Auto-rotate: "+(auto?"on":"off"); };
addEventListener("resize", resize);

buildList(); resize(); select(0); frame();
</script>
"""


if __name__ == "__main__":
    main()
