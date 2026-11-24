'use strict';
const $=id=>document.getElementById(id), clone=x=>JSON.parse(JSON.stringify(x));
const configurations=window.LAB_CONFIGURATIONS, storageKey='erdos99-playground-v1';
let points=[],original=[],source='',selected=-1,history=[],drag=null,view,experiments=[];
function measure(p,tol){
 const distances=p.map(()=>[]),edges=[],violations=[],triangles=[];let min=Infinity,max=0;
 for(let i=0;i<p.length;i++)for(let j=i+1;j<p.length;j++){
  const d=Math.hypot(p[i][0]-p[j][0],p[i][1]-p[j][1]);distances[i][j]=distances[j][i]=d;
  min=Math.min(min,d);max=Math.max(max,d);if(d<1)violations.push([i,j]);if(Math.abs(d-1)<=tol)edges.push([i,j]);
 }
 for(let i=0;i<p.length;i++)for(let j=i+1;j<p.length;j++)for(let k=j+1;k<p.length;k++)if([distances[i][j],distances[i][k],distances[j][k]].every(d=>Math.abs(d-1)<=tol))triangles.push([i,j,k]);
 return {min,max,edges,violations,triangles};
}
function tolerance(){const n=Number($('tol').value);return Number.isFinite(n)&&n>=0&&n<=.1?n:1e-6;}
function fit(){const all=[...points,...original],xs=all.map(p=>p[0]),ys=all.map(p=>p[1]);const cx=(Math.min(...xs)+Math.max(...xs))/2,cy=(Math.min(...ys)+Math.max(...ys))/2;const scale=Math.min(680/Math.max(2,Math.max(...xs)-Math.min(...xs)),480/Math.max(2,Math.max(...ys)-Math.min(...ys)));view={cx,cy,scale};render();}
function screen(p){return [400+(p[0]-view.cx)*view.scale,300-(p[1]-view.cy)*view.scale];}
function position(e){const pt=$('canvas').createSVGPoint();pt.x=e.clientX;pt.y=e.clientY;const s=pt.matrixTransform($('canvas').getScreenCTM().inverse());return [view.cx+(s.x-400)/view.scale,view.cy-(s.y-300)/view.scale];}
function svg(tag,attrs){const el=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v]of Object.entries(attrs))el.setAttribute(k,v);$('canvas').append(el);return el;}
function line(a,b,color,width=2,dash=''){const [x1,y1]=screen(a),[x2,y2]=screen(b);svg('line',{x1,y1,x2,y2,stroke:color,'stroke-width':width,'stroke-dasharray':dash,'pointer-events':'none'});}
function fields(){const p=points[selected];$('selected').textContent=p?`Dot ${selected}`:'Tap or click a dot to select it.';$('x').value=p?p[0]:'';$('y').value=p?p[1]:'';}
function render(){
 if(!view)return;const tol=tolerance(),m=measure(points,tol),o=measure(original,tol);$('canvas').replaceChildren();
 if($('ghost').checked){for(let i=0;i<original.length;i++){const [cx,cy]=screen(original[i]);svg('circle',{cx,cy,r:9,fill:'none',stroke:'#a9b9c5','stroke-dasharray':'3 3'});if(Math.hypot(...points[i].map((v,j)=>v-original[i][j]))>1e-10)line(original[i],points[i],'#b9c4cc',1,'4 4');}}
 for(const tri of m.triangles)svg('polygon',{points:tri.map(i=>screen(points[i]).join(',')).join(' '),fill:'#a171d5','fill-opacity':'.18',stroke:'#8051b0','stroke-width':2,'pointer-events':'none'});
 for(const [i,j]of m.edges)line(points[i],points[j],'#176b91');
 for(const [i,j]of m.violations)line(points[i],points[j],'#d72e42',4);
 const bad=new Set(m.violations.flat());points.forEach((p,i)=>{const [cx,cy]=screen(p);svg('circle',{cx,cy,r:i===selected?10:8,fill:bad.has(i)?'#d72e42':'#18344b',class:'point','data-index':i,tabindex:0,role:'button','aria-label':`Point ${i}: ${p[0].toFixed(4)}, ${p[1].toFixed(4)}`});svg('circle',{cx,cy,r:24,fill:'transparent','data-index':i,class:'touch-target','aria-hidden':'true'});const t=svg('text',{x:cx+12,y:cy-12,'font-size':13,fill:'#18344b','pointer-events':'none'});t.textContent=i;});
 $('metrics').replaceChildren();for(const [label,a,b,count]of [['Closest pair',o.min,m.min],['Furthest pair',o.max,m.max],['Triangles with 3 sides about 1',o.triangles.length,m.triangles.length,true],['Pairs too close (below 1)',o.violations.length,m.violations.length,true]]){const tr=document.createElement('tr');for(const v of [label,count?a:a.toFixed(9),count?b:b.toFixed(9),(b-a>=0?'+':'')+(count?b-a:(b-a).toFixed(9))]){const td=document.createElement('td');td.textContent=v;tr.append(td);}$('metrics').append(tr);}
 $('status').textContent=m.violations.length?`${m.violations.length} pair(s) are too close. Move those dots apart until the red lines disappear.`:'No dots are too close: every distance is at least 1 (computer check).';$('status').style.color=m.violations.length?'#bd2336':'#206c57';$('undo').disabled=!history.length;const quick=$('quick-stats');quick.textContent=`Closest pair: ${m.min.toFixed(4)} · Furthest pair: ${m.max.toFixed(4)} · Triangles: ${m.triangles.length}`;
}
function start(key){source=key;original=clone(configurations[key]);points=clone(original);selected=-1;history=[];$('name').value='';$('notes').value='';fields();fit();}
for(const key of Object.keys(configurations)){const opt=document.createElement('option');opt.value=key;opt.textContent=key.replace(/^n0?/,'Dots: ').replace('_unrestricted',' · best picture found').replace('_margin_',' · triangles opened by ');$('config').append(opt);}
$('config').onchange=()=>start($('config').value);
$('canvas').addEventListener('pointerdown',e=>{if(drag)return;let i=e.target.getAttribute('data-index');if(e.pointerType==='touch'){const at=position(e);let best=Infinity;points.forEach((p,k)=>{const d=Math.hypot(p[0]-at[0],p[1]-at[1]);if(d<best){best=d;i=k;}});const width=$('canvas').getBoundingClientRect().width;if(best*view.scale*width/800>28)return;}if(i===null)return;e.preventDefault();selected=Number(i);history.push(clone(points));drag={id:e.pointerId};$('canvas').setPointerCapture(e.pointerId);fields();render();});
$('canvas').addEventListener('pointermove',e=>{if(!drag||drag.id!==e.pointerId)return;points[selected]=position(e);fields();render();});
function end(e){if(!drag||!e||drag.id===e.pointerId)drag=null;} $('canvas').addEventListener('pointerup',end);$('canvas').addEventListener('pointercancel',end);$('canvas').addEventListener('lostpointercapture',()=>{drag=null;});
$('canvas').addEventListener('keydown',e=>{const i=e.target.getAttribute('data-index');if(i===null||!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown'].includes(e.key))return;e.preventDefault();selected=Number(i);history.push(clone(points));const step=e.shiftKey?.1:.01;points[selected][0]+=e.key==='ArrowLeft'?-step:e.key==='ArrowRight'?step:0;points[selected][1]+=e.key==='ArrowDown'?-step:e.key==='ArrowUp'?step:0;fields();render();$('canvas').querySelector(`[data-index="${selected}"]`).focus();});
$('apply').onclick=()=>{const p=[Number($('x').value),Number($('y').value)];if(selected<0||!p.every(Number.isFinite)||p.some(v=>Math.abs(v)>1e6)||$('x').value===''||$('y').value===''){notice('Select a point and enter finite coordinates between −1,000,000 and 1,000,000.');return;}history.push(clone(points));points[selected]=p;render();};
$('undo').onclick=()=>{if(history.length){points=history.pop();fields();render();}};
$('reset').onclick=()=>{history.push(clone(points));points=clone(original);fields();fit();};$('fit').onclick=fit;$('ghost').onchange=render;$('tol').oninput=render;
function notice(s){$('notice').textContent=s;}
function valid(data){return data&&data.version===1&&typeof data.source==='string'&&configurations[data.source]&&Array.isArray(data.points)&&data.points.length===configurations[data.source].length&&data.points.every(p=>Array.isArray(p)&&p.length===2&&p.every(v=>typeof v==='number'&&Number.isFinite(v)&&Math.abs(v)<=1e6))&&typeof data.name==='string'&&typeof data.notes==='string';}
function capture(){return {version:1,source,points:clone(points),name:$('name').value.trim()||'Untitled experiment',notes:$('notes').value,created:new Date().toISOString(),tolerance:tolerance()};}
function savedList(){ $('saved').replaceChildren();experiments.forEach((e,i)=>{const opt=document.createElement('option');opt.value=i;opt.textContent=e.name;$('saved').append(opt);}); }
try{const items=JSON.parse(localStorage.getItem(storageKey)||'[]');experiments=Array.isArray(items)?items.filter(valid):[];}catch{notice('Browser storage unavailable. Use JSON downloads to keep your work.');}savedList();
$('save').onclick=()=>{const entry=capture(),next=[...experiments,entry];try{localStorage.setItem(storageKey,JSON.stringify(next));experiments=next;savedList();$('saved').value=experiments.length-1;notice('Your idea is saved on this device. Download a backup file to keep a copy or use it on your phone.');}catch{notice('This browser could not save your idea. Download a backup file to keep it.');}};
function restore(e){if(!valid(e))throw Error('Unsupported experiment format or invalid coordinates.');$('config').value=e.source;start(e.source);points=clone(e.points);$('name').value=e.name;$('notes').value=e.notes;if(Number.isFinite(e.tolerance)&&e.tolerance>=0&&e.tolerance<=.1)$('tol').value=e.tolerance;fields();fit();notice('Your idea is open. The grey rings and starting measurements come from its original picture.');}
$('load').onclick=()=>{const e=experiments[Number($('saved').value)];if(e)restore(e);};
$('export').onclick=()=>{const blob=new Blob([JSON.stringify(capture(),null,2)],{type:'application/json'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='erdos99-experiment.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
$('import').onchange=async()=>{try{const f=$('import').files[0];if(!f)return;if(f.size>1000000)throw Error('File is too large.');restore(JSON.parse(await f.text()));}catch(e){notice(`Import failed: ${e.message}`);}finally{$('import').value='';}};
start($('config').value);

const touchMode=matchMedia('(pointer: coarse)');function inputHint(){document.documentElement.classList.toggle('touch-device',touchMode.matches);$('input-hint').textContent=touchMode.matches?'Touch a dot and slide your finger to move it. Use the page edges to scroll.':'Drag a dot with your mouse. Select a dot and use arrow keys for small moves (Shift for larger moves).';}inputHint();touchMode.addEventListener('change',inputHint);
