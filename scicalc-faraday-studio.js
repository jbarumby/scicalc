/* SciCalc — Faraday realistic-looking laboratory.
   Canvas renders an explanatory apparatus; numerical values come from the
   separate simplified physical model and are NOT experimental measurements. */
(function(root){
'use strict';
const old=document.getElementById('fl-canvas'),M=root.SciCalcFaradayModel;
if(!old||!M||document.getElementById('faraday-studio-canvas'))return;
const card=old.closest('.v26-extra-card');
const $=id=>document.getElementById(id);
if(!card||!$('fl-n')||!$('fl-speed'))return;
const studio=document.createElement('div');
studio.className='faraday-studio';
studio.innerHTML=[
 '<header class="faraday-head"><div><div class="faraday-kicker">Física experimental · indução eletromagnética</div>',
 '<h4>Faraday — ímã e solenoide</h4><p>Movimente o ímã e observe a variação do fluxo e a tensão induzida na bobina.</p></div>',
 '<div class="faraday-status" id="faraday-live" data-moving="false"><span class="faraday-status-dot"></span>',
 '<span id="faraday-live-text">APARATO EM REPOUSO</span></div></header>',
 '<div class="faraday-display"><canvas id="faraday-studio-canvas" role="img" aria-label="Ímã com polo N voltado ao solenoide de cobre; linhas de campo e galvanômetro."></canvas>',
 '<div class="faraday-legend"><span><i style="background:#fb7185"></i>Polo N</span>',
 '<span><i style="background:#60a5fa"></i>Polo S</span>',
 '<span><i style="background:#38bdf8"></i>Linhas de campo esquemáticas</span>',
 '<span><i style="background:#4ade80"></i>Campo induzido</span></div></div>',
 '<div class="faraday-data">',
 '<div class="faraday-data-item"><small>Campo axial na bobina</small><strong id="faraday-field">—</strong><span>B(d) (modelo aproximado)</span></div>',
 '<div class="faraday-data-item"><small>Fluxo por espira</small><strong id="faraday-flux">—</strong><span>Φ = B·A, em Wb</span></div>',
 '<div class="faraday-data-item"><small>FEM induzida</small><strong id="faraday-emf">0 V</strong><span id="faraday-emf-note">Sem variação de fluxo</span></div>',
 '<div class="faraday-data-item"><small>Distância ímã-bobina</small><strong id="faraday-distance">—</strong><span id="faraday-distance-note">N espiras · movimento axial</span></div>',
 '</div>',
 '<div class="faraday-explain"><h5>📚 O que está acontecendo?</h5>',
 '<div class="faraday-explain-grid">',
 '<div><b>1. Movimento e campo</b><p id="faraday-step-1">O ímã está parado. A bobina ainda recebe campo magnético.</p></div>',
 '<div><b>2. Fluxo magnético</b><p id="faraday-step-2">Sem movimento relativo, o fluxo permanece constante.</p></div>',
 '<div><b>3. Indução e Lei de Lenz</b><p id="faraday-step-3">Sem variação de fluxo não há FEM induzida.</p></div></div>',
 '<p class="faraday-note">Modelo didático axial simplificado: B(d) = B₀/[1+(d/a)²]^(3/2); ',
 'a é o maior entre o raio equivalente da bobina e 0,045 m. ',
 'Aproximamos o campo como uniforme na área A. ',
 'FEM = −N·dΦ/dt. A escala do desenho não é uma representação em escala real; ',
 'o galvanômetro indica apenas deflexão qualitativa proporcional à FEM, ',
 'pois a resistência do circuito não foi especificada.</p></div>',
 '<div class="faraday-caption">Reprodução visual em câmera lenta (0,55×); velocidade física ajustada pelos controles.</div>'
].join('');
old.insertAdjacentElement('beforebegin',studio);
const canvas=$('faraday-studio-canvas'),ctx=canvas.getContext('2d');
if(!ctx){studio.remove();return;}
card.classList.add('faraday-studio-ready');

let direction='approach',distance=.67,active=false,lastFrame=null,frameId=0,elapsed=0,parameters=null,lastInfo=0;
function read(){
 const get=(id)=>Number(String($(id).value).replace(',','.'));
 return M.parameters(get('fl-n'),get('fl-area'),get('fl-b0'),get('fl-speed'));
}
function sample(){
 if(!parameters){try{parameters=read();}catch(e){return null;}}
 return M.state(parameters,distance,direction,active);
}
function label(str,x,y,color='#e2e8f0',font='700 15px system-ui',align='left'){
 ctx.font=font;ctx.fillStyle=color;ctx.textAlign=align;ctx.fillText(str,x,y);ctx.textAlign='left';
}
function rounded(x,y,w,h,r){
 r=Math.min(r,w/2,h/2);
 ctx.beginPath();ctx.moveTo(x+r,y);ctx.lineTo(x+w-r,y);ctx.quadraticCurveTo(x+w,y,x+w,y+r);
 ctx.lineTo(x+w,y+h-r);ctx.quadraticCurveTo(x+w,y+h,x+w-r,y+h);
 ctx.lineTo(x+r,y+h);ctx.quadraticCurveTo(x,y+h,x,y+h-r);
 ctx.lineTo(x,y+r);ctx.quadraticCurveTo(x,y,x+r,y);ctx.closePath();
}
function arrow(x1,y1,x2,y2,color,width=4){
 const angle=Math.atan2(y2-y1,x2-x1);
 ctx.strokeStyle=color;ctx.fillStyle=color;ctx.lineWidth=width;ctx.lineCap='round';
 ctx.beginPath();ctx.moveTo(x1,y1);ctx.lineTo(x2,y2);ctx.stroke();
 ctx.beginPath();ctx.moveTo(x2,y2);
 ctx.lineTo(x2-12*Math.cos(angle-.5),y2-12*Math.sin(angle-.5));
 ctx.lineTo(x2-12*Math.cos(angle+.5),y2-12*Math.sin(angle+.5));ctx.closePath();ctx.fill();
}
function background(){
 const grd=ctx.createLinearGradient(0,0,1000,440);
 grd.addColorStop(0,'#0b1930');grd.addColorStop(.7,'#122849');grd.addColorStop(1,'#09172b');
 ctx.fillStyle=grd;ctx.fillRect(0,0,1000,440);
 ctx.strokeStyle='rgba(130,162,204,.09)';ctx.lineWidth=1;
 for(let x=20;x<1000;x+=40){ctx.beginPath();ctx.moveTo(x,0);ctx.lineTo(x,440);ctx.stroke();}
 for(let y=20;y<440;y+=40){ctx.beginPath();ctx.moveTo(0,y);ctx.lineTo(1000,y);ctx.stroke();}
 const glow=ctx.createRadialGradient(650,225,5,650,225,310);
 glow.addColorStop(0,'rgba(8,145,178,.20)');glow.addColorStop(1,'rgba(8,145,178,0)');
 ctx.fillStyle=glow;ctx.fillRect(370,0,550,440);
 ctx.setLineDash([6,9]);ctx.strokeStyle='rgba(160,187,220,.28)';
 ctx.beginPath();ctx.moveTo(35,225);ctx.lineTo(790,225);ctx.stroke();ctx.setLineDash([]);
}
function fieldLines(tip,m){
 const gap=560-tip;
 ctx.lineWidth=2;
 const alpha=Math.min(.78,.15+(.55*m.B/parameters.B0));
 for(const sign of [-1,1]){
  ctx.strokeStyle='rgba(56,189,248,'+alpha+')';
  ctx.beginPath();
  ctx.moveTo(tip,215+sign*22);
  ctx.bezierCurveTo(tip+gap*.35,215+sign*100,560-gap*.16,215+sign*105,565,215+sign*67);
  ctx.stroke();
 }
 if(gap>18)arrow(tip+gap*.20,225,tip+gap*.70,225,'#38bdf8',3);
 label('B do ímã →',390,89,'#7dd3fc','700 15px system-ui');
}
function magnet(px){
 ctx.save();ctx.shadowBlur=19;ctx.shadowColor='rgba(0,0,0,.5)';ctx.shadowOffsetY=8;
 rounded(px,184,153,81,12);ctx.fillStyle='#2563eb';ctx.fill();ctx.restore();
 ctx.save();rounded(px,184,153,81,12);ctx.clip();
 const blue=ctx.createLinearGradient(0,184,0,265);
 blue.addColorStop(0,'#a6d8ff');blue.addColorStop(.32,'#3b82f6');blue.addColorStop(1,'#16377d');
 ctx.fillStyle=blue;ctx.fillRect(px,184,77,81);
 const red=ctx.createLinearGradient(0,184,0,265);
 red.addColorStop(0,'#fecdd3');red.addColorStop(.33,'#ed5372');red.addColorStop(1,'#8e1d39');
 ctx.fillStyle=red;ctx.fillRect(px+77,184,77,81);
 ctx.restore();
 ctx.strokeStyle='rgba(255,255,255,.6)';ctx.lineWidth=2;rounded(px,184,153,81,12);ctx.stroke();
 label('S',px+38,236,'white','800 29px system-ui','center');
 label('N',px+115,236,'white','800 29px system-ui','center');
 label('ÍMÃ',px+76,303,'#c3d8f1','700 13px system-ui','center');
}
function coil(){
 const x=585,top=155, h=142,end=750;
 ctx.fillStyle='rgba(148,163,184,.09)';ctx.strokeStyle='rgba(190,202,222,.25)';
 ctx.lineWidth=2;rounded(x,top,end-x,h,16);ctx.fill();ctx.stroke();
 for(let i=0;i<15;i++){
  const cx=x+i*11;
  ctx.lineWidth=4.5;
  const copper=ctx.createLinearGradient(0,top,0,top+h);
  copper.addColorStop(0,'#fff2c6');copper.addColorStop(.3,'#ffbb5c');copper.addColorStop(.7,'#9b4b19');copper.addColorStop(1,'#eaa74c');
  ctx.strokeStyle=copper;
  ctx.beginPath();ctx.ellipse(cx,226,13,80,0,0,2*Math.PI);ctx.stroke();
 }
 ctx.fillStyle='rgba(42,68,93,.36)';
 ctx.beginPath();ctx.ellipse(585,226,13,79,0,0,2*Math.PI);ctx.fill();
 ctx.strokeStyle='#f2c579';ctx.lineWidth=3;ctx.beginPath();ctx.ellipse(585,226,13,79,0,0,2*Math.PI);ctx.stroke();
 label('SOLENOIDE · COBRE',670,96,'#b6d3ea','700 14px system-ui','center');
}
function wires(){
 ctx.strokeStyle='#bb8e52';ctx.lineWidth=3;ctx.beginPath();
 ctx.moveTo(744,158);ctx.bezierCurveTo(777,105,809,119,844,137);
 ctx.moveTo(744,291);ctx.bezierCurveTo(777,339,829,327,844,308);
 ctx.stroke();
}
function gauge(m){
 const cx=862,cy=215,r=62;
 const steel=ctx.createRadialGradient(cx-10,cy-12,6,cx,cy,75);
 steel.addColorStop(0,'#244562');steel.addColorStop(1,'#071629');
 ctx.fillStyle=steel;ctx.beginPath();ctx.arc(cx,cy,r+8,0,Math.PI*2);ctx.fill();
 ctx.strokeStyle='#8b9fb7';ctx.lineWidth=3;ctx.beginPath();ctx.arc(cx,cy,r+8,0,2*Math.PI);ctx.stroke();
 ctx.strokeStyle='#c0d1e4';ctx.lineWidth=2;ctx.beginPath();ctx.arc(cx,cy,r-11,-2.8,-.35,false);ctx.stroke();
 for(let i=-4;i<=4;i++){
  const a=-Math.PI/2+i*.28;
  const c1=38,c2=46;ctx.beginPath();ctx.moveTo(cx+c1*Math.cos(a),cy+c1*Math.sin(a));
  ctx.lineTo(cx+c2*Math.cos(a),cy+c2*Math.sin(a));ctx.stroke();
 }
 let level=0;
 if(m.moving)level=Math.max(-1,Math.min(1,m.emf/.45));
 const ang=-Math.PI/2+level*1.12;
 arrow(cx,cy,cx+42*Math.cos(ang),cy+42*Math.sin(ang),'#fb7185',3);
 ctx.fillStyle='#e2e8f0';ctx.beginPath();ctx.arc(cx,cy,5,0,2*Math.PI);ctx.fill();
 label('−',808,242,'#d4e5f5','700 16px system-ui');
 label('+',914,242,'#d4e5f5','700 16px system-ui');
 label('GALVANÔMETRO',862,99,'#c1d4ee','700 12px system-ui','center');
 label('deflexão relativa',862,308,'#a8c2df','600 12px system-ui','center');
}
function induced(m){
 if(!m.moving){label('Sem FEM induzida',651,365,'#93a8c2','700 14px system-ui','center');return;}
 const left=m.fieldInduced==='left';
 arrow(left?731:604,226,left?625:729,226,'#4ade80',5);
 label('B induzido '+(left?'←':'→'),660,360,'#6ee7b7','700 15px system-ui','center');
 label('Face '+m.inducedFace,586,328,m.inducedFace==='N'?'#fda4af':'#93c5fd','700 14px system-ui','center');
}
function distanceRuler(tip,m){
 const y=380,x0=Math.max(100,Math.min(561,tip)),x1=582;
 ctx.strokeStyle='#8aa3c2';ctx.lineWidth=1.5;ctx.setLineDash([3,5]);
 ctx.beginPath();ctx.moveTo(x0,y);ctx.lineTo(x1,y);ctx.stroke();ctx.setLineDash([]);
 ctx.beginPath();ctx.moveTo(x0,y-6);ctx.lineTo(x0,y+6);ctx.moveTo(x1,y-6);ctx.lineTo(x1,y+6);ctx.stroke();
 label('d = '+m.distance.toFixed(3).replace('.',',')+' m',(x0+x1)/2,y-11,'#bdd1e8','700 12px system-ui','center');
}
function draw(){
 const rect=canvas.getBoundingClientRect(),dpr=Math.min(2,root.devicePixelRatio||1);
 const cw=Math.round(Math.max(320,rect.width||1000)*dpr);
 const ch=Math.round(Math.max(220,rect.height||440)*dpr);
 if(canvas.width!==cw||canvas.height!==ch){canvas.width=cw;canvas.height=ch;}
 ctx.setTransform(cw/1000,0,0,ch/440,0,0);
 const m=sample();if(!m)return;
 const px=383-(m.distance-.09)/(.71-.09)*286,tip=px+153;
 background();fieldLines(tip,m);wires();coil();magnet(px);induced(m);gauge(m);distanceRuler(tip,m);
 if(m.moving){
  const approaching=direction==='approach';
  label(approaching?'MOVIMENTO →':'MOVIMENTO ←',70,83,approaching?'#86efac':'#fdba74','800 16px system-ui');
  arrow(approaching?305:376,324,approaching?372:309,324,approaching?'#4ade80':'#fb923c',5);
 }else label('ÍMÃ PARADO',70,83,'#b8cce5','700 15px system-ui');
 label('ε = −N dΦ/dt',32,423,'#bcd4ea','700 17px system-ui');
}
const asNumber=(n,digits=3)=>Number(n).toLocaleString('pt-BR',{maximumSignificantDigits:digits});
function update(){
 const m=sample();if(!m)return;
 $('faraday-live').dataset.moving=String(m.moving);
 $('faraday-live-text').textContent=m.moving
   ?direction==='approach'?'APROXIMANDO →':'← AFASTANDO'
   :'APARATO EM REPOUSO';
 $('faraday-field').textContent=asNumber(m.B,5)+' T';
 $('faraday-flux').textContent=asNumber(m.flux,5)+' Wb';
 $('faraday-emf').textContent=asNumber(m.emf,5)+' V';
 $('faraday-emf-note').textContent=m.moving?'Sinal da FEM depende do sentido do movimento':'Sem variação do fluxo';
 $('faraday-distance').textContent=asNumber(m.distance,4)+' m';
 $('faraday-distance-note').textContent=parameters.N+' espiras · '+parameters.speed.toLocaleString('pt-BR')+' m/s';
 $('faraday-step-1').textContent=m.moving?
   direction==='approach'?'O polo N aproxima-se da bobina e o campo axial na bobina aumenta.':
     'O polo N afasta-se da bobina e o campo axial na bobina diminui.':
   'O ímã está parado. Há campo magnético na bobina, mas ele permanece constante.';
 $('faraday-step-2').textContent=m.moving
   ?direction==='approach'?'O fluxo Φ = B·A aumenta no sentido do eixo da bobina.':'O fluxo Φ = B·A diminui ao longo do afastamento.'
   :'Não há variação temporal do fluxo por espira (dΦ/dt = 0).';
 $('faraday-step-3').textContent=m.moving
   ?direction==='approach'?'A FEM é negativa na convenção adotada; a face da bobina voltada ao ímã comporta-se como polo N.':
    'A FEM troca de sinal; a face da bobina voltada ao ímã comporta-se como polo S.'
   :'A FEM induzida vale exatamente zero neste modelo ideal.';
}
function render(){draw();update();}
function stop(){
 active=false;lastFrame=null;
 if(frameId)root.cancelAnimationFrame(frameId);
 frameId=0;render();
}
function reset(){
 stop();distance=.67;direction='approach';elapsed=0;
 try{parameters=read();}catch(e){displayError(e);}
 render();
}
function displayError(e){
 const target=$('fl-theory');if(target){
  target.textContent='⚠️ '+String(e.message||e);
  target.setAttribute('role','alert');
 }
}
function start(mode){
 try{
  parameters=read();
  if(frameId)root.cancelAnimationFrame(frameId);
  direction=mode==='away'?'away':'approach';
  distance=M.start(direction);elapsed=0;lastFrame=null;lastInfo=0;
  active=true;render();
  frameId=root.requestAnimationFrame(tick);
 }catch(e){displayError(e);}
}
function tick(now){
 if(!active)return;
 const dt=lastFrame===null?0:Math.min(.06,Math.max(0,(now-lastFrame)/1000));
 lastFrame=now;
 // Playback at 0.55x wall-clock for a clearer, slower educational animation.
 const step=M.advance(distance,direction,parameters.speed,dt*.55);
 distance=step.distance;elapsed+=dt*.55;
 if(step.finished){stop();return;}
 draw();
 if(now-lastInfo>100){update();lastInfo=now;}
 frameId=root.requestAnimationFrame(tick);
}
root.flStart=start;root.flStop=stop;root.flReset=reset;root.flInit=render;
['fl-n','fl-area','fl-b0','fl-speed'].forEach(id=>$(id).addEventListener('change',()=>{
 if(active)return;
 try{parameters=read();render();}catch(e){displayError(e);}
}));
const tab=document.querySelector('nav.tabs [data-tab="simuladores-v26"]');
if(tab)tab.addEventListener('click',()=>root.requestAnimationFrame(render));
if(root.ResizeObserver){
 const observer=new ResizeObserver(()=>{if(document.getElementById('simuladores-v26')?.classList.contains('active'))draw();});
 observer.observe(canvas);
}
document.addEventListener('visibilitychange',()=>{if(document.hidden&&active)stop();});
try{parameters=read();render();}catch(e){displayError(e);}
})(window);
