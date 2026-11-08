/* SciCalc — Estúdio da Lei de Lenz.
   Renderizador alternativo não destrutivo: o canvas legado permanece no DOM,
   mas é ocultado após a inicialização bem-sucedida desta interface. */
(function (root) {
  'use strict';
  const original = document.getElementById('fl-lenz-canvas');
  const readModel = root.scicalcLenzState;
  if (!original || typeof readModel !== 'function' || document.getElementById('lenz-studio-canvas')) return;
  const card = original.closest('.v26-extra-card');
  const speed = document.getElementById('fl-lenz-speed');
  const modeSelect = document.getElementById('fl-lenz-mode');
  const result = document.getElementById('fl-lenz-result');
  if (!card || !speed || !modeSelect) return;

  const studio = document.createElement('div');
  studio.className = 'lenz-studio';
  studio.innerHTML = [
    '<header class="lenz-studio-header">',
      '<div><div class="lenz-studio-kicker">Laboratório interativo · Eletromagnetismo</div>',
      '<h4>Lei de Lenz — campo e corrente</h4>',
      '<p>Observe como a bobina reage à aproximação e ao afastamento do polo norte.</p></div>',
      '<div class="lenz-status-pill" id="lenz-motion-pill" data-moving="false">',
      '<span class="lenz-status-dot" aria-hidden="true"></span><span id="lenz-motion-text">ÍMÃ PARADO</span></div>',
    '</header>',
    '<div class="lenz-diagram-area">',
      '<canvas id="lenz-studio-canvas" class="lenz-display-canvas" role="img"',
       ' aria-label="Ímã com polo norte voltado para uma bobina; setas mostram os campos e um galvanômetro mostra a resposta induzida."></canvas>',
      '<div class="lenz-diagram-legend">',
        '<span><i class="lenz-swatch" style="background:#fb7185"></i>Polo norte (N)</span>',
        '<span><i class="lenz-swatch" style="background:#60a5fa"></i>Polo sul (S)</span>',
        '<span><i class="lenz-swatch" style="background:#38bdf8"></i>Campo do ímã</span>',
        '<span><i class="lenz-swatch" style="background:#34d399"></i>Campo induzido</span>',
      '</div>',
    '</div>',
    '<div class="lenz-readouts" aria-live="off">',
      '<div class="lenz-readout"><small>Movimento</small><strong id="lenz-r-movement">Parado</strong><span id="lenz-r-flow">Fluxo constante</span></div>',
      '<div class="lenz-readout"><small>Face induzida da bobina</small><strong id="lenz-r-pole">—</strong><span>Voltada para o ímã</span></div>',
      '<div class="lenz-readout"><small>Corrente induzida</small><strong id="lenz-r-current">Nenhuma</strong><span>Vista do lado do ímã</span></div>',
      '<div class="lenz-readout"><small>Resposta induzida relativa</small><strong id="lenz-r-emf">0%</strong><span>Modelo didático (não calibrado)</span></div>',
    '</div>',
    '<div class="lenz-info">',
      '<h5>📚 Entenda o que acontece em três passos</h5>',
      '<div class="lenz-steps">',
        '<div class="lenz-step"><span class="num">1</span><b>Movimento do ímã</b><p id="lenz-step1">Escolha aproximar ou afastar o ímã.</p></div>',
        '<div class="lenz-step"><span class="num">2</span><b>Variação do fluxo</b><p id="lenz-step2">Com o ímã parado, o fluxo não varia.</p></div>',
        '<div class="lenz-step"><span class="num">3</span><b>Resposta da bobina</b><p id="lenz-step3">Sem variação de fluxo, não há corrente induzida.</p></div>',
      '</div>',
      '<p class="lenz-footnote"><b>Lei de Faraday–Lenz:</b> ε = −N · dΦ/dt. ',
      'O sinal negativo expressa a oposição à <em>variação</em> do fluxo, não necessariamente ao campo do ímã. ',
      'As intensidades e a velocidade visual são relativas e qualitativas, não valores medidos de corrente, tensão ou campo.</p>',
    '</div>'
  ].join('');
  original.insertAdjacentElement('beforebegin', studio);
  const canvas = document.getElementById('lenz-studio-canvas');
  const ctx = canvas.getContext('2d');
  if (!ctx) { studio.remove(); return; }
  card.classList.add('lenz-modernized');

  const get = id => document.getElementById(id);
  let direction = modeSelect.value === 'away' ? 'away' : 'approach';
  let p = 0.12, playing = false, lastFrame = null, frameId = 0, lastInfoTime = 0;

  const clamp = (v, min, max) => Math.max(min, Math.min(max, v));
  function currentSpeed() {
    const v = Number(String(speed.value).replace(',', '.'));
    return Number.isFinite(v) ? clamp(v, 0.05, 1) : 0.35;
  }
  function physicalState() {
    return readModel(direction, playing, p, currentSpeed());
  }
  function roundedRect(x,y,w,h,r) {
    const rad=Math.min(r,w/2,h/2);
    ctx.beginPath();
    ctx.moveTo(x+rad,y);ctx.lineTo(x+w-rad,y);
    ctx.quadraticCurveTo(x+w,y,x+w,y+rad);
    ctx.lineTo(x+w,y+h-rad);ctx.quadraticCurveTo(x+w,y+h,x+w-rad,y+h);
    ctx.lineTo(x+rad,y+h);ctx.quadraticCurveTo(x,y+h,x,y+h-rad);
    ctx.lineTo(x,y+rad);ctx.quadraticCurveTo(x,y,x+rad,y);
    ctx.closePath();
  }
  function text(content,x,y,color='#e2e8f0',font='600 15px system-ui',align='left') {
    ctx.fillStyle=color;ctx.font=font;ctx.textAlign=align;
    ctx.fillText(content,x,y);ctx.textAlign='left';
  }
  function arrow(x1,y1,x2,y2,color,width=4) {
    const angle=Math.atan2(y2-y1,x2-x1),head=12;
    ctx.strokeStyle=color;ctx.fillStyle=color;ctx.lineWidth=width;ctx.lineCap='round';
    ctx.beginPath();ctx.moveTo(x1,y1);ctx.lineTo(x2,y2);ctx.stroke();
    ctx.beginPath();ctx.moveTo(x2,y2);
    ctx.lineTo(x2-head*Math.cos(angle-.52),y2-head*Math.sin(angle-.52));
    ctx.lineTo(x2-head*Math.cos(angle+.52),y2-head*Math.sin(angle+.52));
    ctx.closePath();ctx.fill();
  }
  function gradientTop() {
    const bg=ctx.createLinearGradient(0,0,960,410);
    bg.addColorStop(0,'#0a182c');bg.addColorStop(.6,'#101e35');bg.addColorStop(1,'#09172b');
    ctx.fillStyle=bg;ctx.fillRect(0,0,960,410);
    ctx.strokeStyle='rgba(91,134,178,.13)';ctx.lineWidth=1;
    for(let x=30;x<960;x+=40){ctx.beginPath();ctx.moveTo(x,0);ctx.lineTo(x,410);ctx.stroke();}
    for(let y=35;y<410;y+=40){ctx.beginPath();ctx.moveTo(0,y);ctx.lineTo(960,y);ctx.stroke();}
    const glow=ctx.createRadialGradient(655,210,20,650,210,270);
    glow.addColorStop(0,'rgba(14,165,233,.17)');glow.addColorStop(1,'rgba(14,165,233,0)');
    ctx.fillStyle=glow;ctx.fillRect(390,0,490,410);
    ctx.setLineDash([7,9]);ctx.strokeStyle='rgba(148,163,184,.32)';
    ctx.beginPath();ctx.moveTo(45,210);ctx.lineTo(775,210);ctx.stroke();ctx.setLineDash([]);
  }
  function drawFieldLines(m, px) {
    const tip=px+140,coil=585, gap=coil-tip;
    ctx.lineWidth=2.5;
    const alpha=clamp(.17+.55*m.relativeFlux,0,.72);
    ctx.strokeStyle='rgba(56,189,248,'+alpha+')';
    for(let sign of [-1,1]){
      ctx.beginPath();
      ctx.moveTo(tip,204+sign*13);
      ctx.bezierCurveTo(tip+gap*.28,204+sign*65,coil-gap*.12,204+sign*88,coil,204+sign*62);
      ctx.stroke();
    }
    if(gap>20){
      arrow(tip+Math.min(gap*.3,35),205,Math.min(coil-12,tip+Math.max(18,gap*.55)),205,'rgba(56,189,248,.87)',3);
    }
    text('B do ímã  →',358,101,'#7dd3fc','700 15px system-ui');
  }
  function drawMagnet(px) {
    ctx.save();ctx.shadowColor='rgba(0,0,0,.45)';ctx.shadowBlur=18;ctx.shadowOffsetY=9;
    const back=ctx.createLinearGradient(px,150,px,266);
    back.addColorStop(0,'#8fcbff');back.addColorStop(.42,'#2563eb');back.addColorStop(1,'#1d3b83');
    roundedRect(px,165,142,79,12);ctx.fillStyle=back;ctx.fill();ctx.restore();
    ctx.save();roundedRect(px,165,142,79,12);ctx.clip();
    ctx.fillStyle='#e34d68';ctx.fillRect(px+71,165,71,79);
    const gloss=ctx.createLinearGradient(0,165,0,244);
    gloss.addColorStop(0,'rgba(255,255,255,.36)');
    gloss.addColorStop(.52,'rgba(255,255,255,.015)');
    gloss.addColorStop(1,'rgba(0,0,0,.24)');
    ctx.fillStyle=gloss;ctx.fillRect(px,165,142,79);ctx.restore();
    ctx.strokeStyle='rgba(255,255,255,.44)';ctx.lineWidth=1.5;
    roundedRect(px,165,142,79,12);ctx.stroke();
    ctx.strokeStyle='rgba(255,255,255,.55)';ctx.beginPath();
    ctx.moveTo(px+71,167);ctx.lineTo(px+71,244);ctx.stroke();
    text('S',px+36,215,'white','800 28px system-ui','center');
    text('N',px+106,215,'white','800 28px system-ui','center');
    text('ÍMÃ PERMANENTE',px+71,286,'#9fb8d9','700 12px system-ui','center');
  }
  function drawCoil(m) {
    const coilX=600, coilEnd=744, top=150, height=117;
    ctx.fillStyle='rgba(148,163,184,.13)';
    roundedRect(coilX,top,coilEnd-coilX,height,17);ctx.fill();
    ctx.strokeStyle='rgba(148,163,184,.36)';ctx.lineWidth=2;
    roundedRect(coilX,top,coilEnd-coilX,height,17);ctx.stroke();
    // Espiras visíveis; um solenoide com eixo horizontal.
    for(let i=0;i<13;i++){
      const cx=coilX+4+i*11;
      const copper=ctx.createLinearGradient(cx,top,cx,top+height);
      copper.addColorStop(0,'#fce6a4');copper.addColorStop(.38,'#edab4a');
      copper.addColorStop(.78,'#aa6127');copper.addColorStop(1,'#f7c66c');
      ctx.strokeStyle=copper;ctx.lineWidth=5;
      ctx.beginPath();ctx.ellipse(cx,208,13,68,0,0,2*Math.PI);ctx.stroke();
    }
    ctx.fillStyle='rgba(203,213,225,.16)';
    ctx.beginPath();ctx.ellipse(coilX,208,14,68,0,0,2*Math.PI);ctx.fill();
    ctx.strokeStyle='#ebc88e';ctx.lineWidth=3;
    ctx.beginPath();ctx.ellipse(coilX,208,14,68,0,0,2*Math.PI);ctx.stroke();
    text('BOBINA / SOLENOIDE',675,104,'#bcd3ec','700 13px system-ui','center');
    if(m.moving){
      const face=m.inducedPole;
      ctx.fillStyle=face==='N'?'rgba(244,63,94,.30)':'rgba(59,130,246,.30)';
      ctx.beginPath();ctx.ellipse(coilX,208,13,63,0,0,2*Math.PI);ctx.fill();
      text('Face '+face,592,310,face==='N'?'#fda4af':'#93c5fd','700 15px system-ui','center');
      if(m.fieldDirection==='left'){
        arrow(723,208,648,208,'#34d399',5);
        text('B induzido  ←',646,342,'#6ee7b7','700 14px system-ui','center');
      }else{
        arrow(637,208,714,208,'#34d399',5);
        text('B induzido  →',646,342,'#6ee7b7','700 14px system-ui','center');
      }
    }else{
      text('Sem campo induzido',660,342,'#95a9c4','600 13px system-ui','center');
    }
  }
  function currentCircle(m) {
    const cx=844,cy=328,r=21;
    ctx.strokeStyle='#344f70';ctx.lineWidth=2;
    ctx.beginPath();ctx.arc(cx,cy,r+11,0,2*Math.PI);ctx.stroke();
    if(!m.moving){
      text('—',cx,cy+7,'#94a3b8','700 22px system-ui','center');
    }else{
      const clockwise=m.currentDirection==='clockwise';
      const start=clockwise?-2.5:1.2, end=clockwise?1.2:-2.5;
      ctx.strokeStyle='#4ade80';ctx.lineWidth=4;ctx.beginPath();
      ctx.arc(cx,cy,r,start,end,!clockwise);ctx.stroke();
      const tangent=end+(clockwise?Math.PI/2:-Math.PI/2);
      const ex=cx+r*Math.cos(end),ey=cy+r*Math.sin(end);
      const back=11;const a=tangent;
      ctx.fillStyle='#4ade80';ctx.beginPath();ctx.moveTo(ex,ey);
      ctx.lineTo(ex-back*Math.cos(a-.45),ey-back*Math.sin(a-.45));
      ctx.lineTo(ex-back*Math.cos(a+.45),ey-back*Math.sin(a+.45));
      ctx.closePath();ctx.fill();
    }
    text('Vista do ímã',844,379,'#9cb4d1','600 12px system-ui','center');
  }
  function gauge(m){
    const cx=844,cy=184,r=54;
    const rim=ctx.createRadialGradient(cx-18,cy-19,3,cx,cy,r+8);
    rim.addColorStop(0,'#233954');rim.addColorStop(1,'#0c1728');
    ctx.fillStyle=rim;ctx.beginPath();ctx.arc(cx,cy,r+9,0,Math.PI*2);ctx.fill();
    ctx.strokeStyle='#58708e';ctx.lineWidth=3;ctx.beginPath();ctx.arc(cx,cy,r+9,0,Math.PI*2);ctx.stroke();
    ctx.strokeStyle='#a9bbd1';ctx.lineWidth=2;ctx.beginPath();
    ctx.arc(cx,cy,r-8,Math.PI*.85,Math.PI*.15,true);ctx.stroke();
    for(let i=-3;i<=3;i++){
      const a=-Math.PI/2+i*.25;
      const x1=cx+36*Math.cos(a),y1=cy+36*Math.sin(a);
      const x2=cx+43*Math.cos(a),y2=cy+43*Math.sin(a);
      ctx.strokeStyle='#94a3b8';ctx.beginPath();ctx.moveTo(x1,y1);ctx.lineTo(x2,y2);ctx.stroke();
    }
    const angle=-Math.PI/2+(m.needle*1.5);
    arrow(cx,cy,cx+38*Math.cos(angle),cy+38*Math.sin(angle),'#fb7185',3);
    ctx.fillStyle='#f8fafc';ctx.beginPath();ctx.arc(cx,cy,5,0,2*Math.PI);ctx.fill();
    text('−',799,202,'#b9cce5','700 17px system-ui');
    text('+',891,202,'#b9cce5','700 17px system-ui');
    text('MEDIDOR DE CORRENTE',844,107,'#bcd3ec','700 11px system-ui','center');
    text('resposta relativa',844,267,'#9cb4d1','600 12px system-ui','center');
  }
  function draw() {
    const box=canvas.getBoundingClientRect(),dpr=Math.min(root.devicePixelRatio||1,2);
    const w=Math.max(280,box.width||960);
    const h=Math.max(180,box.height||410);
    const cw=Math.round(w*dpr),ch=Math.round(h*dpr);
    if(canvas.width!==cw||canvas.height!==ch){canvas.width=cw;canvas.height=ch;}
    ctx.setTransform(canvas.width/960,0,0,canvas.height/410,0,0);
    const m=physicalState();
    const px=83+p*359;
    gradientTop();
    drawFieldLines(m,px);
    drawCoil(m);
    drawMagnet(px);
    gauge(m);
    currentCircle(m);
    if(m.moving){
      const directionArrow=m.mode==='approach'?'→':'←';
      text('MOVIMENTO '+directionArrow,83,86,m.mode==='approach'?'#86efac':'#fdba74','800 15px system-ui');
      arrow(m.mode==='approach'?355:410,305,m.mode==='approach'?414:351,305,
        m.mode==='approach'?'#4ade80':'#fb923c',5);
    }else{
      text('ÍMÃ PARADO',83,86,'#b2c3d9','700 14px system-ui');
    }
    text('ε = −N · dΦ/dt',30,390,'#d3e1f2','700 17px system-ui');
  }
  function updateInfo() {
    const m=physicalState(),moving=m.moving;
    get('lenz-motion-pill').dataset.moving=String(moving);
    get('lenz-motion-text').textContent=moving?(direction==='approach'?'APROXIMANDO → '+Math.round(100*(p-.12)/.82)+'%':'← AFASTANDO '+Math.round(100*(.88-p)/.82)+'%'):'ÍMÃ PARADO';
    get('lenz-r-movement').textContent=moving?(direction==='approach'?'Aproximando':'Afastando'):'Parado';
    get('lenz-r-flow').textContent=m.fluxTrend==='increases'?'Fluxo aumentando':m.fluxTrend==='decreases'?'Fluxo diminuindo':'Fluxo constante';
    get('lenz-r-pole').textContent=m.inducedPole==='—'?'Nenhuma':m.inducedPole+' induzido';
    get('lenz-r-current').textContent=m.currentDirection==='none'?'Nenhuma':m.currentDirection==='clockwise'?'Horária ↻':'Anti-horária ↺';
    get('lenz-r-emf').textContent=(100*m.response).toLocaleString('pt-BR',{maximumFractionDigits:0})+'%';
    get('lenz-step1').textContent=moving
      ? direction==='approach'?'O polo N do ímã se aproxima da face esquerda da bobina.':'O polo N do ímã se afasta da face esquerda da bobina.'
      :'O ímã está parado em relação à bobina.';
    get('lenz-step2').textContent=moving
      ? direction==='approach'?'O fluxo magnético para a direita através das espiras aumenta.':'O fluxo magnético para a direita através das espiras diminui.'
      :'O fluxo é constante no tempo: dΦ/dt = 0.';
    get('lenz-step3').textContent=moving
      ? direction==='approach'
        ?'A bobina cria um polo N voltado ao ímã; o campo induzido aponta para a esquerda. Corrente anti-horária vista do ímã.'
        :'A bobina cria um polo S voltado ao ímã; o campo induzido aponta para a direita. Corrente horária vista do ímã.'
      :'Sem variação do fluxo, a fem induzida é nula neste modelo ideal.';
    if(result){
      result.replaceChildren();
      const box=document.createElement('div');box.className='result';
      const bold=document.createElement('strong');
      bold.textContent=moving?(direction==='approach'?'Aproximação: ':'Afastamento: '):'Ímã parado: ';
      box.appendChild(bold);
      box.appendChild(document.createTextNode(moving
        ?(direction==='approach'?'o campo da bobina se opõe ao aumento do fluxo.':'o campo da bobina se opõe à diminuição do fluxo.')
        :'sem alteração de fluxo, não há corrente induzida.'));
      result.appendChild(box);
    }
  }
  function render() { draw();updateInfo(); }
  function pause(){
    playing=false;lastFrame=null;
    if(frameId)root.cancelAnimationFrame(frameId);
    frameId=0;
    render();
  }
  function reset(){
    playing=false;lastFrame=null;
    if(frameId)root.cancelAnimationFrame(frameId);
    frameId=0;
    direction='approach';modeSelect.value='approach';p=.12;
    render();
  }
  function play(mode){
    const next=mode==='away'?'away':'approach';
    if(frameId)root.cancelAnimationFrame(frameId);
    // Each click starts a clearly visible full trajectory. AFASTAMENTO
    // always begins next to the coil and travels LEFT, away from it.
    direction=next;
    p=typeof root.scicalcLenzStartPosition==='function'
      ?root.scicalcLenzStartPosition(next):(next==='away'?.88:.12);
    modeSelect.value=next;
    playing=true;lastFrame=null;lastInfoTime=0;
    render();
    frameId=root.requestAnimationFrame(tick);
  }
  function tick(now) {
    if(!playing)return;
    const dt=lastFrame===null?0:Math.min(.05,Math.max(0,(now-lastFrame)/1000));
    lastFrame=now;
    const progress=root.scicalcLenzAdvance
      ?root.scicalcLenzAdvance(p,direction,currentSpeed(),dt)
      :{position:clamp(p+dt*currentSpeed()*.46*(direction==='approach'?1:-1),.06,.94),finished:false};
    p=progress.position;
    if(progress.finished){pause();return;}
    draw();
    if(now-lastInfoTime>90){updateInfo();lastInfoTime=now;}
    frameId=root.requestAnimationFrame(tick);
  }

  root.flStartLenz=play;
  root.flStopLenz=pause;
  root.flResetLenz=reset;
  modeSelect.addEventListener('change',function(){
    pause();direction=modeSelect.value==='away'?'away':'approach';
    p=direction==='approach'?.12:.88;render();
  });
  speed.addEventListener('input',function(){
    const label=get('fl-lenz-value');
    if(label)label.textContent=currentSpeed().toLocaleString('pt-BR',{maximumFractionDigits:2})+'×';
    if(!playing)render();
  });
  document.addEventListener('visibilitychange',function(){if(document.hidden&&playing)pause();});
  const section=document.getElementById('simuladores-v26');
  const tab=document.querySelector('nav.tabs [data-tab="simuladores-v26"]');
  if(tab)tab.addEventListener('click',function(){root.requestAnimationFrame(render);});
  if(typeof ResizeObserver==='function'){
    const observer=new ResizeObserver(function(){if(section?.classList.contains('active'))draw();});
    observer.observe(canvas);
  }else{
    root.addEventListener('resize',function(){if(section?.classList.contains('active'))draw();});
  }
  render();
})(window);
