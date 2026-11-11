/* SciCalc Pro: laboratório avançado independente do motor clássico. */
(function () {
  'use strict';
  const $ = id => document.getElementById(id);
  const SCI_PRO_BUILD = '20261111c'; // Versiona também o Web Worker carregado pelo navegador.
  const nav = document.querySelector('nav.tabs');
  const main = document.querySelector('main');
  if (!nav || !main || $('scicalc-pro')) return;

  const tab = document.createElement('button');
  tab.className = 'tab';
  tab.type = 'button';
  tab.dataset.tab = 'scicalc-pro';
  tab.textContent = '🧠 Cálculo Pro';
  const historyTab = nav.querySelector('[data-tab="historico"]');
  nav.insertBefore(tab, historyTab || null);

  const section = document.createElement('section');
  section.id = 'scicalc-pro';
  section.className = 'tab-content';
  section.innerHTML = [
    '<div class="panel">',
      '<div class="sc-pro-intro"><div><h2>🧠 Laboratório de cálculo avançado</h2>',
      '<p class="sub">SymPy no navegador: integrais, derivadas de ordem superior e limites com análise lateral.</p></div>',
      '<span class="sc-pro-chip">Motor independente • processamento local</span></div>',
      '<div class="sc-pro-grid">',
       '<div class="sc-pro-panel"><h3>Definir o cálculo</h3>',
         '<p class="sc-pro-hint">Escreva expressões com x, potências (^), sin, cos, exp, sqrt, log, erf e constantes pi ou E.</p>',
         '<form id="sc-pro-form">',
         '<label class="sc-pro-label" for="sc-pro-action">Operação</label>',
         '<select id="sc-pro-action"><option value="integral">∫ Integral</option><option value="limit">Limite</option><option value="derivative">Derivada</option></select>',
         '<label class="sc-pro-label" for="sc-pro-expr">Função f(x)</label>',
         '<input id="sc-pro-expr" value="exp(-x^2)" placeholder="Ex.: exp(-x^2)" required maxlength="240" spellcheck="false">',
         '<div id="sc-pro-integral-fields" class="sc-pro-row">',
           '<div><label class="sc-pro-label" for="sc-pro-lower">Limite inferior (opcional)</label><input id="sc-pro-lower" placeholder="0 ou -inf"></div>',
           '<div><label class="sc-pro-label" for="sc-pro-upper">Limite superior (opcional)</label><input id="sc-pro-upper" placeholder="1 ou inf"></div>',
         '</div>',
         '<div id="sc-pro-limit-fields" class="sc-pro-row" hidden>',
           '<div><label class="sc-pro-label" for="sc-pro-target">x tende a</label><input id="sc-pro-target" value="0" placeholder="0, inf, -inf"></div>',
           '<div><label class="sc-pro-label" for="sc-pro-direction">Direção</label><select id="sc-pro-direction"><option value="both">Bilateral</option><option value="left">Esquerda</option><option value="right">Direita</option></select></div>',
         '</div>',
         '<div id="sc-pro-derivative-fields" class="sc-pro-row" hidden>',
           '<div><label class="sc-pro-label" for="sc-pro-order">Ordem</label><select id="sc-pro-order"><option>1</option><option>2</option><option>3</option><option>4</option><option>5</option><option>6</option></select></div>',
           '<div><label class="sc-pro-label" for="sc-pro-point">Avaliar em x (opcional)</label><input id="sc-pro-point" placeholder="Ex.: pi/2"></div>',
         '</div>',
         '<p class="sc-pro-help">Primeira utilização requer internet para carregar o módulo Python/SymPy. O processamento matemático ocorre localmente, sem enviar sua expressão a um servidor.</p>',
         '<div class="sc-pro-buttons"><button class="btn" id="sc-pro-solve" type="submit">Resolver com SymPy</button></div>',
         '</form>',
         '<h3>Experimente</h3><div class="sc-pro-buttons" id="sc-pro-examples">',
           '<button type="button" class="sc-pro-example" data-action="integral" data-expr="exp(-x^2)">∫ exp(−x²) dx</button>',
           '<button type="button" class="sc-pro-example" data-action="integral" data-expr="sin(x^2)" data-lower="0" data-upper="1">∫₀¹ sin(x²) dx</button>',
           '<button type="button" class="sc-pro-example" data-action="limit" data-expr="(exp(x)-1-x)/x^2" data-target="0">Limite 0/0</button>',
           '<button type="button" class="sc-pro-example" data-action="limit" data-expr="1/x" data-target="0">Laterais 1/x</button>',
           '<button type="button" class="sc-pro-example" data-action="derivative" data-expr="x^2*exp(x)" data-order="3">3ª derivada</button>',
         '</div>',
       '</div>',
       '<div class="sc-pro-panel"><h3>Resultado matemático</h3>',
         '<p class="sc-pro-hint">Exibição simbólica, valor aproximado quando disponível e diagnóstico da operação.</p>',
         '<div id="sc-pro-status" class="sc-pro-status" aria-live="polite">Escolha um exemplo e clique em resolver.</div>',
         '<div id="sc-pro-result" class="sc-pro-result" aria-live="polite">A resposta será exibida aqui.</div>',
         '<p class="sc-pro-warning">Resultados aproximados não são certificados por aritmética intervalar. Integrais impróprias, divergências e domínios de funções exigem análise matemática.</p>',
       '</div>',
      '</div>',
    '</div>'
  ].join('');
  const historySection = $('historico');
  main.insertBefore(section, historySection || null);

  function activate() {
    nav.querySelectorAll('.tab').forEach(el => el.classList.remove('active'));
    main.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
    tab.classList.add('active');
    section.classList.add('active');
  }
  tab.addEventListener('click', activate);

  function updateFields() {
    const a = $('sc-pro-action').value;
    $('sc-pro-integral-fields').hidden = a !== 'integral';
    $('sc-pro-limit-fields').hidden = a !== 'limit';
    $('sc-pro-derivative-fields').hidden = a !== 'derivative';
  }
  $('sc-pro-action').addEventListener('change', updateFields);
  updateFields();

  let worker = null, isReady = false, activeRequest = null;
  let nextId = 0, loadTimer = null, solveTimer = null;
  function status(text) { $('sc-pro-status').textContent = text; }
  function setBusy(busy) { $('sc-pro-solve').disabled = busy; $('sc-pro-solve').textContent = busy ? 'Processando…' : 'Resolver com SymPy'; }
  function stopTimers() { clearTimeout(loadTimer); clearTimeout(solveTimer); }
  function resetWorker(reason) {
    stopTimers();
    if (worker) worker.terminate();
    worker = null;
    isReady = false;
    activeRequest = null;
    setBusy(false);
    if (reason) { status('Não foi possível concluir o cálculo.'); $('sc-pro-result').textContent = reason; }
  }
  function renderFormula(parent, expression) {
    const block = document.createElement('div');
    block.className = 'sc-pro-math';
    parent.appendChild(block);
    try {
      if (window.katex && typeof katex.render === 'function') {
        katex.render(String(expression), block, { displayMode: true, throwOnError: true, trust: false, strict: 'warn' });
        return;
      }
    } catch (error) { /* Mostrar texto acessível quando KaTeX falhar. */ }
    block.textContent = String(expression);
  }
  function addText(parent, tag, value) {
    const node = document.createElement(tag);
    node.textContent = value;
    parent.appendChild(node);
    return node;
  }

  let lastMathResult = null;
  function renderVerifiedMathSteps(host, data) {
    host.querySelectorAll('.sc-pro-verified-math-steps').forEach(node => node.remove());
    if($('sc-pro-show-steps')?.checked===false || !data?.ok) return;
    if(!Array.isArray(data.steps) || !data.steps.length) return;
    const details=document.createElement('details');
    details.className='didactic sc-pro-verified-math-steps';
    details.open=true;
    const summary=document.createElement('summary');
    summary.textContent='📚 Resolução matemática — passo a passo';
    const body=document.createElement('div');body.className='didactic-body';
    const intro=document.createElement('p');intro.className='didactic-intro';
    intro.textContent='Cada etapa abaixo corresponde a uma identidade, transformação ou verificação matemática identificada pelo motor simbólico.';
    body.appendChild(intro);
    const steps=document.createElement('div');steps.className='step-list';
    data.steps.forEach(function(item){
      if(!item || typeof item.title!=='string')return;
      const stage=document.createElement('div');stage.className='step';
      addText(stage,'div',item.title).className='step-title';
      if(item.tex)renderFormula(stage,item.tex);
      if(item.note)addText(stage,'div',item.note).className='step-note';
      steps.appendChild(stage);
    });
    body.appendChild(steps);
    details.append(summary,body);host.appendChild(details);
  }
  window.scicalcProToggleSteps = function() {
    const node=$('sc-pro-result');
    if(node&&lastMathResult)renderVerifiedMathSteps(node,lastMathResult);
  };

  function displayResult(data) {
    const root = $('sc-pro-result');
    root.replaceChildren();
    if (!data || !data.ok) {
      status('Expressão não resolvida.');
      addText(root, 'p', data && data.error ? data.error : 'Falha inesperada.');
      return;
    }
    status('Concluído • ' + data.status);
    addText(root, 'span', data.status).className = 'sc-pro-exact';
    addText(root, 'h4', 'Resultado:');
    renderFormula(root, data.result_tex || data.result_text);
    if (data.constant) addText(root, 'p', '+ C (constante de integração)');
    if (data.approx) {
      addText(root, 'h4', 'Aproximação decimal:');
      addText(root, 'strong', '≈ ' + data.approx);
    }
    if (data.point_result_tex) {
      addText(root, 'h4', 'Valor no ponto solicitado:');
      renderFormula(root, data.point_result_tex);
      if (data.point_approx) addText(root, 'p', '≈ ' + data.point_approx);
    }
    if (data.left_tex || data.right_tex) {
      addText(root, 'h4', 'Limites laterais:');
      if (data.left_tex) { addText(root, 'p', 'Pela esquerda'); renderFormula(root, data.left_tex); }
      if (data.right_tex) { addText(root, 'p', 'Pela direita'); renderFormula(root, data.right_tex); }
    }
    if (data.verified) addText(root, 'p', '✓ Conferência simbólica: derivar a primitiva recupera o integrando.');
    if (data.note) addText(root, 'p', data.note).className = 'sc-pro-warning';
    if (data.engine_revision) addText(root, 'small', 'Versão do motor SymPy: ' + data.engine_revision).className = 'sc-pro-engine-revision';
    lastMathResult=data;
    renderVerifiedMathSteps(root,data);
  }
  function dispatch() {
    if (!worker || !isReady || !activeRequest) return;
    worker.postMessage({ id: activeRequest.id, request: activeRequest.request });
    solveTimer = setTimeout(() => resetWorker('O cálculo excedeu 25 segundos. Tente uma expressão mais simples ou um intervalo menor.'), 25000);
    status('Calculando com SymPy…');
  }
  function ensureWorker() {
    if (worker) return;
    try {
      worker = new Worker(new URL('./scicalc-pro-worker.mjs?build=' + SCI_PRO_BUILD, document.baseURI), { type: 'module' });
    } catch (err) {
      resetWorker('Este navegador não conseguiu iniciar o motor avançado. ' + err.message);
      return;
    }
    worker.addEventListener('message', event => {
      const msg = event.data || {};
      if (msg.kind === 'ready') {
        isReady = true;
        clearTimeout(loadTimer);
        dispatch();
      } else if (msg.kind === 'failed') {
        resetWorker('Não foi possível carregar SymPy. Confira sua conexão e tente novamente.');
      } else if (msg.kind === 'result' && activeRequest && msg.id === activeRequest.id) {
        clearTimeout(solveTimer);
        activeRequest = null;
        setBusy(false);
        displayResult(msg.payload);
      }
    });
    worker.addEventListener('error', () => resetWorker('Erro ao carregar o módulo de cálculo; confira sua conexão e a disponibilidade do serviço CDN.'));
    loadTimer = setTimeout(() => resetWorker('Tempo de carregamento excedido. Confira sua conexão ou tente novamente.'), 70000);
  }
  function collect() {
    return {
      action: $('sc-pro-action').value,
      expr: $('sc-pro-expr').value.trim(),
      lower: $('sc-pro-lower').value.trim(),
      upper: $('sc-pro-upper').value.trim(),
      target: $('sc-pro-target').value.trim(),
      direction: $('sc-pro-direction').value,
      order: Number($('sc-pro-order').value),
      point: $('sc-pro-point').value.trim()
    };
  }
  function solve(event) {
    if (event) event.preventDefault();
    if (activeRequest) return;
    const request = collect();
    if (!request.expr) { status('Digite uma função.'); return; }
    activeRequest = { id: ++nextId, request };
    setBusy(true);
    status('Carregando o motor simbólico…');
    $('sc-pro-result').textContent = 'A interface permanece disponível durante o processamento.';
    ensureWorker();
    if (isReady) dispatch();
  }
  $('sc-pro-form').addEventListener('submit', solve);

  $('sc-pro-examples').addEventListener('click', e => {
    const b = e.target.closest('button[data-action]');
    if (!b || activeRequest) return;
    $('sc-pro-action').value = b.dataset.action;
    $('sc-pro-expr').value = b.dataset.expr;
    $('sc-pro-lower').value = b.dataset.lower || '';
    $('sc-pro-upper').value = b.dataset.upper || '';
    $('sc-pro-target').value = b.dataset.target || '0';
    $('sc-pro-order').value = b.dataset.order || '1';
    $('sc-pro-direction').value = 'both';
    updateFields();
    solve();
  });

  function addUpgradeButton(sectionId, sourceIds, type) {
    const original = document.querySelector('#' + sectionId + ' .panel > button.btn');
    if (!original) return;
    const b = document.createElement('button');
    b.type = 'button';
    b.className = 'btn small sc-pro-action';
    b.textContent = '🧠 Resolver com motor avançado';
    original.insertAdjacentElement('afterend', b);
    b.addEventListener('click', () => {
      if (activeRequest) { activate(); return; }
      $('sc-pro-action').value = type;
      $('sc-pro-expr').value = $(sourceIds.expr)?.value || '';
      if (type === 'integral') {
        $('sc-pro-lower').value = $(sourceIds.lower)?.value || '';
        $('sc-pro-upper').value = $(sourceIds.upper)?.value || '';
      }
      if (type === 'limit') {
        $('sc-pro-target').value = $(sourceIds.target)?.value || '0';
        const d = $(sourceIds.direction)?.value || 'auto';
        $('sc-pro-direction').value = d === 'left' || d === 'right' ? d : 'both';
      }
      updateFields();
      activate();
      solve();
    });
  }
  addUpgradeButton('integrais', {expr:'int-f',lower:'int-x1',upper:'int-x2'}, 'integral');
  addUpgradeButton('limites', {expr:'lim-f',target:'lim-a',direction:'lim-dir'}, 'limit');
  addUpgradeButton('derivadas', {expr:'der-f'}, 'derivative');

  /* Physics visual: existing simulator parameters; animation is intentionally separate from formula engines. */
  const sims = $('simuladores-v26');
  if (!sims) return;
  const panel = document.createElement('div');
  panel.className = 'sc-pro-sim-panel';
  panel.innerHTML = [
    '<header><h3>🎬 Estúdio de simulações animadas</h3>',
    '<p>Visualização interativa de MHS e lançamento oblíquo. Altere os parâmetros nos simuladores abaixo e veja o movimento atualizado.</p></header>',
    '<div class="sc-pro-row">',
      '<div><label for="sc-pro-sim-mode">Experimento</label><select id="sc-pro-sim-mode"><option value="mhs">Massa–mola (MHS)</option><option value="projectile">Lançamento oblíquo</option></select></div>',
      '<div><label for="sc-pro-speed">Velocidade da reprodução</label><select id="sc-pro-speed"><option value="0.5">0,5×</option><option value="1" selected>1×</option><option value="2">2×</option></select></div>',
    '</div>',
    '<div class="sc-pro-sim"><canvas id="sc-pro-sim-canvas" role="img" aria-label="Visualização de movimento de massa–mola ou lançamento oblíquo"></canvas></div>',
    '<div class="sc-pro-sim-toolbar">',
      '<button class="btn small" type="button" id="sc-pro-play">▶ Reproduzir</button>',
      '<button class="sc-pro-section-btn" type="button" id="sc-pro-reset">↺ Reiniciar</button>',
      '<label for="sc-pro-scrub">Tempo</label><input id="sc-pro-scrub" type="range" min="0" max="1000" value="0">',
      '<strong id="sc-pro-clock">0,00 s</strong>',
    '</div>',
    '<div id="sc-pro-metrics" class="sc-pro-metrics" aria-live="off"></div>',
    '<p class="sc-pro-help">Modelos ideais: MHS sem amortecimento; lançamento sem resistência do ar. Distâncias e tempos são calculados com os parâmetros dos respectivos formulários.</p>'
  ].join('');
  sims.insertBefore(panel, sims.firstElementChild);

  const canvas = $('sc-pro-sim-canvas'), ctx = canvas.getContext('2d');
  let elapsed = 0, running = false, lastFrame = 0, frameId = 0;
  const format = (v, unit) => (Number.isFinite(v) ? Number(v).toLocaleString('pt-BR',{maximumFractionDigits:2}) : '—') + (unit ? ' ' + unit : '');
  const getNum = (id, fallback) => {
    const raw = $(id)?.value;
    const v = Number(String(raw == null ? fallback : raw).replace(',', '.'));
    return Number.isFinite(v) ? v : fallback;
  };
  function model() {
    if ($('sc-pro-sim-mode').value === 'projectile') {
      const v0 = Math.max(0, getNum('v26-obl-v', 30));
      const angle = getNum('v26-obl-a', 45) * Math.PI / 180;
      const h0 = Math.max(0, getNum('v26-obl-h', 0));
      const g = Math.max(.001, getNum('v26-obl-g', 9.81));
      const vx = v0 * Math.cos(angle), vy = v0 * Math.sin(angle);
      const duration = (vy + Math.sqrt(vy * vy + 2 * g * h0)) / g;
      return { kind:'projectile', vx, vy, g, h0, duration:Math.max(0.01,duration) };
    }
    return { kind:'mhs', A:Math.max(.001,Math.abs(getNum('v26-mhs-a',1))),
      w:Math.max(.001,Math.abs(getNum('v26-mhs-w',2))),
      phi:getNum('v26-mhs-p',0),
      duration:Math.max(.1,getNum('v26-mhs-t',10)) };
  }
  function draw() {
    if (!ctx) return;
    const rect = canvas.getBoundingClientRect();
    const W = Math.max(300, rect.width), H = Math.max(220, rect.height);
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = Math.round(W*dpr); canvas.height = Math.round(H*dpr);
    ctx.setTransform(dpr,0,0,dpr,0,0);
    const m = model(), t = Math.min(elapsed, m.duration);
    ctx.fillStyle = '#0b142a'; ctx.fillRect(0,0,W,H);
    ctx.strokeStyle='#22334e'; ctx.lineWidth=1;
    for(let x=28;x<W;x+=45){ctx.beginPath();ctx.moveTo(x,20);ctx.lineTo(x,H-22);ctx.stroke();}
    for(let y=25;y<H;y+=45){ctx.beginPath();ctx.moveTo(20,y);ctx.lineTo(W-20,y);ctx.stroke();}
    ctx.font = '13px system-ui,sans-serif';
    ctx.fillStyle = '#cbd5e1';
    if(m.kind === 'mhs') {
      const mid = W/2, span = Math.min(W*.30,180);
      const x = m.A * Math.cos(m.w*t + m.phi), v = -m.A*m.w*Math.sin(m.w*t+m.phi);
      const bob = mid + span*x/m.A, wall = 35, cy = H*.37;
      ctx.fillStyle='#7086a8';ctx.fillRect(wall-8,cy-48,9,96);
      ctx.strokeStyle='#a5b4fc';ctx.lineWidth=3;ctx.beginPath();ctx.moveTo(wall,cy);
      const length=bob-wall;
      for(let i=0;i<=24;i++) {
        const sx=wall+(length-15)*i/24;
        const sy=cy+(i===0||i===24?0:((i%2)?-14:14));
        ctx.lineTo(sx,sy);
      }
      ctx.stroke();
      ctx.fillStyle='#22d3ee';ctx.beginPath();ctx.arc(bob,cy,21,0,Math.PI*2);ctx.fill();
      ctx.strokeStyle='#7dd3fc';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(mid,cy-30);ctx.lineTo(mid,H-45);ctx.stroke();
      ctx.fillStyle='#f8fafc';ctx.fillText('x(t) = A·cos(ωt+φ)',24,24);
      const top=H*.65, bottom=H-35, left=35, right=W-25;
      ctx.strokeStyle='#fb7185';ctx.lineWidth=2;ctx.beginPath();
      for(let i=0;i<=180;i++){
        const tt=m.duration*i/180;
        const px=left+(right-left)*i/180;
        const py=(top+bottom)/2-((top-bottom)*-.40)*Math.cos(m.w*tt+m.phi);
        if(i===0)ctx.moveTo(px,py);else ctx.lineTo(px,py);
      }ctx.stroke();
      const nowX=left+(right-left)*t/m.duration;
      const nowY=(top+bottom)/2-((top-bottom)*-.40)*Math.cos(m.w*t+m.phi);
      ctx.fillStyle='#fbbf24';ctx.beginPath();ctx.arc(nowX,nowY,5,0,Math.PI*2);ctx.fill();
      $('sc-pro-metrics').innerHTML='';
      addMetric('Posição',format(x,'m'));addMetric('Velocidade',format(v,'m/s'));addMetric('Período',format(2*Math.PI/m.w,'s'));
    } else {
      const xAt = q => m.vx*q, yAt=q=>Math.max(0,m.h0+m.vy*q-.5*m.g*q*q);
      let xMin=0,xMax=0,yMax=Math.max(1,m.h0);
      for(let i=0;i<=150;i++) { const q=m.duration*i/150;const x=xAt(q);xMin=Math.min(xMin,x);xMax=Math.max(xMax,x);yMax=Math.max(yMax,yAt(q)); }
      const pad=44, xRange=Math.max(.001,xMax-xMin), yRange=Math.max(.001,yMax);
      const px=x=>pad+(W-2*pad)*(x-xMin)/xRange;
      const py=y=>H-pad-(H-2*pad)*y/yRange;
      ctx.strokeStyle='#94a3b8';ctx.lineWidth=2;ctx.beginPath();ctx.moveTo(pad,H-pad);ctx.lineTo(W-pad,H-pad);ctx.stroke();
      ctx.setLineDash([5,7]);ctx.strokeStyle='#506582';ctx.beginPath();
      for(let i=0;i<=150;i++) {const q=m.duration*i/150;i===0?ctx.moveTo(px(xAt(q)),py(yAt(q))):ctx.lineTo(px(xAt(q)),py(yAt(q)));}ctx.stroke();ctx.setLineDash([]);
      ctx.strokeStyle='#38bdf8';ctx.lineWidth=3;ctx.beginPath();
      for(let i=0;i<=100;i++) {const q=t*i/100;i===0?ctx.moveTo(px(xAt(q)),py(yAt(q))):ctx.lineTo(px(xAt(q)),py(yAt(q)));}ctx.stroke();
      ctx.fillStyle='#fbbf24';ctx.beginPath();ctx.arc(px(xAt(t)),py(yAt(t)),9,0,Math.PI*2);ctx.fill();
      ctx.fillStyle='#f8fafc';ctx.fillText('Trajetória sem arrasto do ar',24,24);
      ctx.fillText('x: '+format(xAt(t),'m'),24,H-13);
      $('sc-pro-metrics').innerHTML='';
      addMetric('Distância horizontal',format(xAt(t),'m'));addMetric('Altura',format(yAt(t),'m'));addMetric('Tempo de voo',format(m.duration,'s'));
    }
    $('sc-pro-clock').textContent=format(t,'s');
    $('sc-pro-scrub').value=String(Math.round(1000*t/m.duration));
  }
  function addMetric(label, value) {
    const box=document.createElement('div'), tag=document.createElement('span'), v=document.createElement('strong');
    tag.textContent=label;v.textContent=value;box.append(tag,v);$('sc-pro-metrics').appendChild(box);
  }
  function pause() { running=false;cancelAnimationFrame(frameId);$('sc-pro-play').textContent='▶ Reproduzir'; }
  function reset() { pause();elapsed=0;draw(); }
  function tick(timestamp){
    if(!running)return;
    if(lastFrame){elapsed+=Math.min((timestamp-lastFrame)/1000,0.1)*Number($('sc-pro-speed').value);}
    lastFrame=timestamp;
    if(elapsed>=model().duration){elapsed=model().duration;pause();}
    draw();
    if(running)frameId=requestAnimationFrame(tick);
  }
  $('sc-pro-play').addEventListener('click',()=>{
    if(running){pause();return;}
    if(elapsed>=model().duration) elapsed=0;
    running=true;lastFrame=0;$('sc-pro-play').textContent='⏸ Pausar';
    frameId=requestAnimationFrame(tick);
  });
  $('sc-pro-reset').addEventListener('click',reset);
  $('sc-pro-sim-mode').addEventListener('change',reset);
  $('sc-pro-scrub').addEventListener('input',()=>{
    pause();elapsed=model().duration*Number($('sc-pro-scrub').value)/1000;draw();
  });
  ['v26-mhs-a','v26-mhs-w','v26-mhs-p','v26-mhs-t','v26-obl-v','v26-obl-a','v26-obl-h','v26-obl-g']
    .forEach(id=>$(id)?.addEventListener('input',reset));
  const simTab = nav.querySelector('[data-tab="simuladores-v26"]');
  if(simTab) simTab.addEventListener('click',()=>requestAnimationFrame(draw));
  document.addEventListener('visibilitychange',()=>{if(document.hidden)pause();});
  window.addEventListener('resize',()=>{if(sims.classList.contains('active'))draw();});
  requestAnimationFrame(draw);
})();
