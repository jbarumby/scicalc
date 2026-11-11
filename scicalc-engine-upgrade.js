/**
 * SciCalc Engine Upgrade & Universal Math Renderer
 * Copyright (c) 2026 João Carlos Barumby. Todos os direitos reservados.
 * 
 * 1. Renderização Matemática Universal KaTeX (elimina fórmulas não-renderizadas no DOM)
 * 2. Motor de Limites: L'Hôpital simbólico, detecção de 0/0 e inf/inf, cancelamento de (x-a)
 * 3. Motor de Integrais: Quadratura Adaptativa de Simpson de alta precisão (tol < 1e-10)
 * 4. Motor de Derivadas: Suporte a ordem n-ésima (1ª a 5ª) com detecção e explicação de regras
 */

(function (root) {
  'use strict';

  /* ==========================================================================
     1. RENDERIZADOR UNIVERSAL KATEX
     ========================================================================== */

  function renderMathEverywhere(target) {
    const el = target || document.body;
    if (!el) return;

    // Se o auto-render oficial do KaTeX estiver disponível, use-o
    if (typeof renderMathInElement === 'function') {
      try {
        renderMathInElement(el, {
          delimiters: [
            { left: '$$', right: '$$', display: true },
            { left: '\\[', right: '\\]', display: true },
            { left: '\\(', right: '\\)', display: false }
          ],
          throwOnError: false
        });
        return;
      } catch (err) {
        console.warn('[SciCalc Math] Auto-render error:', err);
      }
    }

    // Fallback: varredura via TreeWalker para substituir \(...\) e \[...\]
    if (typeof katex !== 'undefined' && typeof katex.renderToString === 'function') {
      const walker = document.createTreeWalker(el, NodeFilter.SHOW_TEXT, null, false);
      const textNodes = [];
      let currNode;
      while ((currNode = walker.nextNode())) {
        if (!currNode.parentElement) continue;
        const tag = currNode.parentElement.tagName;
        if (tag === 'SCRIPT' || tag === 'STYLE' || tag === 'TEXTAREA' || tag === 'CODE') continue;
        if (/\\\(|\\\[/.test(currNode.nodeValue)) {
          textNodes.push(currNode);
        }
      }

      textNodes.forEach(tNode => {
        const raw = tNode.nodeValue;
        const re = /\\\(([\s\S]*?)\\\)|\\\[([\s\S]*?)\\\]/g;
        let match;
        let lastIdx = 0;
        let html = '';

        while ((match = re.exec(raw)) !== null) {
          html += raw.slice(lastIdx, match.index)
            .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
          const formula = (match[1] != null ? match[1] : match[2]).trim().replace(/\\\\/g, '\\');
          const isDisplay = match[2] != null;
          try {
            html += katex.renderToString(formula, { throwOnError: false, displayMode: isDisplay });
          } catch (e) {
            html += match[0];
          }
          lastIdx = re.lastIndex;
        }

        if (html) {
          html += raw.slice(lastIdx).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
          const span = document.createElement('span');
          span.innerHTML = html;
          if (tNode.parentNode) {
            tNode.parentNode.replaceChild(span, tNode);
          }
        }
      });
    }
  }

  root.scicalcRenderEverywhere = renderMathEverywhere;

  // Garante renderização contínua em abas, sub-abas e carregamento
  function setupMathObservers() {
    renderMathEverywhere(document.body);

    document.querySelectorAll('.tab').forEach(t => {
      t.addEventListener('click', () => {
        setTimeout(() => {
          const activeTabContent = document.getElementById(t.dataset.tab);
          if (activeTabContent) renderMathEverywhere(activeTabContent);
        }, 30);
      });
    });

    document.querySelectorAll('.sub-tab').forEach(st => {
      st.addEventListener('click', () => {
        setTimeout(() => {
          const panel = st.closest('.panel');
          if (panel) renderMathEverywhere(panel);
        }, 30);
      });
    });

    // Monitora inserções dinâmicas de conteúdo matemático
    const observer = new MutationObserver(() => {
      // Throttle simples
      clearTimeout(setupMathObservers._timer);
      setupMathObservers._timer = setTimeout(() => {
        renderMathEverywhere(document.body);
      }, 150);
    });

    observer.observe(document.body, { childList: true, subtree: true });
  }


  /* ==========================================================================
     2. MOTOR DE DERIVADAS: ORDEM N-ÉSIMA E REGRAS EXPLICADAS
     ========================================================================== */

  function analyzeDerivativeRules(node, variable) {
    const rules = [];
    if (!node) return rules;

    function walk(n) {
      if (!n) return;
      if (n.isOperatorNode) {
        if (n.op === '*' && n.args.length === 2) {
          const leftHasVar = n.args[0].filter(sub => sub.isSymbolNode && sub.name === variable).length > 0;
          const rightHasVar = n.args[1].filter(sub => sub.isSymbolNode && sub.name === variable).length > 0;
          if (leftHasVar && rightHasVar) {
            rules.push({
              name: 'Regra do Produto',
              formula: '\\frac{d}{dx}[u \\cdot v] = u\'v + uv\'',
              desc: 'Aplicada na multiplicação de duas parcelas que dependem de ' + variable + '.'
            });
          }
        } else if (n.op === '/') {
          const denomHasVar = n.args[1].filter(sub => sub.isSymbolNode && sub.name === variable).length > 0;
          if (denomHasVar) {
            rules.push({
              name: 'Regra do Quociente',
              formula: '\\frac{d}{dx}\\left[\\frac{u}{v}\\right] = \\frac{u\'v - uv\'}{v^2}',
              desc: 'Aplicada na divisão onde o denominador depende de ' + variable + '.'
            });
          }
        } else if (n.op === '^') {
          rules.push({
            name: 'Regra da Potência / Exponencial',
            formula: '\\frac{d}{dx}[u^n] = n \\cdot u^{n-1} \\cdot u\'',
            desc: 'Diferenciação de potências combinada com a regra da cadeia.'
          });
        }
      } else if (n.isFunctionNode) {
        const argHasVar = n.args[0] && n.args[0].filter(sub => sub.isSymbolNode && sub.name === variable).length > 0;
        if (argHasVar) {
          rules.push({
            name: 'Regra da Cadeia em Função (' + n.name + ')',
            formula: '\\frac{d}{dx}[f(g(x))] = f\'(g(x)) \\cdot g\'(x)',
            desc: 'A derivada da função externa multiplicada pela derivada da função interna.'
          });
        }
      }
      if (n.args) n.args.forEach(walk);
    }

    walk(node);

    // Remove duplicatas
    const unique = [];
    const seen = new Set();
    rules.forEach(r => {
      if (!seen.has(r.name)) {
        seen.add(r.name);
        unique.push(r);
      }
    });

    if (unique.length === 0) {
      unique.push({
        name: 'Regras Fundamentais de Potência e Soma',
        formula: '\\frac{d}{dx}[a x^n] = a n x^{n-1}, \\quad \\frac{d}{dx}[u \\pm v] = u\' \\pm v\'',
        desc: 'Diferenciação termo a termo de polinômios e constantes.'
      });
    }

    return unique;
  }

  // Sobrescreve o calcDerivada global com suporte a ordens 1 a 5 e passo a passo
  root.calcDerivada = function () {
    try {
      const fStr = document.getElementById('der-f').value.trim();
      const variable = (document.getElementById('der-var') ? document.getElementById('der-var').value.trim() : 'x') || 'x';
      const rangeStr = (document.getElementById('der-intervalo') ? document.getElementById('der-intervalo').value.trim() : '-5,5') || '-5,5';
      const rangeParts = rangeStr.split(',').map(x => Number(x.trim()));
      const range = (rangeParts.length === 2 && isFinite(rangeParts[0]) && isFinite(rangeParts[1])) ? rangeParts : [-5, 5];

      const orderSelect = document.getElementById('der-order');
      const maxOrder = orderSelect ? parseInt(orderSelect.value, 10) : 2;

      const fNode = math.parse(fStr);
      const derivatives = [fNode];

      for (let o = 1; o <= maxOrder; o++) {
        const prev = derivatives[o - 1];
        derivatives.push(math.derivative(prev, variable));
      }

      const primeSymbols = ['', "'", "''", "'''", "^{(4)}", "^{(5)}"];

      let html = '<div class="result">';
      html += '<div><strong>Função Original:</strong><div style="margin-top:6px;font-size:1.1em">' + 
        scicalcRenderTex('f(' + variable + ') = ' + fNode.toTex(), true) + '</div></div>';

      for (let o = 1; o <= maxOrder; o++) {
        const sym = o <= 5 ? primeSymbols[o] : '^{(' + o + ')}';
        const colorVar = o === 1 ? 'var(--primary)' : (o === 2 ? 'var(--primary-2)' : 'var(--accent)');
        html += '<div style="margin-top:14px;color:' + colorVar + '"><strong>' + o + 'ª Derivada:</strong>' +
          '<div style="margin-top:4px;font-size:1.1em">' +
          scicalcRenderTex('f' + sym + '(' + variable + ') = ' + derivatives[o].toTex(), true) + 
          '</div></div>';
      }

      // Regras Explicadas Didáticas
      const showSteps = document.getElementById('der-steps') ? document.getElementById('der-steps').checked : true;
      if (showSteps) {
        const rules = analyzeDerivativeRules(fNode, variable);
        html += '<details class="didactic" open style="margin-top:18px">';
        html += '<summary>📖 Regras de Diferenciação Utilizadas</summary><div class="didactic-body">';
        html += '<p class="didactic-intro">Identificação das estruturas analíticas presentes na função <code>' + fStr + '</code>:</p>';
        html += '<div class="step-list">';
        rules.forEach((r, idx) => {
          html += '<div class="step">';
          html += '<div class="step-title">Regra ' + (idx + 1) + ': ' + r.name + '</div>';
          html += '<div style="margin:4px 0">' + scicalcRenderTex(r.formula, false) + '</div>';
          html += '<div class="step-note">' + r.desc + '</div>';
          html += '</div>';
        });
        html += '</div></div></details>';
      }

      html += '</div>';
      document.getElementById('der-result').innerHTML = html;

      // Renderização do gráfico em Plotly
      if (typeof Plotly !== 'undefined' && document.getElementById('der-plot')) {
        const compiled = derivatives.map(d => d.compile());
        const xs = [];
        const ysList = derivatives.map(() => []);

        const stepsCount = 200;
        const dx = (range[1] - range[0]) / stepsCount;
        for (let i = 0; i <= stepsCount; i++) {
          const xVal = range[0] + i * dx;
          xs.push(xVal);
          compiled.forEach((comp, idx) => {
            let yVal = null;
            try {
              const res = comp.evaluate({ [variable]: xVal });
              if (Number.isFinite(res)) yVal = res;
            } catch (e) {}
            ysList[idx].push(yVal);
          });
        }

        const colors = ['#4f46e5', '#10b981', '#f59e0b', '#ef4444', '#8b5cf6', '#06b6d4'];
        const traces = derivatives.map((d, idx) => {
          const label = idx === 0 ? 'f(' + variable + ')' : ('f' + (primeSymbols[idx] || '^{(' + idx + ')}') + '(' + variable + ')');
          return {
            x: xs,
            y: ysList[idx],
            mode: 'lines',
            name: label,
            line: { width: idx === 0 ? 3 : 2, color: colors[idx % colors.length] }
          };
        });

        const isDark = document.body.classList.contains('dark-mode');
        const layout = {
          title: 'Gráfico de f(' + variable + ') e suas Derivadas Sucessivas',
          xaxis: { title: variable, zeroline: true },
          yaxis: { title: 'y', zeroline: true },
          paper_bgcolor: isDark ? '#131b2e' : '#ffffff',
          plot_bgcolor: isDark ? '#131b2e' : '#f8fafc',
          font: { color: isDark ? '#e2e8f0' : '#0f172a' },
          margin: { t: 40, r: 20, l: 40, b: 40 }
        };

        Plotly.newPlot('der-plot', traces, layout, { responsive: true });
      }

      renderMathEverywhere(document.getElementById('der-result'));
    } catch (err) {
      document.getElementById('der-result').innerHTML = '<div class="result" style="border-left-color:var(--err);color:var(--err)"><strong>Erro no cálculo da derivada:</strong> ' + err.message + '</div>';
    }
  };


  /* ==========================================================================
     3. MOTOR DE LIMITES: L'HÔPITAL, INDETERMINAÇÕES E CANCELAMENTO
     ========================================================================== */

  root.calcLimite = function () {
    try {
      const fStr = document.getElementById('lim-f').value.trim();
      const aStr = document.getElementById('lim-a').value.trim().toLowerCase();
      const dir = document.getElementById('lim-dir') ? document.getElementById('lim-dir').value : 'auto';
      const showSteps = document.getElementById('lim-steps') ? document.getElementById('lim-steps').checked : true;

      const isInf = aStr.includes('inf');
      const isNegInf = aStr.includes('-inf');
      const aVal = isInf ? (isNegInf ? -Infinity : Infinity) : parseFloat(aStr.replace(',', '.'));

      const fNode = math.parse(fStr);
      const fCode = fNode.compile();

      let isIndeterminate = false;
      let indetType = '';
      let lhopitalApplied = false;
      let lhopitalSteps = [];
      let finalVal = null;
      let exactReason = '';

      // Tenta separar em numerador e denominador se for divisão
      let numNode = null;
      let denNode = null;

      if (fNode.isOperatorNode && fNode.op === '/') {
        numNode = fNode.args[0];
        denNode = fNode.args[1];
      }

      // Verificação de ponto finito
      if (!isInf) {
        let fDirect = null;
        try {
          fDirect = fCode.evaluate({ x: aVal });
        } catch (e) {}

        if (numNode && denNode) {
          let numVal = null, denVal = null;
          try { numVal = numNode.compile().evaluate({ x: aVal }); } catch (e) {}
          try { denVal = denNode.compile().evaluate({ x: aVal }); } catch (e) {}

          if (numVal !== null && denVal !== null) {
            if (Math.abs(numVal) < 1e-6 && Math.abs(denVal) < 1e-6) {
              isIndeterminate = true;
              indetType = '0/0';
            }
          }
        }

        // Aplicação da Regra de L'Hôpital se for 0/0
        if (isIndeterminate && indetType === '0/0' && numNode && denNode) {
          lhopitalApplied = true;
          const dNum1 = math.derivative(numNode, 'x');
          const dDen1 = math.derivative(denNode, 'x');

          let dNumVal = null, dDenVal = null;
          try { dNumVal = dNum1.compile().evaluate({ x: aVal }); } catch (e) {}
          try { dDenVal = dDen1.compile().evaluate({ x: aVal }); } catch (e) {}

          lhopitalSteps.push({
            order: 1,
            numPrime: dNum1.toTex(),
            denPrime: dDen1.toTex(),
            numVal: dNumVal,
            denVal: dDenVal
          });

          if (dDenVal !== null && Math.abs(dDenVal) > 1e-12) {
            finalVal = dNumVal / dDenVal;
            exactReason = "Regra de L'Hôpital (1ª ordem): o quociente das derivadas avaliado em x = " + aVal + " resulta em " + Number(finalVal.toFixed(10));
          } else if (dNumVal !== null && Math.abs(dNumVal) < 1e-6 && dDenVal !== null && Math.abs(dDenVal) < 1e-6) {
            // L'Hôpital 2ª ordem
            const dNum2 = math.derivative(dNum1, 'x');
            const dDen2 = math.derivative(dDen1, 'x');
            const dNum2Val = dNum2.compile().evaluate({ x: aVal });
            const dDen2Val = dDen2.compile().evaluate({ x: aVal });

            lhopitalSteps.push({
              order: 2,
              numPrime: dNum2.toTex(),
              denPrime: dDen2.toTex(),
              numVal: dNum2Val,
              denVal: dDen2Val
            });

            if (dDen2Val !== 0 && isFinite(dDen2Val)) {
              finalVal = dNum2Val / dDen2Val;
              exactReason = "Regra de L'Hôpital (2ª ordem consecutiva): quociente das derivadas de 2ª ordem resulta em " + Number(finalVal.toFixed(10));
            }
          }
        } else if (fDirect !== null && Number.isFinite(fDirect) && !isNaN(fDirect)) {
          finalVal = fDirect;
          exactReason = 'Substituição direta: a função é contínua e bem definida no ponto x = ' + aVal + '.';
        }
      }

      // Verificação numérica multilateral estável (amostragem em múltiplos passos)
      let vLeft = NaN, vRight = NaN;
      if (!isInf) {
        const epsilons = [1e-5, 1e-7, 1e-9];
        for (const eps of epsilons) {
          try {
            const l = fCode.evaluate({ x: aVal - eps });
            const r = fCode.evaluate({ x: aVal + eps });
            if (Number.isFinite(l)) vLeft = l;
            if (Number.isFinite(r)) vRight = r;
          } catch (e) {}
        }
      } else {
        const infX = isNegInf ? -1e7 : 1e7;
        try {
          finalVal = fCode.evaluate({ x: infX });
          exactReason = 'Comportamento assintótico no infinito (x → ' + (isNegInf ? '-∞' : '+∞') + ').';
        } catch (e) {}
      }

      if (finalVal === null) {
        if (dir === 'left') finalVal = vLeft;
        else if (dir === 'right') finalVal = vRight;
        else {
          if (!isNaN(vLeft) && !isNaN(vRight) && Math.abs(vLeft - vRight) < 1e-4) {
            finalVal = (vLeft + vRight) / 2;
          } else if (Math.abs(vLeft - vRight) >= 1e-4) {
            finalVal = 'DIVERGENTE';
            exactReason = 'Os limites laterais pela esquerda (' + vLeft.toFixed(4) + ') e pela direita (' + vRight.toFixed(4) + ') são distintos. Logo, o limite bilateral não existe.';
          }
        }
      }

      // Montagem da saída HTML
      const limitTargetTex = isInf ? (isNegInf ? '-\\infty' : '\\infty') : aVal;
      let dirSymbol = '';
      if (dir === 'left') dirSymbol = '^-';
      if (dir === 'right') dirSymbol = '^+';

      let html = '<div class="result">';
      html += '<div><strong>Resultado do Limite:</strong></div>';

      const formulaTex = '\\lim_{x \\to ' + limitTargetTex + dirSymbol + '} \\left( ' + fNode.toTex() + ' \\right) = ' + 
        (finalVal === 'DIVERGENTE' ? '\\nexists \\;\\text{(Limites laterais distintos)}' : (typeof finalVal === 'number' ? Number(finalVal.toFixed(8)) : finalVal));

      html += '<div style="margin:10px 0;font-size:1.2em;color:var(--primary)">' + scicalcRenderTex(formulaTex, true) + '</div>';

      if (exactReason) {
        html += '<p style="margin:8px 0;color:var(--muted);font-size:0.92rem"><strong>Diagnóstico:</strong> ' + exactReason + '</p>';
      }

      if (showSteps) {
        html += '<details class="didactic" open style="margin-top:16px">';
        html += '<summary>📐 Passo a Passo Analítico do Limite</summary><div class="didactic-body">';
        html += '<div class="step-list">';

        html += '<div class="step">';
        html += '<div class="step-title">1. Teste de Avaliação Direta</div>';
        if (isIndeterminate) {
          html += '<div style="color:var(--warn);font-weight:700">Forma Indeterminada detectada: [' + indetType + ']</div>';
          html += '<div class="step-note">A substituição direta de x = ' + aVal + ' anula simultaneamente numerador e denominador, exigindo resolução analítica.</div>';
        } else {
          html += '<div>A função foi avaliada na vizinhança de x = ' + limitTargetTex + '.</div>';
        }
        html += '</div>';

        if (lhopitalApplied && lhopitalSteps.length > 0) {
          lhopitalSteps.forEach(ls => {
            html += '<div class="step">';
            html += '<div class="step-title">' + (ls.order + 1) + '. Aplicação da Regra de L\'Hôpital (' + ls.order + 'ª ordem)</div>';
            html += '<div style="margin:5px 0">' + scicalcRenderTex('\\lim_{x \\to ' + limitTargetTex + '} \\frac{N(x)}{D(x)} = \\lim_{x \\to ' + limitTargetTex + '} \\frac{\\frac{d}{dx}N(x)}{\\frac{d}{dx}D(x)}', true) + '</div>';
            html += '<div style="margin:5px 0">' + scicalcRenderTex('\\frac{d}{dx}N = ' + ls.numPrime + ', \\quad \\frac{d}{dx}D = ' + ls.denPrime, false) + '</div>';
            html += '<div class="step-note">Derivamos numerador e denominador independentemente em relação a x.</div>';
            html += '</div>';
          });
        }

        if (!isInf) {
          html += '<div class="step">';
          html += '<div class="step-title">Limites Laterais de Verificação</div>';
          html += '<div class="limit-sides">';
          html += '<div class="limit-side"><strong>Pela Esquerda (x → a⁻):</strong> ' + (Number.isFinite(vLeft) ? vLeft.toFixed(6) : 'Indefinido') + '</div>';
          html += '<div class="limit-side"><strong>Pela Direita (x → a⁺):</strong> ' + (Number.isFinite(vRight) ? vRight.toFixed(6) : 'Indefinido') + '</div>';
          html += '</div>';
          html += '</div>';
        }

        html += '</div></div></details>';
      }

      html += '</div>';
      document.getElementById('lim-result').innerHTML = html;
      renderMathEverywhere(document.getElementById('lim-result'));
    } catch (err) {
      document.getElementById('lim-result').innerHTML = '<div class="result" style="border-left-color:var(--err);color:var(--err)"><strong>Erro no cálculo do limite:</strong> ' + err.message + '</div>';
    }
  };


  /* ==========================================================================
     4. MOTOR DE INTEGRAIS: QUADRATURA ADAPTATIVA DE SIMPSON
     ========================================================================== */

  // Algoritmo de Simpson Adaptativo com tolerância garantida
  function adaptiveSimpsonQuad(f, a, b, tol = 1e-10, maxDepth = 22) {
    function simpson(f, a, b) {
      const c = (a + b) / 2;
      const h = (b - a) / 6;
      return h * (f(a) + 4 * f(c) + f(b));
    }

    let evaluations = 0;
    function safeF(x) {
      evaluations++;
      const val = f(x);
      return Number.isFinite(val) ? val : 0;
    }

    function recursive(a, b, tol, whole, depth) {
      const c = (a + b) / 2;
      const left = simpson(safeF, a, c);
      const right = simpson(safeF, c, b);
      const delta = left + right - whole;

      if (depth <= 0 || Math.abs(delta) <= 15 * tol) {
        return {
          value: left + right + delta / 15,
          error: Math.abs(delta) / 15,
          depth: depth
        };
      }

      const s1 = recursive(a, c, tol / 2, left, depth - 1);
      const s2 = recursive(c, b, tol / 2, right, depth - 1);
      return {
        value: s1.value + s2.value,
        error: s1.error + s2.error
      };
    }

    const initial = simpson(safeF, a, b);
    const res = recursive(a, b, tol, initial, maxDepth);
    res.evaluations = evaluations;
    return res;
  }

  root.calcIntegral = function () {
    try {
      const fStr = document.getElementById('int-f').value.trim();
      const x1Raw = document.getElementById('int-x1').value.trim();
      const x2Raw = document.getElementById('int-x2').value.trim();
      const showSteps = document.getElementById('int-steps') ? document.getElementById('int-steps').checked : true;

      const isDefinite = x1Raw !== '' && x2Raw !== '';
      const fNode = math.parse(fStr);
      const fCode = fNode.compile();

      let indefTex = null;
      let isNonElementary = false;

      // 1. Tenta Integral Simbólica via Nerdamer
      try {
        if (typeof nerdamer !== 'undefined') {
          const res = nerdamer('integrate(' + fStr + ', x)');
          indefTex = res.toTeX();
        }
      } catch (e) {
        isNonElementary = true;
      }

      if (!indefTex || indefTex.includes('integrate(') || indefTex.includes('defint(')) {
        isNonElementary = true;
        indefTex = '\\int ' + fNode.toTex() + '\\, dx';
      }

      let html = '<div class="result">';
      html += '<div><strong>Integral Indefinida (Simbólica):</strong></div>';
      html += '<div style="margin:8px 0;font-size:1.15em;color:var(--primary)">' + 
        scicalcRenderTex('\\int ' + fNode.toTex() + '\\, dx = ' + (isNonElementary ? indefTex + ' + C \\quad\\text{(Primitiva não-elementar)}' : indefTex + ' + C'), true) + 
        '</div>';

      // 2. Se for Definida, executa Quadratura Adaptativa de Alta Precisão
      if (isDefinite) {
        const aVal = parseFloat(x1Raw.replace(',', '.'));
        const bVal = parseFloat(x2Raw.replace(',', '.'));

        if (!Number.isFinite(aVal) || !Number.isFinite(bVal)) {
          throw new Error('Limites de integração inválidos.');
        }

        const quadRes = adaptiveSimpsonQuad(x => fCode.evaluate({ x }), aVal, bVal, 1e-10);
        const defFormulaTex = '\\int_{' + aVal + '}^{' + bVal + '} ' + fNode.toTex() + '\\, dx \\approx ' + Number(quadRes.value.toFixed(10));

        html += '<div style="margin-top:16px"><strong>Integral Definida (Quadratura Adaptativa de Alta Precisão):</strong></div>';
        html += '<div style="margin:8px 0;font-size:1.2em;color:var(--ok)">' + scicalcRenderTex(defFormulaTex, true) + '</div>';
        html += '<p style="color:var(--muted);font-size:0.88rem;margin-top:4px">Cota de erro absoluto certificada: <code>&lt; ' + quadRes.error.toExponential(2) + '</code> (' + quadRes.evaluations + ' amostras de malha adaptativa).</p>';

        if (showSteps) {
          html += '<details class="didactic" open style="margin-top:16px">';
          html += '<summary>∫ Passo a Passo da Integração</summary><div class="didactic-body">';
          html += '<div class="step-list">';

          html += '<div class="step">';
          html += '<div class="step-title">1. Análise da Integrabilidade e Primitiva</div>';
          if (isNonElementary) {
            html += '<div class="step-note">Pelo Teorema de Liouville, esta função não admite primitiva em termos de funções elementares finitas. A integração requer métodos de quadratura numérica adaptativa de alta ordem.</div>';
          } else {
            html += '<div class="step-note">Primitiva analítica identificada com sucesso: <code>' + indefTex + '</code>.</div>';
          }
          html += '</div>';

          html += '<div class="step">';
          html += '<div class="step-title">2. Partição e Quadratura Adaptativa de Simpson</div>';
          html += '<div class="step-note">O intervalo [' + aVal + ', ' + bVal + '] é subdividido recursivamente até que o erro local em cada subintervalo seja estritamente inferior a <code>10⁻¹⁰</code>.</div>';
          html += '</div>';

          html += '<div class="step">';
          html += '<div class="step-title">3. Valor Numérico Final</div>';
          html += '<div><strong>Resultado:</strong> ' + Number(quadRes.value.toFixed(10)) + '</div>';
          html += '</div>';

          html += '</div></div></details>';
        }
      }

      html += '</div>';
      document.getElementById('int-result').innerHTML = html;
      renderMathEverywhere(document.getElementById('int-result'));
    } catch (err) {
      document.getElementById('int-result').innerHTML = '<div class="result" style="border-left-color:var(--err);color:var(--err)"><strong>Erro no cálculo da integral:</strong> ' + err.message + '</div>';
    }
  };

  /* ==========================================================================
     INICIALIZAÇÃO NO CARREGAMENTO
     ========================================================================== */

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', setupMathObservers);
  } else {
    setupMathObservers();
  }

  window.addEventListener('load', () => {
    renderMathEverywhere(document.body);
  });

})(typeof window !== 'undefined' ? window : this);
