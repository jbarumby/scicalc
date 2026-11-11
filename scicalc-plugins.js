/**
 * SciCalc Plugins & Interface Controller
 * Copyright (c) 2026 João Carlos Barumby. Todos os direitos reservados.
 * 
 * Gerenciador modular de plugins, atalhos de navegação e aprimoramento de UX/UI.
 */
(function (root) {
  'use strict';

  root.SciCalc = root.SciCalc || {};

  // Catálogo inicial de plugins e motores integrados
  const defaultPlugins = [
    {
      id: 'scicalc-pro',
      name: 'Cálculo Pro (SymPy)',
      icon: '🧠',
      category: 'Cálculo Simbólico',
      tabId: 'scicalc-pro',
      description: 'Motor SymPy em WebAssembly para cálculo exato de integrais impróprias, derivadas de ordem superior e limites bilaterais sem envio a servidores.',
      status: 'Ativo'
    },
    {
      id: 'lenz-studio',
      name: 'Lenz Studio Pro',
      icon: '⚡',
      category: 'Eletromagnetismo',
      tabId: 'simuladores-v26',
      targetId: 'lenz-studio-canvas',
      description: 'Simulador interativo em tempo real da Lei de Lenz: vetores de campo induzido, polo relativo e deflexão de galvanômetro.',
      status: 'Ativo'
    },
    {
      id: 'faraday-studio',
      name: 'Faraday Studio',
      icon: '🧲',
      category: 'Eletromagnetismo',
      tabId: 'fisica',
      description: 'Laboratório vetorial da Lei da Indução com taxa de variação de fluxo magnético e cálculo de força eletromotriz.',
      status: 'Ativo'
    },
    {
      id: 'edo-fourier-labs',
      name: 'Laboratório EDO & Fourier',
      icon: '📡',
      category: 'Equações Diferenciais',
      tabId: 'fourier',
      description: 'Decomposição harmônica espectral, coeficientes trigonométricos de Fourier e resolução numérica passo a passo de EDOs lineares.',
      status: 'Ativo'
    },
    {
      id: 'number-theory-plus',
      name: 'Teoria dos Números Plus',
      icon: '🔢',
      category: 'Aritmética & Primos',
      tabId: 'numeros',
      description: 'Crivos rápidos de primalidade, equações diofantinas, congruências modulares e cálculo de totiente de Euler.',
      status: 'Ativo'
    },
    {
      id: 'limit-rational-infinity',
      name: 'Limites Assintóticos',
      icon: '📐',
      category: 'Análise Matemática',
      tabId: 'limites',
      description: 'Avaliação exata de limites laterais e comportamento assintótico racional no infinito com certificado analítico.',
      status: 'Ativo'
    },
    {
      id: 'plotly-3d',
      name: 'Plotly 3D Surface Engine',
      icon: '📊',
      category: 'Visualização Científica',
      tabId: 'grafico3d',
      description: 'Motor acelerado por GPU para gráficos de superfícies paramétricas, curvas no espaço e campos escalares tridimensionais.',
      status: 'Ativo'
    },
    {
      id: 'symbolic-calc',
      name: 'Nerdamer & Math.js Core',
      icon: '∫',
      category: 'Álgebra Simbólica',
      tabId: 'integrais',
      description: 'Núcleo clássico de simplificação algébrica, derivadas parciais e cálculo integral analítico no cliente.',
      status: 'Ativo'
    }
  ];

  root.SciCalc.plugins = defaultPlugins.slice();

  // API para registrar novos plugins externamente
  root.SciCalc.registerPlugin = function (plugin) {
    if (!plugin || !plugin.id) return;
    const existingIndex = root.SciCalc.plugins.findIndex(p => p.id === plugin.id);
    if (existingIndex >= 0) {
      root.SciCalc.plugins[existingIndex] = plugin;
    } else {
      root.SciCalc.plugins.push(plugin);
    }
    updatePluginCountBadge();
    renderPluginCards();
  };

  function showToast(message) {
    let toast = document.getElementById('scicalc-toast');
    if (!toast) {
      toast = document.createElement('div');
      toast.id = 'scicalc-toast';
      toast.className = 'sc-toast';
      document.body.appendChild(toast);
    }
    toast.innerHTML = '<span>📋</span> ' + message;
    toast.classList.add('show');
    clearTimeout(toast._hideTimer);
    toast._hideTimer = setTimeout(() => {
      toast.classList.remove('show');
    }, 2400);
  }

  function activateTab(tabId) {
    const nav = document.querySelector('nav.tabs');
    const main = document.querySelector('main');
    if (!nav || !main) return;

    const targetTab = nav.querySelector(`[data-tab="${tabId}"]`);
    const targetContent = document.getElementById(tabId) || main.querySelector(`.tab-content[id="${tabId}"]`);

    if (targetTab && targetContent) {
      nav.querySelectorAll('.tab').forEach(el => el.classList.remove('active'));
      main.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
      targetTab.classList.add('active');
      targetContent.classList.add('active');
      targetTab.scrollIntoView({ behavior: 'smooth', block: 'nearest', inline: 'center' });
    }
  }

  function updatePluginCountBadge() {
    const badge = document.querySelector('.plugin-count-badge');
    if (badge) {
      badge.textContent = root.SciCalc.plugins.length;
    }
  }

  function renderPluginCards(filterText = '') {
    const grid = document.getElementById('sc-plugins-grid');
    if (!grid) return;

    const term = filterText.toLowerCase().trim();
    const filtered = root.SciCalc.plugins.filter(p => {
      return !term || 
        p.name.toLowerCase().includes(term) || 
        p.category.toLowerCase().includes(term) || 
        p.description.toLowerCase().includes(term);
    });

    if (filtered.length === 0) {
      grid.innerHTML = '<div style="grid-column: 1/-1; padding: 24px; text-align: center; color: var(--muted);">Nenhum plugin ou módulo encontrado para o termo pesquisado.</div>';
      return;
    }

    grid.innerHTML = filtered.map(p => `
      <div class="sc-plugin-card" data-plugin-id="${p.id}">
        <div class="sc-plugin-card-top">
          <div class="sc-plugin-icon">${p.icon}</div>
          <div class="sc-plugin-info">
            <span class="sc-plugin-category">${p.category}</span>
            <h4>${p.name}</h4>
            <p class="sc-plugin-desc">${p.description}</p>
          </div>
        </div>
        <div class="sc-plugin-card-bottom">
          <span class="sc-plugin-status-pill">
            <span class="sc-plugin-status-dot"></span> ${p.status || 'Ativo'}
          </span>
          <button type="button" class="sc-plugin-action-btn" data-action="open-tab" data-tab="${p.tabId}">
            Abrir Módulo →
          </button>
        </div>
      </div>
    `).join('');

    grid.querySelectorAll('[data-action="open-tab"]').forEach(btn => {
      btn.addEventListener('click', () => {
        const dialog = document.getElementById('scicalc-plugins-dialog');
        if (dialog && typeof dialog.close === 'function') dialog.close();
        const tabId = btn.dataset.tab;
        if (tabId) activateTab(tabId);
      });
    });
  }

  function setupPluginsModal() {
    if (document.getElementById('scicalc-plugins-dialog')) return;

    const dialog = document.createElement('dialog');
    dialog.id = 'scicalc-plugins-dialog';
    dialog.className = 'sc-modal-dialog';
    dialog.innerHTML = `
      <div class="sc-modal-header">
        <h3><span>🔌</span> Central de Plugins & Módulos Científicos</h3>
        <button type="button" class="sc-modal-close" id="sc-plugins-close" aria-label="Fechar">✕</button>
      </div>
      <div class="sc-modal-body">
        <div class="sc-plugin-filter-bar">
          <input id="sc-plugin-search" type="search" placeholder="🔍 Filtrar plugins por nome, cálculo ou física..." spellcheck="false">
        </div>
        <div class="sc-plugin-grid" id="sc-plugins-grid"></div>
      </div>
    `;

    document.body.appendChild(dialog);

    const closeBtn = document.getElementById('sc-plugins-close');
    closeBtn.addEventListener('click', () => dialog.close());

    dialog.addEventListener('click', (e) => {
      if (e.target === dialog) dialog.close();
    });

    const searchInput = document.getElementById('sc-plugin-search');
    searchInput.addEventListener('input', (e) => {
      renderPluginCards(e.target.value);
    });

    renderPluginCards();
  }

  function setupHeaderButton() {
    const actions = document.querySelector('header .header-actions');
    if (!actions || document.getElementById('btn-open-plugins')) return;

    const btn = document.createElement('button');
    btn.id = 'btn-open-plugins';
    btn.type = 'button';
    btn.className = 'plugins-toggle';
    btn.innerHTML = `<span>🔌</span> Plugins <span class="plugin-count-badge">${root.SciCalc.plugins.length}</span>`;
    
    btn.addEventListener('click', () => {
      setupPluginsModal();
      const dialog = document.getElementById('scicalc-plugins-dialog');
      if (dialog && typeof dialog.showModal === 'function') {
        renderPluginCards();
        dialog.showModal();
      }
    });

    actions.insertBefore(btn, actions.firstChild);
  }

  // Aprimora as caixas de resultado com botão 'Copiar'
  function attachCopyButtons() {
    const containers = document.querySelectorAll('.result, .sc-pro-result, .nt-plus-result, .sc-lab-result');
    containers.forEach(container => {
      if (container.querySelector('.sc-copy-btn')) return;

      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'sc-copy-btn';
      btn.title = 'Copiar resultado para a área de transferência';
      btn.innerHTML = '📋 Copiar';

      btn.addEventListener('click', (e) => {
        e.stopPropagation();
        const textToCopy = container.innerText.replace(/📋\s*Copiar.*/i, '').trim();
        if (!textToCopy) return;

        navigator.clipboard.writeText(textToCopy).then(() => {
          btn.innerHTML = '✓ Copiado!';
          btn.classList.add('copied');
          showToast('Resultado copiado para a área de transferência!');
          setTimeout(() => {
            btn.innerHTML = '📋 Copiar';
            btn.classList.remove('copied');
          }, 2000);
        }).catch(() => {
          showToast('Não foi possível copiar automaticamente.');
        });
      });

      container.appendChild(btn);
    });
  }

  // Monitora alterações no DOM para anexar botões de cópia dinamicamente
  function observeResults() {
    const observer = new MutationObserver(() => {
      attachCopyButtons();
    });
    observer.observe(document.body, { childList: true, subtree: true });
  }

  // Inicialização no carregamento
  function init() {
    setupHeaderButton();
    setupPluginsModal();
    attachCopyButtons();
    observeResults();
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})(typeof window !== 'undefined' ? window : this);
