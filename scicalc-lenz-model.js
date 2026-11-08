/* SciCalc — Lei de Lenz: modelo qualitativo de fluxo e corrente induzida.
   Convenção: polo N do ímã voltado à face esquerda da bobina.
   Corrente observada do lado do ímã. */
(function (root) {
  'use strict';
  function clamp(value, min, max) {
    return Math.min(max, Math.max(min, Number.isFinite(Number(value)) ? Number(value) : min));
  }
  function computeLenzState(mode, moving, position, speed) {
    const approach = mode !== 'away';
    const p = clamp(position, 0, 1);
    const v = clamp(speed, 0, 1);
    // Modelo qualitativo monotônico (não calibra campo, fluxo nem tensão em SI).
    const gap = 1 + 4 * (1 - p);
    const relativeFlux = 1 / (gap * gap * gap);
    const response = moving ? Math.min(1, (3 * v) / Math.pow(gap, 4)) : 0;
    return {
      mode: approach ? 'approach' : 'away',
      moving: Boolean(moving),
      position: p,
      speed: v,
      relativeFlux,
      response,
      inducedPole: moving ? (approach ? 'N' : 'S') : '—',
      fieldDirection: moving ? (approach ? 'left' : 'right') : 'none',
      currentDirection: moving ? (approach ? 'counterclockwise' : 'clockwise') : 'none',
      fluxTrend: moving ? (approach ? 'increases' : 'decreases') : 'constant',
      needle: moving ? (approach ? -response : response) : 0
    };
  }

  // Normalized, testable movement model used by the animation.
  // p=0 means the magnet is far from the coil; p=1 means near the coil.
  function advanceLenz(position, mode, speed, dt) {
    if(!Number.isFinite(dt) || dt<0)throw Error('Tempo de animação inválido.');
    const direction=mode==='away'?-1:1;
    const p=clamp(position,0.06,0.94);
    const delta=clamp(speed,0.05,1)*Math.min(dt,0.25)*0.46;
    const next=clamp(p+direction*delta,0.06,0.94);
    return {position:next,finished:direction<0?next<=0.06:next>=0.94};
  }
  function startPosition(mode) {
    // Always begin a new departure near the coil, or an approach far away.
    return mode==='away'?0.88:0.12;
  }
  root.scicalcLenzAdvance=advanceLenz;
  root.scicalcLenzStartPosition=startPosition;

  root.scicalcLenzState = computeLenzState;
})(typeof window !== 'undefined' ? window : globalThis);
