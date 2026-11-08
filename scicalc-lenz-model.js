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
  root.scicalcLenzState = computeLenzState;
})(typeof window !== 'undefined' ? window : globalThis);
