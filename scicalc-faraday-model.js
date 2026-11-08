/* SciCalc — modelo simplificado do experimento de Faraday.
   Alinhamento axial, campo aproximadamente uniforme nas espiras.
   B(d)=B0[1+(d/a)^2]^(-3/2); a=max(raio de bobina, 0,045m).
   Fluxo por espira: phi=A*B; ligação de fluxo: N*phi.
   FEM=-d(N*phi)/dt. Sem resistência, corrente NÃO é determinada.
*/
(function(root){
  'use strict';
  const clamp=(x,lo,hi)=>Math.max(lo,Math.min(hi,Number(x)));
  function parameters(N,A,B0,speed){
    N=Number(N);A=Number(A);B0=Number(B0);speed=Number(speed);
    if(!Number.isInteger(N)||N<1||N>5000)throw Error('Use N inteiro entre 1 e 5.000 espiras.');
    if(!Number.isFinite(A)||A<0.0001||A>0.3)throw Error('Área entre 0,0001 e 0,3 m².');
    if(!Number.isFinite(B0)||B0<=0||B0>5)throw Error('Campo máximo entre 0 e 5 T.');
    if(!Number.isFinite(speed)||speed<=0||speed>2)throw Error('Velocidade maior que zero e até 2 m/s.');
    return {N,A,B0,speed,a:Math.max(.045,Math.sqrt(A/Math.PI))};
  }
  function state(params,distance,mode,moving){
    const d=clamp(distance,.09,.71),a=params.a;
    const ratio=d/a;
    const B=params.B0/Math.pow(1+ratio*ratio,1.5);
    const derivative=-3*params.B0*d/(a*a)*Math.pow(1+ratio*ratio,-2.5);
    const v=moving ? (mode==='away'?params.speed:-params.speed):0;
    const emf=-params.N*params.A*derivative*v;
    const phi=params.A*B;
    return {distance:d,B,flux:phi,fluxLinkage:params.N*phi,
      emf,distanceSpeed:v, moving:!!moving,
      fieldInduced:!moving?'none':emf<0?'left':'right',
      inducedFace:!moving?'—':emf<0?'N':'S',
      fluxTrend:!moving?'constant':mode==='away'?'decreases':'increases'};
  }
  function advance(distance,mode,speed,seconds){
    const dt=Number(seconds);if(!Number.isFinite(dt)||dt<0)throw Error('Tempo inválido.');
    const v=mode==='away'?1:-1;
    const d=clamp(distance+v*Number(speed)*dt,.09,.71);
    return {distance:d,finished:mode==='away'?d>=.71:d<=.09};
  }
  function start(mode){return mode==='away'?.11:.67;}
  root.SciCalcFaradayModel={parameters,state,advance,start};
})(typeof window!=='undefined'?window:globalThis);
