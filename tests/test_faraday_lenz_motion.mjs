import {test} from 'node:test';
import assert from 'node:assert/strict';
import '../scicalc-lenz-model.js';
import '../scicalc-faraday-model.js';

const L=globalThis.scicalcLenzAdvance;
const F=globalThis.SciCalcFaradayModel;

test('Lenz away always starts next to coil and visibly travels left',()=>{
  const start=globalThis.scicalcLenzStartPosition('away');
  assert.ok(start>.75);
  const next=L(start,'away',0.35,1);
  assert.ok(next.position<start);
  assert.equal(next.finished,false);
  let position=start;
  for(let k=0;k<300;k++){
    const step=L(position,'away',0.35,.05);
    position=step.position;
    if(step.finished)break;
  }
  assert.equal(position,.06);
});

test('Lenz approach starts away from coil and travels right',()=>{
 const start=globalThis.scicalcLenzStartPosition('approach');
 assert.ok(start<.2);
 assert.ok(L(start,'approach',0.4,1).position>start);
});

test('Faraday moving away increases gap and moving toward decreases it',()=>{
 const away=F.advance(.12,'away',.35,.5);
 const near=F.advance(.65,'approach',.35,.5);
 assert.ok(away.distance>.12);
 assert.ok(near.distance<.65);
 assert.equal(F.start('away'),.11);
 assert.equal(F.start('approach'),.67);
});

test('Faraday field and flux decrease with axial magnet distance',()=>{
 const p=F.parameters(200,.01,.2,.35);
 const near=F.state(p,.12,'approach',true),far=F.state(p,.65,'away',true);
 assert.ok(near.B>far.B);
 assert.ok(near.flux>far.flux);
 assert.ok(near.fluxLinkage>far.fluxLinkage);
});

test('Faraday induction sign and induced face invert when direction changes',()=>{
 const p=F.parameters(200,.01,.2,.35);
 const toward=F.state(p,.2,'approach',true),away=F.state(p,.2,'away',true);
 assert.ok(toward.emf<0);
 assert.ok(away.emf>0);
 assert.equal(toward.inducedFace,'N');
 assert.equal(away.inducedFace,'S');
 assert.equal(toward.fieldInduced,'left');
 assert.equal(away.fieldInduced,'right');
 assert.ok(Math.abs(toward.emf+away.emf)<1e-12);
});

test('Faraday stopped magnet gives precisely zero emf but nonzero static flux',()=>{
 const p=F.parameters(200,.01,.2,.35);
 const s=F.state(p,.2,'approach',false);
 assert.equal(s.emf,0);
 assert.equal(s.distanceSpeed,0);
 assert.equal(s.inducedFace,'—');
 assert.equal(s.fieldInduced,'none');
 assert.ok(s.B>0&&s.flux>0);
});

test('Faraday emf proportional to N and speed under same geometry',()=>{
 const p1=F.parameters(100,.01,.2,.2),p2=F.parameters(200,.01,.2,.4);
 const v1=F.state(p1,.2,'approach',true),v2=F.state(p2,.2,'approach',true);
 assert.ok(Math.abs(v2.emf/v1.emf-4)<1e-10);
});

test('Faraday validates physical parameters and keeps magnet inside apparatus',()=>{
 assert.throws(()=>F.parameters(0,.01,.2,.35),/espiras/);
 assert.throws(()=>F.parameters(200,-1,.2,.35),/Área/);
 assert.throws(()=>F.parameters(200,.01,.2,0),/Velocidade/);
 assert.equal(F.advance(.11,'away',1,100).distance,.71);
 assert.equal(F.advance(.67,'approach',1,100).distance,.09);
});
