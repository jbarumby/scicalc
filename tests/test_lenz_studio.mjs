import { test } from "node:test";
import assert from "node:assert/strict";
import "../scicalc-lenz-model.js";

const model=globalThis.scicalcLenzState;

test("approaching north magnet induces north face and counterclockwise current",()=>{
 const s=model("approach",true,.60,.35);
 assert.equal(s.inducedPole,"N");
 assert.equal(s.fieldDirection,"left");
 assert.equal(s.currentDirection,"counterclockwise");
 assert.equal(s.fluxTrend,"increases");
 assert.ok(s.needle<0);
});

test("retreating north magnet induces south face and clockwise current",()=>{
 const s=model("away",true,.60,.35);
 assert.equal(s.inducedPole,"S");
 assert.equal(s.fieldDirection,"right");
 assert.equal(s.currentDirection,"clockwise");
 assert.equal(s.fluxTrend,"decreases");
 assert.ok(s.needle>0);
});

test("stationary magnet gives no emf, no induced field or current",()=>{
 for(const mode of ["approach","away"]){
  const s=model(mode,false,.80,1);
  assert.equal(s.inducedPole,"—");
  assert.equal(s.currentDirection,"none");
  assert.equal(s.fieldDirection,"none");
  assert.equal(s.needle,0);
  assert.equal(s.response,0);
  assert.equal(s.fluxTrend,"constant");
 }
});

test("relative field flux increases monotonically when distance decreases",()=>{
 const far=model("approach",true,.15,.35);
 const near=model("approach",true,.85,.35);
 assert.ok(near.relativeFlux>far.relativeFlux);
 assert.ok(near.response>far.response);
});

test("higher speed increases relative induced response at fixed location",()=>{
 const slow=model("approach",true,.6,.05);
 const fast=model("approach",true,.6,1);
 assert.ok(fast.response>slow.response);
});

test("bounds prevent invalid flux or relative emf",()=>{
 const a=model("away",true,-3,-9),b=model("approach",true,3,9);
 for(const s of [a,b]){
  assert.ok(s.relativeFlux>=0&&s.relativeFlux<=1);
  assert.ok(s.response>=0&&s.response<=1);
  assert.ok(s.position>=0&&s.position<=1);
 }
});
