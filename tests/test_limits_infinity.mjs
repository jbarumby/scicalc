import { test } from "node:test";
import assert from "node:assert/strict";
import "../scicalc-limit-rational-infinity.js";

const compute = globalThis.scicalcExactRationalInfinity;
const n = v => ({type:"ConstantNode",value:String(v)});
const x = {type:"SymbolNode",name:"x"};
const op = (op,...args) => ({type:"OperatorNode",op,args});
const plus=(a,b)=>op("+",a,b);
const minus=(a,b)=>op("-",a,b);
const times=(a,b)=>op("*",a,b);
const over=(a,b)=>op("/",a,b);
const pow=(a,b)=>op("^",a,n(b));

test("(2*x+1)/x at infinity is exactly 2",()=>{
  const exp=over(plus(times(n(2),x),n(1)),x);
  const result=compute(exp,"inf");
  assert.equal(result.value,"2");
  assert.equal(result.tex,"2");
  assert.equal(result.isExact,true);
});

test("rational polynomial equal degrees gives leading coefficient ratio",()=>{
  const numerator=plus(plus(times(n(3),pow(x,2)),times(n(5),x)),n(1));
  const denominator=plus(minus(times(n(2),pow(x,2)),times(n(7),x)),n(4));
  assert.equal(compute(over(numerator,denominator),"∞").value,"3/2");
});

test("constant plus rational term works without a hard-coded example",()=>{
  assert.equal(compute(plus(n(2),over(n(1),x)),"+infinity").value,"2");
});

test("lower polynomial degree yields zero at both infinities",()=>{
  assert.equal(compute(over(n(1),x),"-inf").value,"0");
  assert.equal(compute(over(n(1),pow(x,2)),"inf").value,"0");
});

test("odd degree difference changes sign at negative infinity",()=>{
  assert.equal(compute(over(pow(x,2),x),"inf").value,"+∞");
  assert.equal(compute(over(pow(x,2),x),"-∞").value,"−∞");
});

test("infinite ratio of same degree with negative sign",()=>{
  assert.equal(compute(over(times(n(-3),x),times(n(2),x)),"inf").value,"-3/2");
});

test("finite limits do not trigger infinity logic",()=>{
  assert.equal(compute(over(n(1),x),"0"),null);
});

test("transcendental expressions do not yield fake exact results",()=>{
  assert.throws(()=>compute({type:"FunctionNode",name:"sin",args:[x]},"inf"));
});
