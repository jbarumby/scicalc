import {test} from 'node:test';
import assert from 'node:assert/strict';
import '../scicalc-science-kernels.js';
const k=globalThis.SciCalcLabs;

test('primes and factors',()=>{
 assert.equal(k.isPrime(97),true);assert.equal(k.isPrime(91),false);
 assert.deepEqual(k.factor(360),[{prime:'2',exp:3},{prime:'3',exp:2},{prime:'5',exp:1}]);
 assert.equal(k.sieve(100).length,25);
 assert.throws(()=>k.bi('9007199254740993'),/módulo máximo/);
});
test('exact gcd, Bezout and CRT',()=>{
 const e=k.egcd(84n,30n);
 assert.equal(e.g,6n);assert.equal(84n*e.x+30n*e.y,6n);
 assert.deepEqual(k.crt('2','3','3','5'),{possible:true,value:'8',period:'15',g:'1'});
 assert.equal(k.crt('0','2','1','4').possible,false);
 assert.equal(k.powmod(3n,6n,7n),1n);
});
test('classical arithmetic functions',()=>{
 assert.deepEqual(k.arithmetic('60'),{
  n:'60',factors:[{prime:'2',exp:2},{prime:'3',exp:1},{prime:'5',exp:1}],
  phi:'16',tau:'12',sigma:'168',mu:'0'
 });
 assert.deepEqual(k.arithmetic('1'),{n:'1',factors:[],phi:'1',tau:'1',sigma:'1',mu:'1'});
});
test('RK4 step-doubling converges for y prime = y',()=>{
 const r=k.odeRefined((x,y)=>y,0,1,1,40);
 assert.ok(Math.abs(r.fine.end-Math.E)<1e-7);
 assert.ok(r.estimatedError<1e-7);
 assert.equal(r.coarse.rows[0].k,1);
});
test('RK4 backward integration',()=>{
 const r=k.odeRefined((x,y)=>y,1,Math.E,0,40);
 assert.ok(Math.abs(r.fine.end-1)<1e-6);
});
test('RK4 refuses singular/nonfinite functions',()=>{
 assert.throws(()=>k.odeRefined(()=>Infinity,0,1,1,50),/não|fora|finito|domínio/i);
});
test('FFT identifies two harmonic tones with Hann window',()=>{
 const r=k.fourier(t=>Math.sin(2*Math.PI*2*t)+.5*Math.cos(2*Math.PI*5*t),0,2,256,'hann');
 assert.equal(r.fs,128);assert.equal(r.df,.5);assert.equal(r.nyquist,64);
 const near=(x,y)=>Math.abs(x-y)<1e-10;
 assert.ok(near(r.amps[4],1));
 assert.ok(near(r.amps[10],.5));
 assert.equal(r.peaks[0].freq,2);
});
test('FFT DC and Nyquist normalization',()=>{
 const constant=k.fourier(()=>3,0,1,64,'rectangular');
 assert.ok(Math.abs(constant.amps[0]-3)<1e-12);
 assert.ok(constant.amps.slice(1).every(v=>v<1e-10));
 const cosine=k.fourier(t=>2*Math.cos(2*Math.PI*8*t),0,1,64,'rectangular');
 assert.ok(Math.abs(cosine.amps[8]-2)<1e-10);
});
test('invalid frequency resolution and sample count rejected',()=>{
 assert.throws(()=>k.fourier(t=>t,0,1,250,'hann'),/potência de 2/);
 assert.throws(()=>k.fourier(()=>NaN,0,1,64,'hann'),/fora do domínio/);
});
