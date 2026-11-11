import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';

const html=readFileSync('ciencia.html','utf8');
const engine=readFileSync('scicalc-engine-upgrade.js','utf8');
const picker=readFileSync('scicalc-experiment-picker.js','utf8');
const didactic=readFileSync('scicalc-didactic-restore.js','utf8');

test('new pedagogical modules load after advanced calculation and engine upgrade',()=>{
 const upgrade=html.indexOf('src="./scicalc-engine-upgrade.js"');
 const pro=html.indexOf('src="./scicalc-pro.js?build=20261111c"');
 const explain=html.indexOf('src="./scicalc-didactic-restore.js"');
 const choice=html.indexOf('src="./scicalc-experiment-picker.js"');
 assert.ok(pro>0&&upgrade>pro&&explain>upgrade&&choice>explain);
 assert.match(html,/scicalc-experiment-picker\.css/);
});

test('indefinite integrals also show mathematical steps when enabled',()=>{
 assert.match(engine,/if\(!isDefinite && showSteps\)/);
 assert.match(engine,/Passo a passo da integral indefinida/);
 assert.match(engine,/Integração por partes/);
 assert.match(engine,/Regra da potência/);
});

test('exact rational infinity solver takes priority over numeric approximation',()=>{
 assert.match(engine,/const previousLimitEngine = root\.calcLimite/);
 assert.match(engine,/root\.scicalcExactRationalInfinity\(math\.parse\(input\),target\)/);
 assert.match(engine,/previousLimitEngine\.apply\(this,arguments\)/);
});

test('numeric adaptive Simpson error is not described as certified',()=>{
 assert.doesNotMatch(engine,/Cota de erro absoluto certificada/);
 assert.match(engine,/Estimativa numérica de erro local \(não certificada\)/);
 assert.doesNotMatch(engine,/Pelo Teorema de Liouville, esta função não admite primitiva/);
});

test('experiment picker hides the selected cards before they are clicked',()=>{
 assert.match(picker,/document\.createElement\('details'\)/);
 assert.match(picker,/details\.forEach\(d=>d\.open=false\)/);
 assert.match(picker,/name:'🌀 Massa–mola amortecida'/);
 assert.match(picker,/key:'faraday'/);
 assert.match(picker,/key:'lenz'/);
 assert.match(picker,/openBox\(item\.box,item\.model,item\.key\)/);
 assert.match(picker,/root\.v26InitLab=function\(\)/);
});

test('Cálculo Pro uses verified symbolic steps returned by Python, not generic outline',()=>{
 const pro=readFileSync('scicalc-pro.js','utf8');
 assert.match(didactic,/input\.type='checkbox';input\.id='sc-pro-show-steps';input\.checked=true/);
 assert.match(didactic,/Mostrar resolução matemática passo a passo/);
 assert.doesNotMatch(didactic,/O SymPy usa simplificações e regras de limite/);
 assert.match(pro,/renderVerifiedMathSteps\(root,data\)/);
 assert.match(pro,/Array\.isArray\(data\.steps\)/);
 assert.match(pro,/item\.tex/);
 assert.match(pro,/sc-pro-verified-math-steps/);
});


test('Cálculo Pro scripts and Python core cannot remain stale in browser cache',()=>{
 const pro=readFileSync('scicalc-pro.js','utf8');
 const worker=readFileSync('scicalc-pro-worker.mjs','utf8');
 const core=readFileSync('scicalc-pro-core.py','utf8');
 assert.match(html,/src="\.\/scicalc-pro\.js\?build=20261111c"/);
 assert.match(pro,/scicalc-pro-worker\.mjs\?build=/);
 assert.match(pro,/Versão do motor SymPy:/);
 assert.match(worker,/scicalc-pro-core\.py\?build=20261111c/);
 assert.match(worker,/cache: "no-store"/);
 assert.match(core,/SCICALC_PRO_ENGINE_REV = '20261111c'/);
 assert.match(core,/integral_tex\+r'\\,\'/);
});
