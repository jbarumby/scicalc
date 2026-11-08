/*
 * SciCalc — limites de funções racionais quando x → ±infinito.
 * Aritmética racional exata (BigInt), sem inferir limite por amostragem.
 * Abrange operações +, −, ×, /, potências inteiras e parênteses.
 * Se a expressão exceder essa classe, preserva o resolvedor existente.
 */
(function(root){
  'use strict';
  const MAX_DEGREE=12, MAX_DIGITS=100;

  function gcd(a,b){
    a=a<0n?-a:a; b=b<0n?-b:b;
    while(b){const r=a%b;a=b;b=r;}return a||1n;
  }
  function Q(n,d=1n){
    if(d===0n)throw Error('divisão por zero');
    if(d<0n){n=-n;d=-d;}
    const g=gcd(n,d);
    n/=g;d/=g;
    if(String(n).length>MAX_DIGITS||String(d).length>MAX_DIGITS)throw Error('coeficientes excessivos');
    return {n,d};
  }
  const zero=()=>Q(0n),one=()=>Q(1n);
  function add(a,b){return Q(a.n*b.d+b.n*a.d,a.d*b.d);}
  function neg(a){return Q(-a.n,a.d);}
  function mul(a,b){return Q(a.n*b.n,a.d*b.d);}
  function div(a,b){return Q(a.n*b.d,a.d*b.n);}
  function readConstant(value){
    let s=String(value).trim();
    if(!/^[+-]?\d+(?:\.\d{1,9})?$/.test(s))return null;
    const negative=s.startsWith('-');
    if(s[0]==='-'||s[0]==='+')s=s.slice(1);
    const parts=s.split('.');
    const d=10n**BigInt((parts[1]||'').length);
    const n=BigInt(parts[0]+(parts[1]||''))*(negative?-1n:1n);
    return Q(n,d);
  }
  function trim(poly){
    while(poly.length>1&&poly[poly.length-1].n===0n)poly.pop();
    return poly;
  }
  function pAdd(a,b){
    const r=[];
    for(let i=0;i<Math.max(a.length,b.length);i++)r.push(add(a[i]||zero(),b[i]||zero()));
    return trim(r);
  }
  function pNeg(a){return a.map(neg);}
  function pMul(a,b){
    if(a.length+b.length-2>MAX_DEGREE)throw Error('grau máximo excedido');
    const r=Array.from({length:a.length+b.length-1},zero);
    a.forEach((v,i)=>b.forEach((w,j)=>{r[i+j]=add(r[i+j],mul(v,w));}));
    return trim(r);
  }
  function polyZero(p){return p.length===1&&p[0].n===0n;}
  function ratio(a,b){if(polyZero(b))throw Error('divisão por polinômio nulo');return {n:a,d:b};}
  const c=v=>ratio([v],[one()]);
  const x=()=>ratio([zero(),one()],[one()]);
  function rAdd(a,b){return ratio(pAdd(pMul(a.n,b.d),pMul(b.n,a.d)),pMul(a.d,b.d));}
  function rNeg(a){return ratio(pNeg(a.n),a.d);}
  function rMul(a,b){return ratio(pMul(a.n,b.n),pMul(a.d,b.d));}
  function rDiv(a,b){return ratio(pMul(a.n,b.d),pMul(a.d,b.n));}
  function rPow(a,e){
    if(!Number.isInteger(e)||Math.abs(e)>MAX_DEGREE)throw Error('potência não suportada');
    if(e<0){a=ratio(a.d,a.n);e=-e;}
    let result=c(one());
    for(let k=0;k<e;k++)result=rMul(result,a);
    return result;
  }
  function parseTree(node){
    if(!node)throw Error('árvore vazia');
    if(node.type==='ParenthesisNode')return parseTree(node.content);
    if(node.type==='ConstantNode'){
      const value=readConstant(node.value);
      if(!value)throw Error('constante não racional');
      return c(value);
    }
    if(node.type==='SymbolNode'&&node.name==='x')return x();
    if(node.type!=='OperatorNode')throw Error('fora da classe racional');
    const args=node.args||[];
    if(args.length===1&&node.op==='-')return rNeg(parseTree(args[0]));
    if(args.length===1&&node.op==='+')return parseTree(args[0]);
    if(args.length!==2)throw Error('operador não suportado');
    const a=parseTree(args[0]),b=parseTree(args[1]);
    if(node.op==='+')return rAdd(a,b);
    if(node.op==='-')return rAdd(a,rNeg(b));
    if(node.op==='*')return rMul(a,b);
    if(node.op==='/')return rDiv(a,b);
    if(node.op==='^'){
      if(b.n.length!==1||b.d.length!==1)throw Error('expoente variável');
      const coeff=div(b.n[0],b.d[0]);
      if(coeff.d!==1n||coeff.n < -12n || coeff.n > 12n)throw Error('expoente fora do intervalo');
      return rPow(a,Number(coeff.n));
    }
    throw Error('operador desconhecido');
  }
  function format(q){return q.d===1n?String(q.n):String(q.n)+'/'+String(q.d);}
  function tex(q){return q.d===1n?String(q.n):'\\frac{'+q.n+'}{'+q.d+'}';}
  function infinityValue(expression,target){
    const isNegative=/^\s*(-inf|-infinity|-∞)\s*$/i.test(String(target));
    const isPositive=/^\s*(\+?inf|\+?infinity|\+?∞)\s*$/i.test(String(target));
    if(!isNegative&&!isPositive)return null;
    const r=parseTree(expression);
    if(polyZero(r.n))return {value:'0',tex:'0',degN:0,degD:0,reason:'Numerador identicamente nulo.',isExact:true};
    const degN=r.n.length-1,degD=r.d.length-1;
    const leading=div(r.n[degN],r.d[degD]);
    const difference=degN-degD;
    if(difference<0){
      return {value:'0',tex:'0',degN,degD,reason:'O grau do denominador é maior que o do numerador.',isExact:true};
    }
    if(difference===0){
      return {value:format(leading),tex:tex(leading),degN,degD,reason:'Os graus são iguais: o limite é a razão dos coeficientes líderes.',isExact:true};
    }
    let positive=leading.n>0n;
    if(isNegative&&difference%2===1)positive=!positive;
    return {value:positive?'+∞':'−∞',tex:positive?'\\infty':'-\\infty',
      degN,degD,reason:'O numerador cresce com maior grau; o sinal depende do coeficiente líder e da direção.',isExact:true};
  }
  root.scicalcExactRationalInfinity=infinityValue;
  if(!root.document)return;
  const previous=root.calcLimite;
  root.calcLimite=function(){
    const field=document.getElementById('lim-f'),point=document.getElementById('lim-a');
    if(!field||!point||!root.math||typeof root.math.parse!=='function'){
      return typeof previous==='function'?previous.apply(this,arguments):undefined;
    }
    let answer,tree;
    try{
      tree=root.math.parse(field.value);
      answer=infinityValue(tree,point.value);
    }catch(error){answer=null;}
    if(!answer)return typeof previous==='function'?previous.apply(this,arguments):undefined;
    const container=document.getElementById('lim-result');
    if(!container)return typeof previous==='function'?previous.apply(this,arguments):undefined;
    container.replaceChildren();
    const result=document.createElement('div');
    result.className='result';
    const formula=document.createElement('div');
    const t=String(point.value).trim().startsWith('-')?'-\\infty':'\\infty';
    const latex='\\lim_{x\\to '+t+'}\\left('+tree.toTex()+'\\right)='+answer.tex;
    try{
      if(root.katex&&typeof root.katex.render==='function'){
        root.katex.render(latex,formula,{displayMode:true,throwOnError:true,trust:false});
      } else formula.textContent=latex;
    }catch(e){formula.textContent=latex;}
    result.appendChild(formula);
    const label=document.createElement('p');
    label.className='step-note';label.style.marginTop='12px';
    label.textContent='Resultado exato por análise do termo dominante. Não é uma aproximação numérica.';
    result.appendChild(label);
    const showSteps=document.getElementById('lim-steps')?.checked!==false;
    if(showSteps){
      const details=document.createElement('details');
      details.className='didactic';details.open=true;
      const summary=document.createElement('summary');
      summary.textContent='📚 Passo a passo — graus dos polinômios';
      const body=document.createElement('div');
      body.className='didactic-body';
      const p1=document.createElement('p');
      p1.textContent='Grau do numerador: '+answer.degN+'; grau do denominador: '+answer.degD+'.';
      const p2=document.createElement('p');
      p2.textContent=answer.reason;
      body.append(p1,p2);details.append(summary,body);result.appendChild(details);
    }
    container.appendChild(result);
    try{if(typeof root.salvarHistorico==='function')root.salvarHistorico('Limite',field.value+'; x→'+point.value,answer.value);}catch(e){}
    return answer;
  };
})(typeof window!=='undefined'?window:globalThis);
