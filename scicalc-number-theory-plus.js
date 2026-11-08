/* SciCalc — Teoria dos Números: conferência independente de todos os módulos.
   Preserva os algoritmos e as demonstrações originais; acrescenta análise validada. */
(function(root){
'use strict';
const k=root.SciCalcLabs;
const section=document.getElementById('numeros');
if(!section||!k||section.dataset.labsNt)return;
section.dataset.labsNt='1';
const $=id=>document.getElementById(id);
const fmt=x=>String(x);
const add=(parent,tag,value,cls)=>{
 const el=document.createElement(tag);if(cls)el.className=cls;
 el.textContent=String(value);parent.appendChild(el);return el;
};
function card(title){
 const box=document.createElement('div');box.className='nt-plus-result';
 add(box,'h4',title);
 return box;
}
function line(box,key,value){
 const row=document.createElement('div');row.className='nt-plus-line';
 add(row,'b',key);add(row,'span',value);box.appendChild(row);
}
function explain(box,title,description){
 const wrap=document.createElement('div');wrap.className='nt-plus-explain';
 add(wrap,'b',title);add(wrap,'p',description);box.appendChild(wrap);
}
function safe(run,container){
 container.replaceChildren();
 try{container.appendChild(run());}catch(e){
  const box=card('Não foi possível analisar a entrada');
  add(box,'p',String(e.message||e),'nt-plus-error');container.appendChild(box);
 }
}
function divisors(n){
 const arr=[];
 for(let d=1n;d*d<=n;d++)if(n%d===0n){arr.push(d);if(d*d!==n)arr.push(n/d);}
 arr.sort((a,b)=>a<b?-1:a>b?1:0);
 return arr;
}
function diagnose(type){
 if(type==='prime'){
  const n=k.bi($('num-teste').value,'n',1000000000n),prime=k.isPrime(n);
  const b=card('Conferência independente de primalidade');
  line(b,'Número analisado',n);
  line(b,'Resultado',n<2n?'Não é primo':prime?'Primo':'Composto');
  explain(b,'Critério de divisibilidade','Basta testar os possíveis divisores primos até a raiz quadrada de n. O procedimento é exato nesta faixa de inteiros.');
  if(n>1n&&!prime){
   const f=k.factor(n);line(b,'Fatoração',f.map(x=>x.prime+(x.exp>1?'^'+x.exp:'')).join(' · '));
  }else if(prime){line(b,'Verificação','Nenhum divisor positivo próprio encontrado.');}
  return b;
 }
 if(type==='sieve'){
  const n=Number($('num-crivo').value);
  const primes=k.sieve(n),b=card('Análise do crivo de Eratóstenes');
  line(b,'Primos encontrados',primes.length);
  line(b,'Densidade até N',(100*primes.length/n).toFixed(2)+'%');
  line(b,'Primeiros primos',primes.slice(0,16).join(', '));
  line(b,'Últimos primos',primes.slice(-10).join(', '));
  explain(b,'Como o crivo funciona','Eliminamos os múltiplos de cada primo p começando em p². Os compostos menores já foram eliminados por um divisor menor.');
  return b;
 }
 if(type==='gcd'||type==='lcm'){
  const va=type==='gcd'?'num-a':'num-c',vb=type==='gcd'?'num-b':'num-d';
  const a=k.bi($(va).value,'a'),b=k.bi($(vb).value,'b'),e=k.egcd(a,b);
  if(a===0n&&b===0n)throw Error('MDC(0,0) não é definido nesta calculadora.');
  const box=card('Algoritmo de Euclides e identidade de Bézout');
  line(box,'MDC',e.g);
  line(box,'Coeficientes de Bézout','x = '+e.x+', y = '+e.y);
  explain(box,'Verificação exata',a+' · ('+e.x+') + '+b+' · ('+e.y+') = '+e.g);
  if(type==='lcm'){
   const product=a*b;const l=(a===0n||b===0n)?0n:(product<0n?-product:product)/e.g;
   line(box,'MMC',l);
   explain(box,'Identidade','mmc(a,b) × mdc(a,b) = |a×b|. Se a ou b for zero, mmc = 0.');
  }
  const short=e.steps.slice(0,7);
  if(short.length)explain(box,'Primeiras divisões euclidianas',short.map(x=>x.a+' = '+x.q+' × '+x.b+' + '+x.remainder).join(' ; '));
  return box;
 }
 if(type==='factor'||type==='arithmetic'){
  const id=type==='factor'?'num-fatorar':'num-func-n', n=k.bi($(id).value,'n',1000000000n);
  const a=k.arithmetic(n);
  const b=card(type==='factor'?'Fatoração e propriedades dos divisores':'Verificação das funções aritméticas');
  line(b,'Decomposição',a.factors.length?a.factors.map(f=>f.prime+(f.exp>1?'^'+f.exp:'')).join(' × '):'1 (produto vazio)');
  const ds=divisors(n);
  line(b,'Divisores positivos (quantidade)',ds.length);
  line(b,'Divisores',ds.length<=48?ds.join(', '):ds.slice(0,24).join(', ')+' … '+ds.slice(-8).join(', '));
  line(b,'φ(n) — números coprimos a n',a.phi);
  line(b,'μ(n) — Möbius',a.mu);
  line(b,'τ(n) — número de divisores',a.tau);
  line(b,'σ(n) — soma de divisores',a.sigma);
  const actualSum=ds.reduce((s,v)=>s+v,0n);
  line(b,'Conferência por enumeração',actualSum.toString()===a.sigma&&ds.length===Number(a.tau)?'τ e σ coincidem exatamente.':'Divergência encontrada: revisar cálculo.');
  explain(b,'Interpretação','φ(n) conta os coprimos; μ(n) detecta fatores quadráticos; τ(n) conta os divisores; σ(n) soma os divisores positivos.');
  return b;
 }
 if(type==='fermat'){
  const p=k.bi($('ferm-p').value,'p',1000000n),a=k.bi($('ferm-a').value,'a');
  const b=card('Conferência modular do Teorema de Fermat');
  if(!k.isPrime(p)){line(b,'Hipótese',p+' não é primo.');explain(b,'Atenção','O Pequeno Teorema de Fermat não se aplica na forma a^(p−1) ≡ 1 a este p.');return b;}
  const coprime=k.gcd(a,p)===1n;
  line(b,'p é primo?', 'Sim');
  line(b,'mdc(a,p)',k.gcd(a,p));
  line(b,'a^(p−1) mod p',k.powmod(a,p-1n,p));
  line(b,'a^p mod p',k.powmod(a,p,p));
  line(b,'a mod p',k.mod(a,p));
  explain(b,'Interpretação',coprime?
    'Como p é primo e mdc(a,p)=1, vale a^(p−1) ≡ 1 (mod p). A igualdade a^p ≡ a também foi verificada.':
    'p divide a: a^(p−1) ≡ 1 não se aplica, mas a^p ≡ a (mod p) continua válida.');
  return b;
 }
 throw Error('Módulo desconhecido.');
}
const buttons=[
 ['num-teste','prime'],['num-crivo','sieve'],['num-a','gcd'],['num-c','lcm'],
 ['num-fatorar','factor'],['num-func-n','arithmetic'],['ferm-p','fermat']
];
for(const [id,type] of buttons){
 const input=$(id);if(!input)continue;
 const cardEl=input.closest('.nt-card')||input.closest('.nt-theorem')||input.parentElement;
 const btn=cardEl.querySelector('button.btn');
 if(!btn)continue;
 const host=document.createElement('div');host.className='nt-plus-host';
 btn.insertAdjacentElement('afterend',host);
 const label=add(host,'div','Conferência independente e interpretação','nt-plus-label');
 label.setAttribute('aria-live','polite');
 const render=()=>safe(()=>diagnose(type),host);
 btn.addEventListener('click',render);
 const extra=document.createElement('button');extra.type='button';extra.className='nt-plus-btn';
 extra.textContent='🔬 Conferência avançada';
 host.before(extra);extra.addEventListener('click',render);
}

const fermPanel=$('sub-num-ferm');
if(fermPanel){
 const advanced=document.createElement('button');
 advanced.type='button';
 advanced.className='nt-plus-btn';
 advanced.textContent='🔬 Conferir o Teorema de Fermat';
 const out=document.createElement('div');
 out.className='nt-plus-host';
 const target=$('ferm-result')||fermPanel;
 target.insertAdjacentElement('afterend',advanced);
 advanced.insertAdjacentElement('afterend',out);
 advanced.addEventListener('click',()=>safe(()=>diagnose('fermat'),out));
}

// Fermat tem um botão principal que abre o iframe; o relatório complementar permanece opcional.
const tabs=section.querySelector('.number-tabs');
const panel=section.querySelector('.number-theory-panel');
if(!tabs||!panel)return;
const newTab=document.createElement('button');
newTab.type='button';newTab.className='sub-tab';newTab.textContent='≡ Congruências e inversos';
newTab.setAttribute('data-sub','num-mod-plus');tabs.appendChild(newTab);
const newPanel=document.createElement('div');
newPanel.className='sub-panel nt-plus-mod-panel';newPanel.id='sub-num-mod-plus';
newPanel.innerHTML=[
 '<div class="nt-plus-mod-intro"><h3>Congruências, inversos e Teorema Chinês dos Restos</h3>',
 '<p>Operações com inteiros exatos, incluindo módulos não coprimos quando o sistema possui solução.</p></div>',
 '<div class="nt-plus-mod-grid">',
   '<form id="nt-plus-inv-form" class="nt-plus-mod-card"><h4>Inverso modular e Bézout</h4>',
     '<label for="nt-plus-inv-a">Inteiro a</label><input id="nt-plus-inv-a" value="17">',
     '<label for="nt-plus-inv-m">Módulo m</label><input id="nt-plus-inv-m" value="43">',
     '<button class="btn" type="submit">Calcular inverso</button><div id="nt-plus-inv-out" aria-live="polite"></div>',
   '</form>',
   '<form id="nt-plus-crt-form" class="nt-plus-mod-card"><h4>Sistema de congruências</h4>',
     '<label for="nt-plus-crt-a">x ≡ a (mod m)</label>',
     '<div class="nt-plus-input-pair"><input id="nt-plus-crt-a" aria-label="Residuo a" value="2"><input id="nt-plus-crt-m" aria-label="Módulo m" value="3"></div>',
     '<label for="nt-plus-crt-b">x ≡ b (mod n)</label>',
     '<div class="nt-plus-input-pair"><input id="nt-plus-crt-b" aria-label="Residuo b" value="3"><input id="nt-plus-crt-n" aria-label="Módulo n" value="5"></div>',
     '<button class="btn" type="submit">Resolver sistema (CRT)</button><div id="nt-plus-crt-out" aria-live="polite"></div>',
   '</form>',
 '</div>'
].join('');
panel.appendChild(newPanel);
newTab.addEventListener('click',function(){
 tabs.querySelectorAll('.sub-tab').forEach(b=>b.classList.remove('active'));
 panel.querySelectorAll('.sub-panel').forEach(p=>p.classList.remove('active'));
 newTab.classList.add('active');newPanel.classList.add('active');
});
tabs.querySelectorAll('.sub-tab:not([data-sub="num-mod-plus"])').forEach(b=>b.addEventListener('click',()=>newTab.classList.remove('active')));
$('nt-plus-inv-form').addEventListener('submit',function(e){
 e.preventDefault();safe(()=>{
  const a=k.bi($('nt-plus-inv-a').value,'a'),m=k.bi($('nt-plus-inv-m').value,'m',1000000n);
  if(m<=1n)throw Error('O módulo deve ser maior que 1.');
  const ext=k.egcd(a,m),box=card('Inverso módulo '+m);
  line(box,'mdc(a,m)',ext.g);
  line(box,'Identidade de Bézout',a+'×('+ext.x+') + '+m+'×('+ext.y+') = '+ext.g);
  if(ext.g===1n){const inv=k.mod(ext.x,m);line(box,'Inverso de a',inv);line(box,'Conferência',a+'×'+inv+' ≡ 1 (mod '+m+')');}
  else explain(box,'Inverso inexistente','O inverso existe se, e somente se, mdc(a,m)=1.');
  return box;
 },$('nt-plus-inv-out'));
});
$('nt-plus-crt-form').addEventListener('submit',function(e){
 e.preventDefault();safe(()=>{
  const a=$('nt-plus-crt-a').value,b=$('nt-plus-crt-b').value,
    m=$('nt-plus-crt-m').value,n=$('nt-plus-crt-n').value;
  const result=k.crt(a,m,b,n),box=card('Teorema Chinês dos Restos');
  line(box,'mdc(m,n)',result.g);
  if(result.possible){line(box,'Solução','x ≡ '+result.value+' (mod '+result.period+')');explain(box,'Verificação','A solução satisfaz simultaneamente as duas congruências. Todas as soluções diferem de um múltiplo do MMC(m,n).');}
  else explain(box,'Sistema incompatível','A diferença entre os resíduos não é divisível por mdc(m,n). Não existe solução simultânea.');
  return box;
 },$('nt-plus-crt-out'));
});
})(window);
