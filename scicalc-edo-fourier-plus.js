/* SciCalc - EDO + Fourier, interface didática sobre núcleos testáveis. */
(function(root){
'use strict';
const K=root.SciCalcLabs,$=id=>document.getElementById(id);
if(!K||!$('edo')||!$('fourier'))return;
function num(id){const s=$(id).value.trim().replace(',','.');const n=Number(s);if(!s||!Number.isFinite(n))throw Error('Valor inválido em '+id);return n;}
function fmt(n){if(!Number.isFinite(n))return '—';return n===0?'0':(Math.abs(n)>1e5||Math.abs(n)<1e-4?n.toExponential(5):Number(n.toPrecision(9)).toLocaleString('pt-BR',{maximumFractionDigits:9}));}
function node(parent,tag,content,className){const x=document.createElement(tag);if(className)x.className=className;x.textContent=String(content||'');parent.appendChild(x);return x;}
function row(parent,label,value){const b=node(parent,'div','','sc-lab-metric');node(b,'small',label);node(b,'strong',value);}
function step(parent,title,value){const el=node(parent,'div','','sc-lab-step');node(el,'b',title);node(el,'p',value);}
function error(id,e){const out=$(id);out.replaceChildren();node(out,'div','⚠ '+String(e.message||e),'sc-lab-error').setAttribute('role','alert');}
function evaluator(s,vars){
 s=String(s).trim();
 if(!s||s.length>220||!/^[\dA-Za-z +*\/.^()\-%,]*$/.test(s))throw Error('Expressão não permitida: use funções e operadores matemáticos comuns.');
 const allowed=new Set(['sin','cos','tan','asin','acos','atan','exp','log','ln','sqrt','abs','sinh','cosh','tanh','pi','e','tau','sec','csc','cot','floor','ceil','round',...vars]);
 for(const name of s.match(/[A-Za-z][A-Za-z0-9_]*/g)||[])if(!allowed.has(name))throw Error('Identificador desconhecido: '+name);
 if(!root.math)throw Error('Motor Math.js indisponível.');
 const f=root.math.parse(s).compile();
 return (...args)=>{
  const scope={};vars.forEach((name,i)=>{scope[name]=args[i];});
  const v=f.evaluate(scope);
  if(typeof v!=='number'||!Number.isFinite(v))throw Error('Resultado não real/finito; revise a expressão ou domínio.');
  return v;
 };
}
function plot(id,traces,title,xt,yt){
 if(!root.Plotly)throw Error('Plotly não carregou.');
 return root.Plotly.react(id,traces,{
  title,xaxis:{title:xt},yaxis:{title:yt},
  font:{color:document.body.classList.contains('dark-mode')?'#e2e8f0':'#25314e'},
  paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)',
  margin:{t:48,r:12,l:55,b:50},legend:{orientation:'h',y:1.13}
 },{responsive:true,displaylogo:false});
}
function makeExamples(id,entries,fields){
 const wrapper=node($(id).querySelector('.panel'),'div','','sc-lab-options');
 const label=node(wrapper,'label','Experimente um exemplo:');
 const sel=node(wrapper,'select','');sel.setAttribute('aria-label','Exemplos');
 node(sel,'option','Selecione…').value='';
 Object.keys(entries).forEach(key=>{const option=node(sel,'option',key);option.value=key;});
 sel.addEventListener('change',()=>{
  const values=entries[sel.value];if(!values)return;
  fields.forEach((field,i)=>$(field).value=values[i]);
 });
 const go=$(id).querySelector('button.btn[onclick]');
 if(go)go.insertAdjacentElement('beforebegin',wrapper);
}
makeExamples('edo',{
 'Exponencial — y’ = y':['y','0','1','3','50'],
 'Decaimento — y’ = −2y+x':['-2*y+x','0','1','5','100'],
 'Logística — y’ = y(1−y)':['y*(1-y)','0','0.2','8','100'],
 'Resposta forçada — y’ = cos(x)−y':['cos(x)-y','0','0','12','160']
},['edo-f','edo-x0','edo-y0','edo-xf','edo-n']);
root.calcEDO=function(){
 try{
  const f=evaluator($('edo-f').value,['x','y']);
  const x0=num('edo-x0'),y0=num('edo-y0'),xf=num('edo-xf'),N=num('edo-n');
  const r=K.odeRefined(f,x0,y0,xf,N),a=r.coarse,b=r.fine;
  const out=$('edo-result');out.replaceChildren();
  const box=node(out,'div','','sc-lab-result');
  node(box,'h3','EDO — Runge–Kutta de 4ª ordem com refinamento');
  const cards=node(box,'div','','sc-lab-metrics');
  row(cards,'Resultado com 2N passos',fmt(b.end));
  row(cards,'Estimativa de erro final',fmt(r.estimatedError));
  row(cards,'Resultado com N passos',fmt(a.end));
  row(cards,'Passo h/2',fmt(b.h));
  node(box,'p','Estimativa |y₂N−yN|/15 baseada na convergência de ordem 4 do RK4. Não é uma cota rigorosa, principalmente perto de singularidades ou instabilidade.','sc-lab-note');
  if($('edo-steps')?.checked!==false){
   const details=node(box,'details','','sc-lab-details');details.open=true;
   node(details,'summary','📚 Resolução passo a passo');
   step(details,'1. Problema de valor inicial','y’ = '+$('edo-f').value+'; y('+fmt(x0)+') = '+fmt(y0)+'.');
   step(details,'2. Malha inicial','N = '+N+' e h = (xf−x₀)/N = '+fmt(a.h)+'.');
   const first=a.rows[0];
   step(details,'3. Quatro inclinações do RK4','k₁='+fmt(first.k1)+'; k₂='+fmt(first.k2)+'; k₃='+fmt(first.k3)+'; k₄='+fmt(first.k4)+'. y₁='+fmt(first.value)+'.');
   step(details,'4. Refinamento e comparação','Com 2N = '+2*N+', obtemos '+fmt(b.end)+'; diferença /15 ≈ '+fmt(r.estimatedError)+'.');
   const wrap=node(details,'div','','sc-lab-tablewrap'),table=node(wrap,'table','','sc-lab-table');
   const header=node(table,'tr');['i','xᵢ','yᵢ','k₁','k₂','k₃','k₄','yᵢ₊₁'].forEach(t=>node(header,'th',t));
   a.rows.forEach(z=>{const tr=node(table,'tr');[z.k,z.x,z.y,z.k1,z.k2,z.k3,z.k4,z.value].forEach(v=>node(tr,'td',fmt(v)));});
  }
  plot('edo-plot',[
   {x:a.xs,y:a.ys,mode:'lines',name:'N passos',line:{color:'#8b5cf6',width:2,dash:'dot'}},
   {x:b.xs,y:b.ys,mode:'lines',name:'2N passos',line:{color:'#06b6d4',width:3}}
  ],'Comparação do RK4: N e 2N passos','x','y(x)');
  if(root.salvarHistorico)root.salvarHistorico('EDO RK4 refinado',"y'="+$('edo-f').value,'y(xf)≈'+b.end);
 }catch(e){error('edo-result',e);}
};
makeExamples('fourier',{
 'Seno de 2 Hz':['sin(2*pi*2*t)','0','2','256'],
 '2 Hz + 5 Hz':['sin(2*pi*2*t)+0.5*cos(2*pi*5*t)','0','2','256'],
 'Portadora de 12 Hz':['2*cos(2*pi*12*t)','0','2','512'],
 'Decaimento exponencial':['exp(-3*t)','0','5','512']
},['four-f','four-tmin','four-tmax','four-n']);
const opts=node($('fourier').querySelector('.panel'),'div','','sc-lab-options');
const label=node(opts,'label','Janela espectral');
const win=node(opts,'select');win.id='four-plus-window';
node(win,'option','Hann (reduz vazamento)').value='hann';
node(win,'option','Retangular').value='rectangular';
const button=$('fourier').querySelector('button.btn[onclick*="calcFourier"]');if(button)button.insertAdjacentElement('beforebegin',opts);
const timeGraph=node($('fourier').querySelector('.panel'),'div','','plot');timeGraph.id='four-plus-time-plot';
$('four-plot').insertAdjacentElement('beforebegin',timeGraph);
root.calcFourier=function(){
 try{
  const signal=evaluator($('four-f').value,['t']);
  const spectrum=K.fourier(signal,num('four-tmin'),num('four-tmax'),num('four-n'),win.value);
  const s=spectrum,out=$('four-result');out.replaceChildren();
  const box=node(out,'div','','sc-lab-result');
  node(box,'h3','Fourier — FFT e espectro unilateral');
  const cards=node(box,'div','','sc-lab-metrics');
  row(cards,'Amostragem fs',fmt(s.fs)+' Hz');row(cards,'Resolução Δf',fmt(s.df)+' Hz');
  row(cards,'Frequência de Nyquist',fmt(s.nyquist)+' Hz');row(cards,'Janela',win.value==='hann'?'Hann':'Retangular');
  node(box,'p','Picos correspondem aos bins da FFT, não a frequências estimadas de forma contínua. As amplitudes foram corrigidas pelo ganho coerente da janela.','sc-lab-note');
  const peak=node(box,'div','','sc-lab-peaks');node(peak,'h4','Principais componentes');
  if(!s.peaks.length)node(peak,'p','Nenhum pico interno detectado.');
  s.peaks.slice(0,5).forEach((p,i)=>step(peak,(i+1)+'. '+fmt(p.freq)+' Hz = '+fmt(p.omega)+' rad/s','Amplitude unilateral ≈ '+fmt(p.amp)));
  if($('four-steps')?.checked!==false){
   const details=node(box,'details','','sc-lab-details');details.open=true;
   node(details,'summary','📚 Passo a passo da transformada');
   step(details,'1. Amostrar','N = '+s.n+', Δt = '+fmt(s.dt)+' s, intervalo sem incluir o extremo final.');
   step(details,'2. Escala de frequências','fs = '+fmt(s.fs)+' Hz; Nyquist = '+fmt(s.nyquist)+' Hz; espaçamento Δf = '+fmt(s.df)+' Hz.');
   step(details,'3. Janela',win.value==='hann'?'w[n] = 0,5−0,5cos(2πn/N), reduz vazamento espectral, mas alarga picos.':'w[n] = 1, sem atenuação dos extremos; pode ocorrer maior vazamento.');
   step(details,'4. FFT','X[k] = Σ f[n]w[n]exp(−2πikn/N). Algoritmo radix-2 O(N log N).');
   step(details,'5. Amplitude e unidades','A[k] = 2|X[k]|/Σw[n], exceto DC e Nyquist, que não são duplicados. Frequência em Hz; ω = 2πf em rad/s.');
   step(details,'6. Interpretação','Frequências acima de Nyquist sofrem aliasing e picos fora de bins exatos apresentam vazamento.');
  }
  plot('four-plus-time-plot',[{x:s.xs,y:s.ys,mode:'lines',name:'Sinal',line:{color:'#8b5cf6',width:2.5}}],
   'Sinal no tempo','t (s)','f(t)');
  plot('four-plot',[{x:s.freqs,y:s.amps,mode:'lines',fill:'tozeroy',name:'FFT',line:{color:'#06b6d4',width:2.5}}],
   'Espectro de amplitudes','Frequência (Hz)','Amplitude');
  if(root.salvarHistorico)root.salvarHistorico('Fourier FFT',$('four-f').value,'N='+s.n+', fs='+s.fs);
 }catch(e){error('four-result',e);}
};
})(window);
