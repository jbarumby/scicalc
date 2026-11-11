/* SciCalc — restauro do passo a passo no Cálculo Pro e abertura das explicações.
 * Acrescenta raciocínio pedagógico visível sem afirmar que o CAS forneceu
 * uma derivação intermediária quando retornou apenas a resposta final.
 */
(function(root){
'use strict';
const $=id=>document.getElementById(id);
const form=$('sc-pro-form'),result=$('sc-pro-result');
if(!form||!result)return;
const submit=$('sc-pro-solve');
const toggle=document.createElement('label');
toggle.className='sc-pro-step-toggle didactic-toggle';
const input=document.createElement('input');
input.type='checkbox';input.id='sc-pro-show-steps';input.checked=true;
toggle.append(input,document.createTextNode(' Mostrar passo a passo do cálculo'));
submit?.parentNode?.insertAdjacentElement('beforebegin',toggle);

function step(parent,title,description,latex){
 const card=document.createElement('div');card.className='step';
 const heading=document.createElement('div');heading.className='step-title';heading.textContent=title;
 const note=document.createElement('div');note.className='step-note';note.textContent=description;
 card.append(heading);
 if(latex){
  const math=document.createElement('div');math.className='sc-pro-math';
  if(root.katex&&typeof root.katex.render==='function'){
   try{root.katex.render(latex,math,{displayMode:true,throwOnError:true,trust:false});}
   catch(e){math.textContent=latex;}
  }else math.textContent=latex;
  card.appendChild(math);
 }
 card.appendChild(note);
 parent.appendChild(card);
}
function currentResultIsValid(){
 return Boolean(result.querySelector('h4')&&$('sc-pro-status')?.textContent?.startsWith('Concluído'));
}
function refresh(){
 result.querySelectorAll('.sc-pro-restored-steps').forEach(el=>el.remove());
 if(!input.checked||!currentResultIsValid())return;
 const action=$('sc-pro-action').value;
 const expression=String($('sc-pro-expr').value||'').trim();
 if(!expression)return;
 const stages=document.createElement('details');
 stages.className='didactic sc-pro-restored-steps';stages.open=true;
 const summary=document.createElement('summary');summary.textContent='📚 Resolução passo a passo — método e conferência';
 const body=document.createElement('div');body.className='didactic-body';
 const intro=document.createElement('p');intro.className='didactic-intro';
 intro.textContent='Estas são as etapas justificadas pelo tipo de operação. Quando o motor não devolve uma derivação detalhada, isso é indicado explicitamente.';
 const list=document.createElement('div');list.className='step-list';
 body.append(intro,list);
 const normalized=expression.replace(/\s+/g,'').toLowerCase();
 const firstFormula=result.querySelector('.sc-pro-math');
 // A fórmula final já foi renderizada acima; não supomos que o CAS tenha
 // devolvido os passos internos ou que uma aproximação seja certificada.
 if(action==='integral'){
  const lower=$('sc-pro-lower')?.value.trim(),upper=$('sc-pro-upper')?.value.trim();
  const definite=Boolean(lower&&upper);
  step(list,'1. Identificar o integrando','Expressão fornecida: '+expression+(definite?'. Limites: '+lower+' até '+upper+'.':'. A integral é indefinida.'));
  if(normalized==='x*exp(x)'||normalized==='exp(x)*x'){
   step(list,'2. Integração por partes','Escolhemos u=x e dv=eˣ dx; assim du=dx e v=eˣ.',
    '\\int x e^x\\,dx=x e^x-\\int e^x\\,dx');
  }else if(normalized==='exp(-x^2)'&&!definite){
   step(list,'2. Identificar função especial','Não há primitiva elementar; a função erro (erf) descreve a integral.',
    '\\int e^{-x^2}\\,dx=\\frac{\\sqrt{\\pi}}{2}\\operatorname{erf}(x)+C');
  }else if(/^x\^\d+$/.test(normalized)){
   const n=Number(normalized.slice(2));
   if(n<=50)step(list,'2. Regra da potência','Aplicamos a regra ∫xⁿ dx = xⁿ⁺¹/(n+1)+C, para n ≠ −1.',
     '\\int x^{'+n+'}\\,dx=\\frac{x^{'+(n+1)+'}}{'+(n+1)+'}+C');
   else step(list,'2. Método simbólico','O motor tentou obter uma primitiva aplicando identidades e regras de integração.');
  }else{
   step(list,'2. Escolha de método','O motor tenta simplificações, integrações elementares e funções especiais. Não foi disponibilizada uma sequência interna detalhada para este caso.');
  }
  step(list,'3. Determinar o resultado',
   definite?'O SymPy avalia a expressão com os limites indicados, analisando convergência quando aplicável.':'O SymPy tenta encontrar uma primitiva. O resultado simbólico aparece acima.');
  if(result.textContent.includes('Conferência simbólica')){
   step(list,'4. Verificar a primitiva','A derivada da primitiva recuperou simbolicamente o integrando.');
  }else{
   step(list,'4. Conferência ou ressalva','Não foi apresentada verificação simbólica completa; confira domínio e singularidades quando houver.');
  }
 }else if(action==='derivative'){
  const order=Number($('sc-pro-order')?.value||1);
  step(list,'1. Identificar função e ordem','Função '+expression+'. Derivar '+order+(order===1?' vez.':' vezes.'));
  if(normalized==='x^2*exp(x)'&&order===1){
   step(list,'2. Regra do produto','Diferenciamos x² e eˣ e somamos os dois produtos.',
    '\\frac{d}{dx}(x^2e^x)=2xe^x+x^2e^x');
  }else{
   step(list,'2. Aplicar a regra apropriada','O SymPy usa regras de potência, produto, quociente e cadeia. Para derivadas de ordem superior repete a diferenciação.');
  }
  step(list,'3. Apresentar a derivada','A expressão simbólica calculada aparece no resultado acima.');
  if($('sc-pro-point')?.value.trim())step(list,'4. Avaliar no ponto','Substituímos o ponto solicitado na derivada obtida; veja o valor acima.');
 }else if(action==='limit'){
  const target=$('sc-pro-target')?.value.trim()||'0';
  const d=$('sc-pro-direction')?.value||'both';
  step(list,'1. Identificar a tendência','Estudamos '+expression+' quando x tende a '+target+'. Direção '+(d==='left'?'pela esquerda':d==='right'?'pela direita':'bilateral')+'.');
  step(list,'2. Analisar o comportamento','O SymPy usa simplificações e regras de limite. Se necessário, compara o comportamento por ambos os lados.');
  if(result.textContent.includes('Limites laterais'))step(list,'3. Comparar os lados','Os limites laterais estão apresentados acima. Para o limite bilateral existir, eles devem coincidir.');
  step(list,'4. Concluir','A resposta simbólica está acima; ela não deve ser substituída por valores de amostragem numérica.');
 }
 if(result.querySelector('strong')&&result.textContent.includes('Aproximação decimal')){
  step(list,'Observação sobre precisão','O valor decimal exibido é uma aproximação. Não equivale a um certificado rigoroso de erro.');
 }
 stages.append(summary,body);result.appendChild(stages);
}
let scheduled=false;
function schedule(){
 if(scheduled)return;scheduled=true;
 queueMicrotask(()=>{scheduled=false;refresh();});
}
const observer=new MutationObserver(mutations=>{
 // Recalcula somente quando o motor atualiza o conteúdo principal.
 if(mutations.every(m=>{
  const nodes=[...m.addedNodes,...m.removedNodes];
  return nodes.length>0&&nodes.every(n=>n.classList?.contains('sc-pro-restored-steps'));
 }))return;
 schedule();
});
observer.observe(result,{childList:true});
input.addEventListener('change',schedule);
root.scicalcRestoreProSteps=refresh;
// Preserve the previously available details in the classic educational forms.
const observed=[
 ['der-result','der-steps'],['int-result','int-steps'],
 ['lim-result','lim-steps'],['par-result','par-steps'],
 ['edo-result','edo-steps'],['four-result','four-steps']
];
observed.forEach(([id,checkId])=>{
 const target=$(id),checked=$(checkId);if(!target)return;
 const ensure=()=>{
  if(checked&&!checked.checked)return;
  target.querySelectorAll('details.didactic,details.sc-lab-details').forEach(details=>details.open=true);
 };
 new MutationObserver(ensure).observe(target,{childList:true});
});
})(window);
