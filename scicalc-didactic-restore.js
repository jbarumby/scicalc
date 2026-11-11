/* SciCalc — controles do passo a passo.
   Cálculo Pro mostra etapas verificadas enviadas pelo motor SymPy.
   Os calculadores clássicos mantêm seus detalhes didáticos acessíveis. */
(function(root){
'use strict';
const $=id=>document.getElementById(id);
const form=$('sc-pro-form'),submit=$('sc-pro-solve');
if(form&&submit&&!$('sc-pro-show-steps')){
 const toggle=document.createElement('label');
 toggle.className='sc-pro-step-toggle didactic-toggle';
 const input=document.createElement('input');
 input.type='checkbox';input.id='sc-pro-show-steps';input.checked=true;
 toggle.append(input,document.createTextNode(' Mostrar resolução matemática passo a passo'));
 submit.parentNode.insertAdjacentElement('beforebegin',toggle);
 input.addEventListener('change',()=>{
   if(typeof root.scicalcProToggleSteps==='function')root.scicalcProToggleSteps();
 });
}
const observed=[
 ['der-result','der-steps'],['int-result','int-steps'],
 ['lim-result','lim-steps'],['par-result','par-steps'],
 ['edo-result','edo-steps'],['four-result','four-steps']
];
observed.forEach(([id,checkId])=>{
 const target=$(id),checked=$(checkId);if(!target)return;
 const ensure=()=>{
  if(checked&&!checked.checked)return;
  target.querySelectorAll('details.didactic,details.sc-lab-details')
    .forEach(details=>{if(!details.open)details.open=true;});
 };
 new MutationObserver(ensure).observe(target,{childList:true});
});
})(window);
