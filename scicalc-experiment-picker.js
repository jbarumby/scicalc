/* SciCalc — catálogo de simuladores por clique.
 * Não remove modelos antigos e não troca as funções físicas existentes.
 * Monta <details> acessíveis a partir dos cards já presentes na página.
 */
(function(root){
'use strict';
const section=document.getElementById('simuladores-v26');
if(!section||section.dataset.experimentPickerReady)return;
const $=id=>document.getElementById(id);
const details=[];
function fold(element,title,description,key){
 if(!element||!element.parentNode)return null;
 const wrapper=document.createElement('details');
 wrapper.className='sc-experiment';
 wrapper.dataset.experiment=key;
 const summary=document.createElement('summary');
 summary.className='sc-experiment-summary';
 const left=document.createElement('span');
 left.className='sc-experiment-summary-text';
 const head=document.createElement('strong');head.textContent=title;
 const desc=document.createElement('small');desc.textContent=description;
 left.append(head,desc);
 const marker=document.createElement('span');
 marker.className='sc-experiment-marker';marker.setAttribute('aria-hidden','true');
 marker.textContent='Abrir experimento';
 summary.append(left,marker);
 const body=document.createElement('div');body.className='sc-experiment-body';
 element.parentNode.insertBefore(wrapper,element);
 wrapper.append(summary,body);
 body.appendChild(element);
 details.push(wrapper);
 return wrapper;
}
const blocks=Array.from(section.children).filter(el=>el.classList.contains('panel'));
const firstPanel=blocks[0],university=blocks[1],emPanel=blocks[2];
const firstGrid=firstPanel?.querySelector('.v26-extra-grid');
const mhs=firstGrid?.querySelector('.v26-extra-card');
const projectile=mhs?.nextElementSibling;
const magneticCards=emPanel?.querySelectorAll('.v26-extra-card')||[];
const animatedStudio=section.querySelector('.sc-pro-sim-panel');

const sim=fold(animatedStudio,
 'Massa–mola animada e lançamento',
 'Animação interativa · reproduzir, pausar e controlar tempo','animation');
const mhsDetails=fold(mhs,'Oscilador harmônico simples (MHS)',
 'Equações de posição, velocidade e gráfico','mhs');
const projectileDetails=fold(projectile,'Lançamento oblíquo',
 'Trajetória, alcance e altura máxima','projectile');
const lab=fold(university,'Laboratório de modelos físicos',
 'Plano inclinado, massa–mola amortecida, circuito RC e onda','lab');
const faraday=fold(magneticCards[0],'Indução eletromagnética de Faraday',
 'Ímã, solenoide, fluxo e força eletromotriz','faraday');
const lenz=fold(magneticCards[1],'Lei de Lenz',
 'Aproximação, afastamento e inversão da corrente','lenz');
if(!lab||!faraday||!lenz)return;
section.dataset.experimentPickerReady='1';
const picker=document.createElement('div');
picker.className='sc-experiment-picker';
picker.setAttribute('role','group');
picker.setAttribute('aria-label','Escolha um experimento para visualizar');
const heading=document.createElement('h3');
heading.textContent='Escolha um experimento para visualizar';
const guide=document.createElement('p');
guide.textContent='As simulações começam recolhidas. Clique no nome para abrir; apenas o experimento escolhido fica visível.';
picker.append(heading,guide);
const choices=document.createElement('div');
choices.className='sc-experiment-choices';
picker.appendChild(choices);
const models=[
 {key:'animation',name:'🎬 Massa–mola animada',subtitle:'Movimento ao vivo',box:sim},
 {key:'mhs',name:'〰️ Oscilador MHS',subtitle:'Equações e gráficos',box:mhsDetails},
 {key:'spring',name:'🌀 Massa–mola amortecida',subtitle:'Massa, mola e atrito',box:lab,model:'spring'},
 {key:'inclined',name:'📐 Plano inclinado',subtitle:'Atrito e aceleração',box:lab,model:'inclined'},
 {key:'rc',name:'🔌 Circuito RC',subtitle:'Carga de capacitor',box:lab,model:'rc'},
 {key:'wave',name:'🌊 Ondas',subtitle:'Propagação harmônica',box:lab,model:'wave'},
 {key:'projectile',name:'🏹 Lançamento',subtitle:'Trajetória balística',box:projectileDetails},
 {key:'faraday',name:'🧲 Faraday',subtitle:'Indução e fluxo',box:faraday},
 {key:'lenz',name:'↔️ Lei de Lenz',subtitle:'Corrente induzida',box:lenz}
].filter(m=>m.box);
const buttons=new Map();
function selectActive(key){
 buttons.forEach((b,id)=>{const isSelected=id===key;b.classList.toggle('active',isSelected);b.setAttribute('aria-pressed',String(isSelected));});
}
function refreshPlots(){
 if(!root.Plotly||!root.Plotly.Plots)return;
 const ids=['v26-mhs-plot','v26-obl-plot','v26-lab-plot'];
 ids.forEach(id=>{
   const el=$(id);
   if(el&&el.offsetParent!==null&&el.data){
     try{root.Plotly.Plots.resize(el);}catch(e){}
   }
 });
}
function openBox(box,model,key){
 if(!box)return;
 // Mark open directly on click so the model initializes on a visible canvas.
 details.forEach(d=>{if(d!==box)d.open=false;});
 box.open=true;
 selectActive(key);
 if(box===lab){
   const select=$('v26-lab-model');
   if(select&&model)select.value=model;
   if(typeof root.v26InitLab==='function'){
     try{root.v26InitLab();}catch(e){console.warn('Lab init:',e);}
   }
 } else if(box===faraday&&typeof root.flInit==='function'){
   try{root.flInit();}catch(e){}
 }
 root.requestAnimationFrame(()=>{
   root.dispatchEvent(new Event('resize'));refreshPlots();
   // Programmatically opening a hidden experiment gives a focused landmark.
   box.scrollIntoView({behavior:'smooth',block:'nearest'});
 });
}
choices.addEventListener('click',event=>{
 const button=event.target.closest('button[data-choice]');
 if(!button)return;
 const item=models.find(m=>m.key===button.dataset.choice);
 if(item)openBox(item.box,item.model,item.key);
});
models.forEach(item=>{
 const b=document.createElement('button');
 b.type='button';b.className='sc-experiment-choice';
 b.dataset.choice=item.key;
 b.setAttribute('aria-pressed','false');
 const label=document.createElement('strong');label.textContent=item.name;
 const sub=document.createElement('small');sub.textContent=item.subtitle;
 b.append(label,sub);choices.appendChild(b);buttons.set(item.key,b);
});
firstPanel?.parentNode?.insertBefore(picker,firstPanel);
const oldInit=root.v26InitLab;
root.v26InitLab=function(){
 if(!lab.open)return;
 if(typeof oldInit==='function')return oldInit.apply(this,arguments);
};
function stopAnimationIfHidden(which){
 if(which===faraday&&typeof root.flStop==='function')root.flStop();
 if(which===lenz&&typeof root.flStopLenz==='function')root.flStopLenz();
 if(which===sim){
   const p=$('sc-pro-play');
   if(p&&p.textContent.includes('Pausar'))p.click();
 }
}
details.forEach(d=>d.addEventListener('toggle',()=>{
 if(d.open){
   details.forEach(other=>{if(other!==d)other.open=false;});
   const key=d===lab?$('v26-lab-model')?.value:d.dataset.experiment;
   selectActive(key);
   if(d===lab&&typeof root.v26InitLab==='function')root.v26InitLab();
   if(d===faraday&&typeof root.flInit==='function')root.flInit();
   root.requestAnimationFrame(()=>{root.dispatchEvent(new Event('resize'));refreshPlots();});
 }else{
   stopAnimationIfHidden(d);
   if(!details.some(other=>other.open))selectActive(null);
 }
}));
const select=$('v26-lab-model');
if(select)select.addEventListener('change',()=>{if(lab.open)selectActive(select.value);});
// Enforce collapsed by default and don't autoplay any model.
details.forEach(d=>d.open=false);
})(window);
