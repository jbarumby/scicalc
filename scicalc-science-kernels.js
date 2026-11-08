/* SciCalc Labs v1: mathematical kernels, independent of the UI.
   No arbitrary JavaScript evaluation. Errors are reported, never silently replaced by zero. */
(function(root){
'use strict';
const abs=n=>n<0n?-n:n;
function bi(value,label='inteiro',max=1000000000000n){
 const s=String(value).trim();
 if(!/^[+-]?\d{1,16}$/.test(s))throw Error(label+': informe um inteiro decimal.');
 const n=BigInt(s);
 if(abs(n)>max)throw Error(label+': módulo máximo de '+max.toString()+' neste laboratório.');
 return n;
}
function gcd(a,b){a=abs(a);b=abs(b);while(b){const r=a%b;a=b;b=r;}return a;}
function egcd(a,b){
 let oldR=a,r=b,oldS=1n,s=0n,oldT=0n,t=1n,steps=[];
 while(r!==0n){
  const q=oldR/r;
  if(steps.length<60)steps.push({a:oldR.toString(),b:r.toString(),q:q.toString(),remainder:(oldR-q*r).toString()});
  [oldR,r]=[r,oldR-q*r];[oldS,s]=[s,oldS-q*s];[oldT,t]=[t,oldT-q*t];
 }
 if(oldR<0n){oldR=-oldR;oldS=-oldS;oldT=-oldT;}
 return {g:oldR,x:oldS,y:oldT,steps};
}
function mod(a,n){return (a%n+n)%n;}
function powmod(a,b,m){
 if(m<=0n||b<0n)throw Error('Módulo positivo e expoente não negativo são necessários.');
 let result=1n%m,base=mod(a,m);
 while(b){if(b&1n)result=result*base%m;base=base*base%m;b>>=1n;}
 return result;
}
function factor(input){
 let n=bi(input,'n',1000000000n);
 if(n<1n)throw Error('Fatoração: informe n ≥ 1.');
 const out=[];let k=2n;
 while(k*k<=n){
  if(n%k===0n){let power=0;do{power++;n/=k;}while(n%k===0n);out.push({prime:k.toString(),exp:power});}
  k=k===2n?3n:k+2n;
 }
 if(n>1n)out.push({prime:n.toString(),exp:1});
 return out;
}
function isPrime(n){n=bi(n,'n',1000000000n);if(n<2n)return false;if(n%2n===0n)return n===2n;
 for(let d=3n;d*d<=n;d+=2n)if(n%d===0n)return false;
 return true;
}
function sieve(n){
 if(!Number.isInteger(n)||n<2||n>1000000)throw Error('Crivo: use 2 ≤ N ≤ 1.000.000.');
 const marked=new Uint8Array(n+1), primes=[];
 for(let p=2;p<=n;p++)if(!marked[p]){primes.push(p);if(p*p<=n)for(let k=p*p;k<=n;k+=p)marked[k]=1;}
 return primes;
}
function arithmetic(value){
 const n=bi(value,'n',1000000000n);if(n<1n)throw Error('Use n ≥ 1.');
 const factors=factor(n);
 let phi=n,tau=1n,sigma=1n,mu=1n;
 for(const item of factors){
  const p=BigInt(item.prime),e=BigInt(item.exp);
  phi=phi/p*(p-1n);tau*=e+1n;sigma*=(p**(e+1n)-1n)/(p-1n);
  if(e>1n)mu=0n;else mu=-mu;
 }
 return {n:n.toString(),factors,phi:phi.toString(),tau:tau.toString(),sigma:sigma.toString(),mu:mu.toString()};
}
function crt(a,m,b,n){
 a=bi(a,'a');b=bi(b,'b');m=bi(m,'m',1000000n);n=bi(n,'n',1000000n);
 if(m<=0n||n<=0n)throw Error('Os módulos devem ser positivos.');
 const e=egcd(m,n),delta=b-a;
 if(delta%e.g!==0n)return {possible:false,g:e.g.toString()};
 const multiplier=mod((delta/e.g)*e.x,n/e.g);
 const period=(m/e.g)*n;
 return {possible:true,value:mod(a+m*multiplier,period).toString(),period:period.toString(),g:e.g.toString()};
}
const finite=(v,label)=>{if(typeof v!=='number'||!Number.isFinite(v))throw Error(label+' fora do domínio real ou divergente.');return v;};
function rk4(f,x,y,h){
 const k1=finite(f(x,y),'k1');
 const k2=finite(f(x+h/2,y+h*k1/2),'k2');
 const k3=finite(f(x+h/2,y+h*k2/2),'k3');
 const k4=finite(f(x+h,y+h*k3),'k4');
 const result=finite(y+h*(k1+2*k2+2*k3+k4)/6,'Próximo y');
 return {value:result,k1,k2,k3,k4};
}
function ode(f,x0,y0,xf,n){
 for(const [v,label] of [[x0,'x inicial'],[y0,'y inicial'],[xf,'x final']])finite(v,label);
 if(typeof f!=='function'||!Number.isInteger(n)||n<2||n>4000||x0===xf)throw Error('EDO: 2 ≤ N ≤ 4.000 e x final diferente de x inicial.');
 const h=(xf-x0)/n, xs=[x0],ys=[y0],rows=[];let y=y0;
 for(let k=0;k<n;k++){
  const x=x0+k*h,step=rk4(f,x,y,h);
  if(k<8||k===n-1)rows.push({k:k+1,x,y,...step});
  y=step.value;xs.push(x0+(k+1)*h);ys.push(y);
 }
 return {xs,ys,rows,h,end:y};
}
function odeRefined(f,x0,y0,xf,n){
 if(n>2000)throw Error('EDO: escolha N ≤ 2.000 para estimar erro com 2N passos.');
 const coarse=ode(f,x0,y0,xf,n),fine=ode(f,x0,y0,xf,n*2);
 // RK4: diferença entre N e 2N, dividida por (2^4-1) para estimar erro da solução refinada.
 const estimatedError=Math.abs(coarse.end-fine.end)/15;
 const diffs=coarse.ys.map((v,k)=>Math.abs(v-fine.ys[k*2])/15);
 return {coarse,fine,estimatedError,diffs};
}
function fft(re,im){
 const n=re.length;
 if(n<4||n>16384||(n&(n-1))!==0||im.length!==n)throw Error('FFT: N deve ser potência de 2 entre 4 e 16.384.');
 for(let i=1,j=0;i<n;i++){
  let bit=n>>1;
  for(;j&bit;bit>>=1)j^=bit;
  j^=bit;
  if(i<j){[re[i],re[j]]=[re[j],re[i]];[im[i],im[j]]=[im[j],im[i]];}
 }
 for(let len=2;len<=n;len*=2){
  const angle=-2*Math.PI/len,wr0=Math.cos(angle),wi0=Math.sin(angle);
  for(let start=0;start<n;start+=len){
   let wr=1,wi=0;
   for(let j=0;j<len/2;j++){
    const i0=start+j,i1=i0+len/2;
    const vr=re[i1]*wr-im[i1]*wi,vi=re[i1]*wi+im[i1]*wr;
    re[i1]=re[i0]-vr;im[i1]=im[i0]-vi;
    re[i0]+=vr;im[i0]+=vi;
    [wr,wi]=[wr*wr0-wi*wi0,wr*wi0+wi*wr0];
   }
  }
 }
 return {re,im};
}
function fourier(f,t0,t1,n,windowName='hann'){
 if(typeof f!=='function')throw Error('Forneça a função do sinal.');
 if(!Number.isInteger(n)||n<16||n>4096||(n&(n-1))!==0)throw Error('Fourier: N deve ser potência de 2 entre 16 e 4096.');
 finite(t0,'Tempo inicial');finite(t1,'Tempo final');
 if(!(t1>t0))throw Error('O tempo final deve ser maior que o inicial.');
 if(!['hann','rectangular'].includes(windowName))throw Error('Janela inválida.');
 const dt=(t1-t0)/n,fs=1/dt,df=1/(t1-t0),xs=[],ys=[],re=[],im=new Array(n).fill(0);
 let weightSum=0;
 for(let i=0;i<n;i++){
  const t=t0+i*dt,value=finite(f(t),'Sinal na amostra '+i);
  const w=windowName==='hann'?.5-.5*Math.cos(2*Math.PI*i/n):1;
  xs.push(t);ys.push(value);re.push(value*w);weightSum+=w;
 }
 fft(re,im);
 const freqs=[],amps=[],phases=[];
 for(let k=0;k<=n/2;k++){
  const amplitude=Math.hypot(re[k],im[k])/(weightSum)*(k===0||k===n/2?1:2);
  freqs.push(k*df);amps.push(amplitude);phases.push(Math.atan2(im[k],re[k]));
 }
 const peaks=[];
 for(let i=1;i<amps.length-1;i++){
  if(amps[i]>=amps[i-1]&&amps[i]>amps[i+1])peaks.push({index:i,freq:freqs[i],omega:2*Math.PI*freqs[i],amp:amps[i]});
 }
 peaks.sort((a,b)=>b.amp-a.amp);
 return {xs,ys,freqs,amps,phases,peaks:peaks.slice(0,6),dt,fs,df,nyquist:fs/2,n,windowName,mean:ys.reduce((a,b)=>a+b,0)/n};
}
root.SciCalcLabs={bi,gcd,egcd,mod,powmod,factor,isPrime,sieve,arithmetic,crt,ode,odeRefined,fft,fourier};
})(typeof window!=='undefined'?window:globalThis);
