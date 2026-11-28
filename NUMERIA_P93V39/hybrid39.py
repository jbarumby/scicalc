from fractions import Fraction as F
from math import lcm, factorial
import numpy as np, json
import opt38 as opt
import schur38 as sch

def need(x,s):
 if not x: raise ArithmeticError(s)
T=4500; B=F(2); A=F(1387,1000)
data=np.loadtxt('deficit_2000000.csv',delimiter=',',skiprows=1,dtype=np.int64)
z=opt.components(data,T); l,r,m,d=z.T
S=2**20; dist=np.maximum(0,np.maximum(l[:,None]-r[None,:],l[None,:]-r[:,None]));cap=(A.numerator*S+A.denominator-1)//A.denominator
ker=np.minimum(cap,(16*S+np.maximum(dist,1)-1)//np.maximum(dist,1));ker[dist==0]=cap;su=l[:,None]+l[None,:];ker+=(16*S+su-1)//su
K=ker.astype(float)*m[None,:]/(S*16e6)*106/333
# choose proposals in float only, then exactly replay the weighted Schur row bounds.
def propose(ks):
 U=ks/(ks-1);mat=K*U[None,:];w=np.ones(len(z))
 for _ in range(32):
  w=mat@w;w/=w.max();w=np.maximum(w,1/4096)
 wi=np.maximum(1,np.rint(w*4096)).astype(np.int64)
 return float(np.max((mat@wi)/wi)),wi

def exact(ks,wi):
 den=lcm(*[int(k)-1 for k in np.unique(ks)]);v=[int(mi)*int(w)*int(k)*(den//(int(k)-1)) for mi,w,k in zip(m,wi,ks)]
 need(int(ker.max())*sum(v)<2**63,'overflow');rows=ker@np.array(v,dtype=np.int64)
 c=max(F(int(x)*106,den*S*16000000*333*int(w)) for x,w in zip(rows,wi));dd=max(F(int(di)*int(k),10**6*(int(k)-1)) for di,k in zip(d,ks));return c,dd

# objective for initial block: weighted moment integral of k_i*density/t^2.
scoreweight=m.astype(float)/(l/16.)**2
candidates=[]
for top in [14,16,20,24,32,40,48]:
 for low in [2,3,4,6,8,10,12,14]:
  if low>top:continue
  for threshold in [0.25,0.5,0.75,1.,1.5]:
   ks=np.where(m/16e6>=threshold,top,low).astype(np.int64);ce,wi=propose(ks)
   if ce>=1.97:continue
   candidates.append((float(scoreweight@ks),ce,top,low,threshold,ks,wi))
# compare a uniform and adaptive candidate; pick a rational feasible lambda with reserve >= .01
out=[]
for _,_,top,low,threshold,ks,wi in sorted(candidates,key=lambda x:x[0]):
 try:c,dd=exact(ks,wi)
 except ArithmeticError:continue
 choices=[]
 for i in range(80,100):
  lam=F(i,100);s=lam*dd/((1-lam)*c);f1=c/lam*s/(1+s)
  if B-f1>=F(1,100):choices.append((lam,s,f1))
 if not choices:continue
 lam,s,f1=choices[0]
 out.append({'top':top,'low':low,'threshold':threshold,'c':str(c),'d':str(dd),'lambda':str(lam),'s':str(s),'f1':str(f1),'reserve':float(B-f1),'localized_reserve':float(B-c),'ks':ks.tolist(),'weights':wi.tolist(),'score':float(scoreweight@ks),'uniform24_score':float(scoreweight@np.full(len(z),24))});break
need(bool(out),'no candidate')
# exact Schur gate for the hybrid: off-diagonal penalty <= d'/(1-lambda)*sqrt(g); tail deviation from A0 <=50g+eps.
row=out[0];c=F(row['c']);dd=F(row['d']);lam=F(row['lambda']);s=F(row['s']);f1=F(row['f1']);eps=sch.EPS;h=sch.H;target=sch.M
# optimize rho and N in rational replay; float proposals unnecessary. candidate rho 5/4, 6/5, 7/6, 9/8.
for rho in [F(5,4),F(6,5),F(7,6),F(9,8),F(11,10)]:
 sch.T=T;sch.RHO=rho;sch.E=sch.expup(A*T*(rho-1/rho));coef=50+dd/(1-lam)
 for N in range(3150,3501,5):
  gg=sch.g(N,0);dev=50*gg+eps;b2=2*coef**2*gg+2*eps**2;res=B-f1-dev-target
  if res>0 and (h-target)*res>b2:
   need(F(48,factorial(2*N))<eps and F(72,factorial(2*N)**2)<eps,'poles');row.setdefault('schur',[]).append({'rho':str(rho),'N':N,'g_even':float(gg),'slack_even':float((h-target)*res-b2),'both_parities':True});break
json.dump(row,open('hybrid39.json','w'),indent=2)
print({k:v for k,v in row.items() if k not in ['ks','weights']})
