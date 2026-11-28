import numpy as np,json
from fractions import Fraction as F
from pathlib import Path
S=2**20
A=F(1387,1000)

def components(data,T):
 z=data[data[:,0]>=16*T]; parts=[]
 for l,r,d in z.tolist():
  if parts and parts[-1][1]==l:parts[-1][1]=r;parts[-1][2]+=d;parts[-1][3]=max(parts[-1][3],d)
  else:parts.append([l,r,d,d])
 return np.array(parts,dtype=np.int64)

def cert(z):
 if not len(z):return F(0),F(0),[]
 l,r,m,d=z.T; dist=np.maximum(0,np.maximum(l[:,None]-r[None,:],l[None,:]-r[:,None]));cap=(A.numerator*S+A.denominator-1)//A.denominator
 kernel=np.minimum(cap,(16*S+np.maximum(dist,1)-1)//np.maximum(dist,1));kernel[dist==0]=cap
 su=l[:,None]+l[None,:];kernel+=(16*S+su-1)//su
 # kernel upper in 2^-20; propose weights numerically then exact rational replay
 K=kernel.astype(float)*m[None,:]/(S*16e6)*106/333
 weights=np.ones(len(z));
 for _ in range(40):
  weights=K@weights;weights/=weights.max();weights=np.maximum(weights,1/4096)
 wi=np.maximum(1,np.rint(weights*4096)).astype(np.int64)
 maxprod=int(kernel.max())*sum(int(x)*int(y) for x,y in zip(m,wi))
 if maxprod>=2**63: raise ValueError('overflow')
 rows=kernel@(m*wi)
 c=max(F(int(v)*106,S*16000000*333*int(w)) for v,w in zip(rows,wi));return c,F(int(max(d)),1000000),wi.tolist()

def main():
 results=[]
 for b in [1500000,1625000,1750000,1875000,2000000]:
  data=np.loadtxt(f'deficit_{b}.csv',delimiter=',',skiprows=1,dtype=np.int64)
  for T in [1600,3100,4500,6400,10000]:
   z=components(data,T);c,d,w=cert(z);beta=F(b,1000000)
   # optimized harmonic at H=I: no improvement over min(c,d) at endpoint;
   # select lambda with positive reserve and best slope tradeoff
   choices=[]
   for i in range(1,100):
    lam=F(i,100);s=lam*d/((1-lam)*c);f1=c/lam*s/(1+s)
    if f1<beta:choices.append((f1,s,lam))
   # prefer at least 1/10 reserve, minimum slope d/(1-lambda)
   feasible=[x for x in choices if x[0]<=beta-F(1,10)]
   pick=min(feasible,key=lambda x:x[1]) if feasible else (min(choices,key=lambda x:x[0]) if choices else None)
   k=next((k for k in range(2,100) if c*k/(k-1)<beta),None)
   row={'beta':str(beta),'T':T,'components':len(z),'c':str(c),'c_display':float(c),'d':str(d),'d_display':float(d),'k_uniform':k,'eta_uniform':float(c*k/(k-1)) if k else None,'harmonic':{'lambda':str(pick[2]),'s':str(pick[1]),'f1':float(pick[0]),'reserve':float(beta-pick[0])} if pick else None,'weights':w}
   results.append(row);print({k:v for k,v in row.items() if k!='weights'},flush=True)
 with open('opt38.json','w') as f:json.dump(results,f,indent=2)
if __name__=='__main__':main()
