"""Certified sufficient prefix via signed endpoint moment Gram, conditional on V38 symbol bounds."""
import prefix39 as p
from fractions import Fraction as F
from math import factorial
import json, argparse
from pathlib import Path
import numpy as np
from mpmath import iv

def moments(T,ks,n):
 data=np.loadtxt('deficit_2000000.csv',delimiter=',',skiprows=1,dtype=np.int64);data=data[data[:,0]>=16*T]
 scale=2**512; M={};group=-1;prev=-1;cells=[]
 for l,r,d in data.tolist():
  if l!=prev:group+=1
  prev=r;cells.append((l,r,d*ks[group]))
 if group+1!=len(ks):raise ArithmeticError('components mismatch')
 for j in range(2,2*n+1,2):
  e=j-1;lo=hi=0;factor=16**e*scale
  for l,r,d in cells:
   lp=l**e;rp=r**e;num=d*factor*(rp-lp);den=10**6*e*lp*rp;q,rem=divmod(num,den);lo+=q;hi+=q+bool(rem)
  M[j]=iv.mpf([p.I(F(lo,scale)).a,p.I(F(hi,scale)).b])
 return M

def correction(par,ks):
 ns=list(range(par,2*p.N,2));n=max(ns);T=4500;M=moments(T,ks,n+1);logT=iv.log(T)
 # 4/pi <=4/3, with positive weight; retain all signed cross terms in same derivative parity.
 J={j:p.I(F(4,3))*((logT+8)/((j-1)*T**(j-1))+p.I(F(1,(j-1)**2*T**(j-1)))+M[j]) for j in M}
 ds=[[(-1)**(j//2)*F(factorial(k+j),2**j*factorial(j)*factorial(k-j))/p.a**j if k>=j else F(0) for j in range(n+1)] for k in ns]
 mat=[[iv.mpf(0) for _ in ns] for _ in ns]
 for i,ni in enumerate(ns):
  for k,nk in enumerate(ns[:i+1]):
   # group equal powers before evaluating intervals to preserve cancellation
   vals={}
   for j in range(ni+1):
    for l in range(j%2,nk+1,2):vals[j+l+2]=vals.get(j+l+2,F(0))+ds[i][j]*ds[k][l]
   mat[i][k]=mat[k][i]=sum((p.I(v)*J[e] for e,v in vals.items()),iv.mpf(0))
 return mat

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--N',type=int,default=8);ap.add_argument('--width',default='6931/5000');ap.add_argument('--adaptive',action='store_true');args=ap.parse_args();p.N=args.N;p.a=F(args.width);p.AA=p.I(p.a);p.TAIL=p.tail();
 import opt38
 z=opt38.components(np.loadtxt('deficit_2000000.csv',delimiter=',',skiprows=1,dtype=np.int64),4500)
 if args.adaptive:ks=json.load(open('hybrid39.json'))['ks'];tag='adaptive_binet'
 else:ks=[24]*len(z);tag='uniform24_binet'
 results=[]
 for par in [0,1]:
  A=p.matrix(par);D=correction(par,ks);R=[[A[i][j]-D[i][j] for j in range(p.N)] for i in range(p.N)];N=p.N
  for m in range(1,N+1):
   for h in [F(1,10**60),F(1,10**10)]:
    B=[r[:] for r in R]
    for i in range(N):B[i][i]-=p.I(h*2*p.a/(2*(2*i+par)+1))
    p.N=m;pv,ok=p.certify(B,h);p.N=N
    results.append({'p':par,'N':m,'h':str(h),'positive':ok,'pivots':len(pv),'last_pivot':pv[-1] if pv else None})
  print(par,[(r['N'],r['positive'],r['pivots']) for r in results if r['p']==par and r['h']==str(F(1,10**60))],flush=True)
  # record interval entries for independent replay
  json.dump({'A':[[str(x) for x in r] for r in A],'D':[[str(x) for x in r] for r in D]},open(f'matrices_{tag}_{args.N}_p{par}.json','w'),indent=2)
 json.dump({'width':str(p.a),'N':args.N,'tag':tag,'results':results,'scope':'conditional on V38 deficit bounds, and analytic explicit formula; prefixes only'},open(f'moments_{tag}_{args.N}.json','w'),indent=2)
if __name__=='__main__':main()
