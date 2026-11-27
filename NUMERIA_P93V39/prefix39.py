"""Interval certificate for full form and stronger finite sufficient lower matrices."""
from fractions import Fraction as F
from math import factorial
import mpmath as mp
from mpmath import iv
import sympy as sp
import json
mp.mp.dps=120;iv.dps=110
a=F(6931,5000); u,z=sp.symbols('u z'); N=8; L=40

def I(x):
 x=F(x);return iv.mpf(x.numerator)/x.denominator

def poly_eval(coeff,x):
 out=iv.mpf(0)
 for c in reversed(coeff):out=out*x+I(c)
 return out

def coeff(poly,var):
 p=sp.Poly(poly,var);return [F(p.nth(k)) for k in range(p.degree()+1)] if not p.is_zero else [F(0)]

def positive(x):return x.a>0

def gauss():
 nodes,_=mp.gauss_quadrature(L,'legendre');p=coeff(sp.legendre(L,u),u);dp=[(k+1)*p[k+1] for k in range(L)]; result=[];prev=F(-1)
 for node in nodes:
  center=F(str(node));lo=center-F(1,10**85);hi=center+F(1,10**85)
  if not(prev<lo<hi<1):raise ValueError('node ordering')
  fl=poly_eval(p,I(lo));fh=poly_eval(p,I(hi));
  if not((fl.a>0 and fh.b<0) or (fl.b<0 and fh.a>0)):raise ValueError('node sign')
  X=iv.mpf([I(lo).a,I(hi).b]);W=2/((1-X*X)*poly_eval(dp,X)**2)
  if not positive(W):raise ValueError('weight')
  result.append((X,W));prev=hi
 return result

G=gauss(); AA=I(a)
# Composite GL24 on eight panels in [0,1], E4: sup |u|<6/5, |Im(2au)|<1, |Re(2au)|< 10/3.
# Uniform error <= 16 M 4^-48, loose bound for quadrature/series tail.
def quad(fun,M):
 out=iv.mpf(0)
 for j in range(8):
  for X,W in G:out+=W*fun((I(F(2*j+1))+X)/16)/16
 err=I(F(16)*M/F(4**(2*L)));return out+iv.mpf([-err.b,err.b])

from math import comb
from functools import lru_cache
@lru_cache(None)
def leg(n):return [(-1)**(n-j)*comb(n,j)*comb(n+j,j) for j in range(n+1)]
def corr(n,m):
 aa,bb=leg(n),leg(m); groups={}
 for i,ai in enumerate(aa):
  for k in range(i+1):
   for j,bj in enumerate(bb):
    key=(i-k,k+j+1);groups[key]=groups.get(key,0)+ai*bj*comb(i,k)
 out=[F(0)]*(n+m+2)
 for (v,w),c in groups.items():
  if c:
   f=F(c,w)*2*a
   for j in range(w+1):out[v+j]+=f*comb(w,j)*(-1)**j
 return out

def pow_norm(p):return sum(abs(x)*F(6,5)**i for i,x in enumerate(p))

def tail():
 out=sum(2*iv.exp(-I((2*k+F(1,2))*2*a))/I(2*k+F(1,2)) for k in range(64))
 rest=2*iv.exp(-I((128+F(1,2))*2*a))/I(128+F(1,2))/(1-iv.exp(-I(4*a)))
 return out+iv.mpf([0,rest.b])

TAIL=tail();const=-iv.euler-iv.pi/2-3*iv.log(2)-iv.log(iv.pi)
q=[2,3,4,5,7,8,9,11,13]; primes=[2,3,2,5,7,2,3,11,13]

def pole(n):
 p=coeff(sp.legendre(n,2*u-1),u)
 return 2*AA*quad(lambda x:iv.exp(AA*(2*x-1)/2)*poly_eval(p,x), 4*pow_norm(p))

def matrix(par):
 ns=list(range(par,2*N,2));ps=[pole(n) for n in ns];A=[[None]*N for _ in range(N)]
 for i,n in enumerate(ns):
  for j in range(i+1):
   m=ns[j];C=corr(n,m);C0=F(2)*a/(2*n+1) if n==m else F(0)
   if C[0]!=C0:raise ValueError('C(0)')
   D=[-c for c in C[1:]]
   def fun(x):
    v=2*AA*x;return poly_eval(D,x)*v*iv.exp(v/2)/((iv.exp(v)-iv.exp(-v))/2)
   B=I(C0)*(const+TAIL)+quad(fun,12*pow_norm(D))
   for qq,pp in zip(q,primes):B-=2*iv.log(pp)/iv.sqrt(qq)*poly_eval(C,iv.log(qq)/(2*AA))
   B+=2*(-1)**par*ps[i]*ps[j]
   A[i][j]=A[j][i]=B
 return A

def certify(A,h):
 # Normalize by mass using diagonal shift h*norm^2. h refers to orthonormal basis.
 B=[row[:] for row in A];piv=[]
 for i in range(N):
  v=B[i][i];
  if not positive(v):return piv,False
  piv.append(str(v.a));
  for j in range(i+1,N):
   for k in range(j,N):B[k][j]=B[j][k]=B[k][j]-B[k][i]*B[j][i]/v
 return piv,True

def correction(par,mode):
 # In each parity |F_f(t)| <= 2 sum |f^(k)(a)| t^(-k-1).
 # Cauchy-Schwarz with weights 16^-k gives <= (64/15) sum 16^k |f^(k)(a)|^2 t^(-2k-2).
 T=PARAM_T; d=PARAM_D; lam=PARAM_LAM; k=PARAM_K
 extra=d/(1-lam) if mode=='harmonic' else k*d
 C=25+extra; logT=iv.log(T); ns=list(range(par,2*N,2)); mat=[[iv.mpf(0) for _ in ns] for _ in ns]
 for j in range(max(ns)+1):
  ds=[F(factorial(n+j),2**j*factorial(j)*factorial(n-j))/a**j if n>=j else F(0) for n in ns]
  m=2*j+1
  coef=I(F(64*16**j,45*T**m)) *((2*logT+I(C))/m+I(F(2,m*m)))
  for ii in range(N):
   for jj in range(N):mat[ii][jj]+=coef*I(ds[ii]*ds[jj])
 return mat

PARAM_T=6400;PARAM_D=F(949191,250000);PARAM_LAM=F(43,50);PARAM_K=8

def main():
 global PARAM_T,PARAM_D,PARAM_LAM,PARAM_K,N,a,AA,TAIL
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--N',type=int,default=8);ap.add_argument('--width',default='6931/5000');ap.add_argument('--k',type=int,default=24);args=ap.parse_args()
 N=args.N;a=F(args.width);AA=I(a);TAIL=tail();PARAM_T=4500;PARAM_D=F(3870563,1000000);PARAM_LAM=F(19,20);PARAM_K=args.k
 results=[]
 for p in [0,1]:
  A=matrix(p)
  for mode in ['full','localized']:
   R=[r[:] for r in A] if mode=='full' else [[A[i][j]-d for j,d in enumerate(row)] for i,row in enumerate(correction(p,mode))]
   for m in range(1,N+1):
    for h in [F(1,10**60),F(1,10**10)]:
     B=[r[:] for r in R]
     for i in range(N):B[i][i]-=I(h*2*a/(2*(2*i+p)+1))
     old=N;N=m;pv,ok=certify(B,h);N=old
     results.append({'p':p,'N':m,'mode':mode,'h':str(h),'positive':ok,'pivots':len(pv),'last_pivot':pv[-1] if pv else None})
   print(p,mode,[(r['N'],r['positive'],r['pivots']) for r in results if r['p']==p and r['mode']==mode and r['h']==str(F(1,10**60))],flush=True)
 out={'width':str(a),'T':PARAM_T,'k':PARAM_K,'L':L,'results':results,'scope':'hybrid compression equals localized for N<=600; conditional on V38 symbol and analytic bounds; no full head'}
 json.dump(out,open(f'prefix39_N{N}_k{PARAM_K}.json','w'),indent=2)
if __name__=='__main__':main()
