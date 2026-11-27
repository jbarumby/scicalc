from fractions import Fraction as F
from math import factorial
import json
from pathlib import Path
A=F(1387,1000);T=6400;BETA=F(2);H=F(1,10**60);M=F(9,10)*H;RHO=F(5,4);EPS=F(1,2**2000)
GRID=2**160
def need(ok,msg):
 if not ok:raise ArithmeticError(msg)
def up(x):
 v=x*GRID;return F(-(-v.numerator//v.denominator),GRID)
def expup(x):
 y=x/8192;need(0<=y<=1,'exp range');s=t=F(1)
 for j in range(1,41):t*=y/j;s+=t
 s+=t*y/41/(1-y/42);s=up(s)
 for _ in range(13):s=up(s*s)
 return s
E=expup(A*T*(RHO-1/RHO))
def g(N,p):
 return (8*A*T/3)*E*RHO**(-2*(2*N+p))/((1-1/RHO)**2*(1-RHO**(-4)))
def logiv(x):
 x=F(x);k=0
 while x>=2:x/=2;k+=1
 def ser(z):
  t=(z-1)/(z+1);s=sum((2*t**(2*j+1)/(2*j+1) for j in range(110)),F(0));return s,s+2*t**221/(221*(1-t*t))
 a,b=ser(x);c,d=ser(F(2));return a+k*c,b+k*d

def main():
 global T,RHO,E
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--T',type=int,default=6400);ap.add_argument('--lam',default=None);ap.add_argument('--k',type=int,default=None);args=ap.parse_args();T=args.T
 row=next(r for r in json.load(open('opt38.json')) if r['beta']=='2' and r['T']==T);c=F(row['c']);d=F(row['d']);lam=F(args.lam or row['harmonic']['lambda']);s=lam*d/((1-lam)*c);f1=c/lam*s/(1+s);k=args.k or row['k_uniform'];eta=k*c/(k-1)
 results=[]
 # choose rho by approximate optimizer for candidate only; all accepted gates rational
 RHO=F(5,4);E=expup(A*T*(RHO-1/RHO))
 for mode,pen,coef in [('harmonic',f1,50+d/(1-lam)),(f'localized_M600_k{k}',eta,F(50))]:
  for N in range(int(A*T/2)//25*25, int(A*T/2)+1001,25):
   rs=[]
   for p in [0,1]:
    gg=g(N,p);b2=2*coef**2*gg+2*EPS**2;dev=50*gg+EPS;res=BETA-pen-dev-M;slack=(H-M)*res-b2
    rs.append({'parity':p,'g':float(gg),'b2':float(b2),'reserve':float(res),'slack':float(slack),'positive':bool(res>0 and slack>0)})
   if all(r['positive'] for r in rs):
    need(F(48,factorial(2*N))<EPS and F(72,factorial(2*N)**2)<EPS,'pole bound');results.append({'mode':mode,'N':N,'rho':str(RHO),'sectors':rs});break
  else:raise ValueError('no N')
 # a16 strictly above log16/2; round upward on a 10^-75 grid.
 lo,hi=logiv(16);scale=10**75;a16=F((hi*scale/2).numerator//(hi*scale/2).denominator+1,scale)
 delta=2*a16-lo;need(2*a16>hi and a16<A and delta<F(3,10**75),'crossing width')
 omega=4000000;lofar,_=logiv(F(omega*7,44));need(lofar-F(1,omega)-F(9944,1000)-2>1,'high reserve')
 _,log2=logiv(2);need(log2/4<F(1,5),'kappa');loss=F(omega,6)*delta;need(loss<M-F(89999999,10**8)*H,'loss')
 out={'scope':'conditional on two complete certified heads at a16 and symbol/Binet inputs; no actual crossing', 'T':T,'k':k,'head_h':str(H),'target_m':str(M),'c':str(c),'d':str(d),'lambda':str(lam),'f1':str(f1),'eta':str(eta),'schur':results,'a16':str(a16),'a16_decimal':str(a16.numerator*scale//a16.denominator//scale)+'.'+str(a16.numerator*scale//a16.denominator%scale).zfill(75),'delta_upper':float(delta),'loss_upper':float(loss),'omega':omega,'conditional_Q_lower':str(F(89999999,10**8)*H),'complete_heads_certified':False}
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
