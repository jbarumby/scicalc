"""Exact rational Schur gates; same conditional hypotheses as V38."""
from fractions import Fraction as F
from math import factorial
import schur38 as s
import json
r=json.load(open('hybrid39.json'));C=F(r['c']);D=F(r['d']);lam=F(r['lambda']);f1=F(r['f1']);a=s.A;s.T=4500;h=s.H;m=s.M;eps=s.EPS
rows=[]
for theta in [F(0),F(1,2),F(1)]:
 pen=(1-theta)*C+theta*f1;coef=50+theta*D/(1-lam)
 best=None
 for rho in [F(i,1000) for i in range(1180,1321,5)]:
  s.RHO=rho;s.E=s.expup(a*s.T*(rho-1/rho))
  # Monotone g permits binary search for exact smallest N at each rho in this range.
  def gate(N):
   g=s.g(N,0);reserve=2-pen-50*g-eps-m;slack=(h-m)*reserve-2*coef**2*g-2*eps**2
   return reserve>0 and slack>0,g,slack
  lo,hi=3000,3600
  if not gate(hi)[0]:continue
  while hi-lo>1:
   mid=(lo+hi)//2
   if gate(mid)[0]:hi=mid
   else:lo=mid
  ok,g,slack=gate(hi)
  if best is None or hi<best['N']:
   if not(F(48,factorial(2*hi))<eps and F(72,factorial(2*hi)**2)<eps):raise ArithmeticError('poles')
   best={'theta':str(theta),'rho':str(rho),'N':hi,'g_even':float(g),'slack_even':float(slack),'reserve_limit':float(2-pen),'both_parities':True,'conditional_on_full_heads':True}
 rows.append(best)
print(json.dumps(rows,indent=2))
