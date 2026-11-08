#!/usr/bin/env python3
"""ES8-P90: independent short Python-FLINT API smoke test for interval Riemann
and cancellation-free exp quotient. NOT a full Weil matrix certificate."""
from fractions import Fraction as F
from flint import arb,fmpq,ctx
import flint
ctx.prec=256
def q(x):
    if not isinstance(x,F):x=F(x)
    return arb(fmpq(x.numerator,x.denominator))
def cell(a,b):
    return arb(fmpq(((a+b)/2).numerator,((a+b)/2).denominator),
               fmpq(((b-a)/2).numerator,((b-a)/2).denominator))
def eps_bound(x):
    lo=x.lower().fmpq();hi=x.upper().fmpq()
    return F(int(lo.p),int(lo.q)),F(int(hi.p),int(hi.q))
def rimann(f,a,b,N):
    total=arb(0)
    for k in range(N):
        lo=a+(b-a)*F(k,N); hi=a+(b-a)*F(k+1,N)
        total+=q(hi-lo)*f(cell(lo,hi))
    return total
c=cell(F(1,3),F(2,3))
lo,hi=eps_bound(c)
assert lo<=F(1,3) and hi>=F(2,3)
assert ctx.prec==256
I=rimann(lambda x:x*x,F(0),F(1),2048)
lo,hi=eps_bound(I)
assert lo<=F(1,3)<=hi
assert hi-lo<F(1,100)
z=cell(F(0),F(1,10**6))
x=(-z/4).exp()
stable=-x*(1+x+x*x)/(1+x+x*x+x*x*x)
lo,hi=eps_bound(stable)
assert lo<=F(-3,4)<=hi
h=q(F(2,5))
def bump(X):
    if bool(X.lower()>=h or X.upper()<=-h):return arb(0)
    if not bool(X.lower()>-h and X.upper()<h):
        return cell(F(0),F(1))
    u=X/h
    den=1-u*u
    if not bool(den.lower()>0):return cell(F(0),F(1))
    return (1-1/den).exp()
R0=rimann(lambda X:bump(X)*bump(X),F(-2,5),F(2,5),4096)
lo,hi=eps_bound(R0)
assert lo<F(3933523,10000000)<hi
assert hi-lo<F(1,100)
print('SMOKE_PASS',flint.__version__,ctx.prec,'R0',str(R0),'quadratic_integral',str(I),'base0',str(stable))
