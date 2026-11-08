#!/usr/bin/env python3
"""ES8-P83: independent exact rational series/identities vs 256-bit FLINT balls."""
import json,sys
from fractions import Fraction as Q
from pathlib import Path
from flint import arb,acb,ctx
import flint
ctx.prec=256
assert ctx.prec==256
checks=[]
def endpoints(a):
    m,e=a.mid().man_exp(); r,f=a.rad().man_exp()
    mid=Q(int(m))*Q(2)**int(e); rad=Q(int(r))*Q(2)**int(f)
    return mid-rad,mid+rad
def overlap(x,y):return max(x[0],y[0])<=min(x[1],y[1])
def verify(label,actual,oracle):
    if not overlap(endpoints(actual),oracle):
        raise AssertionError('DISJOINT: '+label)
    checks.append(label)
def dfact(p):
    d=1
    for k in range(1,2*p+2,2):d*=k
    return d
def bessel_ball(p,x,modified=False,K=90):
    x=Q(x);s=Q(0);term=Q(1);sign=1 if modified else -1
    for k in range(K):
        s+=term
        term*=sign*x*x/Q(2*(k+1)*(2*p+2*k+3))
    ratio=x*x/Q(2*(K+1)*(2*p+2*K+3))
    if ratio>=1:raise ArithmeticError('bessel geometric ratio >=1')
    tail=abs(term)/(1-ratio)
    fac=x**p/dfact(p)
    return fac*(s-tail),fac*(s+tail)
for p in (0,1,12,21,25):
    for x in (Q(0),Q(1,10),Q(1),Q(3),Q(10)):
        z=arb(x.numerator)/x.denominator
        j=z**p/dfact(p)*(-z*z/4).hypgeom_0f1(arb(p)+arb(3)/2)
        i=z**p/dfact(p)*(z*z/4).hypgeom_0f1(arb(p)+arb(3)/2)
        verify(f'j({p},{x})',j,bessel_ball(p,x))
        verify(f'i({p},{x})',i,bessel_ball(p,x,True))
def digamma_increment(r,N=2048):
    r=Q(r);u2=r*r/4;s=Q(0)
    for k in range(N):
        q=Q(4*k+1,4)
        s+=u2/(q*(q*q+u2))
    return s,s+u2/(2*(Q(N)-Q(3,4))**2)
psi0=acb(arb(1)/4).digamma().real
for r in (0,1,5,15,30):
    s=acb(0,arb(r)/2); b=acb(arb(1)/4)
    d=((b+s).digamma()+(b-s).digamma()).real/2-psi0
    verify('psi_delta_'+str(r),d,digamma_increment(r))
for n in (16,32,64,110):
    weights=arb(0)
    for k in range(n):
        x,w=arb.legendre_p_root(n,k,weight=True)
        lo,hi=endpoints(x);wl,wh=endpoints(w)
        if not (-1<lo<=hi<1 and wl>0):raise AssertionError('bad GL node/weight')
        weights+=w
        if k in (0,n//4,n//2,n-1):
            prev,curr=arb(1),x
            for d in range(2,n+1):
                newer=((2*d-1)*x*curr-(d-1)*prev)/d
                prev,curr=curr,newer
            if not overlap(endpoints(curr),(Q(0),Q(0))):
                raise AssertionError('P_n(x) lacks zero')
            deriv=n*(x*curr-prev)/(x*x-1)
            expected=2/((1-x*x)*deriv*deriv)
            verify(f'weight_{n}_{k}',w,endpoints(expected))
    verify('weight_sum_'+str(n),weights,(Q(2),Q(2)))
out={'status':'FLINT_COMPONENT_CHECKS_PASSED','bits':ctx.prec,
     'python_flint':flint.__version__,'checks':checks,
     'count':len(checks),'scope':'rational Bessel and digamma oracles plus independent Legendre weight equations, not an end-to-end RH certificate'}
dest=Path(sys.argv[1] if len(sys.argv)>1 else 'ES8P83_COMPONENT_AUDIT.json')
if dest.exists():raise FileExistsError(dest)
dest.write_text(json.dumps(out,indent=2)+'\n')
print('PASS',len(checks),'independent component checks; 256 bits')
