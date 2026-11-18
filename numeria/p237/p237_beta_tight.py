#!/usr/bin/env python3
"""Second, tighter independent rational beta(504) enclosure.

Machin's pi identity, alternating rational arctan bounds, a rational
atanh-series with geometric tail, and integer sqrt brackets. All acceptance
checks are exact Fraction comparisons; mpmath/float is never needed.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from fractions import Fraction as F


def atan_inv_bounds(q,terms=72):
    if q<=1 or terms<=0 or terms%2!=0:raise ValueError('invalid alternating arctan setup')
    termsum=F()
    for k in range(terms):
        term=F((-1)**k,(2*k+1)*q**(2*k+1))
        termsum+=term
    # N even: first neglected term is positive, so true atan ∈ [S_N,S_N+next].
    next_term=F(1,(2*terms+1)*q**(2*terms+1))
    return termsum,termsum+next_term


def pi_machin():
    al,au=atan_inv_bounds(5);bl,bu=atan_inv_bounds(239)
    lo=16*al-4*bu;hi=16*au-4*bl
    if not F(333,106)<lo<hi<F(355,113):raise ArithmeticError('pi gate')
    return lo,hi


def atanh_interval(z, terms=110):
    if not F(0)<=z<=F(1,3):raise ValueError('invalid atanh z')
    s=F();power=z
    for k in range(terms):
        s+=2*power/(2*k+1)
        power*=z*z
    tail=2*power/(F(2*terms+1)*(1-z*z))
    return s,s+tail

def log_core(x,terms=110):
    if x<=0:raise ValueError('nonpositive log')
    scale=0
    while x>=2:x/=2;scale+=1
    while x<1:x*=2;scale-=1
    z=(x-1)/(x+1)
    if not F(0)<=z<=F(1,3):raise ValueError('normalization error')
    lo,hi=atanh_interval(z,terms)
    return lo,hi,scale


def log_enclosure(x, log2,terms=110):
    lo,hi,k=log_core(x,terms)
    if k>=0:return lo+k*log2[0],hi+k*log2[1]
    return lo+k*log2[1],hi+k*log2[0]


def inv_sqrt(n,bits=256):
    m=math.isqrt(n<<(2*bits));l=F(m,1<<bits);h=F(m+1,1<<bits)
    if not l*l<=n<h*h:raise ArithmeticError('sqrt gate')
    return 1/h,1/l


def beta_tight():
    pi_lo,pi_hi=pi_machin()
    log2=atanh_interval(F(1,3))
    logs={p:log_enclosure(F(p),log2) for p in (2,3,5)}
    sqrt={p:inv_sqrt(p) for p in (2,3,5)}
    Alo=(2*sqrt[2][0]*logs[2][0]+2*sqrt[3][0]*logs[3][0]
         +logs[2][0]+2*sqrt[5][0]*logs[5][0])
    Ahi=(2*sqrt[2][1]*logs[2][1]+2*sqrt[3][1]*logs[3][1]
         +logs[2][1]+2*sqrt[5][1]*logs[5][1])
    tl=log_enclosure(F(504)/(2*pi_hi),log2)[0]
    tu=log_enclosure(F(504)/(2*pi_lo),log2)[1]
    lo=tl-F(1,504)-Ahi
    hi=tu-F(1,504)-Alo
    if not F(1,1000)<lo<hi<F(2,1000):raise ArithmeticError('beta gate')
    if not hi-lo<F(1,10**70):raise ArithmeticError('beta precision below gate')
    return lo,hi,(pi_lo,pi_hi),(Alo,Ahi)


def decimal_enclosure(lo,hi,digits=65):
    s=10**digits
    a=lo.numerator*s//lo.denominator
    b=-((-hi.numerator*s)//hi.denominator)
    def render(x):
        sign='-' if x<0 else ''
        v=str(abs(x)).zfill(digits+1)
        return sign+v[:-digits]+'.'+v[-digits:]
    return render(a),render(b)


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    lo,hi,pi,A=beta_tight()
    low,high=decimal_enclosure(lo,hi)
    data={'schema':'NUMERIA.ES8P237.BETA_MACHIN_RATIONAL.1',
      'Tsharp':504,'window_L':'161/200','prime_terms':[2,3,4,5],
      'method':'Machin pi with alternating arctan 72 terms; atanh logs 110 terms; integer sqrt 256 bits',
      'beta_lower_decimal_outward_65':low,
      'beta_upper_decimal_outward_65':high,
      'beta_width_proven_lt_1e-70':hi-lo<F(1,10**70),
      'pi_inside_original_rational_bounds':F(333,106)<pi[0] and pi[1]<F(355,113),
      'positive_beta':lo>0,'A_between_4_and_5':F(4)<A[0] and A[1]<F(5),
      'analytic_Weil_Zhu_bridge_proven':False,'RH_proved':False}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(data,indent=2)+'\n')
    print('PASS certified beta between',low,high)

if __name__=='__main__':main()
