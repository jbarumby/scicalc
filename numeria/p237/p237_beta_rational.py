#!/usr/bin/env python3
"""Independent exact rational enclosure for beta(T=504), L=161/200.

Uses pi interval 333/106 < pi < 355/113, rational atanh series
with a geometric remainder, and integer-square-root enclosures.
This DOES NOT prove the full analytic Weil reduction or matrix positivity.
"""
from __future__ import annotations
import json,argparse,math
from fractions import Fraction as Q
from pathlib import Path

PILO=Q(333,106);PIHI=Q(355,113)

def log_series_at_z(z:Q,terms=96):
    if not Q(0)<=z<Q(1):raise ValueError('z outside [0,1)')
    s=Q(0);power=z
    for k in range(terms):
        s+=2*power/(2*k+1)
        power*=z*z
    tail=2*power/(Q(2*terms+1)*(1-z*z))
    return s,s+tail

L2LO,L2HI=log_series_at_z(Q(1,3))

def log_interval(x:Q):
    if x<=0:raise ValueError('log must be positive')
    s=0
    while x>=2:x/=2;s+=1
    while x<1:x*=2;s-=1
    z=(x-1)/(x+1)
    lo,hi=log_series_at_z(z)
    if s>=0:return lo+s*L2LO,hi+s*L2HI
    return lo+s*L2HI,hi+s*L2LO


def inverse_sqrt(n:int,bits=160):
    if n<=0:raise ValueError('nonpositive sqrt')
    floor=math.isqrt(n<<(2*bits))
    low=Q(floor,1<<bits)
    high=Q(floor+1,1<<bits)
    assert low*low<=n and high*high>n
    return 1/high,1/low


def beta_bounds():
    logs={p:log_interval(Q(p)) for p in (2,3,5)}
    sq={p:inverse_sqrt(p) for p in (2,3,5)}
    # A = (2/sqrt2)*log2 + (2/sqrt3)*log3 + log2 + (2/sqrt5)*log5
    A_lo=2*sq[2][0]*logs[2][0]+2*sq[3][0]*logs[3][0]+logs[2][0]+2*sq[5][0]*logs[5][0]
    A_hi=2*sq[2][1]*logs[2][1]+2*sq[3][1]*logs[3][1]+logs[2][1]+2*sq[5][1]*logs[5][1]
    log_T_lo=log_interval(Q(504)/(2*PIHI))[0]
    log_T_hi=log_interval(Q(504)/(2*PILO))[1]
    betalo=log_T_lo-Q(1,504)-A_hi
    betahi=log_T_hi-Q(1,504)-A_lo
    if not 0<betalo<betahi<2:raise ArithmeticError('rational beta positivity not certified')
    return betalo,betahi,A_lo,A_hi


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path)
    a=p.parse_args();lo,hi,Alo,Ahi=beta_bounds()
    # Conservative decimal comparisons, proven using the original Q fractions.
    gates={'beta_gt_1_per_1000':lo>Q(1,1000),
           'beta_lt_2_per_1000':hi<Q(2,1000),
           'A_between_4_and_5':Q(4)<Alo<Ahi<Q(5)}
    if not all(gates.values()):raise ArithmeticError('rational gate failure')
    d={'schema':'NUMERIA.ES8P237.BETA_RATIONAL.1',
       'Tsharp':504,'window_L':'161/200','prime_terms':[2,3,4,5],
       'method':'rational atanh at z<=1/3, 96 terms, exact geometric tail; integer sqrt bounds; pi 333/106..355/113',
       'pi_lower':str(PILO),'pi_upper':str(PIHI),
       'beta_lower_fraction':str(lo),'beta_upper_fraction':str(hi),
       'beta_lower_decimal_diagnostic':format(float(lo),'.12g'),
       'beta_upper_decimal_diagnostic':format(float(hi),'.12g'),
       'A_lower_decimal_diagnostic':format(float(Alo),'.12g'),
       'A_upper_decimal_diagnostic':format(float(Ahi),'.12g'),
       'exact_gates':gates,
       'requires_independent_analytic_bridge':True,'RH_proved':False}
    a.out.write_text(json.dumps(d,indent=2)+'\n')
    print('PASS rational β in (1/1000, 2/1000)',float(lo),float(hi))

if __name__=='__main__':main()
