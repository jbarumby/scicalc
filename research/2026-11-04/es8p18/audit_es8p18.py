#!/usr/bin/env python3
"""ES8-P18: exact rational uniform pole-projection and tail estimates.

Inputs: historic B_N formula and E7 conditional archimedean tail bound.
No assertion about analytic E7 validity, Weil positivity, RH, or uniform gaps.
"""
from fractions import Fraction as F
import json, argparse
from pathlib import Path

L=F(9,20)
A_ODD=F(409,400)
A_EVEN=F(1609,800)
OLD=F(5148881,2560000)
ODD_SERIES_BOUND=F(1,4)
EVEN_SERIES_BOUND=F(25,144)
ODD_SQ=A_ODD**2*ODD_SERIES_BOUND
EVEN_SQ=A_EVEN**2*EVEN_SERIES_BOUND
KAPPA=1+max(ODD_SQ,EVEN_SQ)
KAPPA_SIMPLE=F(171,100)
RAW_DELTA=F(406,91125)

def coeff(k:int, X:F)->F:
    assert k>=3 and X>=9
    p=1 if k%2 else 2
    return F(k,p)*(p*p*X+L*L)/(k*k*X+L*L)

def parity_sums(n:int,X:F)->tuple[F,F]:
    assert n>=3 and X>=9
    odd=sum((coeff(k,X)**2 for k in range(3,n+1,2)),F(0))
    even=sum((coeff(k,X)**2 for k in range(4,n+1,2)),F(0))
    return odd,even

def exact_norm_sq(n:int,X:F)->F:
    u,v=parity_sums(n,X)
    return 1+max(u,v)

def odd_reciprocal_bound(t:int)->bool:
    # 1/(2m+1)^2 <= 1/(4m(m+1)) for m>=1; telescopes to 1/4
    return sum((F(1,(2*m+1)**2) for m in range(1,t+1)),F(0)) <= F(1,4)*(1-F(1,t+1))

def even_reciprocal_bound(t:int)->bool:
    # First terms m=2,3 exact; for m>=4, 1/(4m^2) <= 1/[4m(m-1)].
    if t<2: return True
    vals=sum((F(1,4*m*m) for m in range(2,t+1)),F(0))
    if t==2: return vals==F(1,16)
    telescopic= F(1,16)+F(1,36)+F(1,4)*(F(1,3)-F(1,t))
    return vals<=telescopic<=F(25,144)

def full_validation():
    assert ODD_SQ==F(167281,640000)
    assert EVEN_SQ==F(2588881,3686400)
    assert KAPPA==F(6275281,3686400) < KAPPA_SIMPLE < OLD
    assert KAPPA_SIMPLE*RAW_DELTA==F(3857,506250) < F(1,100)
    assert F(1,100)-KAPPA_SIMPLE*RAW_DELTA==F(2411,1012500)
    # C(k,X) decreases with X for k>p; X=9 is a rational coefficient envelope
    rows=[]
    for n in (3,4,5,6,8,10,16,20,40,80,120):
        for X in (F(9),F(10)):
            odd,even=parity_sums(n,X)
            assert odd<=ODD_SQ and even<=EVEN_SQ
            assert exact_norm_sq(n,X)<=KAPPA
            for k in range(3,n+1):
                a=A_ODD if k%2 else A_EVEN
                assert 0<coeff(k,X)<=a/k
        assert exact_norm_sq(n,F(10))<=exact_norm_sq(n,F(9))
        o,e=parity_sums(n,F(9))
        rows.append({'N':n,'odd_sum_at_X9':str(o),'even_sum_at_X9':str(e),'norm_sq_upper_at_X9':str(1+max(o,e))})
    for t in range(1,150):
        assert odd_reciprocal_bound(t)
        assert even_reciprocal_bound(t)
    return {'schema':'NUMERIA.ES8P18.POLE_PROJECTION.1',
     'scope':'EXACT_RATIONAL_PROOF_CHECKS; INFINITE_ESTIMATE_PAPER_PROOF; E7_TAIL_CONDITIONAL',
     'input_formula':'(k/p)*(p^2*X+L^2)/(k^2*X+L^2), X=pi^2>=9, L=9/20, p=1 odd/p=2 even',
     'odd_inverse_square_sum_upper':'1/4',
     'even_inverse_square_sum_upper':'25/144',
     'odd_coefficient_sq_sum_upper':str(ODD_SQ),
     'even_coefficient_sq_sum_upper':str(EVEN_SQ),
     'improved_norm_squared_upper':str(KAPPA),
     'convenient_norm_squared_bound':str(KAPPA_SIMPLE),
     'historical_norm_squared_bound':str(OLD),
     'ES8P15_raw_tail_at_N1':str(RAW_DELTA),
     'conditional_projected_arch_tail_upper':str(KAPPA_SIMPLE*RAW_DELTA),
     'conditional_projected_arch_tail_budget':'1/100',
     'slack':str(F(1,100)-KAPPA_SIMPLE*RAW_DELTA),
     'operator_extension_statement':'B_infinity x=(sum_odd c_k x_k,sum_even c_k x_k,x), bounded in l2; nonfinite-rank extension; proof in report conditional on exact Fourier coefficients',
     'pivot_correction_truncation':'||C-C_N||op <= (1609/800)/sqrt(N) for N>=3, for pivot correction only',
     'warning':'The full finite-rank truncations of B_infinity do NOT converge in operator norm because the identity tail remains',
     'not_proved':['analytic E7 tail for true Weil integrals','identification of global Weil form','uniform positive spectral gap','RH'],
     'test_dimensions':rows}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',default='RESULTADOS_ES8P18.json')
    args=ap.parse_args()
    result=full_validation()
    Path(args.out).write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+'\n')
    print('ES8-P18 passed exact rational checks: kappa <= 171/100, projected conditional tail < 1/100')

if __name__=='__main__':main()
