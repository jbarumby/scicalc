#!/usr/bin/env python3
"""ES8-P17 — exact rational audit of pole-killing B_N projection.

No floating point is used in certificates; decimal strings are display only.
Analytic premise: pi^2 >= 9; historical formula for B_N.
This module proves no global Weil positivity or archimedean integral enclosure.
"""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import argparse
import json

L = F(9, 20)
COEFF_ODD = F(409, 400)
COEFF_EVEN = F(1609, 800)
ODD_SSQ_MAX = COEFF_ODD**2 / 2
EVEN_SSQ_MAX = COEFF_EVEN**2 / 4
KAPPA_EXACT = 1 + max(ODD_SSQ_MAX, EVEN_SSQ_MAX)
KAPPA_BUDGET = F(21, 10)
DELTA_N1 = F(16, 3645) + F(2, 30375)


def coeff(k: int, x: F) -> F:
    assert k >= 3 and x >= 9
    p = 1 if k % 2 else 2
    return F(k, p) * (p*p*x + L*L)/(k*k*x + L*L)


def alpha(k: int) -> F:
    return COEFF_ODD if k % 2 else COEFF_EVEN


def bound_coefficient(k: int, x: F) -> bool:
    return 0 < coeff(k,x) <= alpha(k)/k


def norm_bound_exact(N: int, x: F) -> F:
    """Exact ||B_N||^2 = 1+max(sum c_odd^2, sum c_even^2).
    B_N is [C_N;I], C_N has two rows with disjoint parity supports.
    """
    if N < 3: raise ValueError('N>=3 required')
    odd = sum((coeff(k,x)**2 for k in range(3,N+1,2)), F())
    even = sum((coeff(k,x)**2 for k in range(4,N+1,2)), F())
    return 1 + max(odd,even)


def projected_gram(N: int, x: F) -> list[list[F]]:
    """Explicit Gram B'B; off-diagonal = two shared-pivot rank-one parts."""
    ks = list(range(3,N+1))
    return [[F(int(k==j)) + (coeff(k,x)*coeff(j,x) if k%2==j%2 else F())
             for j in ks] for k in ks]


def b_matrix(N: int, x: F) -> list[list[F]]:
    arr=[[F() for j in range(N-2)] for i in range(N)]
    for j,k in enumerate(range(3,N+1)):
        p=1 if k%2 else 2
        arr[p-1][j]=-coeff(k,x)
        arr[k-1][j]=1
    return arr


def gram_multiplication(N: int, x: F) -> list[list[F]]:
    b=b_matrix(N,x)
    return [[sum((b[r][i]*b[r][j] for r in range(N)),F()) for j in range(N-2)]
            for i in range(N-2)]


def delta(N: int) -> F:
    if N < 1: raise ValueError('N>=1')
    return F(16,3645*N) + F(2,30375*N**3)


def check_premises():
    assert L*L / F(9) == F(9,400)
    assert COEFF_ODD == 1+L*L/9
    assert COEFF_EVEN == (4+L*L/9)/2
    assert ODD_SSQ_MAX == F(167281,320000)
    assert EVEN_SSQ_MAX == F(2588881,2560000)
    assert KAPPA_EXACT == F(5148881,2560000) < KAPPA_BUDGET
    assert DELTA_N1 == F(406,91125)
    assert KAPPA_BUDGET*DELTA_N1 == F(1421,151875)<F(1,100)
    assert F(1,100)-KAPPA_BUDGET*DELTA_N1 == F(391,607500)


def run():
    check_premises()
    details=[]
    for n in (6,8,10,20,40,80):
        k_lo=norm_bound_exact(n,F(10))
        k_hi=norm_bound_exact(n,F(9))
        assert 1 <= k_lo <= k_hi <= KAPPA_EXACT < KAPPA_BUDGET
        # at x=9 use larger coefficients since c_k is decreasing for x>0 and k>p
        assert all(bound_coefficient(k,F(9)) for k in range(3,n+1))
        if n<=20:
            for x in (F(9),F(10)):
                assert projected_gram(n,x)==gram_multiplication(n,x)
        assert KAPPA_BUDGET*delta(n)<F(1,100)
        details.append({'N':n,'norm_squared_lower_at_X10':str(k_lo),
                        'norm_squared_upper_at_X9':str(k_hi),
                        'upper_decimal_diagnostic':f'{float(k_hi):.12f}',
                        'raw_tail_upper_fraction':str(delta(n)),
                        'projected_tail_budget_fraction':str(KAPPA_BUDGET*delta(n))})
    return {
        'schema':'NUMERIA.ES8P17.PROJECTION_BUDGET.1',
        'status':'EXACT_RATIONAL_CONDITIONAL_ON_E7_TAIL_AND_PI_GE_3',
        'L':'9/20',
        'projection_formula':'B[k-1,k-3]=1; B[p-1,k-3]=-(k/p)*(p^2*X+L^2)/(k^2*X+L^2); p=1 odd, 2 even; X=pi^2',
        'proof':'B^T B=I+u_odd u_odd^T+u_even u_even^T with disjoint supports; Cauchy-Schwarz and telescoping square-sum bounds',
        'odd_amplitude_upper':'409/400',
        'even_amplitude_upper':'1609/800',
        'odd_squared_coefficients_sum_upper':str(ODD_SSQ_MAX),
        'even_squared_coefficients_sum_upper':str(EVEN_SSQ_MAX),
        'norm_squared_upper':str(KAPPA_EXACT),
        'convenient_uniform_norm_squared_upper':str(KAPPA_BUDGET),
        'truncation':'R(N)=200*N^2',
        'raw_tail_upper_at_N1':str(DELTA_N1),
        'uniform_projected_tail_upper':str(KAPPA_BUDGET*DELTA_N1),
        'projected_tail_budget':'1/100',
        'budget_slack':str(F(1,100)-KAPPA_BUDGET*DELTA_N1),
        'limits':['E7 true-integral tail estimate remains an analytic premise',
                  'Bound applies only to archimedean raw tail after projection',
                  'Interpolation and prime terms are NOT covered',
                  'No uniform spectral gap for the central Weil matrices',
                  'No infinite Hilbert-Schmidt proof from finite samples'],
        'dimensions':details,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output', default='RESULTADOS_ES8P17.json')
    args=parser.parse_args()
    obj=run()
    Path(args.output).write_text(json.dumps(obj,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('ES8-P17: rational projection bounds passed; kappa <= 21/10; projected raw tail <1/100 (conditional E7)')

if __name__=='__main__':main()
