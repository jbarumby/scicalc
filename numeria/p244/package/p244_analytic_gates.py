#!/usr/bin/env python3
"""ES8-P244: independent exact-rational checks of elementary input budgets.

This verifies finite rational inequalities, NOT analytic theorems such as the
Weil explicit formula, Binet integral representation, Gauss quadrature or the
Legendre/Bessel integral identity (which must be separately justified).
"""
from fractions import Fraction as F
from math import isqrt
import json
import hashlib
from pathlib import Path

L = F(161, 200)
T = 504


def exp_lower(x: F, n: int) -> F:
    """Finite positive Taylor sum <= exp(x) for x >= 0."""
    if x < 0 or n < 0:
        raise ValueError('nonnegative x and n required')
    term = total = F(1)
    for k in range(1, n+1):
        term *= x / k
        total += term
    return total


def exp_upper(x: F, n: int) -> F:
    """Strict upper bound on exp(x) by a geometric majorant for its tail."""
    if x < 0 or n+2 <= x:
        raise ValueError('invalid exponential upper-bound parameters')
    term=F(1)
    for k in range(1,n+1):
        term *= x/k
    first=term*x/(n+1)
    return exp_lower(x,n)+first/(1-x/F(n+2))


def ceil_sqrt(q: F) -> int:
    if q < 0:
        raise ValueError('sqrt negative')
    c = isqrt(q.numerator // q.denominator)
    while c*c * q.denominator < q.numerator:
        c += 1
    return c


def verify():
    # Strict elementary log bounds, by positive rational Taylor sums.
    for x, threshold, n in ((F(7,10),2,3), (F(11,10),3,5),
                            (F(161,100),5,7), (F(28,5),253,9)):
        if not exp_lower(x,n)>threshold:
            raise ArithmeticError(f'exp lower failed {x}')
    # Bounds used for A_L: log2<.7; log3<1.1; log5<1.61;
    # sqrt2>1.4, sqrt3>1.7, sqrt5>2.2.
    if not (F(7,5)**2<2 and F(17,10)**2<3 and F(11,5)**2<5):
        raise ArithmeticError('sqrt comparisons')
    a_upper=F(2)*F(7,10)/F(7,5) + F(2)*F(11,10)/F(17,10) + F(7,10) + F(2)*F(161,100)/F(11,5)
    if not a_upper<5:raise ArithmeticError('A_L<5')
    # Principal logarithm for z=5/4+iy, 0<=y<=252.
    # |z|<253, Re(z)>=5/4, arg z in [0,pi/2].
    # log |z| < 28/5 since exp(28/5)>253 above. 
    # pi < 22/7 is classical and gives arg z < 11/7.
    if not (F(5,4)**2 + 252**2 < 253**2 and
            F(28,5)**2 + F(11,7)**2 < 6**2):
        raise ArithmeticError('|log z|<6 proof missing')
    m_bound = F(6)+F(2,5)+F(4,75)+F(4)+F(2)+F(5)+F(1)
    if m_bound != F(1384,75) or not m_bound<20:
        raise ArithmeticError('M20 bound')
    # Complex E3 geometrical mapping (1008 half-unit panels).
    real_t = F(3025,6)
    imag_t = F(1,3)
    if real_t != F(2*1007+1,4)+F(5,12) or imag_t != (F(3)-F(1,3))/8:
        raise ArithmeticError('Bernstein ellipse geometry')
    if not (F(1,12)>0 and F(5,12)+F(3025,12)<253):
        raise ArithmeticError('digamma strip')
    if not (F(274,7)<40 and F(274,7)<50):
        raise ArithmeticError('complex symbol sup')
    # Chebyshev k>=120 tail + positivity/exactness of GL60 leads to E=1008*50 a_i a_j/3^120.
    e_factor=F(1008*50,3**120)
    aa=[];df=1
    for n in range(754):
        df*=2*n+1
        r = min(F(2), F(2*406**n, df))
        aa.append(ceil_sqrt(2*L*(2*n+1))*r)
    rows={}
    for parity,label in ((0,'even'),(1,'odd')):
        aa_s=aa[parity::2]
        exact_upper=e_factor*max(aa_s)*sum(aa_s)
        if not (0<exact_upper<F(4,10**47)):
            raise ArithmeticError('GL60 row bound')
        rows[label]={ 'row_error_upper_strict':'4/10^47',
            'max_mode_amplitude_exact':str(max(aa_s)),
            'sum_mode_amplitudes_sha256':hashlib.sha256(str(sum(aa_s)).encode()).hexdigest(),
            'row_error_sci_illustrative':f'{float(exact_upper):.12e}'}
    # Sharper digamma lower bound at all t>=T by the Binet split.
    # low integral <=2/(3t^2), high integral <= (2/3+4/(9t))*16/(9t^3)
    t=F(T)
    high=(F(2,3)+F(4)/(9*t))*F(16)/(9*t**3)
    if not high<F(1)/(3*t*t):
        raise ArithmeticError('Binet high integral')
    # Hence full Re digamma > log(t/2)-3/(2t^2).
    high_freq_floor=F(1,T)-F(3,2*T*T)
    if high_freq_floor != F(335,169344) or high_freq_floor<=0:
        raise ArithmeticError('uniform positivity above T')
    if not (F(5)<exp_lower(2*L,8) and exp_upper(2*L,8)<6):
        raise ArithmeticError('prime-power cutoff 5<exp(2L)<6 not proved')
    return {
      'schema':'NUMERIA.ES8P244.EXACT_ANALYTIC_BUDGET_CHECK.1',
      'parameter_L':'161/200', 'Tsharp':T,
      'imported_classical_lemmas':['Binet representation for psi', 'Legendre plane wave and Rodrigues',
       'Bernstein Chebyshev coefficient bound', 'Gauss-Legendre 60 exactness and positive weights',
       'Weil explicit formula / Zhu one-stroke reduction'],
      'beta_504_sign_is_imported_from_P237_rational_certificate': True,
      'digamma_Binet_complex_log_repair':'proved |log(5/4+iy)|<6, 0<=y<=252, using exp(28/5)>253 and pi<22/7',
      'symbol_M20_rational_bound':str(m_bound),
      'Bernstein_ellipse_rho':3,
      'symbol_complex_bound':'274/7 < 50',
      'GL60_entry_error_factor_exact':str(e_factor),
      'GL60_row_bounds_by_parity':rows,
      'high_frequency_symbol_excess_strict_lower_bound':str(high_freq_floor),
      'L_0p805_window_weil_positivity_proved_independently':False,
      'RH_proved':False}

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();r=verify();args.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n')
    print('P244 rational analytic budgets: PASS')
