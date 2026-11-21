"""NUMERIA ES8-P260: rational Euler/Weierstrass audit of the archimedean tail.

Independent of Arb/MPFR, floating point, and Riemann hypothesis.
Proof: see PROVA_P260.md. This DOES NOT certify Weil positivity.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json

SCALE = 10**30
REFERENCE_FREQUENCY = 504
PRIME_POWERS = ((2, 2), (3, 3), (4, 2), (5, 5))
FLOORS = ((850, F(1, 2)), (1100, F(3, 4)),
          (1500, F(1)), (2000, F(5, 4)))

def require(condition: bool, message: str):
    if not condition:
        raise ArithmeticError(message)

def log_interval(x: F, terms: int = 80):
    """Exact rational bounds: log(x)=2*atanh((x-1)/(x+1))."""
    x = F(x)
    if x <= 1 or terms < 1:
        raise ValueError('require x > 1 and terms >= 1')
    v = (x - 1) / (x + 1)
    lo = 2 * sum((v**(2*k+1) / (2*k+1) for k in range(terms)), F(0))
    tail = 2*v**(2*terms+1) / ((2*terms+1)*(1-v*v))
    return lo, lo+tail

def pi_interval(terms=12):
    """Machin identity, alternating arctan bounds, exact rational."""
    if terms < 1:
        raise ValueError('terms must be positive')
    def atan_bracket(den):
        u = F(1, den)
        lo = sum(((-1)**k*u**(2*k+1)/(2*k+1)
                  for k in range(2*terms)), F(0))
        hi = lo + u**(4*terms+1)/(4*terms+1)
        return lo, hi
    a, b = atan_bracket(5)
    c, d = atan_bracket(239)
    lo, hi = 16*a-4*d, 16*b-4*c
    require(F(314159, 100000)<lo<hi<F(355, 113), 'pi brackets invalid')
    return lo, hi

def euler_gamma_upper(harmonic_n=4096):
    """gamma < H_N-log N; H_N upper from integer ceilings."""
    if harmonic_n < 1 or harmonic_n & (harmonic_n-1):
        raise ValueError('N must be positive power of 2')
    h_up = F(sum((SCALE+k-1)//k for k in range(1, harmonic_n+1)), SCALE)
    log2_lower = log_interval(2)[0]
    return h_up - (harmonic_n.bit_length()-1)*log2_lower

def sqrt_lower(n: int):
    if type(n) is not int or n <= 0:
        raise ValueError('n must be positive integer')
    lo = F(isqrt(n*SCALE*SCALE), SCALE)
    require(lo*lo <= n < (lo+F(1,SCALE))**2, 'root enclosure invalid')
    return lo

def prime_amplitude_upper():
    return sum((2*log_interval(p)[1]/sqrt_lower(q) for q,p in PRIME_POWERS), F(0))

def psi_at_504_lower(count=8192):
    """Lower bound for Re psi(1/4+252i), Euler series with tail."""
    if count < 1:
        raise ValueError('count must be positive')
    b = REFERENCE_FREQUENCY // 2
    partial_scaled = 0
    for n in range(count):
        m = 4*n + 1
        numerator = 16*b*b - 3*m
        denominator = (n+1)*(m*m + 16*b*b)
        partial_scaled += (numerator*SCALE)//denominator
    x = F(4*count+1, 4)
    missing_tail_lower = -F(3,4)*(1/x + 1/(x*x))
    return -euler_gamma_upper() + F(partial_scaled, SCALE) + missing_tail_lower

def archimedean_anchor_lower():
    """Lower bound for Re psi(1/4+252i)-log(pi)-sum prime amplitudes."""
    pi_interval()
    return psi_at_504_lower()-log_interval(F(355,113))[1]-prime_amplitude_upper()

def digamma_increment_lower(frequency: int, count: int = 2048):
    """Positive-term series for Re psi(1/4+it/2)-Re psi(1/4+252i).

    Every omitted term is positive. Individual terms are rounded DOWN
    by exact integer arithmetic. No Binet estimates needed.
    """
    if type(frequency) is not int or frequency <= REFERENCE_FREQUENCY or frequency % 2:
        raise ValueError('frequency must be an even integer > 504')
    if count < 1:
        raise ValueError('count must be positive')
    t0 = REFERENCE_FREQUENCY
    d = frequency*frequency - t0*t0
    s = 0
    for n in range(count):
        m=4*n+1
        num=16*m*d
        den=(m*m+4*t0*t0)*(m*m+4*frequency*frequency)
        s += num*SCALE//den
    return F(s, SCALE)

def conservative_grid_lower(x: F, places=18):
    base=10**places
    return F((x*base).numerator//(x*base).denominator, base)

def audit():
    anchor = archimedean_anchor_lower()
    require(anchor > F(1, 400), 'unproved positive archimedean anchor')
    results=[]
    for threshold,floor in FLOORS:
        increment = digamma_increment_lower(threshold)
        lower = anchor+increment
        require(lower > floor, f'floor {floor} not certified from {threshold}')
        results.append({
           'threshold_frequency_exact':threshold,
           'certified_constant_floor_exact':str(floor),
           'digamma_increment_strict_lower_grid':str(conservative_grid_lower(increment)),
           'four_prime_symbol_lower_grid':str(conservative_grid_lower(lower)),
           'slack_above_floor_lower_grid':str(conservative_grid_lower(lower-floor)),
           'positive_tail_terms_omitted':True,
        })
    g=euler_gamma_upper()
    return {
      'schema':'NUMERIA.ES8P260.AHLFORS_ARCHIMEDEAN_RATIONAL.1',
      'reference_frequency':REFERENCE_FREQUENCY,
      'log_gamma_analytic_basis':'Euler-Weierstrass series; digamma=Gamma_prime/Gamma',
      'pi_bound_method':'Machin identity and arctan alternating series',
      'log_bound_method':'atanh geometric tail exact rational',
      'anchor_lower_grid':str(conservative_grid_lower(anchor)),
      'gamma_upper_grid':str(-conservative_grid_lower(-g)),
      'prime_powers':[[q,p] for q,p in PRIME_POWERS],
      'new_piecewise_floors':results,
      'background_floor_504_to_850':'1/4 from P253 (not rederived here)',
      'prior_floor_850_to_1100':'1/2 from P258 (independently recertified)',
      'new_model_correction_vs_P259':'1/4*(I-G_1100)+1/4*(I-G_1500)+1/4*(I-G_2000)',
      'source_ahlfors':'Complex Analysis, 3rd ed. chapter 4 residues and chapter 5 Gamma/Stirling/zeta (bibliographic guide, not a copied proof)',
      'uses_arb':False,
      'uses_float_sign_tests':False,
      'sector_odd_377_positive_certified':False,
      'previous_negative_witness_repaired':None,
      'Weil_crossing_established':False,
      'Riemann_hypothesis_proved':False,
    }

def main():
    cmd=argparse.ArgumentParser(description=__doc__)
    cmd.add_argument('--out',type=Path,required=True)
    args=cmd.parse_args()
    data=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('P260: rational floors checked',[(x['threshold_frequency_exact'],x['certified_constant_floor_exact']) for x in data['new_piecewise_floors']])

if __name__ == '__main__':main()
