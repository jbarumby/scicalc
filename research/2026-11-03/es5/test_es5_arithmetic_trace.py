"""Regression gates for ES5: strict separation of certified and exploratory claims."""
import math

import mpmath as mp
import pytest

from es5_arithmetic_trace import (
    SCALE, ceildiv, compact_support_check, factorization,
    gaussian_trace_exploration, log_int_bounds, prime_power_base,
    rational_log_bounds, series_intervals, spf_sieve, completed_log_derivative
)


def test_ceildiv_exact_signs():
    assert ceildiv(10, 3) == 4
    assert ceildiv(-10, 3) == -3
    assert ceildiv(0, 7) == 0


def test_log_rational_remainder_nonnegative():
    for numerator,denominator in [(1,1),(3,2),(2,1)]:
        from fractions import Fraction
        lo,hi=rational_log_bounds(Fraction(numerator,denominator))
        assert 0 <= lo <= hi
        assert hi-lo < 1 / (10**40)


def test_log_integer_enclosures_against_independent_high_precision():
    with mp.workdps(80):
        for n in (1,2,3,7,31,997,4096):
            lo,hi=log_int_bounds(n)
            x=mp.log(n)*SCALE
            assert mp.mpf(lo) <= x <= mp.mpf(hi)


def test_prime_power_n_or_not():
    sieve=spf_sieve(120)
    assert [prime_power_base(n,sieve) for n in (1,2,4,8,9,16,18,27,36,49,64,81,97)] == \
           [None,2,2,2,3,2,None,3,None,7,2,3,97]
    assert factorization(72,sieve)=={2:3,3:2}


def test_compact_support_log_two_only():
    out=compact_support_check()
    assert out['L']=='9/20'
    assert 'only prime power n=2' in out['certified']


def test_independent_certificates_overlap_for_sigma_four():
    x=series_intervals(4,512)
    assert x['overlap_certified']
    assert x['prime_powers_count_through_N']>0
    assert int(x['dirichlet_interval_scaled_exact'][0]) < int(x['dirichlet_interval_scaled_exact'][1])
    assert int(x['primes_interval_scaled_exact'][0]) < int(x['primes_interval_scaled_exact'][1])


def test_certificate_enclosures_nest_at_increased_N():
    short=series_intervals(4,128)
    long=series_intervals(4,512)
    assert int(long['primes_interval_scaled_exact'][0]) >= int(short['primes_interval_scaled_exact'][0])
    assert int(long['primes_interval_scaled_exact'][1]) <= int(short['primes_interval_scaled_exact'][1])


def test_independent_zeta_value_within_exact_interval():
    v=series_intervals(4,1024)
    with mp.workdps(80):
        target=-mp.zeta(4,derivative=1)/mp.zeta(4)
        for key in ('primes_interval_scaled_exact','dirichlet_interval_scaled_exact'):
            low,hi=[mp.mpf(z)/SCALE for z in v[key]]
            assert low < target < hi


def test_rigorous_prime_tail_positive_and_shrinking():
    a=series_intervals(4,128)
    b=series_intervals(4,1024)
    assert mp.mpf(a['prime_tail_bound'])>mp.mpf(b['prime_tail_bound'])>0


def test_gaussian_trace_numerical_exploration_only():
    r=gaussian_trace_exploration(sigma=4,N=256,dps=35)
    with mp.workdps(40):
        difference=mp.mpf(r['numerical_abs_difference'])
        tails=(mp.mpf(r['rigorous_analytic_bound_for_missing_prime_tail'])
               +mp.mpf(r['rigorous_analytic_bound_for_missing_t_tail']))
        # Diagnostic uses nonrigorous quadrature, with conservative test tolerance.
        assert difference < tails + mp.mpf('1e-16')
        assert 'NOT interval certified' in r['note']


def test_functional_equation_derivative_numeric():
    with mp.workdps(50):
        s=mp.mpc('2.25','0.75')
        residual=completed_log_derivative(s)+completed_log_derivative(1-s)
        assert abs(residual)<mp.mpf('1e-45')


def test_bad_inputs():
    with pytest.raises(ValueError): series_intervals(1,128)
    with pytest.raises(ValueError): series_intervals(4,2)
    with pytest.raises(ValueError): gaussian_trace_exploration(3,128)
