"""ES8-P gates: rigorous rational envelopes distinguished from numeric diagnostics."""
import pytest
from fractions import Fraction as Q
import sympy as sp
import mpmath as mp
from es8p_projection import (L, constants, atan_bounds, atanh_bounds, prove_pole_projection,
    b_upper, projected_mode_abs, v24_enclosure_entry,projected_v24_error, fmt,
    overlap_closed,sym_overlap,mode,pole_moment,signed_projector,prime_numeric,mode_fourier,cross_parity_phase_diagnostic)


def test_pi_machin_rational_bounds():
    x=constants()
    with mp.workdps(75):
        assert mp.mpf(x['pi_lo'].numerator)/x['pi_lo'].denominator <mp.pi<mp.mpf(x['pi_hi'].numerator)/x['pi_hi'].denominator

def test_log_support_exact():
    c=constants()
    assert c['ln2_hi']<2*L<c['ln3_lo']
    assert c['ln2_lo']<c['ln2_hi'] and c['ln3_lo']<c['ln3_hi']

def test_sqrt_interval_rational_exact():
    c=constants()
    assert c['sqrtL_lo']**2<L<c['sqrtL_hi']**2

def test_invalid_atan_arguments():
    with pytest.raises(ValueError):atan_bounds(Q(2),3)
    with pytest.raises(ValueError):atanh_bounds(Q(1),3)

def test_pole_identically_zero_for_N6_N8_N10():
    for N in (6,8,10):
        d=prove_pole_projection(N)
        assert len(d)==N-2 and all(x['identity_zero_for_both_pole_moments'] for x in d)

def test_pole_cancellation_symbolic_independent():
    X=sp.Symbol('X',positive=True)
    l=sp.Rational(9,20)
    for k in range(3,11):
        p=1 if k%2 else 2
        b=-sp.Rational(k,p)*(p*p*X+l*l)/(k*k*X+l*l)
        assert sp.cancel(k/(k*k*X+l*l)+b*p/(p*p*X+l*l))==0

def test_projector_pivot_parity():
    c=constants(); B=projected_mode_abs(10,c)
    assert len(B)==10 and len(B[0])==8
    for j,k in enumerate(range(3,11)):
        p=1 if k%2 else 2
        assert B[k-1][j]==1 and B[p-1][j]>0
        assert B[2-p][j]==0

def test_projection_coeffs_upper_actual():
    c=constants()
    with mp.workdps(70):
        X=mp.pi**2
        for k in range(3,11):
            p=1 if k%2 else 2
            q=b_upper(k,p,c)
            hi=mp.mpf(q.numerator)/q.denominator
            actual=mp.mpf(k)/p*(p*p*X+mp.mpf(9)**2/400)/(k*k*X+mp.mpf(9)**2/400)
            assert hi>actual

def test_v24_error_symmetry_and_delta_quadratic():
    c=constants()
    x=v24_enclosure_entry(3,9,Q(1,100000),c)
    assert x==v24_enclosure_entry(9,3,Q(1,100000),c)
    y=v24_enclosure_entry(3,9,Q(1,1000000),c)
    assert x>=99*y and x<=101*y

def test_projected_rational_error_monotonic():
    c=constants()
    for N in (6,8,10):
        a=projected_v24_error(N,Q(1,10000),c)
        b=projected_v24_error(N,Q(1,100000),c)
        d=projected_v24_error(N,Q(1,1000000),c)
        assert Q(a['projected_max_row_exact_fraction'])>Q(b['projected_max_row_exact_fraction'])>Q(d['projected_max_row_exact_fraction'])

def test_projection_does_not_imply_positivity():
    # Perturbation bound cannot certify lambda_min without a matching verified matrix margin.
    d=projected_v24_error(10,Q(1,1000000),constants())
    assert Q(d['projected_max_row_exact_fraction'])>0
    assert 'positive_definite' not in d

def test_cross_parity_symmetric_overlap_zero():
    with mp.workdps(60):
        u=mp.log(2)
        for i,j in [(1,2),(1,4),(2,3),(3,4),(9,10)]:
            assert abs(sym_overlap(i,j,u))<mp.mpf('1e-55')

def test_closed_overlap_independent_quadrature():
    with mp.workdps(55):
        u=mp.log(2);l=mp.mpf(9)/20
        for i,j in [(1,1),(2,4),(3,5)]:
            val=mp.quad(lambda x:mode(i,x)*mode(j,x-u),[-l+u,l])
            assert abs(val-overlap_closed(i,j,u))<mp.mpf('1e-48')

def test_pole_moments_numerical():
    with mp.workdps(60):
        B=signed_projector(10)
        for sign in (-1,1):
            v=mp.matrix([pole_moment(k,sign) for k in range(1,11)])
            res=B.T*v
            assert max(abs(x) for x in res)<mp.mpf('1e-55')

def test_prime_matrix_only_numerical():
    v=prime_numeric(6,45)
    assert abs(float(v['sample_prime_entries']['P12_cross_parity']))<1e-35
    assert float(v['projected_prime_max_absolute_entry'])>0
    assert 'NOT interval' in v['numeric_status']

def test_fmt_upper():
    assert fmt(Q(1,3),5)=='0.33334'
    assert fmt(Q(1,2),5)=='0.50000'


def test_cross_parity_archimedean_phase_mismatch():
    with mp.workdps(60):
        f1=mode_fourier(1,mp.mpf(1))
        f2=mode_fourier(2,mp.mpf(1))
        product=f1*mp.conj(f2)
        assert abs(mp.re(product))<mp.mpf('1e-55')
        assert abs(mp.im(product))>mp.mpf('0.1')

def test_phase_gate_reported():
    v=cross_parity_phase_diagnostic(60)
    assert 'POINTWISE_MISMATCH' in v['mathematical_gate']
    assert abs(float(v['unsym_arch_imaginary_integral_0_to_4_numeric']))>0.1
