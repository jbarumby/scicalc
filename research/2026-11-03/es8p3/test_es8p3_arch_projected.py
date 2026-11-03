"""Tests for interval operations, analytic identity, and rigorous finite gates.

mpmath comparisons are diagnostics only; they are NOT the origin of
any certified interval endpoint. Scientific status remains conditional.
"""
import math
from fractions import Fraction as Q

import mpmath as mp
import pytest

import es8p3_arch_projected as e
from es8p2_prime_cert import I, constants, projected_prime, projection_column


def in_interval_mpmath(x, iv):
    def mpq(q): return mp.mpf(q.numerator)/q.denominator
    return mpq(iv.lo) <= x <= mpq(iv.hi)


def test_pi_log_pi_gamma_omega0_independent_precision():
    with mp.workdps(90):
        assert in_interval_mpmath(mp.pi,constants()['pi'])
        assert in_interval_mpmath(mp.log(mp.pi),e.ln_pi())
        assert in_interval_mpmath(mp.euler,e.gamma_rational())
        assert in_interval_mpmath(mp.re(mp.digamma(mp.mpf(1)/4))-mp.log(mp.pi), e.omega_zero())


def test_hurwitz_euler_maclaurin_encloses_11_exponents():
    with mp.workdps(90):
        for s in range(3,16):
            for M in (32,128):
                q=Q(M)+Q(1,4)
                val=mp.zeta(s,mp.mpf(q.numerator)/q.denominator)
                assert in_interval_mpmath(val,e.zeta_power_enclosure(s,q))


def test_exp_rational_bounds_and_conservative_far_tail():
    with mp.workdps(90):
        for x in (Q(0),Q(9,20),Q(9,4),Q(40),Q(101,2)):
            assert in_interval_mpmath(mp.exp(-mp.mpf(x.numerator)/x.denominator),e.exp_minus_interval(x))
    assert e.exp_minus_interval(Q(50)).lo == 0


def test_orthogonal_opposite_parities_are_exact_zero():
    for i in (1,3,5,7,9):
        for j in (2,4,6,8,10):
            assert e.arch_raw(i,j,64)==I(0)


def test_exponential_kernel_formula_independent_one_dimensional_quadrature():
    with mp.workdps(50):
        L=mp.mpf(9)/20;D=2*L
        def kernel_overlap(i,j,q):
            a=2*q; wi=i*mp.pi/D;wj=j*mp.pi/D
            def Ti(x):
                return (2*a*mp.sin(wj*x)+wj*mp.exp(-a*x)-(-1)**j*wj*mp.exp(-a*(D-x)))/(a*a+wj*wj)
            return 2/D*mp.quad(lambda x:mp.sin(wi*x)*Ti(x),[0,D])
        def closed(i,j,q):
            a=2*q;wi=i*mp.pi/D;wj=j*mp.pi/D
            val=2*a/(a*a+wi*wi) if i==j else 0
            if (i-j)%2==0:
                val+=2/L*wi*wj*(1-(-1)**i*mp.exp(-a*D))/((a*a+wi*wi)*(a*a+wj*wj))
            return val
        for i,j,q in [(1,1,mp.mpf('0.25')),(1,3,mp.mpf('1.25')),(2,4,mp.mpf('0.25')),(1,2,mp.mpf('0.25'))]:
            assert abs(kernel_overlap(i,j,q)-closed(i,j,q))<mp.mpf('1e-45')


def test_direct_fourier_spectral_arch_diag_numeric_crosscheck():
    """Numerically compare fully independent digamma quadrature with the
    certified series; only a diagnostic, quadrature not interval-certified."""
    with mp.workdps(50):
        L=mp.mpf(9)/20
        def F(k,r):
            w=k*mp.pi/(2*L)
            return w/mp.sqrt(L)*(mp.exp(-1j*r*L)-(-1)**k*mp.exp(1j*r*L))/(w*w-r*r)
        def g(r):
            omega=mp.re(mp.digamma(mp.mpf(1)/4+1j*r/2))-mp.log(mp.pi)
            return omega*abs(F(1,r))**2/mp.pi
        # First resonance at r=3.4906; integrate to 400, so missing Fourier
        # tail <= explicit V24-style constant; diagnostic threshold 2e-4.
        num=mp.quad(g,[0,1,2,3,4,6,10,20,50,100,200,400])
        iv=e.arch_raw(1,1,128)
        mid=(mp.mpf(iv.lo.numerator)/iv.lo.denominator+mp.mpf(iv.hi.numerator)/iv.hi.denominator)/2
        assert abs(num-mid)<mp.mpf('2e-4')


def test_arch_refinement_nests_and_is_not_prematurely_truncated():
    for i,j in [(1,1),(3,3),(3,5),(8,10),(10,10)]:
        a=e.arch_raw(i,j,64);b=e.arch_raw(i,j,128)
        assert max(a.lo,b.lo)<=min(a.hi,b.hi)
        assert b.width()<=a.width()


def test_prime_projection_inherited_exact_signs_and_parity():
    for N in (6,8,10):
        assert projected_prime(3,4,N)==I(0)
        assert projected_prime(3,3,N).hi<0
        assert projected_prime(4,4,N).lo>0


def test_pole_projection_coefficients_cancel_both_moments():
    with mp.workdps(75):
        L=mp.mpf(9)/20
        def moment(k,sgn):
            w=k*mp.pi/(2*L)
            return w/mp.sqrt(L)*(mp.exp(-sgn*L/2)-(-1)**k*mp.exp(sgn*L/2))/(w*w+mp.mpf('0.25'))
        for k in range(3,11):
            c=projection_column(k,10)
            p=1 if k%2 else 2
            cp=(c[p].lo+c[p].hi)/2
            cpm=mp.mpf(cp.numerator)/cp.denominator
            for sgn in (-1,1):
                assert abs(moment(k,sgn)+cpm*moment(p,sgn))<mp.mpf('1e-37')


def test_all_three_dimensions_shifted_ldl():
    for N in (6,8,10):
        mat=[[e.full_projected(i,j,N,128) for j in range(3,N+1)] for i in range(3,N+1)]
        result=e.ldl_certify(mat,Q(1,20))
        assert result['status']=='CERTIFIED_POSITIVE_DEFINITE'
        assert all(v[0] != '0.000000000000000000' for v in result['pivots'])


def test_strict_negative_control_ldl_rejects_overshift():
    N=6
    mat=[[e.full_projected(i,j,N,128) for j in range(3,N+1)] for i in range(3,N+1)]
    assert e.ldl_certify(mat,Q(2))['status']=='INCONCLUSIVE_OR_NONPOSITIVE'


def test_v24_taper_error_below_certified_margin():
    for N in (6,8,10):
        assert 0 < e.taper_projected_row_bound(N,Q(1,1000000)) < Q(1,100)
    assert e.taper_projected_row_bound(10,Q(1,1000000)) < Q(1,1000)


def test_taper_quadratic_scale_and_positive_shape():
    N=10
    bound1=e.taper_projected_row_bound(N,Q(1,1000000))
    bound2=e.taper_projected_row_bound(N,Q(1,100000))
    assert bound1*100==bound2


def test_full_replay_scientific_gates():
    data=e.run(128)
    assert data['schema']=='NUMERIA.ES8P3.ARCH_PROJECTED_INTERVAL_LDL.1'
    assert data['status'].startswith('CONDITIONAL')
    assert all(x['conditional_tapered_positive'] for x in data['finite_matrices'].values())
    assert all(x['cross_parity_exact_zero'] for x in data['finite_matrices'].values())
    assert any('historical' in z.lower() for z in data['not_proven'])
    assert 'RH' in data['status']


def test_reject_invalid_parameters():
    with pytest.raises(ValueError):e.run(16)
    with pytest.raises(ValueError):e.arch_raw(11,11)
    with pytest.raises(ValueError):e.arch_tail(1,1,16)
    with pytest.raises(ValueError):e.perturb_bound_raw(1,1,Q(0))


def test_projection_numeric_eigenvalue_diagnostic_only():
    import numpy as np
    with mp.workdps(40):
        for N in (6,8,10):
            mat=np.array([[float((e.full_projected(i,j,N,128).lo+e.full_projected(i,j,N,128).hi)/2)
                           for j in range(3,N+1)] for i in range(3,N+1)])
            assert np.min(np.linalg.eigvalsh(mat))>.09

