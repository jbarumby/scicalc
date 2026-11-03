"""NUMERIA ES8-M regression and adversarial tests. Tests of explicit formula are diagnostic unless marked rational exact."""
import json
from fractions import Fraction as Q
from pathlib import Path
import mpmath as mp
import pytest

from es6_weil_gaussian_baseline import (
    SCALE, arctan_alt, pi_bounds, exp_neg_bounds, log_int_bounds,
)
from es8m_moment_audit import (
    prime_derivative_cert, full_moments_0_1_certified,
    moments_of_finite_rational_positive_atoms, hankel, exact_ldl,
    toy_benchmarks, laguerre_minus_half, interval_from_scaled,
    run,
)
ROOT=Path(__file__).resolve().parent
BASE=json.loads((ROOT/'ES7_RESULTADOS_BASELINE.json').read_text())


def test_arctan_5_rational_bound_corrected():
    lo,hi=arctan_alt(Q(1,5),44)
    assert lo<hi
    with mp.workdps(110):
        v=mp.atan(mp.mpf(1)/5)
        assert mp.mpf(lo.numerator)/lo.denominator<v<mp.mpf(hi.numerator)/hi.denominator


def test_arctan_239_requires_sufficient_independent_precision():
    lo,hi=arctan_alt(Q(1,239),44)
    with mp.workdps(260):
        v=mp.atan(mp.mpf(1)/239)
        assert mp.mpf(lo.numerator)/lo.denominator<v<mp.mpf(hi.numerator)/hi.denominator


def test_pi_exact_scaled_enclosure():
    lo,hi=pi_bounds()
    assert hi-lo<=2
    with mp.workdps(90):
        assert mp.mpf(lo)/SCALE < mp.pi < mp.mpf(hi)/SCALE


def test_prime_derivative_enclosures_overlap_at_two_cutoffs():
    c64=prime_derivative_cert(64)
    c80=prime_derivative_cert(80)
    l64,u64=map(int,c64['scaled_exact'])
    l80,u80=map(int,c80['scaled_exact'])
    assert 0<l64<=u64 and 0<l80<=u80
    assert max(l64,l80)<=min(u64,u80)
    assert max(u64-l64,u80-l80)<5000000


def test_prime_derivative_matches_independent_high_precision_diagnostic():
    p=prime_derivative_cert(64)
    lo,hi=map(int,p['scaled_exact'])
    with mp.workdps(80):
        a=mp.mpf(1)/50
        def fn(x):
            from es6_weil_gaussian_baseline import prime_powers_mpmath
            return prime_powers_mpmath(256,x)[0]
        v=mp.diff(fn,a)
        assert mp.mpf(lo)/SCALE < v < mp.mpf(hi)/SCALE


def test_mu0_mu1_positive_rational_intervals():
    p=full_moments_0_1_certified(BASE)
    for k in ('mu0','mu1'):
        a=p[k]['rational_interval']
        lo,hi=(Q(int(x),int(y)) for x,y in a)
        assert 0<lo<hi
        assert p[k]['certified_positive']


def test_mu0_mu1_recheck_with_independent_numeric_expression():
    p=full_moments_0_1_certified(BASE)
    with mp.workdps(60):
        from es6_weil_gaussian_baseline import analytic_parts
        x=analytic_parts(mp.mpf(1)/50,N=128,T=96,derivative=True)
        values={'mu0':x['total'],'mu1':-x['d_total']}
        for k,v in values.items():
            a=p[k]['rational_interval'];lo,hi=[Q(int(x),int(y)) for x,y in a]
            assert mp.mpf(lo.numerator)/lo.denominator < v < mp.mpf(hi.numerator)/hi.denominator


def test_classical_rational_stieltjes_fixture_both_families_pd():
    d=toy_benchmarks()['positive_atoms']['matrices']
    for n in ('0','1','2','3'):
        assert d[n]['H']['state']=='CERTIFIED_PD_RATIONAL_FIXTURE'
        assert d[n]['H_shift']['state']=='CERTIFIED_PD_RATIONAL_FIXTURE'


def test_offaxis_fixture_has_negative_ldl_pivot_exact():
    fixture=toy_benchmarks()['offaxis_fixture']
    assert fixture['H1']['pivots']==['2','-2']
    assert fixture['H1_shift']['pivots']==['2','-4']
    assert fixture['H1']['state']=='CERTIFIED_NOT_PSD'


def test_zero_pivot_cannot_be_called_positive_definite():
    assert exact_ldl([[Q(1),Q(1)],[Q(1),Q(1)]])['state']=='INCONCLUSIVE_ZERO_PIVOT'


def test_stieltjes_shift_is_not_derivative_same_sign():
    # m'_k(a)=-m_{k+1}(a); sign matters for H'_n = -H_shift_n.
    nodes=[Q(1),Q(3)];weights=[Q(2),Q(5)]
    moments=moments_of_finite_rational_positive_atoms(nodes,weights,3)
    assert hankel(moments,1,1)==[[moments[1],moments[2]],[moments[2],moments[3]]]


def test_laguerre_derivative_identity_at_orders_zero_one_two():
    with mp.workdps(60):
        x=mp.mpf('1.7')
        assert laguerre_minus_half(0,x)==1
        assert laguerre_minus_half(1,x)==mp.mpf('.5')-x
        assert mp.almosteq(2*laguerre_minus_half(2,x), x*x-3*x+mp.mpf('.75'))


def test_bad_derivative_certificate_arguments():
    with pytest.raises(ValueError):prime_derivative_cert(32)


def test_conservative_scientific_gates():
    result=run(zeros=12,degree=1,dps=35)
    assert result['gates']['Weil_universal_positivity']=='NOT_PROVED'
    assert result['gates']['H_n_actual_zeta_for_n_ge_1']=='NOT_CERTIFIED'
    assert 'RH proof' in result['scientific_status']
