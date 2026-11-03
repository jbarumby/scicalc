"""ES8-M2: rigorous arithmetic checks, adversarial and integration tests.

Fast tests use the committed complete-run JSON. The Github workflow additionally
recomputes that JSON with M=16384 and checks signs and containment.
"""
import json
from fractions import Fraction as F
from pathlib import Path
import pytest
import mpmath as mp
import es8m2_certify as mod

ROOT=Path(__file__).parent

def expected():
    return json.loads((ROOT/'RESULTADOS_ES8M2.json').read_text(encoding='utf8'))

def iq(raw):return [F(int(n),int(d)) for n,d in raw]

def test_exponential_kernel_constants_are_correct():
    assert mod.A==F(1,50)
    assert mod.C==400
    assert mod.S==10**28 and mod.S2==10**40


def test_laguerre_moment_factors_at_exact_arguments():
    x=mod.interval(F(1))
    assert [mod.laguerre_derivative_factor(k,x)[0] for k in range(4)] == [F(1),F(1,2),F(-5,4),F(23,8)]


def test_laguerre_differentiation_independent_mpmath():
    with mp.workdps(55):
        a=mp.mpf(1)/50
        u=mp.log(2)
        x=u*u/(4*a)
        fun=lambda z:z**(-mp.mpf('0.5'))*mp.exp(-u*u/(4*z))
        for k in range(4):
            numerical=mp.diff(fun,a,k)/fun(a)
            expected=(a**(-k))*sum([
                mp.mpf(1), x-mp.mpf('0.5'), x*x-3*x+mp.mpf('0.75'),
                x**3-mp.mpf('7.5')*x*x+mp.mpf('11.25')*x-mp.mpf('1.875')][k:k+1])
            assert abs(numerical-expected)<mp.mpf('1e-42')


def test_tail_bound_entire_prime_derivatives():
    bound=F(5*625**3*4,91)*F(10,27)**182
    assert 0<bound<F(1,10**70)
    assert mod.laguerre_derivative_factor(3,mod.interval(200))[0]>0


def test_uniform_prime_moment_certificate_positive():
    obj=expected()['prime_derivatives']
    assert obj['N']>=64
    assert len(obj['P_derivatives_rational'])==4
    for p in obj['P_derivatives_rational']:
        lo,hi=iq(p)
        assert lo>0 and hi>=lo


def test_rational_display_is_outward():
    for a in expected()['moments_normalized_by_400_power']:
        lo,hi=iq(a['interval_rational'])
        shown=[F(v) for v in a['display']]
        assert shown[0]<=lo<=hi<=shown[1]


def test_all_four_full_moments_positive():
    for x in expected()['moments_normalized_by_400_power']:
        lo,hi=iq(x['interval_rational'])
        assert 0<lo<=hi


def test_two_complete_hankel_matrices_are_certified_positive():
    obj=expected()
    assert len(obj['hankel_and_stieltjes_2x2'])==2
    for matrix in obj['hankel_and_stieltjes_2x2']:
        assert matrix['status']=='CERTIFIED_POSITIVE_DEFINITE'
        lo,hi=iq(matrix['det_rational'])
        assert lo>0 and hi>lo


def test_recompute_determinants_independently_from_rational_moments():
    vals=[tuple(iq(x['interval_rational'])) for x in expected()['moments_normalized_by_400_power']]
    for i,d in enumerate(expected()['hankel_and_stieltjes_2x2']):
        # Independent interval product, exact fractions: positive moment endpoints.
        low=vals[i][0]*vals[i+2][0]-vals[i+1][1]**2
        high=vals[i][1]*vals[i+2][1]-vals[i+1][0]**2
        assert iq(d['det_rational'])==[low,high]
        assert low>0


def test_diagnostic_matches_enclosures_but_is_not_a_certificate():
    out=expected()
    assert 'floating_point_crosscheck_NOT_CERTIFIED' in out
    for i, val in enumerate(out['floating_point_crosscheck_NOT_CERTIFIED']):
        lo,hi=iq(out['moments_normalized_by_400_power'][i]['interval_rational'])
        assert lo<F(val)<hi


def test_input_csv_bibliography_deduplicated():
    r=mod.classify_refs(ROOT/'sources')
    keys=[x['DOI'].strip().lower() or x['Title'].lower() for x in r]
    assert len(keys)==len(set(keys))==19
    assert any('Connes' in x['Title'] or 'Trace formula in noncommutative geometry' in x['Title'] for x in r)
    assert any('Li coefficients' in x['Title'] for x in r)


def test_invalid_certification_parameters_rejected():
    with pytest.raises(ValueError):mod.prime_moment_intervals(63)
    with pytest.raises(ValueError):mod.arch_moments(M=16)
    with pytest.raises(ValueError):mod.arch_moments(terms=4)
    with pytest.raises(ValueError):mod.laguerre_derivative_factor(4,mod.interval(1))
    with pytest.raises(ValueError):mod.upper_gaussian_integral_coeff(3,mod.A,48)


def test_gaussian_tail_coeff_recurrence():
    v0=mod.upper_gaussian_integral_coeff(0,mod.A,48)
    v2=mod.upper_gaussian_integral_coeff(2,mod.A,48)
    assert v0==F(1,2*mod.A*48)
    assert v2==F(48,2*mod.A)+v0/(2*mod.A)


def test_scale_invariant_of_results():
    r=expected()
    assert r['variable']=='x=z_rho^2/400'
    assert r['a']=='1/50'
    assert r['archimedean']['M']==16384
    assert 'NOT RH proof' in r['status']
