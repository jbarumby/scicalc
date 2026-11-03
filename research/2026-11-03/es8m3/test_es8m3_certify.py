"""ES8-M3 scientific gates. Analytic inequalities and tests are both auditable.
No test here claims the universal Weil positivity or RH.
"""
import math
from fractions import Fraction as Q
import mpmath as mp
import pytest
import es8m3_certify as m


def test_gamma_em_encloses_euler_constant():
    lo,hi=m.gamma_precise()
    with mp.workdps(65):
        assert mp.mpf(lo)<=mp.euler*m.S<=mp.mpf(hi)
        assert mp.mpf(hi-lo)/m.S < mp.mpf('1e-15')


def test_omega0_matches_digamma():
    lo,hi=m.omega0()
    with mp.workdps(65):
        v=(mp.digamma(mp.mpf('0.25'))-mp.log(mp.pi))*m.S
        assert mp.mpf(lo)<=v<=mp.mpf(hi)


@pytest.mark.parametrize('terms',[32,64,96])
def test_omega_euler_maclaurin_encloses_numerical_values(terms):
    with mp.workdps(65):
        for t in (Q(0),Q(1),Q(4),Q(14),Q(32),Q(48)):
            lo,hi=m.omega_interval_em(t,terms)
            val=mp.re(mp.digamma(mp.mpf('0.25')+mp.j*t.numerator/(2*t.denominator)))-mp.log(mp.pi)
            assert mp.mpf(lo)<=val*m.S<=mp.mpf(hi)


def test_omega_em_width_shrinks():
    for t in (Q(4),Q(16),Q(48)):
        x=m.omega_interval_em(t,32)
        y=m.omega_interval_em(t,64)
        assert y[1]-y[0]<x[1]-x[0]
        assert x[0]<=y[0]<=y[1]<=x[1]


def test_omega_derivative_majorants():
    # Bounds from differentiating the absolutely convergent positive digamma series.
    assert m.omega_derivative_bounds()[1:]==[Q(10),Q(36),Q(208),Q(1632)]


def test_fourth_bound_against_numeric_differentiation():
    with mp.workdps(45):
        aa=mp.mpf(1)/50
        for k in (0,1,3,5):
            for t in (0,4,12,30):
                def f(y):
                    return (y*y/400)**k*mp.exp(-aa*y*y)*(mp.re(mp.digamma(mp.mpf(1)/4+1j*y/2))-mp.log(mp.pi))
                # For complex-step derivatives the real-part extension is not
                # holomorphic. Use real differentiation via mp.diff's step method.
                num=abs(mp.diff(f,mp.mpf(t),4))
                bound=m.fourth_bound(Q(t),k)*m.es7.exp_negative_interval(Q(t*t,50))[1]/m.S
                assert num < mp.mpf(bound.numerator)/bound.denominator if bound else num < mp.mpf('1e-60')


def test_laguerre_polynomial_recurrence_exact():
    pol=m.p_coeffs()
    assert len(pol)==6
    assert pol[0]==[Q(1)]
    assert pol[1]==[Q(-1,2),Q(1)]
    assert pol[2]==[Q(3,4),Q(-3),Q(1)]
    assert pol[5][-1]==1
    assert all(sum(abs(c)*Q(200)**(j-k) for j,c in enumerate(p))<2 for k,p in enumerate(pol))


def test_derivative_factor_by_high_precision_difference():
    with mp.workdps(55):
        a=mp.mpf('0.02'); b=mp.mpf('0.12'); x=b/a
        for k in range(6):
            true=mp.diff(lambda z: z**(-mp.mpf('0.5'))*mp.exp(-b/z), a, k)
            polynomial=sum(mp.mpf(c.numerator)/c.denominator*x**j for j,c in enumerate(m.p_coeffs()[k]))
            candidate=a**(-mp.mpf('0.5')-k)*mp.exp(-x)*polynomial
            assert mp.almosteq(true,candidate,rel_eps=mp.mpf('1e-45'))


def test_prime_tail_and_derivatives():
    r=m.prime_derivatives(64)
    assert r['tail_bound_per_derivative']=='1e-64'
    assert len(r['derivatives'])==6
    for item in r['derivatives']:
        lo,hi=m.from_exact(item['interval_rational'])
        assert lo<=hi


def test_prime_rejects_unjustified_truncation():
    with pytest.raises(ValueError):m.prime_derivatives(63)


def test_simpson_rejects_odd_M():
    with pytest.raises(ValueError):m.arch_simpson(M=1031,terms=16)


def test_paper_psd_gates_full_experiment():
    # Recomputes six full analytic Weil moments, not just known critical zeros.
    d=m.run(M=8192,terms=96)
    assert len(d['moments'])==6
    assert len(d['matrix_certificates'])==2
    for row in d['matrix_certificates']:
        assert row['status']=='CERTIFIED_POSITIVE_DEFINITE'
        assert all(m.from_exact(z['interval_rational'])[0]>0 for z in row['leading_principal_minors'])
    reference=['0.0370825834626389645','0.0187053407106170003','0.0095492225342477344',
               '0.0050033221982893095','0.00276629437148112548','0.00169129215687991098']
    with mp.workdps(55):
        for item, ref in zip(d['moments'],reference):
            lo,hi=m.from_exact(item['interval_rational'])
            assert mp.mpf(lo.numerator)/lo.denominator<mp.mpf(ref)<mp.mpf(hi.numerator)/hi.denominator
    assert 'no RH proof' in d['status']


def test_refinement_is_nested():
    # Independent read-back of two fully computed outputs, rather than an
    # unrelated float-only regression.
    import json
    from pathlib import Path
    p=Path(__file__).resolve().parent
    a=json.loads((p/'RESULTADOS_ES8M3.json').read_text())
    b=json.loads((p/'RESULTADOS_ES8M3_REFINED.json').read_text())
    assert a['parameters']['M']==8192 and b['parameters']['M']==16384
    for x,y in zip(a['moments'],b['moments']):
        xl,xu=m.from_exact(x['interval_rational'])
        yl,yu=m.from_exact(y['interval_rational'])
        assert xl<=yl<=yu<=xu
    for mat in b['matrix_certificates']:
        assert mat['status']=='CERTIFIED_POSITIVE_DEFINITE'


def test_determinant_adversarial_negative_gate():
    # Negative 2x2 principal minor must never be recorded as positive.
    fake=[m.interval(i) for i in (1,2,1,4,1,4)]
    _,rs=m.moments_and_matrices({'moments':[{'interval_rational':m.exact(z)} for z in fake]},
                                {'derivatives':[{'interval_rational':m.exact(m.interval(0))} for _ in range(6)]})
    assert rs[0]['status']=='CERTIFIED_NOT_PSD'


def test_no_claim_of_equivalence_to_v24():
    import json
    from pathlib import Path
    p=Path(__file__).resolve().parent
    r=json.loads((p/'RESULTADOS_ES8M3_REFINED.json').read_text())
    assert 'No direct E7/V24 matrix-entry certification' in r['limitations'][1]
    assert 'No proof of Weil positivity for all test functions' in r['limitations'][0]
