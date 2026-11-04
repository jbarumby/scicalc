import json
from fractions import Fraction as Q
from pathlib import Path
import mpmath as mp
import pytest
from es8p7_independent import (
    pi, omega0, omega, amplitude, qij, I, sin_interval,
    sinc_interval, grid, comp, check_admissibility, reference_triage,
    hashfile, width,
)

ROOT=Path(__file__).parent


def ininterval(z,iv):
    return mp.mpf(iv[0].numerator)/iv[0].denominator <= z <= mp.mpf(iv[1].numerator)/iv[1].denominator


def test_pi_interval():
    mp.mp.dps=60
    assert ininterval(mp.pi,pi())


def test_omega_origin():
    mp.mp.dps=60
    z=mp.digamma(mp.mpf('0.25'))-mp.log(mp.pi)
    assert ininterval(z,omega0())


@pytest.mark.parametrize('r',[Q(0),Q(1),Q(4),Q(10),Q(160)])
def test_digamma_nodes(r):
    mp.mp.dps=60
    exact=mp.re(mp.digamma(mp.mpf('0.25')+mp.j*mp.mpf(r.numerator)/(2*r.denominator)))-mp.log(mp.pi)
    assert ininterval(exact,omega(r,256))


@pytest.mark.parametrize('x',[Q(1,7),Q(3,2),Q(27,1),Q(-41,3)])
def test_sinc_positive_negative_and_reduction(x):
    mp.mp.dps=60
    val=mp.sin(mp.mpf(x.numerator)/x.denominator)/(mp.mpf(x.numerator)/x.denominator)
    assert ininterval(val,sinc_interval(I(x)))


@pytest.mark.parametrize('k,r',[(1,Q(0)),(1,Q(4)),(2,Q(4)),(3,Q(10)),(10,Q(160))])
def test_amplitude_against_direct_fourier(k,r):
    mp.mp.dps=65
    L=mp.mpf(9)/20
    w=k*mp.pi/(2*L)
    r_mp=mp.mpf(r.numerator)/r.denominator
    f=lambda t:mp.sin(w*(t+L))/mp.sqrt(L)*mp.exp(1j*r_mp*t)
    val=mp.quad(f,[-L,L])
    a=amplitude(r,k)
    # The Fourier transform can be purely imaginary for even k. The
    # magnitude is |A_k| for either parity; unlike the previous weak test,
    # this checks actual containment with high precision.
    al=mp.mpf(a[0].numerator)/a[0].denominator
    ah=mp.mpf(a[1].numerator)/a[1].denominator
    lower=mp.mpf(0) if al<=0<=ah else min(abs(al),abs(ah))
    upper=max(abs(al),abs(ah))
    assert lower <= abs(val) <= upper


def test_cross_parity_real_entry_exactly_zero():
    for i,j in [(1,2),(2,3),(5,6),(8,9)]:
        assert qij(Q(4),i,j)==I(0)


def test_grid_exact_nodes_and_coverage():
    ss=grid()
    assert len(ss)==15
    assert sorted(set(x[0] for x in ss))==[0,40,100,400,1599]
    assert all(r==Q(c,10)+Q(j,150) for c,j,r in ss)
    assert all(0<=r<=160 for _,_,r in ss)


def test_intersection_negative_control():
    assert comp(I(-1,1),I(0,2))
    assert not comp(I(-2,-1),I(1,2))


def test_admissibility_scope():
    d=check_admissibility()
    assert Q(d['strip_half_height'])>Q(1,2)
    assert 'O((1+|x|)^-4)' in d['uniform_decay']
    assert 'Individual H_ij' in d['scope']


def test_references_canonical():
    rs=reference_triage(ROOT/'REFERENCIAS_ORIGINAIS_16.csv')
    assert len(rs)==8  # 7 user records, plus one external reference
    assert sum(x['source']=='user_supplied_csv_metadata' for x in rs)==7
    assert any(x['doi']=='10.1016/j.amc.2006.08.061' for x in rs)


def test_historical_fixed_sources_pinned():
    assert hashfile(ROOT/'fase4_fixedpoint_sem_arb.py')=='0935bb988f7a7bd4c83b36a4a873a71b1612eda00c77ac9bd6ba7daeef600bda'
    assert hashfile(ROOT/'fase4_nos_sem_arb.py')=='151ec9a0a6c249f443bcbac278e8f004e02d268ef7590f81e78e5280872cbe49'


def test_original_phase4_metadata():
    d=json.loads((ROOT/'RESUMO_FASE4_1600CELULAS_D40.json').read_text())
    assert d['num_pairs']==30 and d['cells']==1600 and d['digits']==40
    assert len(d['records'])==30


def test_result_file_if_present():
    p=ROOT/'RESULTADOS_ES8P7.json'
    if not p.exists():
        pytest.skip('Re-run computation to create results')
    d=json.loads(p.read_text())
    assert d['statistics']['compared']==450
    assert d['statistics']['incompatible']==0
    assert d['status']=='450_SAMPLE_NODAL_INTERSECTIONS_VERIFIED'
    assert len(d['independent_digamma_mp65_not_certified'])==15
    assert all(r['contained_in_rational_interval'] for r in d['independent_digamma_mp65_not_certified'])
