import hashlib, json, math
from fractions import Fraction as Q
from pathlib import Path
import pytest
from es8p8_full import IndependentFixed, FixedRationalNodes, PAIRS, HIST_SHA, full_grid, scan

DIR=Path(__file__).resolve().parent

def overlap(a,bs,b,as_):
    return a[0]*bs <= b[1]*as_ and b[0]*as_ <= a[1]*bs

def test_grid_size_and_pair_count():
    assert len(PAIRS)==30
    assert len(set(PAIRS))==30
    assert all(1<=i<=j<=10 and (i-j)%2==0 for i,j in PAIRS)
    assert sum(1 for _ in full_grid())==25600
    assert Q(1599,10)+Q(15,150)==Q(160)

def test_historical_sha_pinning():
    for name,expect in HIST_SHA.items():
        assert hashlib.sha256((DIR/name).read_bytes()).hexdigest()==expect

def test_pi_and_omega0_independent():
    x=IndependentFixed()
    import mpmath as mp
    mp.mp.dps=48
    pi=mp.pi
    assert mp.mpf(x.pi[0])/x.S<pi<mp.mpf(x.pi[1])/x.S
    y=mp.digamma(mp.mpf('0.25'))-mp.log(mp.pi)
    assert mp.mpf(x.omega_at_zero[0])/x.S < y < mp.mpf(x.omega_at_zero[1])/x.S

def test_sinc_at_zero_and_sine_values():
    import mpmath as mp
    mp.mp.dps=55
    x=IndependentFixed()
    for v in (Q(0),Q(1,100),Q(1,2),Q(1),Q(3),Q(9),Q(19),Q(-43,11)):
        r=x.num(v);n=x.sinc(r);s=x.sin(r)
        z=mp.mpf(v.numerator)/v.denominator
        a=mp.sin(z)/z if z else mp.mpf(1)
        b=mp.sin(z)
        n0=mp.mpf(n[0])/x.S;n1=mp.mpf(n[1])/x.S
        s0=mp.mpf(s[0])/x.S;s1=mp.mpf(s[1])/x.S
        assert n0<=a<=n1
        assert s0<=b<=s1

def test_omega_vs_high_precision_diagnostic():
    import mpmath as mp
    mp.mp.dps=50
    e=IndependentFixed()
    for r in (Q(0),Q(1,150),Q(1),Q(7,2),Q(20),Q(50),Q(160)):
        lo,hi=e.omega(r)
        rr=mp.mpf(r.numerator)/r.denominator
        t=mp.re(mp.digamma(mp.mpf(1)/4+1j*rr/2))-mp.log(mp.pi)
        assert mp.mpf(lo)/e.S <= t <= mp.mpf(hi)/e.S

def test_historical_and_independent_intersect_various_points():
    a=FixedRationalNodes(digits=40,N=64)
    b=IndependentFixed()
    for r in (Q(0),Q(1,150),Q(17,3),Q(20),Q(400,3),Q(160)):
        for i,j in PAIRS:
            assert overlap(a.qij(r,i,j),b.S,b.qij(r,i,j),a.S),(r,i,j)

def test_negative_mutations_detected():
    a=(1,2);b=(5,7);S=10
    assert not overlap(a,S,b,S)
    assert overlap(a,S,(2,3),S)
    assert not overlap((10,12),S,(-12,-10),S)

def test_partial_scan_labeled_partial(tmp_path):
    v=scan(0,1,chunk=1)
    assert v['status']=='PARTIAL_OR_MISMATCH'
    assert v['count_compared']==480
    assert v['count_positions']==16
    assert v['failures_total']==0
    assert len(v['chunks'])==1
    assert len(v['chunks'][0]['sha256_pairwise_intervals'])==64

def test_full_artifact_when_present():
    p=DIR/'RESULTADOS_ES8P8.json'
    if not p.exists():pytest.skip('full result not yet run')
    x=json.loads(p.read_text())
    assert x['scope']['full_expected_comparisons']==768000
    assert x['count_positions']==25600
    assert x['count_compared']==768000
    assert x['failures_total']==0
    assert x['status']=='FULL_GRID_PAIRWISE_INTERVAL_INTERSECTIONS_CHECKED'
    assert sum(c['comparisons'] for c in x['chunks'])==768000
    assert all(c['incompatible']==0 and len(c['sha256_pairwise_intervals'])==64 for c in x['chunks'])
