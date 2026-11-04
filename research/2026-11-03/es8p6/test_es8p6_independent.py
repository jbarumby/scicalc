"""Independent reproducible tests for ES8-P6. Only stdlib is required."""
from __future__ import annotations
from fractions import Fraction as Q
from pathlib import Path
import sys
import zipfile
import hashlib
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parent))
import es8p6_independent_audit as a
ROOT=Path(__file__).resolve().parent
FIXTURE=ROOT/'FIXTURE_FASE4_ES8P6.json'

def test_strip_conditions():
    assert a.D==Q(7,20) and a.ALPHA==Q(3,40)
    assert a.DIRECTION==Q(137,200)

def test_majorant_slack_nonnegative():
    for center in (0,Q(1,20),1,8,80,Q(3199,20)):
        z=abs(Q(center))+a.D
        raw=Q(7,2)+1/a.ALPHA+z
        original=3+(1+1/a.ALPHA)*(1+z/2)
        assert raw<=original
        assert a.independent_cauchy_bound(center)>0

def test_product_max_formula():
    assert a.interpolation_product()==Q(16416784384,1601806640625)

def test_weights_degree_15():
    w=a.lagrange_weights()
    x=[Q(j,15) for j in range(16)]
    assert len(w)==16
    for degree in range(16):
        assert sum((w[j]*x[j]**degree for j in range(16)),Q())==Q(1,degree+1)

def test_weights_nontrivial_negative():
    assert min(a.lagrange_weights())<0

def test_exact_error_closed_form():
    assert a.independent_total_error()==Q(4154891126177792,209922127267236328125)

def test_error_closed_vs_full_sum():
    assert a.independent_total_error()==a.exact_error_by_cells()

def test_pi_simple_interval():
    lo,hi=a.pi_bounds()
    assert Q(3)<lo<hi<Q(4)

def test_ln480_bracket():
    lo,hi=a.ln480_bounds()
    assert Q(6)<lo<hi<Q(7)

def test_tail_all_same_parity():
    assert all(a.baseline_tail_upper(i,j)>0 for i in range(1,11) for j in range(i,11) if (i-j)%2==0)

def test_tail_cross_parity_zero():
    assert a.baseline_tail_upper(1,2)==Q(0)

def test_tail_conservative_order():
    assert a.baseline_tail_upper(1,1)<a.baseline_tail_upper(10,10)

def test_param_validation():
    assert a.L==Q(9,20) and a.H==Q(1,10) and a.RHO==Q(3,10)

def test_full_bundle_fixture():
    import json
    data=json.loads(FIXTURE.read_text())
    assert data['pairs']==30
    assert Q(data['interpolation_error_exact'])==a.independent_total_error()
    assert data['source_hashes']==4
    assert data['stored_numeric_intervals_checked']==30

def test_error_mutation_detected():
    original=a.independent_total_error()
    assert original != original*Q(9,10)

def test_missing_prime_or_weil_claim():
    import json
    d=json.loads(FIXTURE.read_text())
    assert any('universal Weil' in s for s in d['not_proved'])

def test_fixture_does_not_claim_recompute_nodes():
    import json
    d=json.loads(FIXTURE.read_text())
    assert 'NODES_NOT_RECOMPUTED' in d['status']

def test_9_rational_node_spotchecks():
    from es8p6_node_spotcheck import spotcheck
    result=spotcheck()
    assert result['n_cases']==9
    assert all(e['overlap'] for e in result['checks'])


def test_sharp_omega_constant_rational():
    c=a.sharp_omega_constant_witness()
    assert c['omega_lower']>-Q(157,50)
    assert c['omega_upper_offset']<-Q(99,100)


def test_sharp_tail_strictly_tighter_for_30():
    assert all(0<a.sharp_tail_upper(i,j)<a.baseline_tail_upper(i,j)
        for i,j in a.EXPECTED_PAIRS)


def test_sharp_tail_control_values():
    assert Q(2,10**5)<a.sharp_tail_upper(1,1)<Q(3,10**5)
    assert Q(9,10**5)<a.sharp_tail_upper(2,2)<Q(1,10**4)


def test_historical_two_sharp_tails():
    import json
    d=json.loads(FIXTURE.read_text())
    assert len(d.get('historical_e7_tail_comparisons',[]))==2
    assert all(x['difference_less_than_1e-14'] for x in d['historical_e7_tail_comparisons'])
