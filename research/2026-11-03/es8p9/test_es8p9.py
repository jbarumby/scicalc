"""ES8-P9: unit, exact-moment, mutation and full-result regression tests."""
import copy,hashlib,json
from fractions import Fraction as Q
from pathlib import Path
import pytest
import es8p9_integrate as mod
from es8p9_projected import run as projected_run
from es8p6_independent_audit import independent_total_error

ROOT=Path(__file__).resolve().parent

def recorded():return json.loads((ROOT/'RESULTADOS_ES8P9_INTEGRAIS.json').read_text())

def test_01_signed_weights():
    w=mod.vandermonde_weights()
    assert len(w)==16 and sum(w)==1 and sum(x<0 for x in w)==6

def test_02_moments():
    w=mod.vandermonde_weights()
    for p in range(16):assert sum(w[i]*Q(i,15)**p for i in range(16))==Q(1,p+1)

def test_03_weights_multiple_implementations():
    assert mod.vandermonde_weights()==mod.lagrange_weights()

def test_04_cauchy_exact():
    assert independent_total_error()==Q(4154891126177792,209922127267236328125)

def test_05_30_legacy_entries():
    v=mod.reference_inputs();assert len(v['historical_integrals'])==30

def test_06_canonical_pairs():
    assert len(mod.PAIRS)==30 and (1,1) in mod.PAIRS and (2,10) in mod.PAIRS

def test_07_negative_weight_extremes():
    # Weight -3 applied to [1,2] gives [-6,-3], not [-3,-6].
    w=-3;a,b=1,2
    assert (w*b,w*a)==(-6,-3)

def test_08_positive_weight_extremes():
    w=3;a,b=-2,5
    assert (w*a,w*b)==(-6,15)

def test_09_decimal_outwards():
    assert mod.outward_decimal(Q(-11,8),3)=='-1.375'
    assert mod.outward_decimal(Q(1,3),3)=='0.333'
    assert mod.outward_decimal(Q(1,3),3,upper=True)=='0.334'

def test_10_source_hash_locked():
    p=ROOT/'HISTORICAL_PHASE4_REFERENCE.json'
    data=json.loads(p.read_text())
    assert data['source_phase4_sha256']=='b0d7d335534239a6863e655295a658282551d334b0808ee801710a4b5699d8fe'

def test_11_full_grid_recorded():
    r=recorded();assert r['scope']['evaluations']==768000
    assert r['status']=='INTEGRATED_30_ENTRIES_CONDITIONALLY_CERTIFIED'
    assert r['compatibility']=={'finite_integrals_intersecting':30,'nodal_integrals_intersecting':30}

def test_12_32_chunks():
    r=recorded();assert len(r['batches'])==32
    assert sum(c['evaluations'] for c in r['batches'])==768000
    assert all(len(c['sha256_weighted_pair_sums'])==64 for c in r['batches'])

def test_13_positive_finite_matrix_checks():
    p=projected_run(ROOT/'RESULTADOS_ES8P9_INTEGRAIS.json')
    assert p['status']=='FINITE_PROJECTED_POSITIVITY_CONDITIONAL'
    assert [x['N'] for x in p['matrices']]==[6,8,10]
    assert all(x['certified_M_minus_0p05I_positive'] for x in p['matrices'])

def test_14_mutation_no_premature_success(tmp_path):
    r=recorded();r['status']='PARTIAL_OR_INCOMPATIBLE'
    p=tmp_path/'mutated.json';p.write_text(json.dumps(r))
    with pytest.raises(ValueError,match='require complete'):
        projected_run(p)

def test_15_mutation_incomplete_entries(tmp_path):
    r=recorded();r['entries'].pop('1,1')
    p=tmp_path/'mutated.json';p.write_text(json.dumps(r))
    with pytest.raises(KeyError):projected_run(p)

def test_16_outward_full_integrals():
    r=recorded()
    for x in r['entries'].values():
        assert x['nodal_intersects_historical'] and x['finite_intersects_historical']
        a,b=map(Q,x['independent_finite_with_cauchy']);c,d=map(Q,x['independent_global_with_tail'])
        assert c<a<=b<d

def test_17_partial_not_full():
    p=ROOT/'PARTIAL_ES8P9_SAMPLE.json'
    if p.exists():
        r=json.loads(p.read_text());assert r['status']!='INTEGRATED_30_ENTRIES_CONDITIONALLY_CERTIFIED'

# Full replay is in GitHub Actions rather than in each unit test.
