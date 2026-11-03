"""ES8-P5 historical-reconciliation regression and deliberately negative gates."""
import copy
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
import pytest

import es8p5_reconcile as audit

ROOT=Path(__file__).resolve().parent
FIXTURE=ROOT/'FIXTURE_HISTORICO_ES8P5.json'
RESULT=ROOT/'RESULTADOS_ES8P5.json'

@pytest.fixture(scope='module')
def fixture():
    return json.loads(FIXTURE.read_text(encoding='utf8'))

@pytest.fixture(scope='module')
def replay(fixture):
    return audit.analyze_fixture(fixture,64)


def test_sha_of_uploaded_archives_matches_record(fixture):
    p=fixture['provenance']
    assert p['v16_zip_sha256']==audit.V16_SHA
    assert p['e7_zip_sha256']==audit.E7_SHA
    assert p['phase4_zip_sha256']==audit.PHASE4_SHA


def test_full_uploaded_v16_manifest_checked_in_original_run(fixture):
    p=fixture['provenance']
    assert p['v16_manifest_entries']==7
    assert p['all_v16_internal_sha256'] is True
    assert p['e7_manifest_entries']==32
    assert p['phase4_hashed_integral_files']==30
    # This fixture assertion does NOT independently reopen the ZIP in CI.


def test_exact_replay_number_of_archived_raw_entries(fixture):
    assert len(fixture['raw_arch_intervals'])==30
    actual={tuple(r['pair']) for r in fixture['raw_arch_intervals']}
    expected={(i,j) for i in range(1,11) for j in range(i,11) if (i-j)%2==0}
    assert actual==expected


def test_symmetric_prime_source_has_expected_same_parity_pairs(fixture):
    assert len(fixture['prime_intervals'])==30
    assert all((r['pair'][0]-r['pair'][1])%2==0 for r in fixture['prime_intervals'])


def test_eight_projection_columns(fixture):
    assert [r['mode'] for r in fixture['projection_columns']]==list(range(3,11))
    assert all(r['pivot']==(1 if r['mode']%2 else 2) for r in fixture['projection_columns'])


def test_full_matrix_projected_dimensions(fixture):
    assert sum(r['N']==6 for r in fixture['projected_intervals'])==10
    assert sum(r['N']==8 for r in fixture['projected_intervals'])==21
    assert sum(r['N']==10 for r in fixture['projected_intervals'])==36


def test_recomputed_arch_30_of_30(replay):
    assert replay['comparison_counts']['raw_arch']=={'pass':30,'total':30}


def test_analytic_tail_is_crucial(replay):
    assert replay['comparison_counts']['raw_finite_only']=={'pass':2,'total':30}
    assert replay['comparison_counts']['raw_arch']['pass']==30


def test_recomputed_prime_30_and_projection_8(replay):
    assert replay['comparison_counts']['prime']=={'pass':30,'total':30}
    assert replay['comparison_counts']['pole_projection_columns']=={'pass':8,'total':8}


def test_recomputed_projected_67_of_67(replay):
    assert replay['comparison_counts']['projected']=={'pass':67,'total':67}
    assert replay['comparison_counts']['per_N']=={'6':{'pass':10,'total':10},
        '8':{'pass':21,'total':21},'10':{'pass':36,'total':36}}


def test_wrong_arch_interval_rejected(fixture):
    bad=copy.deepcopy(fixture)
    bad['raw_arch_intervals'][0]['historical_full_interval']=['111','112']
    with pytest.raises(AssertionError,match='failed raw_arch'):
        audit.analyze_fixture(bad,64)


def test_wrong_prime_sign_rejected(fixture):
    bad=copy.deepcopy(fixture)
    bad['prime_intervals'][0]['historical_prime_interval']=['111','112']
    with pytest.raises(AssertionError,match='failed prime'):
        audit.analyze_fixture(bad,64)


def test_wrong_projection_ratio_rejected(fixture):
    bad=copy.deepcopy(fixture)
    bad['projection_columns'][0]['historical_coefficient_interval']=['11','12']
    with pytest.raises(AssertionError,match='failed pole_projection_columns'):
        audit.analyze_fixture(bad,64)


def test_wrong_projected_interval_rejected(fixture):
    bad=copy.deepcopy(fixture)
    bad['projected_intervals'][0]['historical_interval']=['111','112']
    with pytest.raises(AssertionError,match='failed projected'):
        audit.analyze_fixture(bad,64)


def test_report_does_not_claim_rh_proof():
    result=json.loads(RESULT.read_text(encoding='utf8'))
    assert result['historical_v16_replay']['replayed_v16_json_byte_identical'] is True
    assert result['historical_v16_replay']['replayed_v16_gates']==12
    assert 'NOT_GLOBAL_WEIL_OR_RH' in result['scientific_status']
    assert 'full_source_binary_checked'==result['current_mode']


def test_exact_fraction_interval_intersections():
    assert audit.intersects(['1/10','1/5'],['1/6','1/4'])
    assert not audit.intersects(['1/10','1/5'],['2','3'])
    assert F('1/3')<F('1/2')


def test_invalid_archived_zip_sha_rejected(tmp_path):
    blob=tmp_path/'mutated-v16.zip'
    blob.write_bytes(b'not the authenticated historical ZIP')
    with pytest.raises(AssertionError,match='SHA256 mismatch'):
        audit.extract_original(blob)
