from __future__ import annotations
from pathlib import Path
import json
import sys
import hashlib
from fractions import Fraction as Q
import zipfile
import pytest

BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE))
import es8p4_reconcile as e
import triage_es8p4 as b


def test_git_blob_two_paths_two_commits():
    o=json.loads((BASE/'GITHUB_TREE_E7_SNAPSHOT.json').read_text())
    g=e.git_object_gate(o)
    assert g['object_identity_by_git'] and len(g['records'])==4


def test_git_mutation_fails_closed():
    o=json.loads((BASE/'GITHUB_TREE_E7_SNAPSHOT.json').read_text())
    o[next(iter(o))]['NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip']['sha']='0'*40
    assert not e.git_object_gate(o)['object_identity_by_git']


def test_git_missing_path_rejected():
    o=json.loads((BASE/'GITHUB_TREE_E7_SNAPSHOT.json').read_text())
    del o[next(iter(o))]['numeria_e8k/v15/E7_original.bin']
    with pytest.raises(ValueError):e.git_object_gate(o)


def test_vendor_same_bytes_as_prior_es8p2():
    assert e.source_gate()['matches']
    assert e.sha256(BASE/'vendor/es8p2_prime_cert.py')==e.PRIOR_SOURCE_SHA['es8p2_prime_cert.py']


def test_double_pole_moment_cancellation():
    assert all(e.pole_exact_gate(n)['all_cancel'] for n in (6,8,10))


def test_matrix_replay_shifted_ldl_positive():
    r=e.finite_replay(64)
    assert len(r)==3
    for x in r.values():
        assert x['ldl_status']=='CERTIFIED_POSITIVE_DEFINITE'
        assert x['cross_parity_exact_zero']
        assert x['strict_taper_less_than_shift']


def test_finite_gate_does_not_pretend_historical_identity():
    r=e.audit(M=64)
    assert r['verified_v17_model'] is True
    assert r['full_historical_e7_v16_equivalence_established'] is False
    assert r['optional_archive_check']['checked'] is False
    assert 'NOT_ESTABLISHED' in r['scientific_status']


def test_strict_mode_rejects_without_original():
    with pytest.raises(RuntimeError,match='NOT ESTABLISHED'):
        e.audit(M=64,strict=True)


def test_wrong_zip_is_rejected(tmp_path):
    p=tmp_path/'impostor.zip'
    with zipfile.ZipFile(p,'w') as z:z.writestr('manifest.json','{}')
    x=e.optional_e7_gate(p)
    assert x['checked'] is True
    assert x['source_binary_sha256_validated'] is False
    assert x['historical_entry_reconciliation'] is False


def test_missing_zip_file_fails():
    with pytest.raises(FileNotFoundError):e.optional_e7_gate(BASE/'not-a-file.zip')


def test_invalid_cutoffs():
    for M in (0,31,2049):
        with pytest.raises(ValueError):e.audit(M=M)


def test_csv_triage_counts_and_top_sources():
    t=b.triage([BASE/'References_LeapSpace (14).csv',BASE/'References_LeapSpace (15).csv'])
    assert t['input_records']==31 and t['unique_records']==28 and t['duplicate_excess']==3
    direct=[r for r in t['records'] if r['level'].startswith('ALTA')]
    assert {r['doi'] for r in direct}=={'10.1016/S0377-0427(00)00336-8','10.1016/j.jmaa.2017.12.045'}


def test_duplicate_doi_not_double_counted():
    t=b.triage([BASE/'References_LeapSpace (14).csv',BASE/'References_LeapSpace (15).csv'])
    assert next(r for r in t['records'] if r['doi']=='10.1016/j.jmaa.2005.01.020')['occurrences']==3
    assert next(r for r in t['records'] if r['doi']=='10.1016/j.jnt.2013.03.004')['occurrences']==2


def test_reconciliation_json_replay_matches_saved():
    data=e.audit(M=64)
    committed=json.loads((BASE/'RESULTADOS_ES8P4.json').read_text())
    assert data==committed


def test_nonorthonormal_basis_projected_shift_is_not_orthonormal_eigenvalue():
    # B columns are not orthonormal; both the raw basis-coordinate eigenvalues
    # and change-of-basis to Gram-normalized metric are distinct claims.
    from es8p4_reconcile import projection_column
    c=projection_column(3,6)
    assert set(c)=={1,3}
    assert c[1].hi<0
    # Our strict interpretation is positivity of B^T M B minus .05 I
    # in the B-coordinates, not a uniform Hilbert-space spectral gap.
    assert e.audit(M=64)['full_historical_e7_v16_equivalence_established'] is False
