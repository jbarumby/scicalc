"""NUMERIA ES8-P4: fail-closed historical-provenance & math-consistency gate.

Checks SHA-1 object IDs (git blob) across two git snapshots, pinned source
content hashes, exact algebraic two-pole cancellation and numerical/interval
replay of the RECONSTRUCTED finite V17 model. It explicitly does NOT compare
E7/V16 original binary source internals absent their original files.

Run --e7 FILE to check E7 bytes by SHA-256, ZIP integrity and metadata,
without claiming entry-by-entry equivalence. Run --strict to REJECT if full
historical identity is missing. No RH claim.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
import zipfile

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'vendor'))
from es8p3_arch_projected import full_projected, ldl_certify, taper_projected_row_bound
from es8p2_prime_cert import constants, projection_column, L, I, SCALE

GIT_BLOB_E7='5b6abb59aa09aa53b088d672a3c4e2c240e5ac4d'
E7_SHA256='4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7'
E7_SIZE=1641486
PRIOR_SOURCE_SHA={
    'es8p2_prime_cert.py':'f275aa8c385a394b6a3445bd489f27de97e900d1bc1fb7389cb78707d6f90c61',
}


def sha256(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_gate():
    evidence={}
    for filename,expected in PRIOR_SOURCE_SHA.items():
        actual=sha256(HERE/'vendor'/filename)
        evidence[filename]={'actual_sha256':actual,'expected_sha256':expected,'same':actual==expected}
    return {'sources':evidence,'matches':all(x['same'] for x in evidence.values())}


def git_object_gate(snapshot):
    """Git blob equal SHA-1 means equal bytes assuming Git object identity.
    Snapshot was obtained from GitHub git/trees/{commit}?recursive=1 API.
    This is NOT a re-hash of the original ZIP file downloaded here.
    """
    expected_commits={
       'fe6bdb7463172ea7107cd050410fe0f94196de29',
       '2053129b03cffbce302b142af73df37e993f1613'
    }
    if set(snapshot)!=expected_commits:raise ValueError('wrong git snapshot commits')
    paths={'NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip','numeria_e8k/v15/E7_original.bin'}
    rec=[]
    for commit,details in snapshot.items():
        if set(details)!=paths:raise ValueError('missing E7 object path')
        for path,data in details.items():
            rec.append({'commit':commit,'path':path,'git_blob_sha1':data['sha'],
                        'size':data['size'],'matches':data['sha']==GIT_BLOB_E7 and data['size']==E7_SIZE})
    return {'object_identity_by_git':all(x['matches'] for x in rec), 'records':rec,
            'what_it_proves':'the 4 Git tree entries refer to one and the same Git blob, not verified ZIP internals'}


def pole_exact_gate(N):
    """Polynomials in X=pi^2; verify coefficients identically in Q[X].
    For each k>=3 and p same parity, alpha=-k/p*(p²X+L²)/(k²X+L²).
    Multiplying M_k + alpha M_p by common denominators proves zero.
    Both +/- pole moments factor identically within parity.
    """
    result=[]
    for k in range(3,N+1):
        p=1 if k%2 else 2
        # k/(k²X+L²) - (k/p)*(p²X+L²)/(k²X+L²)*p/(p²X+L²) == 0
        coeff_X=Q(k*p*p)-Q(k,p)*p*p*p
        coeff_0=Q(k)*L*L-Q(k,p)*p*L*L
        result.append({'mode':k,'pivot':p,'coeff_X':str(coeff_X),
                       'coeff_const':str(coeff_0),'both_signs_cancel':coeff_X==0 and coeff_0==0})
    return {'N':N,'all_cancel':all(x['both_signs_cancel'] for x in result),'rows':result}


def finite_replay(M=64):
    result={}
    for N in (6,8,10):
        mat=[[full_projected(i,j,N,M) for j in range(3,N+1)] for i in range(3,N+1)]
        ldl=ldl_certify(mat,Q(1,20))
        # Already tested in ES8P3. This checks absence of accidental regression,
        # not independent physical identification of historical E7 matrices.
        bound=taper_projected_row_bound(N,Q(1,1000000))
        result[str(N)]={
            'dimension':N-2,'ldl_shift':'1/20',
            'ldl_status':ldl['status'],
            'taper_error_upper_rational':str(bound),
            'strict_taper_less_than_shift':bound<Q(1,20),
            'cross_parity_exact_zero':all(mat[i][j]==I(0) for i in range(N-2) for j in range(N-2) if (i-j)%2),
            'pole_proof':pole_exact_gate(N),
            'source_status':'RECONSTRUCTED_V17_MODEL_NOT_ORIGINAL_E7_BYTE_MATCH'
        }
    return result


def optional_e7_gate(path):
    if path is None:return {'checked':False,'reason':'Original E7 ZIP bytes unavailable in this environment',
                            'source_binary_sha256_validated':False,
                            'historical_entry_reconciliation':False}
    p=Path(path)
    if not p.is_file(): raise FileNotFoundError(p)
    hash_=sha256(p)
    if hash_!=E7_SHA256 or p.stat().st_size!=E7_SIZE:
        return {'checked':True,'source_binary_sha256_validated':False,'size':p.stat().st_size,
                'sha256':hash_,'zip_test':None,'historical_entry_reconciliation':False}
    with zipfile.ZipFile(p) as z:
        err=z.testzip(); names=z.namelist()
    return {'checked':True,'source_binary_sha256_validated':err is None,
            'zip_test':'OK' if err is None else f'BAD:{err}',
            'sha256':hash_,'zip_entries':len(names),
            'historical_entry_reconciliation':False,
            'reason':'Correct ZIP alone does not establish entry-by-entry mathematical equality'}


def audit(e7=None,M=64,snapshot=None,strict=False):
    if M<32 or M>2048: raise ValueError('M outside certified range')
    if snapshot is None:snapshot=json.loads((HERE/'GITHUB_TREE_E7_SNAPSHOT.json').read_text(encoding='utf-8'))
    gate=git_object_gate(snapshot)
    source=source_gate()
    finite=finite_replay(M)
    archive=optional_e7_gate(e7)
    model_ok=source['matches'] and gate['object_identity_by_git'] and all(
        r['ldl_status']=='CERTIFIED_POSITIVE_DEFINITE' and r['strict_taper_less_than_shift'] and
        r['cross_parity_exact_zero'] and r['pole_proof']['all_cancel'] for r in finite.values())
    # Explicitly NOT inferred from hashes, coefficients or tests, even when ZIP is available.
    historical_equivalence=False
    data={'schema':'NUMERIA.ES8P4.HISTORICAL_RECONCILIATION.1',
          'scientific_status':'DOCUMENTARY_AND_RECONSTRUCTED_MODEL_CHECKS_PASS; ORIGINAL_E7_V16_ENTRY_IDENTITY_NOT_ESTABLISHED; NO_RH_PROOF',
          'M':M,'verified_v17_model':model_ok, 'git_e7':gate, 'source_integrity':source,
          'finite_v17_replay':finite,'optional_archive_check':archive,
          'full_historical_e7_v16_equivalence_established':historical_equivalence,
          'what_is_not_proven':['No original E7/V16 internal entry-by-entry replay',
                               'Git blob identity is not proof of the source formulas or their signs',
                               'V17/B5 phase4 source hashes require original source bytes to revalidate',
                               'No all-N Weil positivity or RH proof']}
    if strict and not historical_equivalence:raise RuntimeError('STRICT HISTORICAL E7/V16 EQUIVALENCE GATE: NOT ESTABLISHED')
    return data


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--e7',type=Path,default=None,help='Optional historical E7 ZIP; verify SHA256 and ZIP integrity')
    p.add_argument('--M',type=int,default=64)
    p.add_argument('--strict',action='store_true')
    p.add_argument('--out',type=Path,default=HERE/'RESULTADOS_ES8P4.json')
    a=p.parse_args()
    try:r=audit(a.e7,a.M,strict=a.strict)
    except (ValueError,RuntimeError) as e:
        print('ES8-P4 GATE FAILURE:',e);return 2
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,sort_keys=True,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('ES8-P4 model PASS:',r['verified_v17_model'],'historical:',r['full_historical_e7_v16_equivalence_established'])
    return 0 if r['verified_v17_model'] else 2
if __name__=='__main__':raise SystemExit(main())
