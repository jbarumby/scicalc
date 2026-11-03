#!/usr/bin/env python3
"""NUMERIA ES8-P5: independent old ZIP / ES8-P3 interval reconciliation.

Never infer a proof of the entire Weil criterion or RH from 30/30 and 67/67
interval intersections. The E7, V16, and ES8-P3 enclosure lemmas and the
Guinand-Weil specialization remain mathematical premises to audit separately.

Two modes:
  --v16-upload <the exact original V16 ZIP>: full source+archive SHA replay.
  --fixture <the certified-by-local-extraction JSON>: source snapshot replay
     in GitHub Actions WITHOUT claiming original ZIP binary authentication.
"""
import argparse
import hashlib
import importlib.util
import io
import json
import sys
import tempfile
import zipfile
from fractions import Fraction as F
from pathlib import Path

V16_SHA = 'aee5bfedf977f7a0a77d66949753faaab673bb9fa18b01b4425160e464936826'
E7_SHA = '4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7'
PHASE4_SHA = 'b0d7d335534239a6863e655295a658282551d334b0808ee801710a4b5699d8fe'
ROOT = Path(__file__).resolve().parent
P3_SOURCE = ROOT / 'es8p3_arch_projected.py'
P2_SOURCE = ROOT / 'es8p2_prime_cert.py'


def sha(data): return hashlib.sha256(data).hexdigest()


def ensure_sources():
    if not P3_SOURCE.is_file() or not P2_SOURCE.is_file():
        raise FileNotFoundError('Supply the ES8-P3 and ES8-P2 sources beside this auditor')
    if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
    import es8p3_arch_projected as new
    import es8p2_prime_cert as p2
    return new, p2


def intersects(a,b):
    return F(a[0]) <= F(b[1]) and F(b[0]) <= F(a[1])


def narrow_interval(v):return [str(v.lo),str(v.hi)]


def counts(records,key='overlap'):
    return {'pass':sum(bool(x[key]) for x in records),'total':len(records)}


def independent_matches(data, M=128):
    """Recompute P3 rather than accepting stored new-value fields in fixture."""
    new,p2=ensure_sources()
    raw=[];primes=[];columns=[];projected=[]
    for row in data['raw_arch_intervals']:
        i,j=row['pair']; ref=row['historical_full_interval']
        a=new.arch_raw(i,j,M)
        raw.append({'pair':[i,j], 'overlap':intersects(ref,narrow_interval(a)),
                    'historical_full_interval':ref,
                    'current_recomputed_interval':narrow_interval(a),
                    'tail_radius':row['tail_radius'],
                    'finite_only_overlap':intersects(row['historical_finite_interval'],narrow_interval(a))})
    for row in data['prime_intervals']:
        i,j=row['pair']; ref=row['historical_prime_interval']; b=p2.raw_prime(i,j)
        primes.append({'pair':[i,j], 'overlap':intersects(ref,narrow_interval(b))})
    for row in data['projection_columns']:
        k,p=row['mode'],row['pivot'];newcoef=p2.projection_column(k,10)[p]
        columns.append({'mode':k,'pivot':p,'overlap':intersects(row['historical_coefficient_interval'],narrow_interval(newcoef))})
    for row in data['projected_intervals']:
        N=row['N']; i,j=row['pair'];ref=row['historical_interval']
        a=new.full_projected(i,j,N,M)
        projected.append({'N':N,'pair':[i,j],'overlap':intersects(ref,narrow_interval(a))})
    result={'M':M,'raw_arch':counts(raw),'prime':counts(primes),
            'pole_projection_columns':counts(columns),'projected':counts(projected),
            'per_N':{str(N):counts([x for x in projected if x['N']==N]) for N in [6,8,10]},
            'raw_finite_only':counts(raw,'finite_only_overlap'),
            'raw_detail':raw,'projected_failures':[x for x in projected if not x['overlap']]}
    expected={'raw_arch':30,'prime':30,'pole_projection_columns':8,'projected':67}
    for k,n in expected.items():
        if result[k]!={'pass':n,'total':n}:
            raise AssertionError(f'failed {k}: {result[k]} (expected {n}/{n})')
    if result['per_N']!={'6':{'pass':10,'total':10},'8':{'pass':21,'total':21},'10':{'pass':36,'total':36}}:
        raise AssertionError('Unexpected projected dimensions')
    return result


def zip_integrity(z):
    bad=z.testzip()
    if bad is not None:raise AssertionError('ZIP CRC failed at '+bad)
    return True


def load_v16_module(source_bytes,folder):
    path=Path(folder)/'v16_original.py';path.write_bytes(source_bytes)
    spec=importlib.util.spec_from_file_location('v16_original_replay',str(path))
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def extract_original(v16_path, extra_e7_paths=(), nested_e7_path=None):
    b=Path(v16_path).read_bytes()
    if sha(b)!=V16_SHA:raise AssertionError('Original V16 outer ZIP SHA256 mismatch')
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        zip_integrity(z)
        hashes={}
        lines=z.read('SHA256SUMS.txt').decode().splitlines()
        for line in lines:
            digest,name=line.split('  ',1)
            valid=sha(z.read(name))==digest
            hashes[name]=valid
            if not valid:raise AssertionError('Bad V16 internal SHA256: '+name)
        e7=z.read('NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip')
        p4=z.read('NUMERIA_FASE4_30_ENTRADAS_AUDITORIA_PROJECAO_2026-10-24.zip')
        if sha(e7)!=E7_SHA or sha(p4)!=PHASE4_SHA:
            raise AssertionError('E7 or Phase4 ZIP SHA256 mismatch')
        v16source=z.read('v16/AUDITAR_E8K_V16_FORMULA.py')
        expected_v16_result=z.read('v16/RESULTADO_E8K_V16.json')
    for path in extra_e7_paths:
        if sha(Path(path).read_bytes())!=E7_SHA:raise AssertionError('Nonidentical E7 archive '+str(path))
    if nested_e7_path:
        with zipfile.ZipFile(nested_e7_path) as n:
            zip_integrity(n)
            name='NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip'
            if sha(n.read(name))!=E7_SHA:raise AssertionError('Nested E7 differs')
    with zipfile.ZipFile(io.BytesIO(e7)) as ez:
        zip_integrity(ez)
        manifest=json.loads(ez.read('MANIFESTO_SHA256.json'))['sha256']
        failures=[name for name,digest in manifest.items() if sha(ez.read(name))!=digest]
        if failures:raise AssertionError('E7 manifest hash failed: '+repr(failures))
    with zipfile.ZipFile(io.BytesIO(p4)) as pz:
        zip_integrity(pz)
        summ=json.loads(pz.read('RESUMO_FASE4_1600CELULAS_D40.json'))
        if len(summ['records'])!=30:raise AssertionError('Wrong 30-record Phase4 count')
        for r in summ['records']:
            name=f"FASE4_G{r['i']}{r['j']}_1600CELULAS_D40.json"
            if sha(pz.read(name))!=r['sha256']:
                raise AssertionError('Phase4 individual source mismatch '+name)
    return {'v16_zip_sha256':sha(b),'e7_zip_sha256':sha(e7),'phase4_zip_sha256':sha(p4),
            'v16_manifest_entries':len(hashes),'all_v16_internal_sha256':all(hashes.values()),
            'e7_manifest_entries':len(manifest),'phase4_hashed_integral_files':30,
            'extra_identical_e7_copies':len(extra_e7_paths),'nested_identical_e7':bool(nested_e7_path)},v16source,e7,p4,expected_v16_result


def from_full_archive(v16_path, extra_e7_paths=(), nested_e7_path=None):
    integrity,source,e7,p4,expected = extract_original(v16_path,extra_e7_paths,nested_e7_path)
    with tempfile.TemporaryDirectory() as t:
        temp=Path(t)
        (temp/'e7.zip').write_bytes(e7)
        (temp/'phase4.zip').write_bytes(p4)
        v16=load_v16_module(source,temp)
        replay=v16.run(temp/'e7.zip',temp/'phase4.zip')
        text=json.dumps(replay,ensure_ascii=False,indent=2)+'\n'
        bytes_equal=text.encode()==expected
        if not bytes_equal:raise AssertionError('V16 replay does not equal historical JSON bytewise')
        with zipfile.ZipFile(temp/'phase4.zip') as phase:
            G,hashes,summary=v16.phase4_data(phase)
        P={}
        for i in range(1,11):
            for j in range(i,11):
                P[i,j]=v16.Iv(0) if (i-j)%2 else (-2*v16.LN2/v16.SQ2)*v16.overlap(i,j)
        B=[[v16.Iv(0) for _ in range(8)]for _ in range(10)]
        for k in range(3,11):
            p=1 if k%2 else 2
            B[k-1][k-3]=v16.Iv(1)
            B[p-1][k-3]=v16.residue_ratio(k,p)
        raw=[]
        for (i,j),g in sorted(G.items()):
            tail=(4*i*j)*v16.PI/(v16.LF*v16.LF*v16.LF)*(8+v16.LOG480+F(1,3))/(3*160**3)
            tail_bound=F(tail.b,v16.S)
            full=g+v16.Iv(-tail_bound,tail_bound)
            lo,hi=g.bounds();l2,h2=full.bounds()
            raw.append({'pair':[i,j], 'historical_finite_interval':[str(lo),str(hi)],
                        'tail_radius':str(tail_bound),'historical_full_interval':[str(l2),str(h2)]})
        primes=[{'pair':[i,j],'historical_prime_interval':list(map(str,P[i,j].bounds()))} for i,j in sorted(P) if (i-j)%2==0]
        poles=[{'mode':k,'pivot':(1 if k%2 else 2),
                'historical_coefficient_interval':list(map(str,v16.residue_ratio(k,1 if k%2 else 2).bounds()))}
               for k in range(3,11)]
        projected=[]
        for N in [6,8,10]:
            A=v16.build(N,G,P,B)
            for c in range(N-2):
                for d in range(c,N-2):
                    projected.append({'N':N,'pair':[c+3,d+3],
                                      'historical_interval':list(map(str,A[c][d].bounds()))})
        fixture={'provenance':integrity,'raw_arch_intervals':raw,'prime_intervals':primes,
                 'projection_columns':poles,'projected_intervals':projected,
                 'source':'exact rational fields extracted from the uploaded V16 ZIP with all SHA256 checks'}
        return fixture,{'replayed_v16_json_byte_identical':bytes_equal,
                        'replayed_v16_gates':sum(replay['gates'].values()),
                        'replayed_v16_gates_total':len(replay['gates']),
                        'historical_pivots':{str(x['N']):x['min_shifted_pivot_approx'] for x in replay['projected_matrices']},
                        **integrity}


def analyze_fixture(data,M=128,original_replay=None):
    check=independent_matches(data,M)
    return {'schema':'NUMERIA.ES8P5.HISTORICAL_RECONCILIATION.1',
            'scientific_status':'CONDITIONAL_HISTORICAL_INTERVAL_COMPATIBILITY_NOT_GLOBAL_WEIL_OR_RH',
            'source_provenance':data['provenance'],
            'current_mode':'full_source_binary_checked' if original_replay else 'text_fixture_replay_only',
            'historical_v16_replay':original_replay,
            'comparison_counts':{key:check[key] for key in ['raw_arch','raw_finite_only','prime','pole_projection_columns','projected','per_N']},
            'raw_details':check['raw_detail'],
            'caveats':['Interval intersection is compatibility, NOT an identity proof of true integrals',
                 'Phase4 Taylor/Cauchy majorant and E7 tail lemma remain analytic premises',
                 'ES8-P3 Hurwitz Euler-Maclaurin enclosures are separately conditional on stated lemmas',
                 'Guinand-Weil explicit-formula normalization is an external mathematical premise',
                 'Only N=6,8,10 in projected subspace; no universal positivity and no Riemann Hypothesis proof']}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--v16-upload',type=Path)
    p.add_argument('--extra-e7',type=Path,action='append',default=[])
    p.add_argument('--nested-e7',type=Path)
    p.add_argument('--fixture',type=Path)
    p.add_argument('--out',type=Path,default=Path('RESULTADOS_ES8P5.json'))
    p.add_argument('--fixture-out',type=Path)
    p.add_argument('--M',type=int,default=128)
    args=p.parse_args()
    if bool(args.v16_upload)==bool(args.fixture):p.error('select one of --v16-upload or --fixture')
    if args.v16_upload:
        fixture,replay=from_full_archive(args.v16_upload,args.extra_e7,args.nested_e7)
    else:
        fixture=json.loads(args.fixture.read_text(encoding='utf8'))
        replay=None
    result=analyze_fixture(fixture,args.M,replay)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf8')
    if args.fixture_out:
        args.fixture_out.parent.mkdir(parents=True,exist_ok=True)
        args.fixture_out.write_text(json.dumps(fixture,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf8')
    print('ES8P5',result['current_mode'],'| raw:',result['comparison_counts']['raw_arch'],
          'prime:',result['comparison_counts']['prime'],
          'projection B:',result['comparison_counts']['pole_projection_columns'],
          'projected:',result['comparison_counts']['projected'])
    print('NO UNIVERSAL WEIL POSITIVITY; NO RH PROOF')

if __name__=='__main__':main()
