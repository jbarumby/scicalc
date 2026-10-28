#!/usr/bin/env python3
# coding: utf-8
"""NUMERIA E8-K V12: conferidor de evidencias SEM executar codigo B5.

Fonte: ZIP de release V12 contendo o ZIP E7. A conferencia compara bytes,
JSONL, fracoes racionais e saidas E7 diretamente com o relatorio V12.
Nao e prova independente da corretude matematica do integrando, das cotas
analiticas de quadratura/cauda ou da Hipotese de Riemann.

Uso: python CONFERIR_E8K_V12_V2.py --pacote E8K_V12_GITHUB_COMPLETO.zip \
 --out CONFERENCIA_E8K_V12_V2.json --selftest
Retorno 0: gates aprovados; 2: falha; 3: erro de integracao.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import io
import json
import re
import sys
import zipfile
from types import SimpleNamespace
from fractions import Fraction as F
from pathlib import Path

EXPECTED_SCHEMA='NUMERIA.E8K.V12.FULL_POINT_REPLAY.1'
TICK=F(1,10**42)
EXPECTED_E7_SHA='4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7'
EXPECTED_B5_SHA='a1753d3010e34f29371a9706cd6bb05b3b0c8491ca3cafb44bf98b7cd53d3df9'
NAMES={'G11':'PONTOS_G11_B5.jsonl','G22':'PONTOS_G22_B5.jsonl'}
B3_NAMES={'G11':'proveniencia_b3/RESULTADO_B3_G11_WINDOWS.json', 'G22':'proveniencia_b3/RESULTADO_G12_G22_WINDOWS.json'}

class FolderReader:
    """Lê somente arquivos dentro da raiz local (modo GitHub Actions)."""
    def __init__(self,root):self.root=Path(root).resolve()
    def read(self,name):
        target=(self.root/name).resolve()
        if not target.is_relative_to(self.root):raise ValueError('caminho fora da raiz: '+name)
        return target.read_bytes()
    def infolist(self):
        return [SimpleNamespace(filename=str(p.relative_to(self.root)).replace('\\','/'),is_dir=lambda: False)
                for p in self.root.rglob('*') if p.is_file()]

def sha(b):return hashlib.sha256(b).hexdigest()
def f(x):return F(x)
def dec_err(x):
    """Metade della unita dell'ultima casa esplicita, inclusive expoente cientifico."""
    t=str(x).lower().strip()
    m=re.fullmatch(r'([+-]?\d+)(?:\.(\d+))?(?:e([+-]?\d+))?',t)
    if not m:raise ValueError('decimal invalido '+repr(x))
    exponent=int(m.group(3) or 0)-len(m.group(2) or '')
    return F(10**exponent,2) if exponent>=0 else F(1,2*10**(-exponent))

def assert_true(cond,msg):
    if not cond: raise AssertionError(msg)

def check_sha_manifest(pack):
    text=pack.read('SHA256SUMS.txt').decode('utf-8')
    names={x.filename for x in pack.infolist() if not x.is_dir()}
    listed=set(); checked=0
    for line in text.splitlines():
        if not line.strip() or line.startswith('#'):continue
        m=re.fullmatch(r'([0-9a-f]{64})  ([^\r\n]+)',line)
        assert_true(m is not None,'SHA256SUMS malformado')
        digest,name=m.groups()
        assert_true(name in names and name not in listed,'entrada ausente/duplicada SHA256SUMS '+name)
        listed.add(name)
        assert_true(sha(pack.read(name))==digest,'hash divergente '+name)
        checked+=1
    assert_true(checked>=8,'poucas entradas no SHA256SUMS')
    assert_true('RESULTADO_E8K_V12.json' in listed,'resultado nao listado no manifesto')
    assert_true('NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip' in listed,'E7 nao listado')
    return {'files_checked':checked,'sha256sums_present':True}

def check_e7_integrity(pack):
    zbytes=pack.read('NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip')
    assert_true(sha(zbytes)==EXPECTED_E7_SHA,'hash ZIP E7 divergente')
    e7=zipfile.ZipFile(io.BytesIO(zbytes))
    manifest=json.loads(e7.read('MANIFESTO_SHA256.json'))
    mapping=manifest['sha256']
    checked=0
    for name,digest in mapping.items():
        assert_true(sha(e7.read(name))==digest,'hash manifesto E7 diverge: '+name)
        checked+=1
    b5= json.loads(e7.read('b5/RESULTADO_B5.json'))
    assert_true(sha(e7.read('fontes/AUDITOR_PONTOS_B5.py'))==EXPECTED_B5_SHA,'SHA fonte B5 divergente')
    assert_true(b5['source_sha256']==EXPECTED_B5_SHA,'fonte B5 reportada divergente')
    assert_true(b5['points']==15754,'B5 total diferente')
    return e7,b5,{'sha256_e7':sha(zbytes),'source_b5_sha256':EXPECTED_B5_SHA,'manifest_entries_checked':checked}

def check_b3_provenance(pack,b5,doc):
    result={}
    for item in b5['interval_results']:
        name=item['component']; path=B3_NAMES[name]
        data=pack.read(path)
        real_sha=sha(data)
        assert_true(real_sha==item['source_input_sha256'],'B3 original nao corresponde entrada B5: '+name)
        assert_true(doc['entries'][name]['original_b3_input_sha256_recovered_and_verified'] is True, 'V12 nao afirma B3: '+name)
        result[name]={'path':path,'sha256':real_sha,'b5_input_sha256':item['source_input_sha256']}
    return result

def check_report_structure(doc,b5,e7):
    assert_true(doc.get('schema')==EXPECTED_SCHEMA,'schema divergente')
    assert_true(doc.get('all_point_exports_reproduced') is True,'all_point_exports_reproduced falso')
    assert_true(doc.get('status')=='ALL_POINTS_REPRODUCED_B5_INPUT_HASHES_MATCHED','status diferente')
    assert_true(doc['source_integrity']['e7_zip_sha256']==EXPECTED_E7_SHA,'hash E7 no JSON V12 falso')
    assert_true(doc['source_integrity']['b5_source_sha256']==EXPECTED_B5_SHA,'hash fonte B5 no JSON V12 falso')
    assert_true(doc.get('provenance_gaps')==[],'provenance_gaps nao vazio')
    assert_true(all(v is True for v in doc['negative_checker_controls'].values()),'controle negativo nao passou')
    assert_true(len(doc['negative_checker_controls'])>=4,'faltam controles negativos')
    counts=doc['point_counts']; assert_true(counts=={'G11':9524,'G22':6230,'total':15754},'contagem global incorreta')
    assert_true(doc['exact_integrand_enclosures_replayed']==15754,'replay de integrands incorreto')
    assert_true(doc['exact_weighted_enclosures_replayed']==15754,'replay de pesos incorreto')
    assert_true(doc['quadrature_and_tail_formulas_arithmetically_recomposed']==2,'recomposicoes divergentes')
    for en in b5['interval_results']:
        name=en['component']; entry=doc['entries'][name]; bands=entry['bands']; base=en['bands']
        assert_true(len(bands)==len(base)==entry['bands_reproduced'],'numero de bandas divergente: '+name)
        assert_true(entry['points_regenerated_and_compared']==en['points']==counts[name],'contagem distinta B5/V12: '+name)
        assert_true(entry['integrands_reproduced']==counts[name] and entry['weighted_products_reproduced']==counts[name],'contagens de operacoes: '+name)
        assert_true(sum(x['count'] for x in bands)==counts[name],'soma das bandas incorreta: '+name)
        assert_true(f(bands[0]['from'])==0 and f(bands[-1]['to'])==160,'limites de bandas incorretos: '+name)
        for i,(band,orig) in enumerate(zip(bands,base)):
            assert_true(band['sum_match'] is True,'banda sem sum_match: '+name)
            assert_true(f(band['from'])==f(orig['from']) and f(band['to'])==f(orig['to']), 'faixa divergente do B5: '+name)
            assert_true(band['count']==orig['cells'],'contagem de faixa divergente: '+name)
            assert_true(f(band['from'])<f(band['to']),'banda vazia/invertida: '+name)
            if i: assert_true(f(bands[i-1]['to'])==f(band['from']),'bandas nao contiguas: '+name)
        assert_true(f(entry['finite_sum_lower'])==f(en['certified_finite_midpoint_sum'][0]) and f(entry['finite_sum_upper'])==f(en['certified_finite_midpoint_sum'][1]),'soma finita divergente: '+name)
        assert_true(f(entry['finite_sum_lower'])<=f(entry['finite_sum_upper']),'intervalo finito invertido')
    return {'points':15754,'bands':11,'negative_controls':len(doc['negative_checker_controls'])}

def _check_decimal_display(exact,decimal_string,name):
    assert_true(abs(exact-f(decimal_string))<=dec_err(decimal_string),f'decimal incorreto em {name}')

def check_jsonl_bands(e7,b5,doc):
    results={}
    for en in b5['interval_results']:
        name=en['component']; entry=doc['entries'][name]
        path='b5/'+NAMES[name]
        blob=e7.read(path)
        expected_file=next(x for x in b5['point_files'] if x['file']==NAMES[name])
        assert_true(sha(blob)==expected_file['sha256']==entry['jsonl_sha256'],'hash JSONL diverge: '+name)
        lines=blob.splitlines()
        assert_true(len(lines)==en['points']==expected_file['count'],'total de linhas JSONL diverge: '+name)
        idx=0;total_lo=F(0);total_hi=F(0);max_point_w=F(0)
        for band in en['bands']:
            n=int(band['cells']); sumlo=F(0);sumhi=F(0)
            a=f(band['from']); b=f(band['to']);step=(b-a)/n
            for j,line in enumerate(lines[idx:idx+n]):
                obj=json.loads(line)
                assert_true(f(obj['r'])==a+(2*j+1)*step/2,'malha JSONL incorreta '+name)
                w0,w1=(f(x) for x in obj['weighted']); i0,i1=(f(x) for x in obj['integrand'])
                assert_true(i0<=i1 and w0<=w1,'extremos invertidos '+name)
                # B5 multiplica intervalos quantizados; igualdade exata com
                # integrand*step nao e esperada. Exigir inclusao dirigida e
                # limite independente de excesso por quantizacao de um tick.
                assert_true(w0<=i0*step and i1*step<=w1,
                            'produto nao contem integrand*step: '+name)
                slack=(abs(step)+max(abs(i0),abs(i1))+3)*TICK
                assert_true(i0*step-w0<=slack and w1-i1*step<=slack,
                            'excesso de arredondamento supera cota: '+name)
                max_point_w=max(max_point_w,i1-i0)
                sumlo+=w0;sumhi+=w1
            assert_true(sumlo==f(band['sum_interval'][0]) and sumhi==f(band['sum_interval'][1]),'soma de banda incorreta: '+name)
            total_lo+=sumlo;total_hi+=sumhi;idx+=n
        assert_true(idx==len(lines),'linhas nao consumidas: '+name)
        assert_true(total_lo==f(en['certified_finite_midpoint_sum'][0]) and total_hi==f(en['certified_finite_midpoint_sum'][1]),'soma total distinta: '+name)
        _check_decimal_display(total_hi-total_lo,entry['finite_sum_width_approx'],name+' width')
        _check_decimal_display(max_point_w,entry['max_point_interval_width_approx'],name+' max point width')
        results[name]={'sha256':sha(blob),'points':len(lines),'finite_sum_interval':[str(total_lo),str(total_hi)],'bands_checked':len(en['bands'])}
    return results

def check_e7_composite(e7,doc,b5):
    official=json.loads(e7.read('RESULTADO_INTEGRAIS_E7.json'))
    results={}
    for en in b5['interval_results']:
        name=en['component'];entry=doc['entries'][name]
        o=next(x for x in official['results'] if x['entry']==name)
        finite_lo=f(entry['finite_sum_lower']);finite_hi=f(entry['finite_sum_upper'])
        assert_true([finite_lo,finite_hi]==[f(x) for x in o['sum_b5']],'soma finita E7 diverge: '+name)
        q=f(o['quadrature_error_b3']); t=f(o['tail_proposition6']); assert_true(q>=0 and t>=0,'erro negativo')
        composite=[finite_lo-q-t,finite_hi+q+t]
        assert_true(composite==[f(x) for x in o['combined_interval']],'intervalo E7 composto diverge: '+name)
        _check_decimal_display(q,entry['quadrature_bound_declared_approx'],name+' quad')
        _check_decimal_display(t,entry['tail_bound_declared_approx'],name+' tail')
        for i,s in enumerate(entry['e7_composite_interval_approx']):
            _check_decimal_display(composite[i],s,name+' composite '+str(i))
        results[name]={'quadrature_bound_exact':str(q),'tail_bound_exact':str(t),'composite_lower_exact':str(composite[0]),'composite_upper_exact':str(composite[1]),'analytical_bounds_certified':False}
    return results

def check_b3_b5_non_nesting(b5,doc):
    metrics={}
    for en in b5['interval_results']:
        name=en['component'];entry=doc['entries'][name]
        a,b=[f(x) for x in en['certified_finite_midpoint_sum']]
        c,d=[f(x) for x in en['b3_midpoint_interval']]
        nested=(c<=a and b<=d); intersects=(max(a,c)<=min(b,d))
        assert_true(entry['b5_enclosed_in_b3'] is nested,'aninhamento B3/B5 contraditorio: '+name)
        assert_true(entry['b5_intersects_b3'] is intersects,'intersecao B3/B5 contraditoria: '+name)
        assert_true(entry['b5_width_at_most_twice_b3'] is ((b-a)<=2*(d-c)),'largura B3/B5 contraditoria: '+name)
        metrics[name]={'b5_inside_b3':nested,'intersects':intersects,
            'lower_endpoint_difference_exact':str(a-c),
            'upper_endpoint_difference_exact':str(b-d),
            'max_endpoint_displacement_exact':str(max(abs(a-c),abs(b-d))),
            'overlap_length_exact':str(max(F(0),min(b,d)-max(a,c))),
            'comparison_scope':'total finite midpoint intervals only; NOT fraction of individual points'}
    return metrics

def check_all(pack, run_selftest=False):
    gates={};metrics={}; e7=None;b5=None;doc=None
    def run(name,func):
        try:
            metrics[name]=func();gates[name]={'status':'PASS','detail':'verificado com artefatos originais'}
        except Exception as exc:
            gates[name]={'status':'FAIL','detail':str(exc)[:500]}
    run('sha256sums_integrity',lambda:check_sha_manifest(pack))
    def load():
        nonlocal e7,b5,doc
        doc=json.loads(pack.read('RESULTADO_E8K_V12.json'))
        e7,b5,metrics0=check_e7_integrity(pack)
        return metrics0
    run('e7_manifest_and_b5_integrity',load)
    if e7 is not None and doc is not None:
        run('b3_provenance',lambda:check_b3_provenance(pack,b5,doc))
        run('report_structure',lambda:check_report_structure(doc,b5,e7))
        run('jsonl_bands_exact_recompute',lambda:check_jsonl_bands(e7,b5,doc))
        run('e7_composite_exact_recompute',lambda:check_e7_composite(e7,doc,b5))
        run('b3_b5_non_nesting_quantified',lambda:check_b3_b5_non_nesting(b5,doc))
    else:
        for name in ['b3_provenance','report_structure','jsonl_bands_exact_recompute','e7_composite_exact_recompute','b3_b5_non_nesting_quantified']:
            gates[name]={'status':'FAIL','detail':'nao executado: dependencia de integridade indisponivel'}
    if run_selftest:
        # Mutations evaluated against a passing baseline of the SPECIFIC check.
        tests=[]
        if doc is not None and b5 is not None and e7 is not None:
            cases=[
                ('band_count', 'report_structure',lambda d:d['entries']['G11']['bands'][0].__setitem__('count',d['entries']['G11']['bands'][0]['count']+1),lambda d:check_report_structure(d,b5,e7)),
                ('band_gap', 'report_structure',lambda d:d['entries']['G22']['bands'][1].__setitem__('from','13'),lambda d:check_report_structure(d,b5,e7)),
                ('finite_endpoint', 'report_structure',lambda d:d['entries']['G22'].__setitem__('finite_sum_upper','0'),lambda d:check_report_structure(d,b5,e7)),
                ('jsonl_hash', 'jsonl_bands_exact_recompute',lambda d:d['entries']['G11'].__setitem__('jsonl_sha256','0'*64),lambda d:check_jsonl_bands(e7,b5,d)),
                ('composite_display','e7_composite_exact_recompute',lambda d:d['entries']['G11']['e7_composite_interval_approx'].__setitem__(0,'-1.5'),lambda d:check_e7_composite(e7,d,b5)),
                ('b3_nesting_claim','b3_b5_non_nesting_quantified',lambda d:d['entries']['G22'].__setitem__('b5_enclosed_in_b3',True),lambda d:check_b3_b5_non_nesting(b5,d)),
            ]
            for label,target,mut,checker in cases:
                if gates[target]['status']!='PASS':tests.append({'mutation':label,'detected':False,'reason':'baseline gate nao PASS'});continue
                changed=copy.deepcopy(doc);mut(changed)
                try:checker(changed); detected=False
                except Exception:detected=True
                tests.append({'mutation':label,'detected':detected})
        else:tests.append({'mutation':'not_run','detected':False})
        metrics['selftest']=tests
        gates['selftest']={'status':'PASS' if tests and all(x['detected'] for x in tests) else 'FAIL','detail':f'{sum(x["detected"] for x in tests)}/{len(tests)} mutacoes isoladamente detectadas'}
    passed=bool(gates) and all(x['status']=='PASS' for x in gates.values())
    return {'schema':'NUMERIA.E8K.V12.INDEPENDENT_EVIDENCE_CHECK.2','status':'ALL_AVAILABLE_EVIDENCE_CHECKS_PASSED' if passed else 'EVIDENCE_CHECK_FAILURES', 'all_gates_passed':passed,'gates':gates,'metrics':metrics,'not_proved':['corretude transcendental de cada integrando','validade analitica das cotas de quadratura e cauda','formula explicita de Weil global','Hipotese de Riemann']}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--pacote',type=Path,help='E8K_V12_GITHUB_COMPLETO.zip')
    g.add_argument('--root',type=Path,help='pasta raiz extraida (uso no GitHub)')
    ap.add_argument('--out',type=Path,default=Path('CONFERENCIA_E8K_V12_V2.json'))
    ap.add_argument('--selftest',action='store_true')
    args=ap.parse_args()
    try:
        if args.pacote:
            with zipfile.ZipFile(args.pacote) as z:report=check_all(z,args.selftest)
        else:
            report=check_all(FolderReader(args.root),args.selftest)
    except Exception as exc:
        print('ERRO DE INTEGRACAO:',repr(exc),file=sys.stderr)
        return 3
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    for name,g in report['gates'].items():print(f"[{g['status']}] {name}: {g['detail']}")
    print('STATUS:',report['status'],'RELATORIO:',args.out)
    return 0 if report['all_gates_passed'] else 2

if __name__=='__main__':sys.exit(main())
