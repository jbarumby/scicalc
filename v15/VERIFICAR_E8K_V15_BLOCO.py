#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NUMERIA E8-K V15: certificação RACIONAL CONDICIONAL de sinal da matriz
arquimediana 2x2, com controle de proveniência, testemunhas e mutações isoladas.

Este arquivo NÃO prova a identidade de Weil nem todas as hipóteses analíticas
usadas no E7/V13. Não importa nem executa o B5 como oráculo independente.
"""
from __future__ import annotations
import argparse, copy, hashlib, json, re, sys
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from pathlib import Path
from zipfile import ZipFile

SCHEMA = 'NUMERIA.E8K.V15.ARCHIMEDEAN_2X2_CONDITIONAL.1'
V14_SCHEMA = 'NUMERIA.E8K.V13.FULL_POINT_REPLAY.1'
E7_SHA = '4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7'
B5_SHA = 'a1753d3010e34f29371a9706cd6bb05b3b0c8491ca3cafb44bf98b7cd53d3df9'


def sha(b: bytes) -> str: return hashlib.sha256(b).hexdigest()
def fr(x): return Q(str(x))
def dec(x, digits=24):
    with localcontext() as ctx:
        ctx.prec = digits + 12
        return format(Decimal(x.numerator)/Decimal(x.denominator), f'.{digits}g')


def check_manifest(root: Path):
    p=root/'SHA256SUMS.txt'
    assert p.is_file(), 'SHA256SUMS.txt ausente'
    passed=[]; errors=[]
    for line in p.read_text(encoding='utf-8').splitlines():
        if not line.strip() or line.lstrip().startswith('#'): continue
        match=re.fullmatch(r'([a-fA-F0-9]{64})  (.+)',line)
        if not match:
            errors.append(f'linha manifest invalida: {line[:60]}');continue
        expected, name=match.groups()
        dest=(root/name).resolve()
        if not dest.is_relative_to(root.resolve()):
            errors.append('caminho fora do diretorio');continue
        if not dest.is_file():
            errors.append(f'arquivo ausente: {name}');continue
        if sha(dest.read_bytes())!=expected.lower():
            errors.append(f'hash divergente: {name}');continue
        passed.append(name)
    if len(passed)<20: errors.append(f'menos de 20 arquivos verificados: {len(passed)}')
    return {'ok':not errors, 'n_verified':len(passed),'errors':errors, 'paths':passed}


def audit_math(doc, v13):
    """Somente dados publicados e aritmética exata. Sem valores float."""
    gates={}
    def gate(name,ok,detail): gates[name]={'status':'PASS' if ok else 'FAIL','detail':detail}
    gate('v14_schema_counts',doc.get('schema')==V14_SCHEMA and doc.get('point_counts')=={'G11':9524,'G22':6230,'total':15754} and doc.get('rounding')=='outward',
         'schema V14, 15754 pontos e arredondamento outward declarados')
    gate('v13_analytic_premises_status', v13.get('all_gates_passed') is True and len(v13.get('gates',{}))==10 and all(v is True for v in v13.get('gates',{}).values()),
         'V13 conferiu cotas numéricas, condicionadas aos lemas analíticos enumerados')
    a={}; band_ok=True; comp_ok=True; margins={}; proofs={}
    for k, expected in [('G11',9524),('G22',6230)]:
        ent=doc.get('entries',{}).get(k,{})
        try:
            lo,hi=fr(ent['finite_sum_lower']),fr(ent['finite_sum_upper'])
            q,t=fr(ent['quadrature_bound_declared_exact']),fr(ent['tail_bound_declared_exact'])
            clo,chi=map(fr,ent['e7_composite_interval_exact'])
            bands=ent['bands']
            # Check exact aggregation and location of every band; do not assume B3 per-node values.
            sumlo=sum((fr(b['sum'][0]) for b in bands),Q(0))
            sumhi=sum((fr(b['sum'][1]) for b in bands),Q(0))
            cells=sum(b['count'] for b in bands)
            boundaries=[fr(b['from']) for b in bands]+[fr(bands[-1]['to'])]
            continuous=all(fr(bands[i]['to'])==fr(bands[i+1]['from']) for i in range(len(bands)-1))
            thisband=(sumlo==lo and sumhi==hi and cells==expected and continuous and boundaries[0]==0 and boundaries[-1]==160)
            band_ok &= thisband
            thiscomp=(clo==lo-q-t and chi==hi+q+t and clo<=chi and q>=0 and t>=0)
            comp_ok &= thiscomp
            a[k]={'lo':clo,'hi':chi}
            margins[k]=-chi
            proofs[k]={
              'band_count':len(bands), 'cells':cells, 'all_band_exact_sums_match':thisband,
              'composite_recomposed_exactly':thiscomp,
              'finite_interval':[str(lo),str(hi)], 'quadrature_bound':str(q),
              'tail_bound':str(t),'composite_interval':[str(clo),str(chi)],
              'composite_approx':[dec(clo),dec(chi)],'margin_to_zero_exact':str(-chi),
              'margin_to_zero_approx':dec(-chi),
              'bounds_analytic_status':'CONDITIONAL_ON_E7_V13_ANALYTIC_LEMMAS'}
        except (KeyError,ValueError,ZeroDivisionError,IndexError,TypeError) as exc:
            band_ok=False;comp_ok=False;proofs[k]={'error':repr(exc)}
    gate('v14_exact_band_aggregation',band_ok,'11 bandas contiguas; somas racionais exatas; 15754 celulas')
    gate('v14_exact_composite_intervals',comp_ok,'Gkk=intervalo B5 +- quadratura E2/B3 +- cauda Proposicao 6')
    signs=all(k in margins and margins[k]>0 and a[k]['lo']<=a[k]['hi'] for k in ('G11','G22'))
    gate('strict_negative_diagonals',signs,'Extremos superiores exatos de G11 e G22 sao negativos')
    if signs:
        m11,m22=margins['G11'],margins['G22']
        det_margin=m11*m22
        # Even if G12 is not exactly zero, any symmetric |G12|<=1/10 suffices.
        delta=Q(1,10)
        robust=(delta*delta<det_margin)
        gate('offdiag_robustness_0p1',robust,'Se |G12|<=1/10, (G11*G22-G12^2)>0 por menores principais')
        proofs['matrix']={
            'G12_phase':'cos((1-2)*pi/2)=0 (identidade analitica assumida do E7)',
            'diagonal_negative':'G11.upper<0; G22.upper<0',
            'lambda_max_upper_exact':str(-min(m11,m22)),
            'lambda_max_upper_approx':dec(-min(m11,m22)),
            'determinant_lower_if_g12_zero_exact':str(det_margin),
            'determinant_lower_if_g12_zero_approx':dec(det_margin),
            'tested_offdiagonal_error_radius_exact':str(delta),
            'determinant_lower_even_if_abs_g12_le_1_10_exact':str(det_margin-delta*delta),
            'determinant_lower_even_if_abs_g12_le_1_10_approx':dec(det_margin-delta*delta),
            'interpretation':'NEGATIVE_DEFINITE_ARCHIMEDEAN_2X2_CONDITIONAL'}
    else:
        gate('offdiag_robustness_0p1',False,'Sem diagonais estritamente negativas, robustez nao estabelecida')
    return gates,proofs


def selftest(base, v13):
    """Cada mutação deve provocar a falha de um gate especificamente escolhido;
    nenhuma falha pré-existente serve como detecção espúria."""
    candidates=[]
    def test(name, mutate, target):
        d=copy.deepcopy(base); vv=copy.deepcopy(v13)
        try:
            mutate(d,vv)
            gates,_=audit_math(d,vv)
            detected=gates[target]['status']=='FAIL'
        except Exception: detected=False
        candidates.append({'mutation':name,'expected_gate':target,'detected':detected})
    test('upper_G11_positive',lambda d,v: d['entries']['G11']['e7_composite_interval_exact'].__setitem__(1,'1/10'),'v14_exact_composite_intervals')
    test('upper_G22_positive',lambda d,v: d['entries']['G22']['e7_composite_interval_exact'].__setitem__(1,'1/10'),'v14_exact_composite_intervals')
    test('band_sum_tampered',lambda d,v: d['entries']['G11']['bands'][0]['sum'].__setitem__(0,'0'),'v14_exact_band_aggregation')
    test('band_end_shifted',lambda d,v: d['entries']['G22']['bands'][0].__setitem__('to','13'),'v14_exact_band_aggregation')
    test('quadrature_tampered',lambda d,v: d['entries']['G11'].__setitem__('quadrature_bound_declared_exact','1/10'),'v14_exact_composite_intervals')
    test('v13_false_analytic_gate',lambda d,v: v['gates'].__setitem__(next(iter(v['gates'])),False),'v13_analytic_premises_status')
    test('changed_total_count',lambda d,v: d['point_counts'].__setitem__('total',15753),'v14_schema_counts')
    return candidates


def run(root, do_selftest=True):
    manifest=check_manifest(root)
    d=json.loads((root/'RESULTADO_E8K_V14.json').read_text(encoding='utf-8'))
    v13=json.loads((root/'RESULTADO_E8K_V13.json').read_text(encoding='utf-8'))
    zip_path=root/'NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip'
    integrity={'e7_zip_sha256':sha(zip_path.read_bytes()), 'b5_source_sha256':None}
    with ZipFile(zip_path) as z:
        integrity['b5_source_sha256']=sha(z.read('fontes/AUDITOR_PONTOS_B5.py'))
    integrity_ok=(integrity['e7_zip_sha256']==E7_SHA and integrity['b5_source_sha256']==B5_SHA and
       d.get('source_integrity',{}).get('e7_zip_sha256')==E7_SHA and
       d.get('source_integrity',{}).get('b5_source_sha256')==B5_SHA)
    g,proofs=audit_math(d,v13)
    g['manifest_all_files_intact']={'status':'PASS' if manifest['ok'] else 'FAIL','detail':f"{manifest['n_verified']} arquivos comparados; {len(manifest['errors'])} problemas"}
    g['e7_b5_exact_sha256']={'status':'PASS' if integrity_ok else 'FAIL', 'detail':'SHA256 do ZIP E7, fonte B5 e declaracoes V14'}
    for k in ('G11','G22'):
        ent=d['entries'][k]
        path=root/('b5/PONTOS_G11_B5.jsonl' if k=='G11' else 'b5/PONTOS_G22_B5.jsonl')
        ok=sha(path.read_bytes())==ent['jsonl_sha256']
        g[f'{k}_jsonl_sha256']={'status':'PASS' if ok else 'FAIL', 'detail':'hash JSONL reproduzivel; confere o arquivo, nao os teoremas'}
    muts=selftest(d,v13) if do_selftest else []
    if do_selftest:
        g['negative_controls_isolated']={'status':'PASS' if all(x['detected'] for x in muts) and len(muts)==7 else 'FAIL',
          'detail':f"{sum(x['detected'] for x in muts)}/{len(muts)} mutacoes isoladas detectadas"}
    ok=all(x['status']=='PASS' for x in g.values())
    return {'schema':SCHEMA,'phase':'E8K_V15','status':'CONDITIONAL_2X2_ARCHIMEDEAN_NEGATIVE_DEFINITE' if ok else 'FAIL_INCOMPLETE',
      'all_gates_passed':ok,'gates':dict(sorted(g.items())), 'exact_witnesses':proofs,
      'negative_controls':muts,'input_sha256':integrity,
      'manifest':{'verified_file_count':manifest['n_verified'],'problems':manifest['errors']},
      'theorem_dependency':'A interpretacao como integral real requer identidades e cotas analiticas E7/V13, incluindo normalizacao e paridade do F_k; nao demonstradas inteiramente pelo programa.',
      'not_certified':['Lean 4 nao compilado neste ambiente','Fórmula explícita de Weil global / termos primos / polos','RH ou positividade em dimensoes arbitrarias','Cotas analiticas integralmente mecanizadas','Autenticidade de hashes sem assinatura externa']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('.'))
    parser.add_argument('--out',type=Path,default=Path('v15/RESULTADO_E8K_V15.json'))
    parser.add_argument('--no-selftest',action='store_true')
    args=parser.parse_args()
    try: report=run(args.root.resolve(),not args.no_selftest)
    except Exception as exc:
        print('ERRO INTEGRACAO:', repr(exc),file=sys.stderr);return 3
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')
    print('V15:',report['status']); print('Gates:',len(report['gates']),'PASS:',sum(v['status']=='PASS' for v in report['gates'].values()))
    print('G11/G22 margem:',*[f"{k}:{report['exact_witnesses'][k].get('margin_to_zero_approx')}" for k in ('G11','G22')])
    return 0 if report['all_gates_passed'] else 2

if __name__=='__main__':raise SystemExit(main())
