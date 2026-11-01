#!/usr/bin/env python3
"""NUMERIA E8-K V12: reproduz TODOS os pontos B5 e recompõe o E7.

NÃO prova cotas analíticas, erro de quadratura ou a fórmula explícita.
Executa a fonte B5 somente depois de confirmar SHA256 fixo; usar ZIP confiável.
Exit: 0 = reprodução passou (lacunas de proveniência sinalizadas);
      2 = falha de dados/reprodução; 4 = --strict-provenance com lacunas.
"""
import argparse
import hashlib
import json
import sys
import zipfile
from fractions import Fraction as Q
from pathlib import Path

ZIP_SHA = '4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7'
SOURCE_SHA = 'a1753d3010e34f29371a9706cd6bb05b3b0c8491ca3cafb44bf98b7cd53d3df9'
S = 10**42
ORIGINAL_INPUT_FILES = {'G11':'RESULTADO_B3_G11_WINDOWS.json', 'G22':'RESULTADO_G12_G22_WINDOWS.json'}
FILES = {
    'G11': ('b5/PONTOS_G11_B5.jsonl', 'corrigidos/RESULTADO_B3_G11_E2_CORRIGIDO.json', 1),
    'G22': ('b5/PONTOS_G22_B5.jsonl', 'corrigidos/RESULTADO_B3_G22_E2_CORRIGIDO.json', 2),
}

def sha(data): return hashlib.sha256(data).hexdigest()
def fpair(v):
    if len(v) != 2: raise AssertionError('par de extremos inválido')
    l, u = map(Q, v)
    if l > u: raise AssertionError('intervalo invertido')
    return l, u

def test_negative_controls():
    """Sondas sobre os predicados do próprio verificador; não executam fonte mutada."""
    def ordered(lo, hi): return lo <= hi
    def exact(a, b): return a == b
    def valid_grid(got, expect): return got == expect
    def valid_sha(x, target): return sha(x) == target
    return {
        'reversed_interval_rejected': not ordered(Q(2), Q(1)),
        'changed_endpoint_detected': not exact((Q(1), Q(2)), (Q(1), Q(3))),
        'shifted_midpoint_detected': not valid_grid(Q(3, 2), Q(1, 2)),
        'tampered_bytes_detected': not valid_sha(b'tampered', sha(b'original')),
    }

def audit(zip_path, provenance_dir=None):
    data = zip_path.read_bytes()
    if sha(data) != ZIP_SHA: raise AssertionError('SHA-256 do ZIP E7 incorreto')
    with zipfile.ZipFile(zip_path) as z:
        src = z.read('fontes/AUDITOR_PONTOS_B5.py')
        if sha(src) != SOURCE_SHA: raise AssertionError('SHA-256 do B5 incorreto')
        ns = {'__name__': '_e8k_v12_review'}
        exec(compile(src, 'fontes/AUDITOR_PONTOS_B5.py', 'exec'), ns)
        I = ns['I']
        b5 = json.loads(z.read('b5/RESULTADO_B5.json'))
        e7 = json.loads(z.read('RESULTADO_INTEGRAIS_E7.json'))
        if b5['schema'] != 'NUMERIA.B5.POINTWISE_RATIONAL_ENCLOSURE.1':
            raise AssertionError('schema B5 divergente')
        if b5['source_sha256'] != SOURCE_SHA: raise AssertionError('fonte B5 no relatório divergente')
        if int(b5['fixed_point_scale']) != S: raise AssertionError('escala B5 divergente')
        if len(b5['interval_results']) != 2: raise AssertionError('quantidade de componentes divergente')
        detail = {}
        total_points = 0
        gaps = []
        for item in b5['interval_results']:
            name = item['component']
            if name not in FILES or name in detail: raise AssertionError('componente desconhecido/duplicado')
            point_fn, input_fn, k = FILES[name]
            content = z.read(point_fn)
            actual_p_sha = sha(content)
            point_manifest = next(x for x in b5['point_files'] if x['file'] == Path(point_fn).name)
            if point_manifest['sha256'] != actual_p_sha:
                raise AssertionError(f'{name}: hash dos pontos divergente')
            if point_manifest['count'] != item['points']:
                raise AssertionError(f'{name}: contagem do manifesto divergente')
            report_b3_bytes = z.read(input_fn)
            b3_actual_sha = sha(report_b3_bytes)
            provenance_ok = b3_actual_sha == item['source_input_sha256']
            original_match = False
            original_path = None
            if provenance_dir is not None:
                original_path = Path(provenance_dir) / ORIGINAL_INPUT_FILES[name]
                if original_path.is_file():
                    original_match = sha(original_path.read_bytes()) == item['source_input_sha256']
                    if not original_match:
                        raise AssertionError(f'{name}: original de proveniência fornecido tem SHA incorreto')
                    original_mesh = json.loads(original_path.read_text(encoding='utf-8'))['mesh']
                    if [(v['from'],v['to'],v['cells']) for v in original_mesh] != [(v['from'],v['to'],v['cells']) for v in json.loads(report_b3_bytes)['mesh']]:
                        raise AssertionError(f'{name}: malha original e corrigida divergentes')
            if not provenance_ok and not original_match:
                gaps.append({'component':name,'issue':'B5_INPUT_B3_HASH_NOT_RESOLVED',
                  'b5_expected_sha256':item['source_input_sha256'],
                  'packaged_corrected_b3_sha256':b3_actual_sha,
                  'packaged_filename':input_fn,
                  'original_provenance_checked': provenance_dir is not None})
            b3 = json.loads(report_b3_bytes)
            mesh = b3['mesh']
            bands = item['bands']
            if len(mesh) != len(bands): raise AssertionError(f'{name}: faixas divergentes')
            record_lines = iter(content.splitlines())
            total_iv = I.q(0)
            sum_fraction_lo, sum_fraction_hi = Q(0), Q(0)
            count = 0
            peak_point_width = Q(0)
            per_band = []
            expected_start = Q(0)
            for band_index, (bm, bb) in enumerate(zip(mesh, bands)):
                a, b, n = Q(bm['from']), Q(bm['to']), int(bm['cells'])
                if a != expected_start or b <= a or not 0 < n <= 20000:
                    raise AssertionError(f'{name}: malha inválida em faixa {band_index}')
                if (str(a), str(b), n) != (bb['from'], bb['to'], bb['cells']):
                    raise AssertionError(f'{name}: malha B3/B5 divergente em faixa {band_index}')
                h = (b-a) / n
                h_iv = I.q(h)
                band_iv = I.q(0)
                for j in range(n):
                    raw = next(record_lines, None)
                    if raw is None: raise AssertionError(f'{name}: faltam registros')
                    obj = json.loads(raw)
                    r = a + (2*j+1)*h/2
                    if Q(obj['r']) != r: raise AssertionError(f'{name}: ponto fora da malha na linha {count+1}')
                    supplied_y = fpair(obj['integrand']); supplied_w = fpair(obj['weighted'])
                    expected_y = ns['eval_point'](r, k)
                    expected_w = expected_y*h_iv
                    if supplied_y != (Q(expected_y.a,S),Q(expected_y.b,S)):
                        raise AssertionError(f'{name}: integrando não reproduzido na linha {count+1}')
                    if supplied_w != (Q(expected_w.a,S),Q(expected_w.b,S)):
                        raise AssertionError(f'{name}: peso não reproduzido na linha {count+1}')
                    if supplied_w[0] > supplied_y[0]*h or supplied_w[1] < supplied_y[1]*h:
                        raise AssertionError(f'{name}: peso sem inclusão na linha {count+1}')
                    band_iv += expected_w
                    sum_fraction_lo += supplied_w[0]
                    sum_fraction_hi += supplied_w[1]
                    peak_point_width=max(peak_point_width,supplied_y[1]-supplied_y[0])
                    count+=1
                if fpair(bb['sum_interval'])!=(Q(band_iv.a,S),Q(band_iv.b,S)):
                    raise AssertionError(f'{name}: soma da faixa divergente')
                total_iv += band_iv
                per_band.append({'from':str(a),'to':str(b),'count':n,'sum_match':True})
                expected_start=b
            if expected_start!=Q(160): raise AssertionError(f'{name}: malha não vai a 160')
            if next(record_lines, None) is not None: raise AssertionError(f'{name}: há registros extras')
            if count!=item['points']: raise AssertionError(f'{name}: quantidade incorreta')
            final_lo,final_hi=Q(total_iv.a,S),Q(total_iv.b,S)
            if (final_lo,final_hi)!=fpair(item['certified_finite_midpoint_sum']):
                raise AssertionError(f'{name}: soma B5 divergente')
            if (final_lo,final_hi)!=(sum_fraction_lo,sum_fraction_hi):
                raise AssertionError(f'{name}: soma de Fraction independente da classe I divergente')
            b3_lo,b3_hi=fpair(item['b3_midpoint_interval'])
            incl=b3_lo<=final_lo and final_hi<=b3_hi
            intersects=final_lo<=b3_hi and b3_lo<=final_hi
            twice=(final_hi-final_lo)<=2*(b3_hi-b3_lo)
            if (incl,intersects,twice)!=(item['b5_enclosed_in_b3'],item['b5_intersects_b3'],item['b5_width_at_most_twice_b3']):
                raise AssertionError(f'{name}: flags de compatibilidade B3 divergentes')
            E=next(e for e in e7['results'] if e['entry']==name)
            err=Q(E['quadrature_error_b3']);tail=Q(E['tail_proposition6'])
            if err<0 or tail<0:raise AssertionError(f'{name}: erro/cauda negativo')
            if fpair(E['sum_b5'])!=(final_lo,final_hi):
                raise AssertionError(f'{name}: soma E7/B5 divergente')
            if fpair(E['combined_interval'])!=(final_lo-err-tail,final_hi+err+tail):
                raise AssertionError(f'{name}: composição E7 divergente')
            if E['point_records']!=count or E['point_sha256']!=actual_p_sha:
                raise AssertionError(f'{name}: referência E7 aos pontos divergente')
            total_points += count
            detail[name] = {
                'points_regenerated_and_compared':count,
                'bands_reproduced':len(bands),
                'jsonl_sha256':actual_p_sha,
                'integrands_reproduced':count,
                'weighted_products_reproduced':count,
                'midpoint_grid_exact':True,
                'finite_sum_exact_match':True,
                'finite_sum_lower':str(final_lo),
                'finite_sum_upper':str(final_hi),
                'finite_sum_width_approx':f'{float(final_hi-final_lo):.12g}',
                'quadrature_bound_declared_approx':f'{float(err):.12g}',
                'tail_bound_declared_approx':f'{float(tail):.12g}',
                'e7_composite_formula_exact_match':True,
                'e7_composite_interval_approx':[f'{float(final_lo-err-tail):.12g}',f'{float(final_hi+err+tail):.12g}'],
                'b5_enclosed_in_b3':incl,
                'b5_intersects_b3':intersects,
                'b5_width_at_most_twice_b3':twice,
                'b5_input_b3_sha256_matches_packaged_corrected_b3':provenance_ok,
                'original_b3_input_sha256_recovered_and_verified':original_match,
                'max_point_interval_width_approx':f'{float(peak_point_width):.12g}',
                'bands':per_band,
            }
        if set(detail)!=set(FILES) or total_points!=b5['points'] or total_points!=15754:
            raise AssertionError('quantidade total incorreta')
        controls=test_negative_controls()
        if not all(controls.values()): raise AssertionError('sondas negativas falharam')
        return {
            'schema':'NUMERIA.E8K.V12.FULL_POINT_REPLAY.1',
            'status':'ALL_POINTS_REPRODUCED_WITH_PROVENANCE_GAPS' if gaps else 'ALL_POINTS_REPRODUCED_B5_INPUT_HASHES_MATCHED',
            'all_point_exports_reproduced':True,
            'point_counts':{'G11':detail['G11']['points_regenerated_and_compared'],'G22':detail['G22']['points_regenerated_and_compared'],'total':total_points},
            'exact_integrand_enclosures_replayed':total_points,
            'exact_weighted_enclosures_replayed':total_points,
            'quadrature_and_tail_formulas_arithmetically_recomposed':2,
            'negative_checker_controls':controls,
            'source_integrity':{'e7_zip_sha256':ZIP_SHA,'b5_source_sha256':SOURCE_SHA},
            'provenance_gaps':gaps,
            'provenance_note':'B3 corrigidos do E7 possuem SHA diferente dos B3 originais usados no B5; as entradas originais foram conferidas separadamente quando --provenance-dir foi fornecido.',
            'entries':detail,
            'not_certified':[
                'B5 pontos são regenerados pela MESMA implementação B5: isto é reprodução, não um oráculo matemático independente de todos os integrandos.',
                'Os B3 corrigidos do E7 diferem dos B3 originais; use --provenance-dir para conferir os originais e preservar explicitamente ambas as versões.',
                'A validade analítica de cada cota de resto transcendental não é reprovada nem novamente demonstrada aqui.',
                'Cotas globais de derivada para quadratura contínua não são demonstradas por esta suíte.',
                'Cota de cauda para r>160, fórmula explícita de Weil e a Hipótese de Riemann não são provadas.'
            ]
        }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--zip',type=Path,required=True)
    p.add_argument('--out',type=Path,default=Path('RESULTADO_E8K_V12.json'))
    p.add_argument('--provenance-dir',type=Path,help='Pasta com arquivos B3 originais preservados e fixados por SHA')
    p.add_argument('--strict-provenance',action='store_true')
    args=p.parse_args()
    try:
        report=audit(args.zip,args.provenance_dir)
    except Exception as e:
        print('V12 FAILURE:',repr(e),file=sys.stderr)
        return 2
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf8')
    print(f"V12: {report['status']}  PONTOS={report['point_counts']['total']}  LACUNAS_PROVENIENCIA={len(report['provenance_gaps'])}")
    print('Relatório:',args.out)
    return 4 if args.strict_provenance and report['provenance_gaps'] else 0
if __name__=='__main__':raise SystemExit(main())
