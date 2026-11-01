#!/usr/bin/env python3
"""E8-K V14: replay de testemunhas V13+ e oraculo racional independente estratificado.

NÃO executa a fonte B5: lê seus JSONL e confere todos os pontos, somas por banda,
SHA, proveniência, pesos e intervalos compostos. Amostra 128 integrandos usando
séries racionais independentes para pi, digamma, log, seno, cos e sinc.
Não prova os majorantes analíticos globais de quadratura/cauda nem RH.
"""
import argparse,hashlib,json,zipfile
from fractions import Fraction as F
from pathlib import Path
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
import ORACULO_RACIONAL_V14 as oracle
ZIP_NAME='NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip'
ZIP_SHA='4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7'
B5_SRC='a1753d3010e34f29371a9706cd6bb05b3b0c8491ca3cafb44bf98b7cd53d3df9'
B3_NAMES={'G11':'RESULTADO_B3_G11_WINDOWS.json','G22':'RESULTADO_G12_G22_WINDOWS.json'}
SCHEMA='NUMERIA.E8K.V13.FULL_POINT_REPLAY.1'

def sha(b):return hashlib.sha256(b).hexdigest()
def get_iv(x):
    a,b=map(F,x)
    if a>b:raise AssertionError('intervalo invertido')
    return a,b

def outward(x,digits=45,rounding=ROUND_FLOOR):
    x=F(x)
    with localcontext() as ctx:
        ctx.prec=digits+60
        val=Decimal(x.numerator)/Decimal(x.denominator)
        return format(val.quantize(Decimal(1).scaleb(-digits),rounding=rounding),'f')

def pair_display(lo,hi):
    return [outward(lo,rounding=ROUND_FLOOR),outward(hi,rounding=ROUND_CEILING)]

def metrics(a,b):
    """Somente UM par de intervalos de SOMA FINITA. Nunca percentuais nodais."""
    al,ah=a;bl,bh=b
    overlap=max(F(0),min(ah,bh)-max(al,bl))
    gap=max(F(0),bl-ah,al-bh)
    ab=bl<=al and ah<=bh
    ba=al<=bl and bh<=ah
    return {'scope':'aggregate_finite_sum_only','comparison_count':1,
      'nested_fraction':1.0 if (ab or ba) else 0.0,
      'b5_in_b3_fraction':1.0 if ab else 0.0,
      'b3_in_b5_fraction':1.0 if ba else 0.0,
      'disjoint_count':1 if gap>0 else 0,
      'max_gap':str(gap),'overlap_width':str(overlap),
      'b5_minus_b3_lower':str(al-bl),'b5_minus_b3_upper':str(ah-bh),
      'width_ratio_max':str(max((ah-al)/(bh-bl),(bh-bl)/(ah-al))) if ah>al and bh>bl else None,
      'per_point_metrics':'NOT_AVAILABLE_NO_B3_NODE_INTERVALS'}

def choose_indices(n, band_i):
    want=12 if band_i<7 else 11  # 7*12 + 4*11 = 128
    return sorted({((2*j+1)*n)//(2*want) for j in range(want)})

def build(root):
    zbytes=(root/ZIP_NAME).read_bytes()
    if sha(zbytes)!=ZIP_SHA:raise AssertionError('SHA E7 nao coincide')
    with zipfile.ZipFile(root/ZIP_NAME) as z:
        manifest=json.loads(z.read('MANIFESTO_SHA256.json'))['sha256']
        for fname,expected in manifest.items():
            if sha(z.read(fname))!=expected:raise AssertionError('manifesto E7 incorreto: '+fname)
        if sha(z.read('fontes/AUDITOR_PONTOS_B5.py'))!=B5_SRC:raise AssertionError('SHA fonte B5')
        b5=json.loads(z.read('b5/RESULTADO_B5.json'))
        e7=json.loads(z.read('RESULTADO_INTEGRAIS_E7.json'))
        if b5['source_sha256']!=B5_SRC or int(b5['fixed_point_scale'])!=10**42:raise AssertionError('escala/fonte B5')
        const=oracle.make_constants()
        entries={};sample_details=[];total=0;band_idx=0
        for item in b5['interval_results']:
            comp=item['component'];k=1 if comp=='G11' else 2
            raw=(root/'b5'/f'PONTOS_{comp}_B5.jsonl').read_bytes()
            if raw!=z.read('b5/'+f'PONTOS_{comp}_B5.jsonl'):raise AssertionError('JSONL raiz != ZIP')
            point_manifest=next(a for a in b5['point_files'] if a['file']==f'PONTOS_{comp}_B5.jsonl')
            if sha(raw)!=point_manifest['sha256']:raise AssertionError('hash JSONL errado')
            b3raw=(root/'proveniencia_b3'/B3_NAMES[comp]).read_bytes()
            if sha(b3raw)!=item['source_input_sha256']:raise AssertionError('SHA B3 original nao coincide')
            b3cor=json.loads(z.read('corrigidos/'+('RESULTADO_B3_G11_E2_CORRIGIDO.json' if k==1 else 'RESULTADO_B3_G22_E2_CORRIGIDO.json')))
            if [(v['from'],v['to'],v['cells']) for v in json.loads(b3raw)['mesh']]!=[(v['from'],v['to'],v['cells']) for v in b3cor['mesh']]:
                raise AssertionError('malha B3 original/corrigida nao coincide')
            lines=raw.splitlines()
            if len(lines)!=item['points']:raise AssertionError('numero de pontos nao coincide')
            bands=[];offset=0;all_lo=F(0);all_hi=F(0);prev=F(0);maxpt=F(0)
            for bb in item['bands']:
                a,b,n=F(bb['from']),F(bb['to']),int(bb['cells'])
                if a!=prev or b<=a:raise AssertionError('banda nao contigua')
                step=(b-a)/n;blo=F(0);bhi=F(0)
                picked=choose_indices(n,band_idx)
                if len(picked)!=(12 if band_idx<7 else 11):raise AssertionError('amostragem nao unica')
                for j in range(n):
                    line=json.loads(lines[offset+j]);r=a+(2*j+1)*step/2
                    if F(line['r'])!=r:raise AssertionError('ponto medio errado')
                    i1,i2=get_iv(line['integrand']);w1,w2=get_iv(line['weighted'])
                    if w1>i1*step or w2<i2*step:raise AssertionError('peso nao encerra intervalo')
                    if maxpt<i2-i1:maxpt=i2-i1
                    blo+=w1;bhi+=w2
                    if j in picked:
                        independent,om=oracle.evaluate_independent(r,k,const)
                        passed=i1<=independent[0] and independent[1]<=i2
                        sample_details.append({'component':comp,'band':[str(a),str(b)],'global_index':offset+j,
                            'r':str(r),'exported_integrand':[str(i1),str(i2)],
                            'independent_rational_integrand':[str(independent[0]),str(independent[1])],
                            'independent_rational_omega':[str(om[0]),str(om[1])],
                            'encloses_independent':passed})
                        if not passed:raise AssertionError(f'oraculo fora do B5 {comp}/{offset+j}')
                if get_iv(bb['sum_interval'])!=(blo,bhi):raise AssertionError('soma por banda divergente')
                bands.append({'from':str(a),'to':str(b),'count':n,'cells':n,'sum':[str(blo),str(bhi)],
                    'sum_match':True,'rounding':'exact_rational'})
                all_lo+=blo;all_hi+=bhi;offset+=n;prev=b;band_idx+=1
            if prev!=160 or offset!=len(lines):raise AssertionError('bandas nao cobrem 0..160')
            if (all_lo,all_hi)!=get_iv(item['certified_finite_midpoint_sum']):raise AssertionError('soma global divergente')
            ev=next(x for x in e7['results'] if x['entry']==comp)
            q=F(ev['quadrature_error_b3']);tail=F(ev['tail_proposition6'])
            if (all_lo,all_hi)!=get_iv(ev['sum_b5']):raise AssertionError('soma E7 nao coincide')
            if get_iv(ev['combined_interval'])!=(all_lo-q-tail,all_hi+q+tail):raise AssertionError('composicao E7 incorreta')
            if ev['point_sha256']!=sha(raw):raise AssertionError('SHA pontos E7')
            b3interval=get_iv(item['b3_midpoint_interval']);met=metrics((all_lo,all_hi),b3interval)
            if bool(met['b5_in_b3_fraction'])!=item['b5_enclosed_in_b3'] or bool(met['disjoint_count']==0)!=item['b5_intersects_b3']:
                raise AssertionError('metrica B3 B5 nao coincide')
            entries[comp]={
                'points_regenerated_and_compared':offset,'integrands_reproduced':offset,'weighted_products_reproduced':offset,
                'bands_reproduced':len(bands),'bands':bands,'jsonl_sha256':sha(raw),
                'original_b3_input_sha256_recovered_and_verified':True,
                'original_b3_input_sha256':sha(b3raw),'original_b3_input_bytes':len(b3raw),
                'provenance':{'original_b3_input_sha256':sha(b3raw),'corrected_b3_same_mesh':True},
                'b5_enclosed_in_b3':item['b5_enclosed_in_b3'],
                'b5_intersects_b3':item['b5_intersects_b3'],
                'b5_b3_metrics':met,
                'finite_sum_lower':str(all_lo),'finite_sum_upper':str(all_hi),
                'finite_sum_interval_outward':pair_display(all_lo,all_hi),
                'finite_sum_width_approx':outward(all_hi-all_lo,digits=45,rounding=ROUND_CEILING),
                'max_point_interval_width_approx':outward(maxpt,digits=45,rounding=ROUND_CEILING),
                'quadrature_bound_declared_exact':str(q), 'tail_bound_declared_exact':str(tail),
                'quadrature_bound_declared_approx':outward(q,digits=45,rounding=ROUND_CEILING),
                'tail_bound_declared_approx':outward(tail,digits=45,rounding=ROUND_CEILING),
                'e7_composite_interval_exact':[str(all_lo-q-tail),str(all_hi+q+tail)],
                'e7_composite_interval_approx':pair_display(all_lo-q-tail,all_hi+q+tail),
                'e7_composite_formula_exact_match':True,
                'band_sum_rounding':'exact_rational_no_interband_rounding'
            }
            total+=offset
        if total!=15754 or len(sample_details)!=128:raise AssertionError('contagem global/amostra nao coincide')
        details={'schema':'NUMERIA.E8K.V14.ORACLE_DETAILS.1',
          'method':'Machin/pi; Euler-Maclaurin N128 gamma; atanh/log; Taylor trig sinc; 128-term digamma series and integral tail',
          'scope':'128 pointwise B5 exported integrands; exact rational interval containment',
          'samples':sample_details}
        detail_bytes=(json.dumps(details,ensure_ascii=False,sort_keys=True,indent=2)+'\n').encode()
        (root/'oraculo'/'AMOSTRA_128.json').write_bytes(detail_bytes)
        replay_path=root/'REPLAY_E8K_V12_INTEGRAL.json'
        if not replay_path.is_file(): raise AssertionError('replay V12 integral e obrigatorio: execute VERIFICAR_E8K_V12.py antes')
        replay=json.loads(replay_path.read_text(encoding='utf8'))
        if not (replay.get('all_point_exports_reproduced') is True and replay.get('point_counts',{}).get('total')==total and not replay.get('provenance_gaps')):
            raise AssertionError('replay V12 nao esta aprovado')
        for comp in ('G11','G22'):
            if replay['entries'][comp]['finite_sum_lower']!=entries[comp]['finite_sum_lower'] or replay['entries'][comp]['finite_sum_upper']!=entries[comp]['finite_sum_upper']:
                raise AssertionError('replay V12 difere da soma V14')
        report={
          'schema':SCHEMA,'phase':'E8K_V14','status':'FULL_EXPORTED_POINT_REPLAY_AND_128_INDEPENDENT_RATIONAL_ENCLOSURES_PASSED',
          'rounding':'outward','source_integrity':{'e7_zip_sha256':ZIP_SHA,'b5_source_sha256':B5_SRC,
              'b5_result_sha256':sha((root/'b5/RESULTADO_B5.json').read_bytes()),
              'e7_manifest_sha256':sha((root/'MANIFESTO_SHA256.json').read_bytes())},
          'full_replay_v12':{'status':replay['status'],'sha256':sha(replay_path.read_bytes()),
              'path':'REPLAY_E8K_V12_INTEGRAL.json','method':'reexecuta fonte B5 sha-fixada; distinto do oraculo V14'},
          'point_counts':{'G11':len((root/'b5/PONTOS_G11_B5.jsonl').read_bytes().splitlines()),
                          'G22':len((root/'b5/PONTOS_G22_B5.jsonl').read_bytes().splitlines()),'total':total},
          'entries':entries,'exact_integrand_enclosures_replayed':total,'exact_weighted_enclosures_replayed':total,
          'quadrature_and_tail_formulas_arithmetically_recomposed':2,'all_point_exports_reproduced':True,
          'independent_oracle_sample':{'implementation':'Independent rational interval oracle V14: Machin, Euler-Maclaurin, digamma tail, Taylor trigonometric functions',
             'sample_size':len(sample_details),'selection':'stratified per B5 band: 12 for first seven; 11 for other four',
             'mismatches':sum(1 for e in sample_details if not e['encloses_independent']),
             'details_sha256':sha(detail_bytes),'details_path':'oraculo/AMOSTRA_128.json',
             'qualification':'independent formulas and implementation; relies on stated classical analytic tail theorems'},
          'negative_checker_controls':{'execution':'CONFERIR_E8K_V14.py --selftest; aguardando relatorio de gates'},'provenance_gaps':[],
          'not_certified':['Global analytic validity of B3 midpoint second derivative bounds',
              'Independent formal proof of Weil integral tail R>160',
              'Formal Lean compilation and complete B5 fixed-point correspondence',
              'Third-party authenticity of unsigned hashes','Weil explicit formula globally','Riemann Hypothesis'],
          'warning':'hashes sem assinatura são testes de integridade relativa, não autenticação por terceiros'
        }
        return report

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'));a=p.parse_args()
    root=a.root.resolve()
    doc=build(root)
    path=root/'RESULTADO_E8K_V14.json'
    path.write_text(json.dumps(doc,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':doc['status'],'points':doc['point_counts']['total'],
         'oracle_samples':doc['independent_oracle_sample']['sample_size'],
         'b3_b5_scope':'aggregate_only'},ensure_ascii=False))
if __name__=='__main__':main()
