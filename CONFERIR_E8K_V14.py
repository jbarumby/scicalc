#!/usr/bin/env python3
"""Auditor independente de testemunhas da E8-K V14, fail-closed.

Lê fonte E7 de SHA fixo, JSONL e B3 originais; recalcula bandas, pesos,
metadados, métricas de intervalos de soma, decimais para fora e amostra
racional independente de 128 pontos. Não certifica teoremas globais de Weil.
"""
import argparse,copy,hashlib,json,zipfile
from fractions import Fraction as F
from pathlib import Path
import ORACULO_RACIONAL_V14 as ref
from GERAR_E8K_V14 import metrics,choose_indices,sha,get_iv,ZIP_NAME,ZIP_SHA,B5_SRC,SCHEMA,pair_display
NAMES={'G11':'RESULTADO_B3_G11_WINDOWS.json','G22':'RESULTADO_G12_G22_WINDOWS.json'}

def req(x,msg):
    if not x:raise AssertionError(msg)

def manifest(root):
    f=root/'SHA256SUMS.txt'
    seen={}
    for ln in f.read_text(encoding='utf8').splitlines():
        if not ln.strip() or ln.startswith('#'):continue
        a=ln.split('  ',1);req(len(a)==2 and len(a[0])==64,'SHA256SUMS formato')
        digest,name=a
        req(name not in seen,'SHA duplicado '+name)
        target=(root/name).resolve()
        req(target.is_relative_to(root) and target.is_file() and not target.is_symlink(),'caminho invalido '+name)
        req(sha(target.read_bytes())==digest,'hash diferente '+name)
        seen[name]=digest
    required=['RESULTADO_E8K_V14.json','REPLAY_E8K_V12_INTEGRAL.json',
       'RESULTADO_E8K_V13.json',ZIP_NAME,'MANIFESTO_SHA256.json',
       'b5/PONTOS_G11_B5.jsonl','b5/PONTOS_G22_B5.jsonl','b5/RESULTADO_B5.json',
       'proveniencia_b3/'+NAMES['G11'],'proveniencia_b3/'+NAMES['G22'],
       'oraculo/AMOSTRA_128.json','GERAR_E8K_V14.py','CONFERIR_E8K_V14.py','ORACULO_RACIONAL_V14.py']
    req(all(x in seen for x in required),'SHA256SUMS nao cobre arquivos essenciais')
    return {'checked':len(seen)}

def sources(root,doc):
    req(doc['schema']==SCHEMA and doc['phase']=='E8K_V14','schema/phase')
    req(doc['rounding']=='outward','rounding nao outward')
    zbytes=(root/ZIP_NAME).read_bytes();req(sha(zbytes)==ZIP_SHA==doc['source_integrity']['e7_zip_sha256'],'ZIP E7 SHA')
    with zipfile.ZipFile(root/ZIP_NAME) as z:
        req(sha(z.read('fontes/AUDITOR_PONTOS_B5.py'))==B5_SRC==doc['source_integrity']['b5_source_sha256'],'SHA B5')
        mapped=json.loads(z.read('MANIFESTO_SHA256.json'))['sha256']
        for fn,digest in mapped.items():req(sha(z.read(fn))==digest,'manifesto E7 '+fn)
        req((root/'MANIFESTO_SHA256.json').read_bytes()==z.read('MANIFESTO_SHA256.json'),'manifesto espelho')
        req((root/'b5/RESULTADO_B5.json').read_bytes()==z.read('b5/RESULTADO_B5.json'),'B5 espelho')
        b5=json.loads(z.read('b5/RESULTADO_B5.json'))
        for it in b5['interval_results']:
            n=it['component'];orig=(root/'proveniencia_b3'/NAMES[n]).read_bytes()
            req(sha(orig)==it['source_input_sha256']==doc['entries'][n]['original_b3_input_sha256'],'B3 hash original '+n)
            req(doc['entries'][n]['original_b3_input_bytes']==len(orig),'B3 bytes '+n)
            req(doc['entries'][n]['original_b3_input_sha256_recovered_and_verified'] is True,'B3 afirmacao '+n)
            req(json.loads(orig)['mesh']==json.loads(z.read('corrigidos/'+('RESULTADO_B3_G11_E2_CORRIGIDO.json' if n=='G11' else 'RESULTADO_B3_G22_E2_CORRIGIDO.json')))['mesh'],'malha B3 '+n)
        return {'manifest_entries':len(mapped),'b3_original_sha':{x['component']:x['source_input_sha256'] for x in b5['interval_results']}}

def bands(root,doc):
    with zipfile.ZipFile(root/ZIP_NAME) as z:
        b5=json.loads(z.read('b5/RESULTADO_B5.json'));e7=json.loads(z.read('RESULTADO_INTEGRAIS_E7.json'))
        totals={};bandnum=0
        for it in b5['interval_results']:
            n=it['component'];ent=doc['entries'][n];path='b5/PONTOS_'+n+'_B5.jsonl'
            raw=(root/path).read_bytes();req(raw==z.read(path),'JSONL ZIP raiz '+n)
            req(sha(raw)==ent['jsonl_sha256'],'JSONL SHA '+n)
            lines=raw.splitlines();ofs=0;slo=F(0);shi=F(0);prev=F(0)
            req(len(ent['bands'])==len(it['bands']),'num bandas '+n)
            for index,(a,b) in enumerate(zip(it['bands'],ent['bands'])):
                lo,hi=F(a['from']),F(a['to']);count=int(a['cells']);h=(hi-lo)/count
                req(lo==prev and hi>lo,'bandas nao contiguas')
                req((F(b['from']),F(b['to']),b['count'])==(lo,hi,count),'metadata banda')
                sumlo=F(0);sumhi=F(0)
                for j in range(count):
                    obj=json.loads(lines[ofs+j]);r=lo+(2*j+1)*h/2
                    req(F(obj['r'])==r,'ponto errado '+n)
                    al,au=get_iv(obj['integrand']);wl,wu=get_iv(obj['weighted'])
                    req(wl<=al*h and au*h<=wu,'peso nao cobre integrando '+n)
                    sumlo+=wl;sumhi+=wu
                req([str(sumlo),str(sumhi)]==b['sum'],'testemunha soma banda '+n)
                req((sumlo,sumhi)==get_iv(a['sum_interval']),'banda nao coincide B5 '+n)
                req(b['sum_match'] is True,'sum_match flag falsa')
                slo+=sumlo;shi+=sumhi;ofs+=count;prev=hi;bandnum+=1
            req(prev==160 and ofs==len(lines)==it['points'],'pontos e faixa inteira '+n)
            req((slo,shi)==get_iv(it['certified_finite_midpoint_sum']),'soma B5 '+n)
            req((str(slo),str(shi))==(ent['finite_sum_lower'],ent['finite_sum_upper']),'soma V14 '+n)
            e=next(x for x in e7['results'] if x['entry']==n)
            q=F(e['quadrature_error_b3']);t=F(e['tail_proposition6'])
            req(get_iv(e['combined_interval'])==(slo-q-t,shi+q+t),'composicao E7')
            req(ent['e7_composite_interval_exact']==[str(slo-q-t),str(shi+q+t)],'composicao V14')
            req((ent['quadrature_bound_declared_exact'],ent['tail_bound_declared_exact'])==(str(q),str(t)),'erro Q tail')
            totals[n]=ofs
        req(totals=={'G11':9524,'G22':6230} and doc['point_counts']=={'G11':9524,'G22':6230,'total':15754},'pontos totais')
        req(bandnum==11,'11 bandas')
        req(doc['all_point_exports_reproduced'] is True and doc['exact_integrand_enclosures_replayed']==15754,'replay flags')
    return {'points':15754,'bands':11,'composite':2}

def aggregate_metrics(root,doc):
    with zipfile.ZipFile(root/ZIP_NAME) as z:
        b5=json.loads(z.read('b5/RESULTADO_B5.json'))
    for it in b5['interval_results']:
        n=it['component'];m=doc['entries'][n].get('b5_b3_metrics')
        req(isinstance(m,dict),'metrics inexistente '+n)
        req(m.get('scope')=='aggregate_finite_sum_only' and m.get('comparison_count')==1,'escopo agregado ausente')
        for fld in ['nested_fraction','b5_in_b3_fraction','b3_in_b5_fraction']:
            req(type(m.get(fld)) in (float,int) and 0<=m[fld]<=1,'fracao invalida '+fld)
        exp=metrics(get_iv(it['certified_finite_midpoint_sum']),get_iv(it['b3_midpoint_interval']))
        req(m==exp,'desvios B3/B5 divergentes '+n)
        req(m['disjoint_count']==0,'somas finitas B3/B5 disjuntas; investigar '+n)
    return {'scope':'aggregate_finite_sum_only','individual_node_metrics':'not_available'}

def rounding(doc):
    for name,ent in doc['entries'].items():
        lo,hi=F(ent['finite_sum_lower']),F(ent['finite_sum_upper'])
        req(ent['finite_sum_interval_outward']==pair_display(lo,hi),'exibicao soma nao e outward '+name)
        q,t=F(ent['quadrature_bound_declared_exact']),F(ent['tail_bound_declared_exact'])
        c1,c2=get_iv(ent['e7_composite_interval_exact'])
        req((c1,c2)==(lo-q-t,hi+q+t),'limite composto')
        req(ent['e7_composite_interval_approx']==pair_display(c1,c2),'composto decimal nao e outward')
        req(F(ent['finite_sum_width_approx'])>=hi-lo,'largura decimal nao e outward')
        req(F(ent['quadrature_bound_declared_approx'])>=q and F(ent['tail_bound_declared_approx'])>=t,'erro decimal nao e outward')
    return {'components':2,'rounding':'outward'}

def sample(root,doc,details_override=None):
    raw=(root/'oraculo/AMOSTRA_128.json').read_bytes()
    got=doc['independent_oracle_sample']
    req(sha(raw)==got['details_sha256'],'hash testemunhas oraculo')
    details=json.loads(raw) if details_override is None else details_override
    rows=details['samples'];req(len(rows)==128 and got['sample_size']==128 and got['mismatches']==0,'amostra incompleta/mismatch')
    req(all(x['encloses_independent'] is True for x in rows),'oraculo mismatch flag')
    seen=set();band_counts={};const=ref.make_constants()
    for entry in rows:
        name=entry['component'];idx=entry['global_index'];key=(name,idx)
        req(key not in seen,'duplicata amostra');seen.add(key)
        band=(name,*entry['band']);band_counts[band]=band_counts.get(band,0)+1
        ref_int,om=ref.evaluate_independent(F(entry['r']),1 if name=='G11' else 2,const)
        req(entry['independent_rational_integrand']==[str(ref_int[0]),str(ref_int[1])],'oraculo numero adulterado')
        req(entry['independent_rational_omega']==[str(om[0]),str(om[1])],'oraculo omega adulterado')
        req(get_iv(entry['exported_integrand'])[0]<=ref_int[0] and ref_int[1]<=get_iv(entry['exported_integrand'])[1],'oraculo fora de B5')
        line=json.loads((root/'b5'/('PONTOS_'+name+'_B5.jsonl')).read_bytes().splitlines()[idx])
        req(F(line['r'])==F(entry['r']) and list(map(str,map(F,line['integrand'])))==entry['exported_integrand'],'testemunha nao corresponde linha')
    req(len(band_counts)==11 and sorted(band_counts.values())==[11]*4+[12]*7,'estratificacao alterada')
    # Exige os indices deterministas de cada faixa, e nao apenas a contagem.
    expected={}
    with zipfile.ZipFile(root/ZIP_NAME) as z:
        original=json.loads(z.read('b5/RESULTADO_B5.json'))
    ordinal=0
    for item in original['interval_results']:
        comp=item['component'];offset=0
        for band in item['bands']:
            k=(comp,band['from'],band['to'])
            expected[k]={offset+j for j in choose_indices(int(band['cells']),ordinal)}
            ordinal+=1;offset+=int(band['cells'])
    actual={k:set() for k in expected}
    for e in rows:actual[(e['component'],*e['band'])].add(e['global_index'])
    req(actual==expected,'indices deterministicamente estratificados divergem')
    return {'sample_size':128,'bands_covered':11,'mismatches':0}

def replay_v12(root,doc):
    raw=(root/'REPLAY_E8K_V12_INTEGRAL.json').read_bytes()
    x=json.loads(raw)
    req(sha(raw)==doc['full_replay_v12']['sha256'],'replay SHA')
    req(x['status']=='ALL_POINTS_REPRODUCED_B5_INPUT_HASHES_MATCHED' and x['provenance_gaps']==[],'replay V12 nao aprovado')
    req(x['point_counts']==doc['point_counts'],'replay contagem V12')
    for k in ('G11','G22'):
        req(x['entries'][k]['finite_sum_lower']==doc['entries'][k]['finite_sum_lower'] and
            x['entries'][k]['finite_sum_upper']==doc['entries'][k]['finite_sum_upper'],'replay V12 difere V14')
    return {'verified':True,'method':'B5 runtime exact replay'}

def analytic_v13(root,doc):
    x=json.loads((root/'RESULTADO_E8K_V13.json').read_text(encoding='utf8'))
    req(x.get('all_gates_passed') is True and all(x['gates'].values()),'cotas V13 falharam')
    return {'passed_under_analytic_premises':True,'not_formal_theorem_proof':True}

def gates(root,doc,do_selftest=False):
    checks={'integrity_sha256sums':lambda:manifest(root),
            'e7_b5_b3_source_integrity':lambda:sources(root,doc),
            'full_15754_point_band_reconciliation':lambda:bands(root,doc),
            'b3_b5_aggregate_metrics':lambda:aggregate_metrics(root,doc),
            'outward_decimal_rounding':lambda:rounding(doc),
            'independent_oracle_128':lambda:sample(root,doc),
            'v12_full_runtime_replay':lambda:replay_v12(root,doc),
            'v13_analytic_premises':lambda:analytic_v13(root,doc)}
    outcomes={}
    for k,fn in checks.items():
        try:outcomes[k]={'status':'PASS','witness':fn()}
        except Exception as e:outcomes[k]={'status':'FAIL','error':repr(e)[:200]}
    mutations=[]
    if do_selftest:
        suite=[
           ('b3_original_hash','e7_b5_b3_source_integrity',lambda d:d['entries']['G11'].__setitem__('original_b3_input_sha256','0'*64),'src'),
           ('band_sum','full_15754_point_band_reconciliation',lambda d:d['entries']['G22']['bands'][0]['sum'].__setitem__(0,'-999/1'),'bands'),
           ('metric_missing','b3_b5_aggregate_metrics',lambda d:d['entries']['G11'].pop('b5_b3_metrics'),'metrics'),
           ('fraction_outside_domain','b3_b5_aggregate_metrics',lambda d:d['entries']['G22']['b5_b3_metrics'].__setitem__('nested_fraction',2.1),'metrics'),
           ('oracle_mismatch','independent_oracle_128',lambda d:d['independent_oracle_sample'].__setitem__('mismatches',1),'oracle')]
        fnmap={'src':lambda d:sources(root,d),'bands':lambda d:bands(root,d),
            'metrics':lambda d:aggregate_metrics(root,d),'oracle':lambda d:sample(root,d)}
        for name,target,mut,key in suite:
            if outcomes.get(target,{}).get('status')!='PASS':
                mutations.append({'name':name,'detected':False,'reason':'baseline gate did not PASS'});continue
            changed=copy.deepcopy(doc);mut(changed)
            try:fnmap[key](changed);det=False
            except Exception:det=True
            mutations.append({'name':name,'detected':det,'target_gate':target})
        outcomes['negative_mutations']={'status':'PASS' if all(x['detected'] for x in mutations) else 'FAIL',
             'witness':mutations}
    return outcomes

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path('.'))
    ap.add_argument('--out',type=Path,default=Path('CONFERENCIA_E8K_V14.json'))
    ap.add_argument('--selftest',action='store_true');a=ap.parse_args()
    root=a.root.resolve();doc=json.loads((root/'RESULTADO_E8K_V14.json').read_text(encoding='utf8'))
    out=gates(root,doc,a.selftest)
    passed=all(x['status']=='PASS' for x in out.values())
    report={'schema':'NUMERIA.E8K.V14.EVIDENCE_CONFERENCE.1', 'status':'PASSED' if passed else 'FAIL',
        'all_gates_passed':passed,'gates':out,'passed':sum(x['status']=='PASS' for x in out.values()),
        'failed':sum(x['status']=='FAIL' for x in out.values()),
        'not_certified':doc['not_certified']}
    outfile=a.out if a.out.is_absolute() else root/a.out
    outfile.parent.mkdir(parents=True,exist_ok=True)
    outfile.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf8')
    print(f'V14 gates {report["passed"]}/{len(out)} PASS; {report["failed"]} FAIL')
    for n,x in out.items():print(' ',n,x['status'],x.get('error',''))
    raise SystemExit(0 if passed else 2)
if __name__=='__main__':main()
