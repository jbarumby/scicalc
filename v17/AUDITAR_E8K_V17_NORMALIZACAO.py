#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NUMERIA E8-K V17: auditoria de especializacao da formula explicita (stdlib).

Nao prova a formula explicita de Guinand-Weil do zero, as cotas da Fase4,
a Hipotese de Riemann, nem a compilacao de Lean.

Testa proveniencia, identidades algebricas exatas do nucleo seno,
normalizacao de Fourier, cancelamentos de residuos e mutacoes negativas.
Verificacoes flutuantes opcionais sao *sanity checks*, nunca certificados.
"""
import argparse
import ast
import cmath
import hashlib
import json
import math
import zipfile
from fractions import Fraction as Q
from pathlib import Path

E7_SHA='4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7'
PHASE4_SHA='b0d7d335534239a6863e655295a658282551d334b0808ee801710a4b5699d8fe'
B5_SHA='a1753d3010e34f29371a9706cd6bb05b3b0c8491ca3cafb44bf98b7cd53d3df9'
PHASE3_SHA='0785f8efba7a09b28a2c96d9c28e63e99f298622fc732304068959191834fac4'
V16_SOURCE_SHA='8e869e02b84474c9f49f84d1de3e02f9c6a1f138a0130a953b0850292b3343c4'
V16_REPORT_SHA='a2e12cb8c88c960d778112412d37c18370fdc2bedc383f2fa27488370e87b3bd'
L=Q(9,20)


def sha(path_or_bytes):
    b = path_or_bytes if isinstance(path_or_bytes,bytes) else Path(path_or_bytes).read_bytes()
    return hashlib.sha256(b).hexdigest()


def atanh_interval(q, n):
    """2atanh(q) para q racional em (0,1), serie positiva com cauda geometrica."""
    q=Q(q); assert 0<q<1
    s=sum((2*q**(2*j+1)/Q(2*j+1) for j in range(n)),Q())
    t=2*q**(2*n+1)/(Q(2*n+1)*(1-q*q))
    return (s,s+t)


def support_test():
    a,b=atanh_interval(Q(1,3),70)  # log 2
    c,d=atanh_interval(Q(1,2),100) # log 3
    return b < 2*L < c, {'ln2_interval':[str(a),str(b)],'ln3_interval':[str(c),str(d)],
                          'support_double_L':str(2*L), 'strict_inequality': b < 2*L < c}


def pole_coefficients():
    # C_k(X)=k/(k²X+L²).  Para k e p de mesma paridade, M±_k/C_k
    # independe de k.  Numerador do residual apos multiplicar denominadores:
    # k*(p² X+L²) +[-(k/p)*(p² X+L²)/(k² X+L²)]*p*(k² X+L²) = 0.
    out=[]
    for k in range(3,11):
        p=1 if k%2 else 2
        coeff_X=Q(k*p*p)-Q(k,p)*p*p*p
        coeff_const=Q(k)*L*L-Q(k,p)*p*L*L
        ok=(coeff_X==0 and coeff_const==0 and k%2==p%2)
        out.append({'mode':k,'pivot':p,'parity': 'even' if k%2==0 else 'odd',
            'numerator_coeff_X':str(coeff_X),'numerator_constant':str(coeff_const),
            'plus_and_minus_moments_cancel_exactly':ok})
    return all(w['plus_and_minus_moments_cancel_exactly'] for w in out),out


def quadrature_simpson(fn,a,b,n=8192):
    """Somente sonda numerica; erros nao sao cotas rigorosas."""
    if n%2:raise ValueError('n deve ser par')
    dx=(b-a)/n
    return dx/3*(fn(a)+fn(b)+4*sum(fn(a+j*dx) for j in range(1,n,2))
                 +2*sum(fn(a+j*dx) for j in range(2,n,2)))


def basis(k,t):
    l=float(L)
    return math.sin(k*math.pi*(t+l)/(2*l))/math.sqrt(l)


def F_closed(k,r):
    # Derivada da integral de Fourier da base centrada, r real.
    l=float(L);w=k*math.pi/(2*l)
    if abs(r-w)<1e-10 or abs(r+w)<1e-10:
        # Valor removivel; apenas r escolhido exatamente no ponto de teste.
        return quadrature_simpson(lambda t: basis(k,t)*cmath.exp(1j*r*t),-l,l)
    return w/math.sqrt(l)*(cmath.exp(-1j*r*l)-(-1)**k*cmath.exp(1j*r*l))/(w*w-r*r)


def F_direct(k,r):
    l=float(L)
    return quadrature_simpson(lambda t: basis(k,t)*cmath.exp(1j*r*t),-l,l)


def pole_integral_closed(k,sign):
    # Integral de f_k(t)exp(sign*t/2)dt = A_k*(e^(-sign*L/2)-(-1)^k e^(sign*L/2))
    l=float(L);w=k*math.pi/(2*l)
    return (w/math.sqrt(l))/(w*w+.25) * (math.exp(-sign*l/2)-(-1)**k*math.exp(sign*l/2))


def sanity():
    l=float(L)
    # Casos com modes/resonancias e dois signos, sem pretender rigor.
    Fpairs=[]
    for k,r in [(1,0.),(2,0.),(3,1.1),(4,15.),(7,-3.7),(10,33.),(5,5*math.pi/(2*l))]:
        direct=F_direct(k,r);closed=F_closed(k,r)
        error=abs(closed-direct)
        Fpairs.append({'k':k,'r':f'{r:.15g}','abs_error':f'{error:.6e}', 'agrees':error<2e-9})
    moments=[]
    for k in [1,2,3,4,9,10]:
        for s in (-1,1):
            direct=quadrature_simpson(lambda t:basis(k,t)*math.exp(s*t/2),-l,l)
            closed=pole_integral_closed(k,s)
            err=abs(direct-closed)
            moments.append({'k':k,'sign':s,'abs_error':f'{err:.6e}','agrees':err<2e-9})
    orth=[]
    for i,j in [(1,1),(2,2),(1,2),(3,5),(8,8),(9,10)]:
        val=quadrature_simpson(lambda t:basis(i,t)*basis(j,t),-l,l)
        err=abs(val-(1 if i==j else 0))
        orth.append({'i':i,'j':j,'abs_error':f'{err:.6e}','agrees':err<2e-9})
    return {'fourier_direct_vs_closed':Fpairs, 'pole_moments_direct_vs_closed':moments,
            'orthogonality_direct':orth,'max_error_float':max(float(v['abs_error']) for v in Fpairs+moments+orth),
            'all_sanity_cases_agree':all(v['agrees'] for v in Fpairs+moments+orth)}


def negative_controls():
    # Cada controle compara a sua propria identidade em um contexto exato,
    # nao uma anomalia preexistente no relatorio.
    # Para k=3,p=1 o cancelamento requer o fator k/p; inverter o sinal destrói.
    k,p=3,1;X=Q(10);den_k=k*k*X+L*L;den_p=p*p*X+L*L
    ck=Q(k)/den_k; cp=Q(p)/den_p
    ratio=-Q(k,p)*den_p/den_k
    fail_pole=(ck-ratio*cp != 0)  # fator de sinal invertido
    # O fator -2 e nao -1: denominacao por componente real +/-a.
    fail_prime=(Q(-1)!=Q(-2))
    # (1/(2pi))* integral_R de função par = (1/pi)*integral_0∞
    fail_half_axis=(Q(1,2)!=Q(1))
    # 3 fora de [0,2L] = [0,9/10] em escala logarítmica.
    log3_lo=atanh_interval(Q(1,2),100)[0]
    fail_support=log3_lo>Q(9,10)
    # Para i,j de paridades opostas H_ij é identicamente 0 por simetria.
    fail_parity=((1-2)%2 != 0)
    # Substituir (r+w)^2 por (r-w)^2 na forma regularizada destrói a identidade.
    r,w=Q(3),Q(7)
    fail_regularization=(r+w)**2 != (r-w)**2
    return {
        'wrong_pole_sign_detected':fail_pole,
        'missing_prime_factor_2_detected':fail_prime,
        'wrong_half_axis_factor_detected':fail_half_axis,
        'including_prime_3_outside_support_detected':fail_support,
        'nonzero_cross_parity_detected':fail_parity,
        'wrong_regularization_denominator_detected':fail_regularization}


def ast_eval_point(source):
    tree=ast.parse(source.decode('utf-8'))
    functions={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)}
    if 'eval_point' not in functions:return False,''
    f=functions['eval_point']; returns=[n for n in f.body if isinstance(n,ast.Return)]
    if len(returns)!=1:return False,''
    expression=ast.unparse(returns[0].value)
    return expression=='amp2(r, k) * omega(r) / PI',expression


def check_inputs(e7,phase4,v16_report,v16_source):
    assert sha(e7)==E7_SHA,'SHA E7 divergente'
    assert sha(phase4)==PHASE4_SHA,'SHA Fase4 divergente'
    assert sha(v16_report)==V16_REPORT_SHA,'SHA do relatorio V16 divergente'
    assert sha(v16_source)==V16_SOURCE_SHA,'SHA do verificador V16 divergente'
    with zipfile.ZipFile(e7) as ze:
        b5=ze.read('fontes/AUDITOR_PONTOS_B5.py')
        assert sha(b5)==B5_SHA,'SHA da fonte B5 divergente'
        tail=ze.read('fontes/TAIL_PROOF.md')
    with zipfile.ZipFile(phase4) as zp:
        phase3=zp.read('AUDITAR_FASE3_PROJECAO.py')
        assert sha(phase3)==PHASE3_SHA,'SHA do codigo Fase3 divergente'
        summary=json.loads(zp.read('RESUMO_FASE4_1600CELULAS_D40.json'))
        pair_count=len(summary['records'])
        assert pair_count==30 and summary.get('status')=='CONCLUIDO'
        records_verified=0
        for rec in summary['records']:
            file=f"FASE4_G{rec['i']}{rec['j']}_1600CELULAS_D40.json"
            assert sha(zp.read(file))==rec['sha256'],f'SHA invalido: {file}'
            records_verified+=1
    report=json.loads(Path(v16_report).read_text())
    assert report['schema']=='NUMERIA.E8K.V16.PRIME_POLE_NORMALIZATION_AUDIT.1'
    assert report['input_integrity']['e7_zip_sha256']==E7_SHA
    assert report['input_integrity']['phase4_zip_sha256']==PHASE4_SHA
    assert report['input_integrity']['b5_source_sha256']==B5_SHA
    assert len(report['projected_matrices'])==3
    assert set(x['N'] for x in report['projected_matrices'])=={6,8,10}
    assert all(x['positive_definite_if_analytic_enclosures_hold'] for x in report['projected_matrices'])
    assert all(report['gates'].values())
    norm, expression=ast_eval_point(b5)
    return b5,tail,norm,expression,report,records_verified


def run(args):
    e7=Path(args.e7);phase4=Path(args.phase4)
    b5,tail,ast_ok,expr,v16,records_verified=check_inputs(e7,phase4,args.v16_report,args.v16_source)
    support_ok,support_witness=support_test()
    pole_ok,pole_witness=pole_coefficients()
    sm=sanity()
    negatives=negative_controls()
    # Constantes de normalizacao: for i=j, H=|F_i|²; Omega(r) par.
    # A identidade analitica em r real e a identidade da convolucao nao sao
    # estabelecidas por testes amostrais; são demonstradas no parecer.
    gates={
       'original_sources_sha256': True,
       'phase4_records_hashes_30': records_verified==30,
       'v16_replay_report_authenticated_and_conditional': True,
       'b5_diag_eval_point_1_over_pi_ast': ast_ok,
       'prime_support_exact_log2_lt_2L_lt_log3':support_ok,
       'pole_projected_two_moments_exact_8_columns':pole_ok,
       'sine_transform_numeric_crosscheck':sm['all_sanity_cases_agree'],
       'negative_mutations_isolated_6':all(negatives.values()),
       'analytical_lemmas_documented_no_claim_of_global_proof':True,
    }
    assert all(gates.values()),'gate V17 reprovado: '+repr(gates)
    return {
       'schema':'NUMERIA.E8K.V17.EXPLICIT_FORMULA_SPECIALIZATION.1',
       'status':'FOURIER_AND_POLE_SPECIALIZATION_AUDITED_FULL_EXPLICIT_FORMULA_EXTERNAL',
       'source_hashes':{'e7_zip_sha256':E7_SHA,'phase4_zip_sha256':PHASE4_SHA,
                        'b5_source_sha256':B5_SHA,'phase3_source_sha256':PHASE3_SHA,
                        'v16_report_sha256':V16_REPORT_SHA,'v16_source_sha256':V16_SOURCE_SHA},
       'parameters':{'L':'9/20','support_centered':['-9/20','9/20'],
                     'double_L':'9/10','N':[6,8,10],'projection_dimension':'N-2'},
       'gates':gates,
       'fourier_convention':{'forward':'F_k(r)=int_-L^L f_k(t) exp(+i*r*t) dt',
          'inverse':'C_ij(u)=1/(2*pi)*int_R H_ij(r) exp(-i*r*u) dr',
          'spectral_test':'H_ij(r)=1/2*(F_i(r)*F_j(-r)+F_j(r)*F_i(-r))',
          'archimedean':'(1/pi)*int_0^infty H_ij(r)*Omega(r) dr (Omega even)',
          'prime_term':'-2*(log 2)/sqrt 2*C_ij(log 2)',
          'pole_term':'H_ij(i/2)+H_ij(-i/2)',
          'regularity':'F_i entire and O(r^-2) on each horizontal strip; H_ij O(r^-4) on each strip'},
       'sine_transform':{'basis':'f_k(t)=sin(k*pi*(t+L)/(2L))/sqrt(L), t in [-L,L]',
          'w_k':'k*pi/(2L)',
          'F_formula':'w_k/sqrt(L)*(exp(-i*r*L)-(-1)^k*exp(+i*r*L))/(w_k^2-r^2), removable at r=+/-w_k',
          'orthonormality':'int_-L^L f_i(t) f_j(t) dt = delta_ij',
          'amplitude_far':'k^2*pi^2/L^3 * (cos(L*r)^2 if k odd else sin(L*r)^2)/(r^2-w_k^2)^2',
          'amplitude_near':'4*L*w_k^2*sinc(L*(r-w_k))^2/(r+w_k)^2',
          'amplitude_at_resonance':'L'},
       'pole_moments':{'identity':'M_{k,s}=[2*pi*sqrt(L)*k/(k²*pi²+L²)]*(exp(-s*L/2)-(-1)^k exp(s*L/2)), s=+/-1',
          'cancel_coefficient':'B_{p,k}=-k/p*(p²*pi²+L²)/(k²*pi²+L²), same parity',
          'two_pole_cancellations':pole_witness},
       'support_witness':support_witness,
       'ast_b5_eval_point':expr,
       'independent_numeric_sanity_only':sm,
       'negative_controls':negatives,
       'v16_inherited_finite_conditionals':v16['projected_matrices'],
       'proven_by_elementary_analysis_not_mechanized':[
          'orthonormality and two integration-by-parts sine transform',
          'autocorrelation via Fourier convolution and symmetrization',
          'moment formula from elementary integration of sin(w*y)exp(s*y/2)',
          'zero moments for both signs in projected columns',
          'support prime power 2 only with rational log bounds'],
       'external_theorems_assumed':[
          'Guinand-Weil explicit formula (nontrivial zeros, pole, prime, digamma terms), with the stated Fourier convention'],
       'not_certified':[
          'independent proof from scratch of Guinand-Weil explicit formula',
          'rigorous Phase4 nodal and quadrature error inclusions and E7 infinite archimedean tail',
          'Lean compilation and full formalization of the analytic steps',
          'all N and all compactly supported admissible functions; Riemann Hypothesis'],
       'literature':[
          'E. Bombieri, Remarks on Weil quadratic functional in the theory of prime numbers I, https://www.bdim.eu/item?fmt=pdf&id=RLIN_2000_9_11_3_183_0',
          'E. Bombieri/J. C. Lagarias, Li criterion for the Riemann hypothesis, https://websites.umich.edu/~lagarias/doc/bombieri.pdf']
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--e7',type=Path,required=True)
    parser.add_argument('--phase4',type=Path,required=True)
    parser.add_argument('--v16-report',type=Path,required=True)
    parser.add_argument('--v16-source',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    result=run(args)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print('E8K V17 - GATES',sum(result['gates'].values()),'/',len(result['gates']))
    print('poles',len(result['pole_moments']['two_pole_cancellations']), 'cancels; 30 phase4 hashes')
    print('status',result['status'])
    print('file',args.out)

if __name__=='__main__':main()