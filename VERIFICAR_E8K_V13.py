#!/usr/bin/env python3
"""V13 — auditoria independente das cotas B3/Prop6 do E7 por frações exatas.

Limites: a identidade de digama, a cota de derivada de Fourier e o teorema do
ponto médio são premissas analíticas explícitas. Os valores B3 de Omega(0)
são uma premissa de intervalo, não são recomputados neste arquivo. A execução
certifica algebricamente que as cotas NUMÉRICAS publicadas dominam as cotas
recalculadas a partir dessas premissas. Não demonstra Weil ou RH.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from math import factorial, isqrt
from pathlib import Path
from zipfile import ZipFile

ZIP_SHA = '4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7'
ROOT = 'NUMERIA.E8K.V13.ANALYTIC_BOUNDS_AUDIT.1'
R = Q(160)
L = Q(9,20)
M = 128

def sha(b): return hashlib.sha256(b).hexdigest()
def num(obj): return Q(obj['numerator']) / Q(obj['denominator'])
def atanh_interval(q, n=150):
    assert Q(0) <= q < 1
    q=Q(q); q2=q*q; p=q; acc=Q(0)
    for j in range(n): acc += 2*p/Q(2*j+1); p*=q2
    tail=2*p/Q(2*n+1)/(1-q2)
    return acc,acc+tail

def log_pos_bounds(x):
    x=Q(x);assert x>0
    k=0
    while x>=2: x/=2;k+=1
    while x<1: x*=2;k-=1
    t=(x-1)/(x+1)
    a,b=atanh_interval(t)
    c,d=atanh_interval(Q(1,3))
    if k>=0: return a+k*c,b+k*d
    return a+k*d,b+k*c

def atan_enclosure(x,n):
    x=Q(x);tot=Q(0)
    for j in range(n): tot+=(-1)**j*x**(2*j+1)/Q(2*j+1)
    rest=x**(2*n+1)/Q(2*n+1)
    return (tot,tot+rest) if n%2==0 else (tot-rest,tot)

def pi_bounds():
    a=atan_enclosure(Q(1,5),85); b=atan_enclosure(Q(1,239),50)
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]

def sqrt_bounds(x, prec=65):
    t=10**prec;x=Q(x)
    v=isqrt(x.numerator*t*t//x.denominator)
    return Q(v,t),Q(v+1,t)

def interval_omega(r,omega0):
    """Independente das rotinas Python B3 e B5; usa como premissa Omega0 B3."""
    r=Q(r);b2=r*r/4
    partial=sum((b2/(u*(u*u+b2)) for u in (Q(j)+Q(1,4) for j in range(M))),Q(0))
    u=Q(M)+Q(1,4)
    q=b2/(u*u)
    li,ls=log_pos_bounds(1+q)
    first=b2/(u*(u*u+b2))
    return omega0[0]+partial+li/2,omega0[1]+partial+ls/2+first

def upper_derivs(l,h):
    d1=d2=Q(0)
    for j in range(M):
        u=Q(j)+Q(1,4);u2=u*u
        d1+=min(u*h/(2*(u2+l*l/4)**2),Q(1,3)/u2)
        d2+=min(u*max(abs(u2-3*l*l/4),abs(u2-3*h*h/4))/(2*(u2+l*l/4)**3),Q(1,2*u**3))
    u=Q(M)+Q(1,4)
    return d1+(Q(1,u*u)+Q(1,u))/3,d2+Q(1,2*u**3)+Q(1,4*u*u)

def m2_bound(k,l,h,pi_lo,pi_hi,sqrt_lo,sqrt_hi,omega0):
    low_l,high_l=interval_omega(l,omega0)
    low_h,high_h=interval_omega(h,omega0)
    # Omega crescente em r >= 0; |Omega| majorizado pelos valores extremos.
    omega_abs=max(abs(low_l),abs(high_l),abs(low_h),abs(high_h))
    d1,d2=upper_derivs(l,h)
    switch=6*k
    if h<=switch:
        a0=2*sqrt_hi
        a1=sqrt_hi*L
        a2=Q(2,3)*sqrt_hi*L**2
    else:
        if l<switch: raise AssertionError('faixa cruza a transição do polo')
        w2=k*k*pi_hi*pi_hi/(4*L*L)
        D=l*l-w2
        if D<=0: raise AssertionError('denominador fora do dominio')
        C=k*pi_hi/(L*sqrt_lo)
        a0=C/D
        a1=C*(L/D+2*h/D**2)
        a2=C*(L*L/D+4*L*h/D**2+2/D**2+8*h*h/D**3)
    f2=((2*a1*a1+2*a0*a2)*omega_abs+4*a0*a1*d1+a0*a0*d2)/pi_lo
    return f2

def tail_prop6(k,pi_hi,lnr_hi):
    w=k*pi_hi/(2*L)
    assert R>=1 and 2*w<=R
    den=1-(w/R)**2
    return k*k*pi_hi/(L**3*den**2)*(lnr_hi+Q(314,100)+Q(1,3))/(3*R**3)

def one(z,k,tag,pi_lo,pi_hi,sqrt_lo,sqrt_hi,omega0,lnr_hi):
    result=[]
    b3_total=Q(0)
    for band in z['mesh']:
        a,b=Q(band['from']),Q(band['to']);n=int(band['cells'])
        f2_obj=band.get('f2_absolute_bound') or band.get('f2_bound')
        err_obj=band.get('midpoint_error_bound') or band.get('error_bound')
        f2_decl=num(f2_obj);err_decl=num(err_obj)
        computed=m2_bound(k,a,b,pi_lo,pi_hi,sqrt_lo,sqrt_hi,omega0)
        # Ao comparar f2, usamos a MESMA prova analitica, mas constantes
        # e avaliador racionais independentes (sem executar a fonte B3).
        f2_ok=computed<=f2_decl
        midpoint_req=computed*(b-a)**3/(24*n*n)
        err_ok=midpoint_req<=err_decl
        rounded_f2_midpoint_check=(f2_decl*(b-a)**3/(24*n*n)<=err_decl)
        result.append({'faixa':[str(a),str(b)],'cells':n,'f2_dominates_independent_rederivation':f2_ok,
                       'midpoint_formula_dominates':err_ok,'rounded_f2_midpoint_check':rounded_f2_midpoint_check,'published_f2_approx':f2_obj['decimal_approx'],
                       'independent_f2_approx':f'{float(computed):.16g}',
                       'margin_f2_approx':f'{float(f2_decl-computed):.5e}',
                       'published_midpoint_error_approx':err_obj['decimal_approx']})
        b3_total+=err_decl
    if tag=='G11': total=num(z['quadrature_absolute_error_upper'])
    else: total=num(z['G22_quadrature_error_upper'])
    t=tail_prop6(k,pi_hi,lnr_hi)
    return {'bands':result,'bands_total':len(result),'all_band_m2_passed':all(x['f2_dominates_independent_rederivation'] for x in result),
            'all_midpoint_roundings_passed':all(x['midpoint_formula_dominates'] for x in result),
            'sum_published_error_le_total':b3_total<=total,
            'published_quadrature_total_approx':f'{float(total):.16g}',
            'independent_proposition6_tail_upper_approx':f'{float(t):.16g}',
            'independent_proposition6_tail_upper_exact':str(t),
            'tail_rational_upper':t}

def audit(zip_path):
    b=Path(zip_path).read_bytes()
    if sha(b)!=ZIP_SHA: raise AssertionError('Fonte E7 nao confere')
    with ZipFile(zip_path) as z:
        b11=json.loads(z.read('corrigidos/RESULTADO_B3_G11_E2_CORRIGIDO.json'))
        b22=json.loads(z.read('corrigidos/RESULTADO_B3_G22_E2_CORRIGIDO.json'))
        e7=json.loads(z.read('RESULTADO_INTEGRAIS_E7.json'))
        o=b11['omega0_enclosure']
        omega0=(num(o[0]),num(o[1]))
        pi_lo,pi_hi=pi_bounds(); sqrt_lo,sqrt_hi=sqrt_bounds(L)
        _,lnr_hi=log_pos_bounds(R)
        entries={}
        for name,k,doc in [('G11',1,b11),('G22',2,b22)]:
            out=one(doc,k,name,pi_lo,pi_hi,sqrt_lo,sqrt_hi,omega0,lnr_hi)
            orig=next(x for x in e7['results'] if x['entry']==name)
            out['e7_quadrature_match_b3']=Q(orig['quadrature_error_b3'])==(
              num(doc['quadrature_absolute_error_upper']) if k==1 else num(doc['G22_quadrature_error_upper']))
            out['e7_tail_ge_independent_prop6']=Q(orig['tail_proposition6'])>=out['tail_rational_upper']
            out['e7_tail_declared_approx']=f"{float(Q(orig['tail_proposition6'])):.16g}"
            del out['tail_rational_upper']
            entries[name]=out
    flags={k:v for name,entry in entries.items() for k,v in {
        name+'.m2':entry['all_band_m2_passed'],name+'.midpoint_rounding':entry['all_midpoint_roundings_passed'],
        name+'.error_sum':entry['sum_published_error_le_total'],name+'.error_e7_match':entry['e7_quadrature_match_b3'],
        name+'.tail_prop6':entry['e7_tail_ge_independent_prop6']}.items()}
    return {'schema':ROOT, 'status':'NUMERIC_BOUNDS_REDERIVED_FROM_STATED_ANALYTIC_PREMISES' if all(flags.values()) else 'INCONSISTENCY_FOUND',
      'all_gates_passed':all(flags.values()),'gates':flags,'entries':entries,
      'inputs_sha256':{'e7_zip':ZIP_SHA},
      'analytic_premises_not_formalized':['Identidade da digama e inclusao de Omega(0) publicada no B3',
      'Majorantes de Omega,Omega_prime,Omega_double_prime sobre cada faixa',
      'Estimativas A_k,A_k_prime,A_k_double_prime via Fourier e formula trigonometrica',
      'Teorema de erro da regra composta do ponto medio',
      'Lema global |Omega(r)| <= log(r) + 3.14 para r>=1',
      'Proposicao 6 da cauda e desigualdades de monotonicidade envolvidas'],
      'not_certified':['Valores nodais B5 nao reavaliados aqui','Formula explicita de Weil global','Hipotese de Riemann']}

def main():
    p=argparse.ArgumentParser();p.add_argument('--zip',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    try: doc=audit(a.zip)
    except Exception as exc:
        print('V13 FALHA:',repr(exc));return 3
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(doc,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')
    print('V13:',doc['status'],'\n',doc['gates'])
    return 0 if doc['all_gates_passed'] else 2
if __name__=='__main__':raise SystemExit(main())
