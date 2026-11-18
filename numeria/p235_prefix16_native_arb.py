#!/usr/bin/env python3
"""P235 independent interval GL60 16x16 four-window probe; no RH claim.
The E_GL formula is conditional on classical Binet/Bernstein/Bessel theorems.
"""
import argparse, hashlib, json, os, time
from fractions import Fraction
from pathlib import Path
from flint import arb, acb, arb_mat, ctx

def ball(x):
    m,e=x.mid().man_exp()
    r,k=x.rad().upper().man_exp()
    return [[str(m),str(e)],[str(r),str(k)]]

def run(sector,Lnum,Lden,T=504,N=16,bits=900):
    if sector not in ('even','odd') or (Lnum,Lden) not in ((4,5),(161,200)) or T!=504 or N!=16 or bits<900:
        raise ValueError('only P235 specified panels and precision')
    ctx.prec=bits; t0=time.monotonic()
    pi=arb.pi();L=arb(Lnum)/Lden; two=arb(2)
    l2=two.log();l3=arb(3).log();l4=two*l2
    prime=[(2*l2/two.sqrt(),l2),
           (2*l3/arb(3).sqrt(),l3),
           (l2,l4)]
    if (Lnum,Lden)==(161,200):
        prime.append((two*arb(5).log()/arb(5).sqrt(),arb(5).log()))
    AL=sum((a for a,v in prime),arb(0))
    beta=(arb(T)/(two*pi)).log()-arb(1)/T-AL
    if not (arb(0)<beta.lower() and beta.upper()<arb(2)):
        raise ArithmeticError('positive beta <2 prerequisite fails')
    roots=[arb.legendre_p_root(60,k,weight=True) for k in range(60)]
    if not all(roots[i][0].lower()>roots[i+1][0].upper() for i in range(59)):
        raise ArithmeticError('GL60 root order failed')
    if not sum((w for x,w in roots),arb(0)).contains(2):
        raise ArithmeticError('GL60 weight sum does not enclose 2')
    modes=[2*i+int(sector=='odd') for i in range(N)]
    norms=[(-1)**i*(two*L*(2*n+1)).sqrt() for i,n in enumerate(modes)]
    z0=L/2
    pvec=[]
    for n in modes:
        pvec.append(two*((2*n+1)*L/2).sqrt()*(pi/(two*z0)).sqrt()
                    * z0.bessel_i(arb(n)+arb(1)/2))
    acc=arb_mat(N,N)
    quarter=arb(1)/4
    for panel in range(2*T):
        ctr=arb(2*panel+1)/4
        V=[];WV=[]
        for root,wgl in roots:
            t=ctr+quarter*root
            x=L*t
            if x.lower()<=0:
                raise ArithmeticError('nonpositive Bessel argument')
            pre=(pi/(2*x)).sqrt()
            u=[norms[i]*pre*x.bessel_j(arb(n)+arb(1)/2)
               for i,n in enumerate(modes)]
            dig=acb(arb(1)/4,t/2).digamma().real
            psi=dig-pi.log()-sum((a*(t*ln).cos() for a,ln in prime),arb(0))
            w=quarter*wgl*(psi-beta)/pi
            V.extend(u);WV.extend([w*z for z in u])
        J=arb_mat(60,N,V)
        WJ=arb_mat(60,N,WV)
        acc += J.transpose()*WJ
        if panel%168==0:
            print(sector,Lnum,Lden,'panel',panel,'of',2*T,
                  'seconds',round(time.monotonic()-t0,1),flush=True)
    E=Fraction(2*T*12*2100,3**120)
    Eupper=(arb(E.numerator)/arb(E.denominator)).upper()
    matrix=[[None]*N for _ in range(N)]
    sample_C=[]
    for i in range(N):
        for j in range(i,N):
            C=(acc[i,j]+acc[j,i])/2
            # Analytic integral = quadrature + error of absolute size <= E.
            # Add each error just once, AFTER all 1008 panels.
            R=C + (beta if i==j else arb(0))
            pole=2*pvec[i]*pvec[j]
            R=R+pole if sector=='even' else R-pole
            R=R+arb(0,Eupper)
            if not R.is_finite():raise ArithmeticError('nonfinite matrix entry')
            matrix[i][j]=R;matrix[j][i]=R
            if len(sample_C)<6:
                sample_C.append([i,j,ball(C)])
    # Interval LDL: every realized symmetric matrix in the printed entry boxes
    # has positive LDL pivots when these pivot lower bounds are positive.
    LL=[[arb(0) for _ in range(N)] for _ in range(N)]
    diag=[]
    pivots=[]
    for k in range(N):
        p=matrix[k][k]
        for a in range(k):
            p-=LL[k][a]*LL[k][a]*diag[a]
        if p.lower()<=0:
            raise ArithmeticError(f'interval LDL inconclusive/nonpositive pivot at {k}: {p}')
        pivots.append(ball(p))
        diag.append(p)
        for i in range(k+1,N):
            t=matrix[i][k]
            for a in range(k):
                t-=LL[i][a]*LL[k][a]*diag[a]
            LL[i][k]=t/p
    upper=[[i,j,ball(matrix[i][j])] for i in range(N) for j in range(i,N)]
    return {'schema':'NUMERIA.ES8P235.ARB_PREFIX16_V1',
      'version':__import__('flint').__version__,'bits':bits,
      'sector':sector,'L':f'{Lnum}/{Lden}','Tsharp':T,'N':N,
      'GL60_source':'independent arb.legendre_p_root(60,k,weight=True)',
      'special_functions':'native spherical Bessel J_{n+1/2}, native complex digamma, native modified Bessel I_{n+1/2}',
      'prime_terms':[2,3,4]+([5] if Lnum==161 else []),
      'beta_ball':ball(beta),'prime_constant_ball':ball(AL),
      'quadrature_error_bound_exact':str(E),'quadrature_error_added_once':True,
      'quadrature_error_analytic_theorems_imported':True,
      'true_WL_reduction_formalized':False,
      'interval_ldl_16x16_all_pivots_positive':True,
      'interval_ldl_pivot_balls':pivots,
      'quadrature_C_sample_balls':sample_C,
      'interval_R_upper_triangle':upper,
      'nodes':2*T*60,
      'elapsed_seconds_diagnostic':round(time.monotonic()-t0,2),
      'github_sha':os.environ.get('GITHUB_SHA','unknown'),
      'Riemann_hypothesis_proved':False}

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--sector',choices=['even','odd'],required=True)
    p.add_argument('--Lnum',type=int,required=True)
    p.add_argument('--Lden',type=int,required=True)
    p.add_argument('--bits',type=int,default=900)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    r=run(a.sector,a.Lnum,a.Lden,bits=a.bits)
    a.out.write_text(json.dumps(r,indent=2)+'\n')
    print('PASS 16x16 interval LDL for',r['sector'],r['L'],
          'SHA256',hashlib.sha256(a.out.read_bytes()).hexdigest(),
          'secs',r['elapsed_seconds_diagnostic'],flush=True)
