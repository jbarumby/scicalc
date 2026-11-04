"""NUMERIA ES8-P — independent projection/prime sanity and conditional V24 error bounds.

Exact rationals enclose constants pi, log2, log3, sqrt(L). Symbolic pole
projection is exact for any X=pi**2; primal kernel uses analytic trig formulas
and is cross-checked against independent numerical quadrature.

NOT a certified evaluation of the full Weil matrix and NOT a proof of RH.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as Q
from math import comb, isqrt
from pathlib import Path
from functools import lru_cache

import mpmath as mp
import sympy as sp

L = Q(9, 20)
SCALE = 10**30

def floor_q(v: Q) -> Q:
    return Q((v.numerator * SCALE) // v.denominator, SCALE)

def ceil_q(v: Q) -> Q:
    return Q(-((-v.numerator * SCALE)//v.denominator), SCALE)

def atan_bounds(x: Q, terms: int = 45) -> tuple[Q,Q]:
    """Alternating arctan series at 0<x<1, next term bounds the error."""
    if not Q(0)<x<Q(1) or terms<2: raise ValueError('atan inputs')
    s=sum(((-1)**k*x**(2*k+1)/Q(2*k+1) for k in range(terms)),Q())
    err=x**(2*terms+1)/Q(2*terms+1)
    return (s-err,s) if terms%2==0 else (s,s+err)

def atanh_bounds(x:Q,terms:int=90)->tuple[Q,Q]:
    """Positive log series with explicit geometric remainder."""
    if not Q(0)<=x<Q(1) or terms<1: raise ValueError('atanh inputs')
    s=sum((2*x**(2*k+1)/Q(2*k+1) for k in range(terms)),Q())
    remainder=2*x**(2*terms+1)/(Q(2*terms+1)*(1-x*x))
    return s,s+remainder

@lru_cache(maxsize=1)
def constants():
    a,b=atan_bounds(Q(1,5),45)
    c,d=atan_bounds(Q(1,239),45)
    pi_lo=floor_q(16*a-4*d)
    pi_hi=ceil_q(16*b-4*c)
    assert Q(3)<pi_lo<pi_hi<Q(4)
    lg2_lo,lg2_hi=atanh_bounds(Q(1,3),90)
    lg3_lo,lg3_hi=atanh_bounds(Q(1,2),90)
    lo2, hi2=floor_q(lg2_lo),ceil_q(lg2_hi)
    lo3, hi3=floor_q(lg3_lo),ceil_q(lg3_hi)
    assert hi2<Q(9,10)<lo3
    lroot=isqrt((L.numerator*SCALE*SCALE)//L.denominator)
    sqlo=Q(lroot,SCALE)
    sqhi=Q(lroot+1,SCALE)
    assert sqlo**2<=L<sqhi**2
    return {'pi_lo':pi_lo,'pi_hi':pi_hi,'ln2_lo':lo2,'ln2_hi':hi2,
            'ln3_lo':lo3,'ln3_hi':hi3,'sqrtL_lo':sqlo,'sqrtL_hi':sqhi}

def pole_coefficient_exact(k:int,p:int,X):
    return -sp.Rational(k,p)*(p*p*X+sp.Rational(9,20)**2)/(k*k*X+sp.Rational(9,20)**2)

def prove_pole_projection(N:int)->list[dict]:
    if N<3: raise ValueError('N>=3')
    X=sp.Symbol('X', positive=True)
    out=[]
    for k in range(3,N+1):
        p=1 if k%2 else 2
        coeff=pole_coefficient_exact(k,p,X)
        ck=sp.Rational(k)/(k*k*X+sp.Rational(9,20)**2)
        cp=sp.Rational(p)/(p*p*X+sp.Rational(9,20)**2)
        res=sp.factor(ck+coeff*cp)
        assert res==0
        out.append({'mode':k,'pivot':p,'identity_zero_for_both_pole_moments':True})
    return out

def b_upper(k:int,p:int,c):
    """Upper bound on abs(B[p,k]), decreasing in X=pi² for k>p."""
    x=c['pi_lo']**2
    return ceil_q(Q(k,p)*(p*p*x+L*L)/(k*k*x+L*L))

def projected_mode_abs(N:int,c):
    B=[[Q(0) for _ in range(N-2)] for _ in range(N)]
    for col,k in enumerate(range(3,N+1)):
        B[k-1][col]=Q(1)
        p=1 if k%2 else 2
        B[p-1][col]=b_upper(k,p,c)
    return B

def v24_enclosure_entry(i:int,j:int,delta:Q,c)->Q:
    """Replays the explicit V24 component envelope with outward constants."""
    pi_hi=c['pi_hi']; pi_lo=c['pi_lo']; sqlo=c['sqrtL_lo'];sqhi=c['sqrtL_hi']
    def par(k):
        w=Q(k)*pi_hi/(2*L)
        lip=w/sqlo
        curvature=w*w/sqlo
        e=2*delta*delta*lip
        dtaper=2*L*curvature+Q(288)*lip
        return e,dtaper
    ei,di=par(i);ej,dj=par(j)
    Bb=2*sqhi
    q=1/(1-L/2)
    arch=(12*Bb*(ei+ej)+16*(ei*dj+ej*di))/pi_lo
    pole=2*Bb*q*q*(ei+ej)
    prime=2*c['ln2_hi']/sqlo*(ei+ej)
    return ceil_q(arch+pole+prime)

def projected_v24_error(N:int,delta:Q,c)->dict:
    """Bound ||B^T (Delta K) B||_2 <= max row sum of |B|^T E |B|.
    Conditional on the V24 model actually matching the targeted Weil matrix.
    """
    E=[[v24_enclosure_entry(i+1,j+1,delta,c) for j in range(N)] for i in range(N)]
    B=projected_mode_abs(N,c)
    P=[]
    for a in range(N-2):
        row=[]
        for b in range(N-2):
            value=sum((B[i][a]*E[i][j]*B[j][b] for i in range(N) for j in range(N)),Q())
            row.append(ceil_q(value))
        P.append(row)
    row_sums=[ceil_q(sum(row,Q())) for row in P]
    original_rows=[sum(row,Q()) for row in E]
    return {'N':N,'delta':str(delta),'projected_size':N-2,
            'projected_max_row_bound':fmt(max(row_sums),22),
            'projected_max_row_exact_fraction':str(max(row_sums)),
            'unprojected_max_row_bound':fmt(max(original_rows),22),
            'projection_coeffs_upper':[{'k':k,'pivot':1 if k%2 else 2,
                    'abs_coeff_upper':fmt(b_upper(k,1 if k%2 else 2,c),20)} for k in range(3,N+1)]}

def fmt(q:Q,digits=18):
    """Upper rounded positive decimal string, never understated."""
    if q < 0: raise ValueError('format expects nonnegative bounds')
    s=10**digits
    n= -( -q.numerator*s//q.denominator)
    return f'{n//s}.{n%s:0{digits}d}'

def mode(k,x):
    l=mp.mpf(9)/20
    return mp.sin(k*mp.pi*(x+l)/(2*l))/mp.sqrt(l) if -l<=x<=l else mp.mpf(0)

def overlap_closed(i:int,j:int,u):
    """Cij(u) = integral f_i(x) f_j(x-u) dx, 0<=u<=2L."""
    l=mp.mpf(9)/20
    if not 0<=u<=2*l: return mp.mpf(0)
    wi=i*mp.pi/(2*l);wj=j*mp.pi/(2*l)
    lo=u;hi=2*l
    if i==j:
        first=(hi-lo)*mp.cos(wj*u)
    else:
        first=(mp.sin((wi-wj)*hi+wj*u)-mp.sin((wi-wj)*lo+wj*u))/(wi-wj)
    second=(mp.sin((wi+wj)*hi-wj*u)-mp.sin((wi+wj)*lo-wj*u))/(wi+wj)
    return (first-second)/(2*l)

def sym_overlap(i:int,j:int,u):
    return (overlap_closed(i,j,u)+overlap_closed(j,i,u))/2

def pole_moment(k:int,sign:int):
    l=mp.mpf(9)/20
    w=k*mp.pi/(2*l)
    return w/mp.sqrt(l)/(w*w+mp.mpf('0.25'))*(mp.exp(-sign*l/2)-(-1)**k*mp.exp(sign*l/2))

def signed_projector(N):
    l=mp.mpf(9)/20
    X=mp.pi**2
    B=mp.matrix(N,N-2)
    for col,k in enumerate(range(3,N+1)):
        B[k-1,col]=1
        p=1 if k%2 else 2
        B[p-1,col]=-mp.mpf(k)/p*(p*p*X+l*l)/(k*k*X+l*l)
    return B

def mode_fourier(k:int,r):
    """Exact closed Fourier expression away from removable r=±w."""
    l=mp.mpf(9)/20
    w=k*mp.pi/(2*l)
    if abs(r*r-w*w)<mp.mpf('1e-20'):
        return mp.quad(lambda x: mode(k,x)*mp.exp(1j*r*x),[-l,l])
    return (w/mp.sqrt(l))*(mp.exp(-1j*r*l)-(-1)**k*mp.exp(1j*r*l))/(w*w-r*r)

def cross_parity_phase_diagnostic(dps:int=65):
    """Pointwise counterexample: V24 half-line Fi conj Fj vs V17 real sym.
    The nonzero imaginary part is numerical here, while opposite-parity real
    cancellation follows exactly from f_k(-x)=(-1)^(k+1)f_k(x).
    """
    with mp.workdps(dps):
        r=mp.mpf(1)
        prod=mode_fourier(1,r)*mp.conj(mode_fourier(2,r))
        sym=(prod+mp.conj(prod))/2
        omega=mp.re(mp.digamma(mp.mpf(1)/4+1j*r/2))-mp.log(mp.pi)
        assert abs(mp.re(prod))<mp.mpf('1e-55')
        assert abs(mp.im(prod))>mp.mpf('0.1')
        assert abs(sym)<mp.mpf('1e-55')
        def f(t):
            Om=mp.re(mp.digamma(mp.mpf(1)/4+1j*t/2))-mp.log(mp.pi)
            return mp.im(mode_fourier(1,t)*mp.conj(mode_fourier(2,t)))*Om/mp.pi
        bounded_integral=mp.quad(f,[0,1,2,4])
        return {'modes':[1,2], 'frequency':1,
          'halfline_unsym_F1_conjugate_F2':mp.nstr(prod,30),
          'correct_sym_H12':mp.nstr(sym,30),
          'omega_at_1':mp.nstr(omega,25),
          'unsym_arch_imaginary_integral_0_to_4_numeric':mp.nstr(bounded_integral,20),
          'mathematical_gate':'POINTWISE_MISMATCH: on r>0 the V24 unsymmetrized product is not the V17 real symmetrized bilinear kernel for cross parity',
          'status':'pointwise real-part zero exact by parity; nonzero imaginary value is a numerical witness, not a global interval certificate'}

def prime_numeric(N:int,dps:int=65):
    """Numeric sanity for exact support prime 2; not interval certified."""
    with mp.workdps(dps):
        u=mp.log(2)
        K=mp.matrix(N)
        err=mp.mpf(0)
        for i in range(1,N+1):
            for j in range(1,N+1):
                val=sym_overlap(i,j,u)
                if (i-j)%2:
                    err=max(err,abs(val))
                K[i-1,j-1]=-2*u/mp.sqrt(2)*val
        B=signed_projector(N)
        P=B.T*K*B
        poles=[]
        for sign in (-1,1):
            m=mp.matrix([pole_moment(k,sign) for k in range(1,N+1)])
            v=B.T*m
            poles.append(max(abs(v[i]) for i in range(N-2)))
        # Quadrature cross-check analytic overlap, independent integration.
        checks=[]
        for i,j in [(1,1),(1,2),(1,3),(2,4),(3,5),(3,3)]:
            if i>N or j>N:continue
            l=mp.mpf(9)/20
            num=mp.quad(lambda x:mode(i,x)*mode(j,x-u),[-l+u,l])
            d=abs(num-overlap_closed(i,j,u))
            checks.append({'i':i,'j':j,'abs_closed_vs_quad':mp.nstr(d,8)})
        return {'N':N,'prime_support':'ln2<2L<ln3; only prime power 2',
                'prime_sign_and_coefficient':'-2 log(2)/sqrt(2) times symmetric overlap',
                'sample_prime_entries':{'P11':mp.nstr(K[0,0],22),'P13':mp.nstr(K[0,2],22),
                                      'P12_cross_parity':mp.nstr(K[0,1],7)},
                'projected_prime_max_absolute_entry':mp.nstr(max(abs(P[i,j]) for i in range(N-2) for j in range(N-2)),20),
                'projected_prime_trace':mp.nstr(sum(P[i,i] for i in range(N-2)),20),
                'offparity_max_raw':mp.nstr(err,7),
                'two_pole_projected_residuals_numeric':[mp.nstr(x,8) for x in poles],
                'direct_integral_crosschecks':checks,
                'numeric_status':'mpmath high precision diagnostic; NOT interval certificate for prime matrix entries'}

def execute():
    c=constants()
    return {'schema':'NUMERIA.ES8P.MODE_PROJECTION_AND_V24_BOUND.1',
       'scientific_status':'exact pole cancellation; rigorous rational model error bounds CONDITIONAL on V24 normalization and analytic lemmas; numeric prime diagnostic; NOT matrix-positivity certificate or RH proof',
       'test_window':{'L':'9/20','2L':'9/10','log2_lt_2L_lt_log3_certified':True},
       'constant_enclosures_exact_fractions':{k:str(v) for k,v in c.items()},
       'pole_projection_algebra':{str(N):prove_pole_projection(N) for N in (6,8,10)},
       'model_error_projection':[projected_v24_error(N,d,c) for N in (6,8,10) for d in (Q(1,10**4),Q(1,10**5),Q(1,10**6))],
       'prime_matrix_noncertified':[prime_numeric(N) for N in (6,8,10)],
       'cross_parity_archimedean_gate':cross_parity_phase_diagnostic(),
       'normalization_reference':'v17/AUDITAR_E8K_V17_NORMALIZACAO.py; v24/README_V24.md',
       'open_gates':['Fix/prove cross-parity half-line archimedean phase handling: V17 uses Re(Fi conj Fj), V24 writes Fi conj Fj',
          'Identify each full original G+P matrix entry, signs and phases from E7/V16 against the independently defined model',
          'Certify all archimedean cross-integrals F_i F_j* for compactly supported smoothed modes',
          'Certify prime overlap entries (numerical sanity is not interval arithmetic)',
          'Compare projected eigenvalue/LDL minimum to projected smoothing error and quadrature errors',
          'Extend support and basis; no conclusion for all test functions']}

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='RESULTADOS_ES8P.json')
    a=p.parse_args();out=execute();dest=Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')
    print('ES8P: exact pole projection, conditional rational projected V24 perturbation, numeric prime kernel; NO_RH_PROOF')
    for d in out['model_error_projection']:
        if d['N']==10: print('N10',d['delta'],'projected row upper',d['projected_max_row_bound'])
if __name__=='__main__':main()
