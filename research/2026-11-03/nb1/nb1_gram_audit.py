#!/usr/bin/env python3
"""NUMERIA NB1: reproducible, rationally certified UPPER bounds for finite-dimensional Nyman–Beurling approximants.

This certificate is unconditional for each tested fixed rational coefficient vector.
It does NOT certify RH or the global minimum, nor convergence as N -> infinity.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from functools import lru_cache
from math import log, isfinite
import numpy as np

D = 10**26
E = F(1, 10000)

def floor_grid(q): return (q.numerator*D)//q.denominator

def ceil_grid(q): return -((-q.numerator*D)//q.denominator)

def cuts(N, eps=E):
    pts={eps,F(1)}
    for n in range(1,N+1):
        k=1
        while True:
            x=F(1,n*k)
            if x<=eps: break
            if x<1: pts.add(x)
            k+=1
    return sorted(pts)

@lru_cache(maxsize=90000)
def log_rational_enclosure(r):
    """Exact rational atanh series with a rigorous, positive remainder bound (r in (1, 2])."""
    assert 1<r<=2, r
    z=(r-1)/(r+1)
    z2=z*z
    p=z
    partial=F(0)
    k=0
    while True:
        partial += 2*p/F(2*k+1)
        p*=z2
        # Remaining tail from k+1, denominators >=2k+3
        tail=2*p/(F(2*k+3)*(1-z2))
        if tail < F(1,10**30):
            break
        k+=1
        if k>120: raise RuntimeError('unexpected slow series')
    return F(floor_grid(partial),D),F(ceil_grid(partial+tail),D)

def approx_solve(N,pts):
    theta=np.array([1/n for n in range(1,N+1)],dtype=float)
    G=np.outer(theta,theta) # integral 1 to infinity
    b=np.full(N,float(E)/2)
    G+=float(E)*np.full((N,N),0.25)
    np.fill_diagonal(G,G.diagonal()+float(E)/12)
    # approximate unknown contribution x in [0,E] with equidistribution heuristic ONLY for seed
    for u,v in zip(pts[:-1],pts[1:]):
        a=float(u); z=float(v)
        m=(a+z)/2
        k=np.floor(theta/m)
        k=np.where(np.isclose(theta/m,np.round(theta/m),atol=1.e-13,rtol=0), k, k)
        aa=z-a; bb=log(z/a); cc=1/a-1/z
        b+=theta*bb-k*aa
        G+=cc*np.outer(theta,theta)-bb*(np.outer(theta,k)+np.outer(k,theta))+aa*np.outer(k,k)
    G=(G+G.T)/2
    cond=float(np.linalg.cond(G))
    try: c=np.linalg.solve(G,b)
    except np.linalg.LinAlgError: c=np.linalg.lstsq(G,b,rcond=None)[0]
    # use 8-decimal rational coefficients for subsequent exact interval check
    cc=[F(f'{v:.8f}') for v in c]
    if not all(isfinite(float(x)) for x in cc):raise ValueError('nonfinite')
    return cc,cond,float(1-2*np.dot(c,b)+c@G@c)

def certify(N,pts,coef):
    assert len(coef)==N
    B= -sum((coef[n-1]*F(1,n) for n in range(1,N+1)),F(0))
    checked_upper=0
    checked_lower=0
    for u,v in zip(pts[:-1],pts[1:]):
        m=(u+v)/2
        a=F(1)+sum((coef[n-1]*((F(1,n)/m).numerator//(F(1,n)/m).denominator) for n in range(1,N+1)),F(0))
        rr=a*a*(v-u)+B*B*(1/u-1/v)
        Llo,Lhi=log_rational_enclosure(v/u)
        cross=2*a*B
        if cross>=0:
            lower=rr+cross*Llo; upper=rr+cross*Lhi
        else:
            lower=rr+cross*Lhi; upper=rr+cross*Llo
        checked_lower += floor_grid(lower)
        checked_upper += ceil_grid(upper)
    low0=F(0)
    upper0=E*(F(1)+sum(map(abs,coef)))**2
    tail=B*B
    # the certified integrand on (E,1) has nonnegative integral, but interval lower may mildly negative due rounding
    lo=max(F(0),F(checked_lower,D))+tail+low0
    hi=F(checked_upper,D)+tail+upper0
    assert hi>=lo>=0
    return {
      'n':N,'coefficients':[str(x) for x in coef], 'cutoff_epsilon':str(E),
      'breakpoints':len(pts),
      'squared_error_lower_bound':str(lo),'squared_error_upper_bound':str(hi),
      'squared_error_upper_decimal':float(hi),
      'tail_x_greater_than_1_exact':str(tail),
      'head_x_less_than_epsilon_upper_bound':str(upper0),
      'interval_checked_error_upper_bound':str(F(checked_upper,D)),
      'log_enclosures_used':log_rational_enclosure.cache_info().currsize,
      'certification_scope': 'Fixed rational coefficients only; no optimality/global convergence/RH claim'
    }

def run():
    out=[]
    for N in (6,8,10):
        p=cuts(N)
        c,cond,approx=approx_solve(N,p)
        row=certify(N,p,c)
        row['normal_matrix_condition_estimate_float']=cond
        row['approx_minimum_sqr_error_float_unverified']=approx
        out.append(row)
        print(N, 'intervals=',len(p)-1,'cond=',cond,'err_upper=',float(row['squared_error_upper_decimal']), 'coef=',[float(q) for q in c],flush=True)
    payload={'schema':'NUMERIA.NB1.FINITE_RATIONAL_UPPER_BOUNDS.1',
     'status':'CERTIFIED_FINITE_UPPER_BOUNDS_ONLY',
     'space':'L2((0,infinity), dx)',
     'target':'indicator_{0<x<1}',
     'basis':'rho_n(x) = fractional_part(1/(n*x)), 1 <= n <= N',
     'identity_for_x_ge_1':'integral_1^infty (sum c_n rho_n(x))^2 dx = (sum c_n/n)^2',
     'basis_nonzero_at_large_x':'rho_n(x)=1/(n*x) for x>1',
     'log_enclosure':'rational arctanh series with explicit upper bound; every contribution summed with integer outward rounding',
     'proof_limitation':'A sequence of finite upper bounds, without uniform rate or all N, cannot prove RH.',
     'n_results':out}
    data=json.dumps(payload,ensure_ascii=False,indent=2)+'\n'
    from pathlib import Path
    target=Path(__file__).with_name('RESULTADOS_NB1.json')
    target.write_text(data,encoding='utf-8')
    print('wrote',target,'sha256',hashlib.sha256(data.encode()).hexdigest())
if __name__=='__main__': run()
