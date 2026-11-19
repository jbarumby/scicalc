#!/usr/bin/env python3
"""Independent exact-rational audit of the P237 GL60 Bernstein error budget.

Does not prove complex analyticity or the lemma bounding the symbol by 50;
those are imported mathematical hypotheses, separately documented. Given
these bounds on the Bernstein ellipse E_3, verifies the 60-node Chebyshev
quadrature error coefficient and per-row errors over 1008 panels.
"""
from __future__ import annotations
from fractions import Fraction as Q
from math import isqrt
import json, argparse

L=Q(161,200)
N=377
RHO=3
B=406
QGL=60
PANELS=1008

# series Chebyshev coefficients |c_k|<= 2M rho^-k,
# tail from k=2*QGL is 2 M rho^(-2QGL)/(1-1/rho).
# Gauss and exact integral are both bounded linear functionals <=2 ||f||inf.
# error <=4* [2 M rho^(-2QGL)/(1-1/rho)] = 12 M/rho^120.
# g(s)=(Phi(t(s))-beta)*u_i*u_j/(4pi),
# assuming |Phi-beta|<50 and pi>3, M < (50/12) a_i*a_j
# thus per-panel error < 50 a_i a_j/3^120.
COEFF=Q(PANELS*50, RHO**(2*QGL))

def ceil_sqrt_q(a:Q)->int:
    assert a>=0
    # first integer with k²>=a using exact integers
    k=isqrt(a.numerator//a.denominator)
    while Q(k*k)<a:k+=1
    if k and Q((k-1)*(k-1))>=a:raise ArithmeticError('bad square root ceiling')
    return k

def df(n:int)->int:
    p=1
    for j in range(1,2*n+2,2):p*=j
    return p

def amplitudes(par:int)->list[Q]:
    result=[]
    for i in range(N):
        n=2*i+par
        amp=ceil_sqrt_q(2*L*(2*n+1))
        jbound=min(Q(2), Q(2*B**n,df(n)))
        result.append(Q(amp)*jbound)
    return result

def run():
    result={'schema':'NUMERIA.ES8P245C.INDEPENDENT_GL60_BERNSTEIN_BUDGET.1',
      'inputs':{'L':'161/200','Tsharp':504,'GL_nodes':60,'panels':1008,'Bernstein_rho':3,
      'complex_symbol_bound_imported':'|Phi_L(t)-beta| < 50 on every E_3',
      'complex_Bessel_bound_imported':'|j_n(Lt)| < min(2,2*406^n/(2n+1)!!)'},
      'Chebyshev_tail_bound':'3*M/3^120',
      'quadrature_error_bound':'12*M/3^120',
      'per_panel_error_coeff':'50/3^120',
      'all_1008_panels_error_coeff':str(COEFF),
      'max_row_norm_target':'4/10^47',
      'sectors':{},
      'complex_analytic_hypotheses_independently_proved':False,
      'GL60_integrand_evaluations_certified_in_this_script':False,
      'RH_proved':False}
    for sec,par in [('even',0),('odd',1)]:
        aa=amplitudes(par)
        mx=max(aa);summ=sum(aa,Q(0))
        bound=mx*summ*COEFF
        if bound>=Q(4,10**47):raise ArithmeticError(f'{sec} row budget fails')
        result['sectors'][sec]={
           'modes':N,'degree_max':2*(N-1)+par,
           'max_amplitude':str(mx),'sum_amplitudes':str(summ),
           'max_row_error_exact':str(bound),
           'max_row_error_less_than_4e-47':True,
           'max_row_error_approx':f'{float(bound):.12e}'}
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args()
    result=run();open(a.out,'w').write(json.dumps(result,indent=2)+'\n')
    for k,v in result['sectors'].items():print(k, v['max_row_error_approx'])
