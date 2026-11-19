#!/usr/bin/env python3
"""Exact coefficient audit for the Bernstein/GL60 degreewise error budget.

The standard Bessel integral, Legendre bound and Chebyshev/Gauss theorems are
used in the accompanying mathematical proof, NOT implemented or formalized here.
The script verifies all rational constants and max entry / row error budgets.
"""
import argparse,json,math
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path
from p237_symbol_majorant import check_certificate,need

L=F(161,200);N=377;PANELS=1008

def ceil_sqrt(q):
    if q<0:raise ValueError('negative')
    k=math.isqrt(q.numerator//q.denominator)
    if k*k*q.denominator<q.numerator:k+=1
    need((k-1)*(k-1)*q.denominator < q.numerator <=k*k*q.denominator,
         'ceil sqrt invalid')
    return k


@lru_cache(maxsize=1)
def coefficients():
    vals=[];doublefac=1
    for n in range(754):
        doublefac*=2*n+1
        b=min(F(2),F(2*406**n,doublefac))
        normalizer=ceil_sqrt(2*L*(2*n+1))
        vals.append(normalizer*b)
    return tuple(vals)


@lru_cache(maxsize=1)
def build():
    symbol=check_certificate()
    need(symbol['symbol_minus_beta_bound_strict_less']=='274/7','wrong symbol coefficient')
    # The usual Chebyshev coefficient estimate gives for k>=1:
    # |c_k| <= 2*M*3^-k.  Tail k>=120 is 3*M*3^-120.
    # GL60 is exact on degree <=119 and its positive weights sum to 2.
    # Thus |int_{-1}^1 g - GL60(g)| <= 4*3*M*3^-120 = 12*M/3^120.
    # For g(s)=u_i(t(s))*u_j(t(s))*(Psi(t(s))-beta)/(4*pi),
    # pi>3, |Psi-beta|<50, |u_n|<=a_n => M < (50/12)a_i*a_j.
    # So |error on one t-panel| < (50/3^120)*a_i*a_j.
    C=F(PANELS*50,3**120)
    vals=coefficients()
    maxbasis=max(vals)
    need(maxbasis <=84,'P236 max basis exceeded')
    data={}
    for sec,parity in [('even',0),('odd',1)]:
        a=vals[parity::2]
        need(len(a)==N,'wrong parity count')
        s=sum(a,F())
        row_bound=C*max(a)*s
        need(row_bound<F(1,10**46), 'degreewise row error not small')
        diag_last=C*a[-1]*a[-1]
        data[sec]={'modes':N,'max_basis_bound_exact':str(max(a)),
          'sum_basis_bounds_exact':str(s),'row_sum_GL60_error_exact':str(row_bound),
          'row_sum_below_1e-46':row_bound<F(1,10**46),
          'last_mode_degree':752+parity,'last_mode_GL60_diagonal_error_exact':str(diag_last),
          'last_mode_GL60_diagonal_error_below_1e-200':diag_last<F(1,10**200)}
    return {'schema':'NUMERIA.ES8P237.GL60_COMPLEX_ANALYTIC_RATIONAL_GATE.1',
      'rho':3,'Gauss_legendre_nodes':60,'degree_of_exactness':119,
      'chebyshev_truncation_error_coefficient_exact':'3 / 3^120',
      'GL60_integral_minus_quadrature_error_coefficient_exact':'12 / 3^120',
      'rescaled_panel_jacobian':'1/4','reciprocal_pi_conservative_upper':'1/3',
      'symbol_complex_bound_strict_less':symbol['symbol_minus_beta_bound_strict_less'],
      'per_panel_error_coefficient_strict_less':'50 / 3^120',
      'all_panels_error_coefficient_strict_less_exact':str(C),
      'error_added_once_after_summation':True,
      'j_n_complex_bound':'min(2, 2*406^n/(2n+1)!!)',
      'j_n_requires_classical_inputs':['Legendre plane-wave identity',
         'Rodrigues formula with integration by parts',
         '|P_n(u)|<=1 on [-1,1]',
         'e^|Im(Lt)| < 2 for |Im(Lt)|<1/2'],
      'GL60_requires_classical_inputs':['holomorphic Chebyshev coefficient bound on E3',
         'Gauss Legendre exactness degree <=119',
         'Gauss Legendre weights positive summing to 2'],
      'sectors':data,'Leibniz_Gauss_analytic_lemmas_formalized_in_Lean':False,
      'Weil_Zhu_analytic_bridge_proved':False,'LDL_377x377_executed':False,'RH_proved':False}


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();d=build();a.out.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
    print('PASS GL60 degreewise rational budget and both row-sum bounds; classical analytic theorems imported')

if __name__=='__main__':main()
