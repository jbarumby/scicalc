#!/usr/bin/env python3
"""ES8-P10: analytic-error proof witnesses and sharper Phase-4 interpolation.

All certificate outputs use fractions and/or directed original 90-digit
interval arithmetic.  This is a proof-checking aid, not a Lean formalization.

A major improvement is that the *integral* of the absolute degree-16 nodal
polynomial replaces a product of sixteen per-factor maxima.  Existing node
inclusions and Cauchy majorants are kept unchanged.

The digamma-series tail follows from an elementary two-sided coth inequality
proved by positive Taylor coefficients, applied to a positive Laplace kernel.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from decimal import Decimal
from fractions import Fraction as Q
from math import factorial
from pathlib import Path

from es8p6_independent_audit import (
    L, H, RHO, NODES, CELLS, D, ALPHA, DIRECTION,
    independent_total_error, interpolation_product, independent_cauchy_bound,
    sharp_tail_upper,
)
from es8p9_projected import projected_matrix

BASE=Path(__file__).resolve().parent

def newton_polynomial(nodes=16):
    """Exact coefficients of prod_{j=0}^{15}(x-j/15), ascending order."""
    xs=[Q(j,nodes-1) for j in range(nodes)]
    c=[Q(1)]
    for x in xs:
        z=[Q(0)]*(len(c)+1)
        for k,a in enumerate(c):
            z[k]-=a*x
            z[k+1]+=a
        c=z
    return xs,c

def exact_polynomial_integral(coeffs,a,b):
    return sum((c*(b**(i+1)-a**(i+1))/Q(i+1)
                for i,c in enumerate(coeffs)),Q(0))

def absolute_newton_integral():
    """Exact |omega16| integral via sign-constant intervals between its roots."""
    xs,coef=newton_polynomial()
    cells=[exact_polynomial_integral(coef,a,b) for a,b in zip(xs,xs[1:])]
    # Signs alternate because each root is simple, no use of floating point.
    assert all(x != 0 for x in cells)
    assert all(cells[i]*cells[i+1]<0 for i in range(len(cells)-1))
    return sum(map(abs,cells),Q(0)), cells

def original_cell_majorant():
    """Exact old constant = prod_j max(j/15,1-j/15)."""
    return interpolation_product()

def refined_interpolation_error():
    """Exact integrated 16th-derivative/Cauchy error with same M_c as P6.

    On a cell, for real analytic q and x in [0,1],
       |q-P15|(c+H*x) <= M_c*(H/RHO)^16*|omega16(x)|.
    Apply |int|<=int|.|, then sum over cells; M_c is an exact rational
    upper bound from ES8-P6.  The result is strictly smaller than P6.
    """
    a,_=absolute_newton_integral()
    old=independent_total_error()
    new=old*a/original_cell_majorant()
    assert new>0 and new<old
    return new

def coth_coefficients(n):
    """For power t^(2n+1), numerator coefficients × (2n+1)!.

      U(t)=(1+t²/3)sinh(t)-t cosh(t) >= 0.
      V(t)=t cosh(t)-(1+t²/3-t⁴/45)sinh(t) >= 0.

    The coefficient signs establish, for ALL real t>0:
       1 <= t*coth(t) <= 1+t²/3
       t*coth(t) >= 1+t²/3-t⁴/45.
    Since sinh(t)>0, the inequalities follow termwise.
    """
    if n<0: raise ValueError
    # For powers t, t^3, t^5, ...
    # U coefficient = 4n(n-1)/(3(2n+1)!) for n>=1.
    # V coefficient = 16 n(n-1)(n-2)(n+2)/(45(2n+1)!) for n>=2.
    u=Q(0) if n==0 else Q(4*n*(n-1),3*factorial(2*n+1))
    v=Q(0) if n<2 else Q(16*n*(n-1)*(n-2)*(n+2),45*factorial(2*n+1))
    return u,v

def cm_function(u,b2):
    """Positive, completely monotone f(u)=b²/[u(u²+b²)] for u>0."""
    u,b2=Q(u),Q(b2)
    if u<=0 or b2<0:raise ValueError
    return b2/(u*(u*u+b2))

def em_corrections(u,b2):
    """B2=-f'(u)/12 >=0, B4=-f'''(u)/720 >=0 (EXACT)."""
    u,b2=Q(u),Q(b2)
    if u<=0 or b2<0:raise ValueError
    den=u*u+b2
    b2corr=b2*(3*u*u+b2)/(12*u*u*den**2)
    b4corr=b2*(10*u**6+5*b2*u**4+4*b2**2*u**2+b2**3)/(120*u**4*den**4)
    # Independently rearranged historical B4 expression. Direct exact check.
    rearranged=(1/u**4-(u**4-6*u*u*b2+b2*b2)/den**4)/120
    assert b4corr==rearranged and 0<=b4corr and 0<=b2corr
    return b2corr,b4corr

def validate_cauchy_geometry():
    assert (D,RHO,H,ALPHA,DIRECTION)==(Q(7,20),Q(3,10),Q(1,10),Q(3,40),Q(137,200))
    assert 0<ALPHA and D<RHO+H and D<Q(1,2)
    assert 0<L*D<1 and DIRECTION>0
    # The cell supremum majorant has the following strictly positive coefficients.
    M0=independent_cauchy_bound(Q(1,20))
    M1=independent_cauchy_bound(Q(3,20))
    assert 0<M0<M1
    return {'disk_radius':str(RHO),'enlarged_disk_radius':str(D),
            'digamma_real_part_lower':str(ALPHA),
            'fourier_exp_majorant_denominator':str(DIRECTION),
            'sample_cauchy_majorant_cell0':str(M0)}

def new_intervals(ref_path=BASE/'RESULTADOS_ES8P9_INTEGRAIS.json'):
    d=json.loads(ref_path.read_text(encoding='utf-8'))
    assert d['status']=='INTEGRATED_30_ENTRIES_CONDITIONALLY_CERTIFIED'
    assert d['scope']['evaluations']==768000
    assert len(d['entries'])==30
    e_old=independent_total_error();e_new=refined_interpolation_error()
    result={};old_included=0
    for key,v in d['entries'].items():
        lo,hi=(Q(Decimal(a)) for a in v['independent_nodal_0_160'])
        old_lo,old_hi=(Q(Decimal(a)) for a in v['independent_finite_with_cauchy'])
        # Both old and nodal endpoints were rounded separately to 24 places.
        tol=Q(2,10**24)
        assert old_lo<=lo-e_old+tol and hi+e_old<=old_hi+tol
        # Exact rational interval endpoint -> 24-digit decimal *outward*.
        anew=decimal_out(lo-e_new)
        bnew=decimal_out(hi+e_new,upper=True)
        assert Q(Decimal(anew))<=lo-e_new<=hi+e_new<=Q(Decimal(bnew))
        if old_lo<=Q(Decimal(anew)) and Q(Decimal(bnew))<=old_hi:
            old_included+=1
        result[key]=dict(v)
        result[key]['independent_finite_with_cauchy']=[anew,bnew]
        t=sharp_tail_upper(v['i'],v['j'])
        result[key]['independent_global_with_tail']=[decimal_out(Q(Decimal(anew))-t),
                                                     decimal_out(Q(Decimal(bnew))+t,upper=True)]
    assert old_included==30
    return result

def decimal_out(x,precision=24,upper=False):
    x=Q(x);s=10**precision
    n=-((-x.numerator*s)//x.denominator) if upper else x.numerator*s//x.denominator
    sign='-' if n<0 else ''
    a=abs(n)
    return sign+str(a//s)+'.'+str(a%s).zfill(precision)

def audit():
    old=independent_total_error();new=refined_interpolation_error()
    absval,pieces=absolute_newton_integral()
    assert len(pieces)==15 and old/new>Q(6_000_000)
    geo=validate_cauchy_geometry()
    # Coefficient proof of universal coth inequalities; verifying symbolic pattern
    # n=0..64 is a regression control, NOT the justification for all n.
    assert all(min(coth_coefficients(n))>=0 for n in range(65))
    for u in [Q(257,4),Q(513,4),Q(1,4),Q(128)]:
        for b2 in [Q(0),Q(1,4),Q(100),Q(6400)]:
            em_corrections(u,b2)
    outputs=new_intervals()
    mats=[projected_matrix(N,outputs) for N in (6,8,10)]
    assert all(x['certified_M_minus_0p05I_positive'] for x in mats)
    return {
        'schema':'NUMERIA.ES8P10.ANALYTIC_LEMMA_AUDIT.1',
        'status':'DIGAMMA_EM_TAIL_ANALYTIC_PROOF_AND_SHARPER_INTERPOLATION_REPLAY',
        'old_interpolation_error_exact':str(old),
        'old_interpolation_error_outward':decimal_out(old,upper=True),
        'new_interpolation_error_exact':str(new),
        'new_interpolation_error_outward':decimal_out(new,upper=True),
        'error_reduction_factor_exact':str(old/new),
        'error_reduction_factor_floor':(old/new).numerator//(old/new).denominator,
        'old_sup_product_exact':str(original_cell_majorant()),
        'absolute_newton_integral_exact':str(absval),
        'cauchy_disk_geometry':geo,
        'digamma_tail_statement':
          'For u>0,b>=0: integral_u^infty f+f(u)/2-fprime(u)/12+ftripleprime(u)/720 <= sum_{n>=0}f(u+n) <= integral_u^infty f+f(u)/2-fprime(u)/12; f(u)=b^2/(u(u^2+b^2)).',
        'digamma_tail_proof':[
            'f(u)=int_0^infty exp(-u*t)*(1-cos(b*t))dt, positive Laplace kernel',
            'U(t)=(1+t^2/3)sinh(t)-t*cosh(t): coefficient 4n(n-1)/(3(2n+1)!), n>=1',
            'V(t)=t*cosh(t)-(1+t^2/3-t^4/45)sinh(t): coefficient 16n(n-1)(n-2)(n+2)/(45(2n+1)!), n>=2',
            'Therefore t/12-t^3/720 <= 1/(1-exp(-t))-1/t-1/2 <=t/12 for every t>0',
            'Multiply by exp(-u*t)*(1-cos(b*t)) >=0 and integrate by Tonelli',
            'Integral and corrections match historical Euler-Maclaurin digamma expressions; B4 is nonnegative by exact polynomial expansion'
        ],
        'interpolation_proof':[
            'For real q_ij on each cell: q(x)-P15(x)=q^(16)(xi)/16! * product_j(x-x_j), xi real',
            'Cauchy derivative bound |q^(16)| <=16!*M_c/RHO^16 on disks enclosing each real point',
            'Exact absolute polynomial integration: sum_{j=0}^{14}|int_{j/15}^{(j+1)/15} product_{k=0}^{15}(x-k/15)dx|',
            'Replace maximum-product bound with exact absolute integral, retaining M_c and nodal intervals unchanged',
        ],
        'improved_finite_entries_included_in_p9':len(outputs),
        'projected_matrices':mats,
        'assumptions':[
            'Digamma partial-fraction expansion (classical, DLMF 5.7.6)',
            'Holomorphy and Cauchy derivative formula on enlarged cell disks',
            'V17 Fourier normalization and arithmetic prime-2 term inherited from P9',
            'Original ES8-P9 nodal interval endpoints rigorously include exact values',
            'E7 sharp tail bound inherited from P6 and confirmed by rational inequalities',
        ],
        'not_established':[
            'Complete mechanization in Lean/Coq of all special-function identities and endpoint rounding code',
            'An independent theorem of all zero-free regions or universal Weil positivity',
            'Uniform positivity as N increases without bound; Riemann Hypothesis',
        ],
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=BASE/'RESULTADOS_ES8P10.json')
    a=parser.parse_args();d=audit()
    a.out.write_text(json.dumps(d,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(d['status'])
    print('Old/new errors:', d['old_interpolation_error_outward'],d['new_interpolation_error_outward'])
    for m in d['projected_matrices']:
        print('N',m['N'],'positive minus 0.05I:',m['certified_M_minus_0p05I_positive'], 'row delta:',m['perturbation_row_sum_upper_approx'])
if __name__=='__main__':main()
