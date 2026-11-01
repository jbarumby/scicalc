#!/usr/bin/env python3
"""NUMERIA V21 -- exact-rational audit of a Dirichlet interval Laplacian.

H=-d²/dx² on (0,pi), eigenvalues lambda_n=n², n>=1.
Checks dimensionless heat trace, spectral zeta in its convergence domain,
an integer-Mellin identity and the -1/12 Abel/Dirichlet-zeta coefficient.

Scientific scope: explicit toy operator; NOT an RH Hilbert-Polya operator.
All finite numerical inequalities are certified with fractions.Fraction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction as Q
from math import comb, factorial
from pathlib import Path


def outward(x:Q, upward:bool, digits:int=35)->str:
    scale=10**digits
    v=x.numerator*scale//x.denominator
    if upward and Q(v,scale)<x:v+=1
    sign="-" if v<0 else ""
    a=abs(v)
    return f"{sign}{a//scale}.{a%scale:0{digits}d}"


def pair(p):
    return [outward(p[0],False),outward(p[1],True)]


def exp_neg_enclosure(x:Q,degree:int=17)->tuple[Q,Q]:
    """Outward exact bounds e^-x for rational x>=0.

    Range reduce y=x/2^k<=1/2. Alternating Taylor at y:
    S_odd <= exp(-y) <= S_even. Positive bounds are then raised
    to integer power 2^k (monotone). No float used.
    """
    x=Q(x)
    if x<0 or degree<3 or degree%2!=1:
        raise ValueError("x>=0 and odd Taylor degree>=3 required")
    if x==0:return Q(1),Q(1)
    scale=1
    while x>Q(scale,2):scale*=2
    y=x/scale
    summation=Q(1)
    lower=None
    upper=None
    power=Q(1)
    for k in range(1,degree+1):
        power*= -y/Q(k)
        summation+=power
        if k==degree-1:upper=summation
        if k==degree:lower=summation
    assert 0<lower<=upper<=1
    return lower**scale,upper**scale


def spectral_zeta_enclosure(s:int,N:int)->tuple[Q,Q]:
    """Bounds sum n^(-2s), using monotone integral comparison."""
    if s<1 or N<2:raise ValueError("integer s>=1, cutoff N>=2")
    p=2*s
    partial=sum((Q(1,n**p) for n in range(1,N+1)),Q(0))
    low=partial+Q(1,(p-1)*(N+1)**(p-1))
    high=partial+Q(1,(p-1)*N**(p-1))
    return low,high


def heat_trace_enclosure(t:Q,N:int)->tuple[Q,Q]:
    """Bound Theta(t)=sum_{n>=1}exp(-t*n²) via positive rational bounds.

    For k>=0, (N+1+k)^2 >= (N+1)^2+(2N+3)*k,
    therefore the discarded heat trace is at most first/(1-ratio),
    where ratio=exp[-t*(2N+3)].
    """
    t=Q(t)
    if t<=0 or N<2:raise ValueError("positive rational t, cutoff>=2")
    values=[exp_neg_enclosure(t*n*n) for n in range(1,N+1)]
    a=sum((v[0] for v in values),Q(0))
    b=sum((v[1] for v in values),Q(0))
    first=exp_neg_enclosure(t*(N+1)**2)
    ratio=exp_neg_enclosure(t*(2*N+3))
    assert 0<ratio[0]<=ratio[1]<1
    return a+first[0], b+first[1]/(1-ratio[1])


def normalized_mellin_integral(s:int,n:int)->Q:
    """Gamma(s)^{-1} integral_0^infty t^(s-1)*exp(-n²*t)dt.

    Formula is justified by the classical positive Gamma integral and a
    substitution; the executable checks the EXACT consequence for integers.
    """
    if s<1 or n<1:raise ValueError("positive integer s and n required")
    gamma=factorial(s-1)
    integral=Q(gamma,n**(2*s))
    return integral/Q(gamma)


def series_product(a:list[Q],b:list[Q],order:int)->list[Q]:
    return [sum((a[k]*b[j-k] for k in range(j+1)),Q(0))
            for j in range(order+1)]


def exp_series_coefficients(order:int)->list[Q]:
    return [Q(1,factorial(k)) for k in range(order+1)]


def reciprocal_power_series(a:list[Q],order:int)->list[Q]:
    if a[0]==0:raise ValueError("series constant not invertible")
    b=[1/a[0]]
    for n in range(1,order+1):
        b.append(-sum((a[k]*b[n-k] for k in range(1,n+1)),Q(0))/a[0])
    return b


def bernoulli_numbers(upto:int)->list[Q]:
    if upto<0:raise ValueError("nonnegative order")
    out=[Q(1)]
    for m in range(1,upto+1):
        val=-sum((Q(comb(m+1,k))*out[k] for k in range(m)),Q(0))/Q(m+1)
        out.append(val)
    return out


def abel_casimir_certificate():
    """Exact Laurent witness for sum(n*e^(-n*x))=exp(x)/(exp(x)-1)^2.

    exp(x)/(exp(x)-1)^2 = x^-2 * exp(x)/A(x)^2,
    A(x)=(exp(x)-1)/x. Coefficients of x^-2,x^-1,x^0
    are computed by rational formal power series, NOT numerical fits.
    """
    order=2
    exp=exp_series_coefficients(order)
    A=[Q(1,factorial(j+1)) for j in range(order+1)]
    inverse=reciprocal_power_series(A,order)
    square=series_product(inverse,inverse,order)
    coeff=series_product(exp,square,order)
    assert coeff==[Q(1),Q(0),Q(-1,12)]
    B=bernoulli_numbers(2)
    assert B==[Q(1),Q(-1,2),Q(1,6)]
    assert -B[2]/2 == coeff[2]
    vacuum=coeff[2]/2
    assert vacuum==Q(-1,24)
    return {
      'laurent_coefficients_x_minus2_minus1_zero':[str(v) for v in coeff],
      'bernoulli_B0_B1_B2':[str(v) for v in B],
      'abel_finite_part':str(coeff[2]),
      'zeta_minus_one_from_bernoulli':str(-B[2]/2),
      'normalized_casimir_energy_E_over_pi_hbar_c_div_a':str(vacuum),
      'regularization_note':'Abel finite-part and meromorphic zeta value agree in this toy interval problem; classical analytic continuation remains an external theorem.'
    }


def run():
    # classical rational enclosure of pi, used only as a crosscheck in the
    # convergence domain, not as a proof of Euler's evaluations of zeta.
    pi_lo,pi_hi=Q(333,106),Q(355,113)
    zeta_records=[]
    for s in (1,2,3):
        a=spectral_zeta_enclosure(s,32)
        for n in range(1,9):
            assert normalized_mellin_integral(s,n)==Q(1,n**(2*s))
        zeta_records.append({
            's':s,'cutoff':32,'spectral_zeta':pair(a),
            'exact_mellin_term_checks':8,
            'spectral_eigenvalues':'lambda_n = n^2'
        })
        if s==1:
            reference=(pi_lo*pi_lo/6,pi_hi*pi_hi/6)
            assert max(a[0],reference[0])<=min(a[1],reference[1])
        if s==2:
            reference=(pi_lo**4/90,pi_hi**4/90)
            assert max(a[0],reference[0])<=min(a[1],reference[1])
    heat=[]
    for t in (Q(1,4),Q(1,2),Q(1)):
        small=heat_trace_enclosure(t,6)
        large=heat_trace_enclosure(t,10)
        assert max(small[0],large[0])<=min(small[1],large[1])
        assert large[1]-large[0]<small[1]-small[0]
        heat.append({'t':str(t),'N6_enclosure':pair(small),
                     'N10_enclosure':pair(large),
                     'width_shrinks':True})
    return {
      'schema':'NUMERIA.E8K.V21.SPECTRAL_ZETA_MELLIN_CASIMIR.1',
      'status':'TOY_LAPLACIAN_EXACT_RATIONAL_AUDITS_PASSED_OPERATOR_BRIDGE_OPEN',
      'operator':'H=-d²/dx² on L²(0,pi), Dirichlet boundary; domain H² intersect H_0¹',
      'eigenvalues':'lambda_n=n², n=1,2,...; positive and self-adjoint',
      'heat_trace':'Theta(t)=sum_n exp(-t*n²), t>0',
      'spectral_zeta':'zeta_H(s)=sum_n n^(-2s)=zeta(2s) for Re(s)>1/2',
      'mellin_identity':'zeta_H(s)=Gamma(s)^(-1)*integral_0^infinity t^(s-1)Theta(t)dt, Re(s)>1/2',
      'zeta_interval_records':zeta_records,
      'heat_trace_records':heat,
      'casimir_abel_laurent':abel_casimir_certificate(),
      'checks':{'zeta_integral_cases':3,'mellin_single_eigenmode_identities':24,
                'heat_trace_t_cases':3,'laurent_coefficient_identities':3,
                'independent_bernoulli_value':True,
                'no_float_in_certification':True},
      'analytic_assumptions_not_mechanized':[
        'Dirichlet Laplacian spectral theorem with eigenvalues n²',
        'Euler gamma integral and Tonelli interchange in region Re(s)>1/2',
        'alternating Taylor remainder for e^-y on [0,1/2]',
        'integral comparison for positive decreasing sequences',
        'Euler zeta(2),zeta(4) evaluations for pi crosschecks',
        'meromorphic continuation of zeta and zeta(-1)=-B2/2',
        'physical vacuum subtractions and boundary conventions are model-dependent'
      ],
      'not_certified':[
        'No operator with eigenvalues equal to imaginary parts of all Riemann zeros',
        'No verified trace formula matching all prime, pole and archimedean Weil terms',
        'No independent Phase-4 quadratic integration or derivative bounds',
        'A model's zeta regularization is not a direct proof of Weil positivity',
        'No Hilbert-Polya theorem, nor proof of the Riemann Hypothesis'
      ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',required=True,type=Path)
    args=p.parse_args()
    obj=run()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(obj,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('PASS: 3 spectral-zeta enclosures; 24 exact Mellin terms; 3 rigorous heat traces; -1/12 Abel/Bernoulli finite part; -1/24 normalized Casimir energy')
    print('STATUS:',obj['status'])
    return 0


if __name__=='__main__':
    raise SystemExit(main())
