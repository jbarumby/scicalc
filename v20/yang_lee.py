#!/usr/bin/env python3
"""Finite Bose-mode occupation cutoff: zeros of 1+z+...+z^K.

Exact polynomial recurrence and complex-root classification (via algebraic
identity in documentation). Numeric roots, if shown, are illustrations only.
"""
from fractions import Fraction as Q


def polynomial_coeffs(K:int)->list[int]:
    if K<1:raise ValueError('K>=1')
    return [1]*(K+1)


def multiply_z_minus_one(coeff:list[int]):
    # constant term up to degree len(coeff)
    out=[0]*(len(coeff)+1)
    for i,c in enumerate(coeff):out[i]-=c;out[i+1]+=c
    return out


def closest_root_distance_upper(K:int)->Q:
    if K<1:raise ValueError('K>=1')
    # |e^(2 pi i/(K+1)) -1| = 2 sin(pi/(K+1)) < 2 pi/(K+1) < 44/[7(K+1)]
    return Q(44,7*(K+1))


def run():
    seq=[]
    last=None
    for K in [2,3,7,15,31,63,127,255]:
        cs=polynomial_coeffs(K)
        assert multiply_z_minus_one(cs)==[-1]+[0]*K+[1]
        assert sum(cs)==K+1 and all(c>0 for c in cs)
        bound=closest_root_distance_upper(K)
        if last is not None: assert bound<last
        last=bound
        seq.append({'K':K,'z1_value':K+1,'bound_on_closest_root_distance':str(bound),
                    'positive_real_zero':False})
    return {'status':'FINITE_BOSE_MODE_YANG_LEE_ALGEBRA_PASS','cases':seq,
            'exact_polynomial_identity':'(z-1)*sum(z^n,n=0..K)=z^(K+1)-1',
            'zeros':'all (K+1)-th roots of unity except z=1',
            'conditional_physics':'Toy single zero-energy Bose mode with finite K occupation cutoff; not the interacting gas of Wang (2004)',
            'limit_caution':'zero pinch near fugacity z=1 is not an identity between Yang-Lee and Riemann zeros'}
