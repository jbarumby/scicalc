#!/usr/bin/env python3
"""Exact Toeplitz Gram compression for local three-centre moment annihilators.

This is an *algebraic transformation*, not numerical computation of Weil's
functional. q may be rational in unit tests; physically q=exp(d/2)>1.
"""
from fractions import Fraction


def gram_lag_filtered(w: dict, lag: int, c) -> Fraction:
    """For symmetric real translation invariant Gram w_k=w_-k.

    Local filter a_{j}=g_{j-1}-c*g_j+g_{j+1}, c=2*cosh(d/2).
    H_ij = h_{i-j} where
    h_k=w_{k-2}+w_{k+2}-2c*(w_{k-1}+w_{k+1})+(c*c+2)*w_k.
    """
    c=Fraction(c)
    def get(k): return Fraction(w[abs(k)])
    return (get(lag-2)+get(lag+2)
            -2*c*(get(lag-1)+get(lag+1))+(c*c+2)*get(lag))


def explicit_congruence_entry(w: dict, i: int, j: int, c) -> Fraction:
    """Compute D^T W D directly from 9 rational products as a crosscheck."""
    c=Fraction(c)
    cs={-1:Fraction(1),0:-c,1:Fraction(1)}
    out=Fraction(0)
    for a,ca in cs.items():
        for b,cb in cs.items():
            out += ca*cb*Fraction(w[abs((i+a)-(j+b))])
    return out


def verify_toeplitz_congruence(w: dict, n: int, q) -> dict:
    """Compare direct transformed matrix versus compressed filtered lag.

    Input original W is n*n Toeplitz; output H has dimension n-2.
    w contains lags 0,..,n-1. q arbitrary positive rational, c=q+1/q.
    """
    if n < 3:
        raise ValueError('n>=3')
    q=Fraction(q)
    if q<=0 or q==1: raise ValueError('q>0 and q !=1')
    c=q+1/q
    for k in range(n):
        if k not in w: raise ValueError('missing lag')
    pairs=0
    for i in range(1,n-1):
        for j in range(1,n-1):
            p=explicit_congruence_entry(w,i,j,c)
            f=gram_lag_filtered(w,i-j,c)
            if p != f:
                raise AssertionError((i,j,p,f))
            pairs+=1
    return {'dimension':n-2,'checked_entries':pairs,'q':str(q),
            'filter_c':str(c),'toeplitz_lags_needed':n-2,
            'source_lags_needed':n,'status':'EXACT_CONGRUENCE_ALGEBRA'}

if __name__=='__main__':
    n=9
    w={k:Fraction(1,2)**k for k in range(n)}
    print(verify_toeplitz_congruence(w,n,Fraction(3,2)))
