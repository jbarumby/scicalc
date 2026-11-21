"""P260: exact formal-series audit of pole-vector normalization and parity.

Checks the Legendre/Bessel identity underlying +2 even / -2 odd
matrix pole term, CONDITIONAL on the Weil pole functional 2*ell_+*ell_-.
Does not rederive that functional from a contour integral.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import argparse
import json

def require(ok, text):
    if not ok:
        raise ArithmeticError(text)

def legendre_coefficients(n):
    if type(n) is not int or n<0:
        raise ValueError('n must be a nonnegative integer')
    if n==0:return (F(1),)
    p=(F(1),)
    q=(F(0),F(1))
    for m in range(1,n):
        nxt=[F(0)]*(len(q)+1)
        for j,a in enumerate(q):nxt[j+1]+=(2*m+1)*a/F(m+1)
        for j,a in enumerate(p):nxt[j]-=m*a/F(m+1)
        p,q=q,tuple(nxt)
    return q

def legendre_moment(n,k):
    if type(k) is not int or k<0:
        raise ValueError('k must be nonnegative integer')
    return sum((c*F(2,j+k+1) for j,c in enumerate(legendre_coefficients(n))
                if (j+k)%2==0),F(0))

def modified_spherical_bessel_coefficient(n,r):
    """[z^(n+2r)] i_n(z), from I_(n+1/2) Gamma series."""
    if any(type(x) is not int or x<0 for x in (n,r)):
        raise ValueError('indices must be nonnegative integers')
    return F(2**(n+1)*factorial(n+r+1),
             factorial(r)*factorial(2*n+2*r+2))

def audit(max_degree=20,terms_per_degree=12):
    if max_degree<1 or terms_per_degree<1:
        raise ValueError('positive truncation parameters required')
    total=0
    for n in range(max_degree):
        poly=legendre_coefficients(n)
        require(all(c==0 for j,c in enumerate(poly) if j%2!=n%2),
                f'Legendre parity failed degree {n}')
        for r in range(terms_per_degree):
            k=n+2*r
            lhs=legendre_moment(n,k)/factorial(k)
            rhs=2*modified_spherical_bessel_coefficient(n,r)
            require(lhs==rhs,f'Gamma/Bessel coefficient mismatch n={n} r={r}')
            require((-1)**k*lhs==(-1)**n*rhs,
                    f'reflection parity mismatch n={n} r={r}')
            total+=1
    return {
        'schema':'NUMERIA.ES8P260.POLE_LEGENDRE_NORMALIZATION.1',
        'exact_series_coefficients_compared':total,
        'first_n_degrees':max_degree,
        'terms_per_degree':terms_per_degree,
        'legendre_identity':'integral_-1^1 P_n(u) exp(z*u) du = 2*i_n(z)',
        'normalized_pole_vector':'ell_plus(phi_n) = sqrt(2*L*(2*n+1))*i_n(L/2)',
        'parity':'ell_minus(phi_n) = (-1)^n ell_plus(phi_n)',
        'pole_quadratic_assumption':'2*ell_plus(f)*ell_minus(f) for real f',
        'corresponding_even_matrix_factor':'+2',
        'corresponding_odd_matrix_factor':'-2',
        'entire_Weil_explicit_formula_pole_identification_proved_here':False,
        'full_Weil_positivity_proved':False,
        'Riemann_hypothesis_proved':False,
    }

def main():
    cli=argparse.ArgumentParser()
    cli.add_argument('--out',type=Path,required=True)
    args=cli.parse_args()
    x=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('P260 pole normalization conditional',x['exact_series_coefficients_compared'])

if __name__=='__main__':main()
