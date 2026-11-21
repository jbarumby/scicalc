"""Exact algebra of the Cor.1 v2 exponents (math preflight, no certificate)."""
from fractions import Fraction
import mpmath as mp


def coefficient_exponents(sigma):
    s=Fraction(sigma)
    if not Fraction(1,2)<=s<=1:raise ValueError('sigma out of range')
    return 3*(1-s)/(2-s),(7-5*s)/(2-s)


def integrated_first_term(sigma,height):
    """High-precision mpmath diagnostic; not an Arb interval."""
    alpha,p=coefficient_exponents(sigma)
    a=mp.mpf(alpha.numerator)/alpha.denominator
    q=mp.mpf(p.numerator)/p.denominator
    rate=2-a; L=mp.log(height)
    return 2*mp.gammainc(q+1,rate*L,mp.inf)/rate**(q+1)
