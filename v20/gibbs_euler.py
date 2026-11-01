#!/usr/bin/env python3
"""Exact rational Gibbs/Euler and von Mangoldt enclosure for integer beta>2.

Analytic prerequisites (recorded, not formally proved here): integral test,
Euler product for real beta>1, log x <= x for x>0, and product geometric series.
"""
from __future__ import annotations
from fractions import Fraction as Q
from math import isqrt


def primes_upto(n: int) -> list[int]:
    s=[True]*(n+1)
    if n>=0:s[0]=False
    if n>=1:s[1]=False
    for i in range(2,isqrt(n)+1):
        if s[i]:
            for j in range(i*i,n+1,i):s[j]=False
    return [i for i,v in enumerate(s) if v]


def log_interval(x: Q|int, terms: int=20) -> tuple[Q,Q]:
    """Rational directed log via atanh; 1<=y<2 and 0<=q<1/3.

    (2q**(2m+1))/((2m+1)*(1-q*q)) bounds residual from above.
    """
    x=Q(x)
    if x<=0 or terms<2:raise ValueError('log needs positive input, >=2 terms')
    k=0;y=x
    while y>=2:y/=2;k+=1
    while y<1:y*=2;k-=1
    def atanh_enclosure(z):
        s=sum((2*z**(2*j+1)/Q(2*j+1) for j in range(terms)),Q(0))
        remainder=2*z**(2*terms+1)/(Q(2*terms+1)*(1-z*z))
        return s,s+remainder
    lo2,hi2=atanh_enclosure(Q(1,3))
    lo,hi=atanh_enclosure((y-1)/(y+1))
    if k>=0:return lo+k*lo2,hi+k*hi2
    return lo+k*hi2,hi+k*lo2


def dirichlet_zeta(beta:int, N:int) -> tuple[Q,Q]:
    if beta<=1 or N<2:raise ValueError('domain: integer beta>1, N>=2')
    partial=sum((Q(1,n**beta) for n in range(1,N+1)),Q(0))
    return partial+Q(1,(beta-1)*(N+1)**(beta-1)), partial+Q(1,(beta-1)*N**(beta-1))


def euler_zeta(beta:int,P:int,N:int)->tuple[Q,Q]:
    if beta<=1 or P<2:raise ValueError('Euler domain')
    z=Q(1)
    for p in primes_upto(P):z*=Q(p**beta,p**beta-1)
    # Union bound: discarded n divisible by some p>P contributes at most
    # zeta(beta) * sum_{p>P} p^-beta <= U_zeta * P^(1-beta)/(beta-1).
    _,z_upper=dirichlet_zeta(beta,N)
    return z,z+z_upper*Q(1,(beta-1)*P**(beta-1))


def mangoldt_prime_power(n:int)->int|None:
    """p if n is a positive power of prime p, else None."""
    if n<2:return None
    for p in primes_upto(isqrt(n)):
        if n%p==0:
            a=n
            while a%p==0:a//=p
            return p if a==1 else None
    return n


def weighted_interval(terms:list[tuple[Q,tuple[Q,Q]]]):
    low=sum((factor*log_bounds[0] for factor,log_bounds in terms),Q(0))
    high=sum((factor*log_bounds[1] for factor,log_bounds in terms),Q(0))
    return low,high


def mangoldt_dirichlet(beta:int,N:int)->tuple[Q,Q]:
    if beta<=2 or N<2:raise ValueError('integral tail bound needs beta>2')
    coefficients={}
    for n in range(2,N+1):
        p=mangoldt_prime_power(n)
        if p is not None:coefficients[p]=coefficients.get(p,Q(0))+Q(1,n**beta)
    low,high=weighted_interval([(factor,log_interval(p)) for p,factor in coefficients.items()])
    # Lambda(n)<=log n<=n; integral comparison bounds remainder.
    return low,high+Q(1,(beta-2)*N**(beta-2))


def mangoldt_euler(beta:int,P:int)->tuple[Q,Q]:
    if beta<=2 or P<2:raise ValueError('integral tail bound needs beta>2')
    factors=[(Q(1,p**beta-1),log_interval(p)) for p in primes_upto(P)]
    low,high=weighted_interval(factors)
    # For p>P, 1/(p^beta-1) <= 1/[p^beta(1-2^-beta)]
    # and log p<=p: tail <= [(1-2^-beta)*(beta-2)*P^(beta-2)]^-1.
    return low,high+Q(2**beta,2**beta-1)*Q(1,(beta-2)*P**(beta-2))


def overlap(x:tuple[Q,Q],y:tuple[Q,Q])->bool:
    return max(x[0],y[0])<=min(x[1],y[1])


def exact_decimal_bound(x:Q,upward:bool,digits:int=28)->str:
    q=x.numerator*10**digits//x.denominator
    if upward and Q(q,10**digits)<x:q+=1
    sign='-' if q<0 else ''
    q=abs(q)
    return f'{sign}{q//10**digits}.{q%10**digits:0{digits}d}'


def serialize_interval(pair):return [exact_decimal_bound(pair[0],False),exact_decimal_bound(pair[1],True)]


def run():
    results=[]
    for beta in (2,3,4):
        d=dirichlet_zeta(beta,256)
        e=euler_zeta(beta,31,256)
        assert overlap(d,e)
        assert e[0]>0 and d[0]>0
        item={'beta':beta,'dirichlet_interval':serialize_interval(d),'euler_interval':serialize_interval(e),'intersection':True}
        if beta>2:
            a=mangoldt_dirichlet(beta,256)
            b=mangoldt_euler(beta,31)
            assert overlap(a,b)
            assert a[0]>0 and b[0]>0
            item.update(mangoldt_dirichlet_interval=serialize_interval(a),mangoldt_euler_interval=serialize_interval(b),mangoldt_intersection=True)
        results.append(item)
    # Exact algebraic identity of the p-adic coefficient and the k-power truncation.
    for p in (2,3,5,7,11,13,17,19,23,29,31):
        for beta in (2,3,4):
            for K in (1,3,6):
                q=Q(1,p**beta)
                finite=sum((q**k for k in range(1,K+1)),Q(0))
                full=q/(1-q)
                assert full-finite==q**(K+1)/(1-q)
    return {'status':'GIBBS_EULER_EXACT_RATIONAL_ENCLOSURES_PASS','beta_cases':results,
            'finite_prime_power_geometric_identities':11*3*3,
            'domain':'real beta>1 for zeta; beta>2 for simple log-derivative tail',
            'theorems_assumed':['Euler product in Re(s)>1','integral comparison','log x <= x'],
            'note':'Continuation to Re(s)=1/2 does not preserve a probability distribution'}
