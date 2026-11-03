"""NUMERIA–ES5: an unconditional arithmetic gateway to an explicit trace formula.

Two independent EXACT FIXED-POINT interval enclosures for -zeta'(sigma)/zeta(sigma)
for integer sigma>1: (1) prime-power von Mangoldt series, (2) ratio of two
ordinary Dirichlet series, including fully explicit tails. This is a check of
known analytic identities, not a RH result or a global Weil trace formula.

Also performs non-certified mpmath exploratory checks of a Gaussian-smoothed
arithmetic trace on Re(s)=sigma, and of the completed-zeta log derivative.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from functools import lru_cache

import mpmath as mp

SCALE = 10**40


def ceildiv(a: int, b: int) -> int:
    assert b > 0
    return -((-a) // b)


def spf_sieve(N: int) -> list[int]:
    assert N >= 3
    spf = list(range(N + 1))
    for p in range(2, math.isqrt(N) + 1):
        if spf[p] == p:
            for n in range(p * p, N + 1, p):
                if spf[n] == n:
                    spf[n] = p
    return spf


def factorization(n: int, spf: list[int]) -> dict[int, int]:
    f = {}
    while n > 1:
        p = spf[n]
        f[p] = f.get(p, 0) + 1
        n //= p
    return f


def prime_power_base(n: int, spf: list[int]) -> int | None:
    if n <= 1:
        return None
    p = spf[n]
    while n % p == 0:
        n //= p
    return p if n == 1 else None


def rational_log_bounds(x: Fraction, terms: int = 54) -> tuple[Fraction, Fraction]:
    """Rigorous log(x) for 1 <= x <= 2, via atanh and explicit remainder."""
    assert 1 <= x <= 2 and terms >= 1
    y = (x - 1) / (x + 1)
    val = Fraction(0)
    for k in range(terms):
        val += 2 * y**(2*k+1) / (2*k+1)
    rem = 2 * y**(2*terms+1) / ((2*terms+1)*(1-y*y))
    return val, val + rem


@lru_cache(maxsize=None)
def log_int_bounds(n: int) -> tuple[int, int]:
    """[integer floor, integer ceil] of SCALE * ln(n), fully rational."""
    if n == 1:
        return (0, 0)
    assert n > 1
    k = n.bit_length()-1
    lo2, hi2 = rational_log_bounds(Fraction(2))
    lo1, hi1 = rational_log_bounds(Fraction(n,2**k))
    lo = k*lo2 + lo1
    hi = k*hi2 + hi1
    return (lo.numerator*SCALE//lo.denominator, ceildiv(hi.numerator*SCALE,hi.denominator))


def fmt_scaled(n: int, digits: int = 20) -> str:
    with localcontext() as ctx:
        ctx.prec = digits + 30
        return format(Decimal(n)/Decimal(SCALE), f'.{digits}f')


def series_intervals(sigma: int=4, N: int=4096) -> dict:
    """Outward integer enclosures using exact rational operations.

    sigma integer >= 2, N >= 3. C_N=integral_N^inf log(x)x^-sigma dx
    upper-bounds BOTH the derivative-series and prime-series tails.
    """
    if sigma < 2 or N < 3:
        raise ValueError('sigma must be >=2 and N >=3')
    spf = spf_sieve(N)
    logs = {p:log_int_bounds(p) for p in range(2,N+1) if spf[p]==p}
    zlo=zhi=dlo=dhi=plo=phi=0
    count_pp=0
    for n in range(1,N+1):
        denominator = n**sigma
        zlo += SCALE//denominator
        zhi += ceildiv(SCALE,denominator)
        if n > 1:
            f = factorization(n,spf)
            n_lo=sum(exp*logs[p][0] for p,exp in f.items())
            n_hi=sum(exp*logs[p][1] for p,exp in f.items())
            dlo += n_lo//denominator
            dhi += ceildiv(n_hi,denominator)
            p = prime_power_base(n,spf)
            if p is not None:
                count_pp += 1
                plo += logs[p][0]//denominator
                phi += ceildiv(logs[p][1],denominator)
    k=sigma-1
    ln_hi=log_int_bounds(N)[1]
    tail_log = ceildiv(k*ln_hi + SCALE, k*k * N**k)
    tail_z = ceildiv(SCALE,k*N**k)
    zhi_full=zhi+tail_z
    dhi_full=dhi+tail_log
    phi_full=phi+tail_log
    ratio_lo=(dlo*SCALE)//zhi_full
    ratio_hi=ceildiv(dhi_full*SCALE,zlo)
    intersect=max(plo,ratio_lo)<=min(phi_full,ratio_hi)
    return {
        'sigma':sigma,'N':N,'base':'elementary known identity -zeta_prime/zeta=sum von Mangoldt(n)/n^sigma, sigma>1',
        'scale':str(SCALE),'precision':'exact integer/fraction outward bounds, not floating point',
        'prime_powers_count_through_N':count_pp,
        'prime_power_trace_interval':[fmt_scaled(plo),fmt_scaled(phi_full)],
        'independent_dirichlet_log_derivative_interval':[fmt_scaled(ratio_lo),fmt_scaled(ratio_hi)],
        'prime_tail_bound':fmt_scaled(tail_log),
        'zeta_tail_bound':fmt_scaled(tail_z),
        'overlap_certified':intersect,
        'overlap_interval':[fmt_scaled(max(plo,ratio_lo)),fmt_scaled(min(phi_full,ratio_hi))] if intersect else None,
        'primes_interval_scaled_exact':[str(plo),str(phi_full)],
        'dirichlet_interval_scaled_exact':[str(ratio_lo),str(ratio_hi)],
        'zeta_sum_interval_scaled_exact':[str(zlo),str(zhi_full)],
        'minus_zeta_prime_interval_scaled_exact':[str(dlo),str(dhi_full)],
        'analytic_bounds':'sum_{n>N} log(n)/n^sigma <= N^(1-sigma)[log(N)/(sigma-1)+1/(sigma-1)^2]; zeta tail <= N^(1-sigma)/(sigma-1)',
    }


def compact_support_check() -> dict:
    """9/10 = 2L; log(2)<9/10<log(3), rigorously with rational logs."""
    mid=9*SCALE//10
    lo2, hi2 = log_int_bounds(2)
    lo3, hi3 = log_int_bounds(3)
    assert hi2 < mid < lo3
    return {'L':'9/20','twice_L':'9/10',
            'log_2_interval':[fmt_scaled(lo2,25),fmt_scaled(hi2,25)],
            'log_3_interval':[fmt_scaled(lo3,25),fmt_scaled(hi3,25)],
            'certified':'log(2)<2L<log(3); only prime power n=2 falls inside (0,2L)',
            'prime_two_weight':'Lambda(2)/sqrt(2)=log(2)/sqrt(2); sign, pole and projection normalizations must still be checked against V24'}


def completed_log_derivative(s):
    return (1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2
            + mp.zeta(s, derivative=1)/mp.zeta(s))


def gaussian_trace_exploration(sigma:int=4,N:int=4096,dps:int=50) -> dict:
    """Noncertified numerical test; analytic prime truncation and t-tail bounds.

    sqrt(a/pi) int_R exp(-a*t²) exp(i*t*u) [-zeta'/zeta(sigma+it)] dt
       = sum_{n>=2} Lambda(n)/n^sigma exp[-(u-log n)^2/(4*a)].
    Absolute convergence for sigma>1. Here sigma>=4, a=1/4, T=18.
    """
    if sigma < 4 or N < 3:
        raise ValueError('Gaussian experiment restricted to sigma>=4, N>=3')
    with mp.workdps(dps):
        a=mp.mpf(1)/4
        T=mp.mpf(18)
        u=mp.log(2)
        def integrand(t):
            s=sigma+1j*t
            return mp.exp(-a*t*t)*mp.re(mp.exp(1j*t*u)*(-mp.zeta(s,derivative=1)/mp.zeta(s)))
        integ = 2*mp.sqrt(a/mp.pi)*mp.quad(integrand,[0,2,5,10,T])
        spf=spf_sieve(N)
        rhs=mp.fsum([mp.log(p)/n**sigma*mp.exp(-(u-mp.log(n))**2/(4*a))
                     for n in range(2,N+1) if (p:=prime_power_base(n,spf)) is not None])
        # For sigma>=4, -zeta'/zeta(sigma) <= sum_{n>=2} 1/n^3 <= 1/2.
        # Absolute omitted integral bound: sqrt(a/pi)*exp(-a*T^2)/(2*a*T).
        integration_t_tail = mp.sqrt(a/mp.pi)*mp.exp(-a*T*T)/(2*a*T)
        # arithmetic prime tail: Gaussian <=1, Lambda(n)<=log n.
        prime_tail = N**(1-sigma)*(mp.log(N)/(sigma-1)+1/(sigma-1)**2)
        # Check functional equation for completed xi via exact log derivative identity
        s=mp.mpc('2.25','0.75')
        functional_error=abs(completed_log_derivative(s)+completed_log_derivative(1-s))
        return {
            'sigma':sigma,'a':'1/4','u':'log(2)','integral_cutoff_T':'18',
            'LHS_gaussian_logderivative_integral_numeric':mp.nstr(integ,36),
            'RHS_prime_power_gaussian_partial_numeric':mp.nstr(rhs,36),
            'numerical_abs_difference':mp.nstr(abs(integ-rhs),36),
            'rigorous_analytic_bound_for_missing_prime_tail':mp.nstr(prime_tail,20),
            'rigorous_analytic_bound_for_missing_t_tail':mp.nstr(integration_t_tail,20),
            'note':'mpmath quadrature and prime terms are NOT interval certified. Bounds control only truncation tails, not quadrature and roundoff.',
            'completed_xi_logderivative_symmetry_residual_numeric':mp.nstr(functional_error,12),
            'completed_xi_logderivative_identity':'xi_prime/xi(s)=1/s+1/(s-1)-log(pi)/2+digamma(s/2)/2+zeta_prime/zeta(s); xi(s)=xi(1-s)',
            'spectral_limitation':'This arithmetic trace does not identify the ES4 spectrum, prove the global Weil trace, or say anything about zeros off the critical line.'
        }


def run(N:int=4096,sigma:int=4,dps:int=50):
    if sigma<4 or N<3 or dps<30: raise ValueError('N>=3, sigma>=4, dps>=30')
    derived=series_intervals(sigma=sigma,N=N)
    assert derived['overlap_certified']
    smoothed=gaussian_trace_exploration(sigma=sigma,N=N,dps=dps)
    return {'schema':'NUMERIA.ES5.ARITHMETIC_LOGDERIVATIVE_TRACE.1',
            'scientific_status':'known identities and finite interval enclosures in Re(s)>1; numerical Gaussian diagnostic; NOT global Weil formula or RH proof',
            'environment':{'python':platform.python_version(),'mpmath':mp.__version__},
            'finite_interval_certificate':derived,
            'v24_local_prime_support_certificate':compact_support_check(),
            'gaussian_trace_numeric':smoothed,
            'next_gates':['Independently certify gaussian integrals with Arb/FLINT or explicit quadrature remainder',
                          'Prove test-function normalization and pole/archimedean terms match E7/V24',
                          'Establish global Weil explicit trace for admissible tests without RH',
                          'Only then investigate conditional relation to ES4 spectral operator']}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--N',type=int,default=4096)
    p.add_argument('--sigma',type=int,default=4)
    p.add_argument('--dps',type=int,default=50)
    p.add_argument('--out',default='RESULTADOS_ES5.json')
    args=p.parse_args()
    result=run(args.N,args.sigma,args.dps)
    path=Path(args.out)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8')
    print('ES5 | finite interval cert overlap:',result['finite_interval_certificate']['overlap_certified'],
          '| gaussian abs error:',result['gaussian_trace_numeric']['numerical_abs_difference'],
          '| NO_RH_PROOF')

if __name__=='__main__':
    main()
