#!/usr/bin/env python3
"""Independent, prime-free, exact-rational zero-side audit for PR #71.

Claim (using three EXTERNALLY certified zeros, Platt--Trudgian RH up
through H=3e12, Bellotti--Wong 2025 explicit N(T) bound, and Suzuki's formula):
                  Psi(t) > 0 for 18 <= t <= 34.5.
Together with the PR #71 finite certificate for 0<t<=log(1e8),
this gives a continuous positive range 0<t<=34.5 IF the PR certificate
and cited external inputs are accepted.

The script evaluates cosine without floating point: exact-rational Machin
bounds for pi; fixed-point integer interval Horner for Taylor polynomials;
explicit Taylor remainder and exact rational interpolation and tail bounds.
No prime enumeration, mpmath, flint, or external code is used at runtime.

Zero ordinate provenance: the three intervals displayed in the FLINT
acb_dirichlet_zeta_zeros documentation, expanded to +/- 1e-26 each.
Zero count: Bellotti--Wong 2025 (Math. Comp., doi:10.1090/mcom/4133).
Finite RH verification: Platt--Trudgian 2021, H=3e12.
Suzuki 2023 formula (1.3), including off-critical complex ordinates.
"""
from fractions import Fraction as F
from math import factorial
import hashlib
import json
from pathlib import Path

S = 10**80
H = 3*10**12
START_T_MILLI = 18000
END_T_MILLI = 34500
STEP_MILLI = 1
POLY_ORDER = 40
REM_DENOM = 10**50
TAYLOR_REMAINDER_VALID = 4**(2*POLY_ORDER+2)*REM_DENOM < factorial(2*POLY_ORDER+2)

GAMMA_CENTERS = (
    "14.1347251417346937904572519836",
    "21.0220396387715549926284795939",
    "25.0108575801456887632137909926",
)
GAMMA_RAD = F(1,10**26)  # contains the printed FLINT +/- 3.76e-29 etc.


def f_floor(x):
    return x.numerator // x.denominator


def f_ceil(x):
    return -((-x.numerator)//x.denominator)


def i_fraction(lo,hi):
    return (f_floor(lo*S), f_ceil(hi*S))


def i_add(a,b):
    return a[0]+b[0],a[1]+b[1]


def i_mul(a,b):
    c=(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1])
    return min(c)//S,-((-max(c))//S)


def atan_series_interval(q,n_terms=85):
    partial=sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(n_terms)),F(0))
    first_omitted=F((-1)**n_terms,(2*n_terms+1)*q**(2*n_terms+1))
    return (partial + min(F(0),first_omitted),
            partial + max(F(0),first_omitted))


def pi_interval():
    a,b=atan_series_interval(5),atan_series_interval(239)
    lo=16*a[0]-4*b[1]
    hi=16*a[1]-4*b[0]
    if not (3<lo<hi<F(7,2)):
        raise ArithmeticError('failed to enclose pi')
    return i_fraction(lo,hi)


def coefficients():
    return [i_fraction(F((-1)**k,factorial(2*k)),
                       F((-1)**k,factorial(2*k))) for k in range(POLY_ORDER+1)]


def cosine_interval(phase,coeff):
    if not (-4*S<=phase[0] and phase[1]<=4*S):
        raise ArithmeticError('range reduction failed')
    y2=i_mul(phase,phase)
    out=coeff[-1]
    for a in reversed(coeff[:-1]):
        out=i_add(a,i_mul(y2,out))
    # |cos(y)-sum_{k<=40} (-1)^k y^{2k}/(2k)!| <= |y|^82/82!.
    # For |y|<=4, this is < 10^-50, checked without floating point.
    if not TAYLOR_REMAINDER_VALID:
        raise ArithmeticError('Taylor remainder bound unproved')
    r=(S+REM_DENOM-1)//REM_DENOM
    return max(-S,out[0]-r),min(S,out[1]+r)


def zero_intervals():
    result=[]
    for v in GAMMA_CENTERS:
        center=F(v)
        result.append(i_fraction(center-GAMMA_RAD,center+GAMMA_RAD))
    return result


def phase_interval(g,t_milli,k,pi_iv):
    # Since t_milli/1000>0, interval multiplication by positive rational
    low=g[0]*t_milli//1000
    high=-((-g[1]*t_milli)//1000)
    if k>=0:
        return low-2*k*pi_iv[1],high-2*k*pi_iv[0]
    return low-2*k*pi_iv[0],high-2*k*pi_iv[1]


def choose_turns(g,t_milli):
    # Integer-only nearest integer to gamma*t/(2*pi); midpoint pi from Machin.
    # This selection isn't a proof obligation: phase bounds are checked.
    numerator=(g[0]+g[1])*t_milli
    piv=pi_interval_cached
    denominator=2000*(piv[0]+piv[1])
    return (numerator+denominator//2)//denominator


def certify_grid():
    global pi_interval_cached
    pi_interval_cached=pi_interval()
    gs=zero_intervals()
    coeff=coefficients()
    min_lb=None
    min_at=None
    h=hashlib.sha256()
    for t_milli in range(START_T_MILLI,END_T_MILLI+1,STEP_MILLI):
        value_lb=0
        for g in gs:
            k=choose_turns(g,t_milli)
            phase=phase_interval(g,t_milli,k,pi_interval_cached)
            cos_lo,cos_hi=cosine_interval(phase,coeff)
            # since g>0 and 1-cos>=0, denominator upper bound suffices
            term_lo=(2*(S-cos_hi)*S*S)//(g[1]*g[1])
            value_lb+=term_lo
        if min_lb is None or value_lb<min_lb:
            min_lb=value_lb
            min_at=t_milli
        h.update(f"{t_milli},{value_lb}\n".encode('ascii'))
    return (F(min_lb,S),min_at,h.hexdigest(),len(range(START_T_MILLI,END_T_MILLI+1,STEP_MILLI)))


def external_tail_bound():
    # Elementary arithmetic validating the simplifications of the Bellotti--Wong N(T) bound.
    # For y>=H, log y < y/10^10 by monotonicity and log H<30<H/10^10.
    # log H<30 is witnessed by the finite Taylor lower bound to exp(30).
    exp30_lower=sum((F(30**j,factorial(j)) for j in range(90)),F(0))
    if not (exp30_lower>H and F(30)<F(H,10**10)
            and F(8084,1000)/H+F(35,100)/10**10 < F(1,7)):
        raise ArithmeticError('zero-count majorant simplification failed')
    # Consequently Bellotti--Wong -> N(y)<=y*log(y)/(2*pi) for every y>=H.
    # Sum of 1/height^2 above H (positive heights) <= (log H+1)/(pi*H)
    # <31/(3H). For off-critical zeros |Im gamma|<1/2.
    # e^(69/4) <= (1-(69/4)/2048)^-2048 <2^25:
    n=2048
    if not (4*n)**n < 2**25 * (4*n-69)**n:
        raise ArithmeticError('exponential majorant failed')
    return F(2*(1+2**25)*31,3*H)


def certify():
    grid,at,sha,n=certify_grid()
    D=2*(F(1,14)+F(1,21)+F(1,25))
    continuity_loss=D*F(STEP_MILLI,2000)
    spectral_tail=external_tail_bound()
    final_lb=grid-continuity_loss-spectral_tail
    if final_lb<=0:
        raise ArithmeticError(f"not certified: {float(final_lb)}")
    return {
        'scope':'Psi(t)>0 for 18<=t<=34.5; external theorems and certified gamma intervals imported',
        'primeness_or_prime_power_enumeration':False,
        'arithmetic':'exact rationals / integer interval, no floating point evaluations',
        'gamma_source':'FLINT acb_dirichlet_zeta_zeros 128 bit documentation: first three zero intervals',
        'RH_up_to_height':H,
        'zero_count_majorant':'Bellotti-Wong 2025; N(T)<=T log(T)/(2 pi) for T>=3e12',
        'points':n,'intervals':n-1,
        'grid_start_t':'18', 'grid_end_t':'34.5','grid_step':'0.001',
        'grid_minimum_lower':format(float(grid),'.15g'),
        'grid_minimum_at':str(F(at,1000)),
        'continuity_loss_upper':format(float(continuity_loss),'.15g'),
        'high_zero_tail_upper':format(float(spectral_tail),'.15g'),
        'uniform_lower_bound':format(float(final_lb),'.15g'),
        'strictly_positive':True,
        'grid_sha256':sha,
        'gamma_centers':GAMMA_CENTERS,
        'gamma_radius':'1e-26',
        'proof_numerators':{
            'grid_lower_numerator':str(grid.numerator),'grid_lower_denominator':str(grid.denominator),
            'uniform_lower_numerator':str(final_lb.numerator),'uniform_lower_denominator':str(final_lb.denominator),
        },
    }


if __name__=='__main__':
    record=certify()
    destination=Path(__file__).resolve().with_name('NUMERIA_P71_ZERO_SIDE_AUDIT.json')
    with open(destination,'w',encoding='utf-8') as f:
        json.dump(record,f,indent=2,ensure_ascii=False)
        f.write('\n')
    print(json.dumps({k:v for k,v in record.items() if k!='proof_numerators'},indent=2))
    print('audit written:',destination)
