"""NUMERIA–ES6: unconditional Gaussian Guinand–Weil formula, numerical diagnostic.

Two distinct claims:
(1) Exact rational fixed-point enclosures for the GAUSSIAN PRIME-POWER TERM
    at a=1/50, N=64, with an explicit upper bound on its missing tail.
(2) High-precision, NON-CERTIFIED numeric comparison of the full explicit formula
    to a partial list of known critical-line zeta zeros, including d/da.

The full zero-side identity sums ALL nontrivial zeros (including possible
zeros off the critical line); the mpmath zetazero list is NOT a certification
of the complete zero sum. Neither ES4 spectral identification nor RH is proved.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import platform
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import mpmath as mp

SCALE = 10**40

def ceildiv(x:int,y:int)->int:
    return -((-x)//y)

def arctan_alt(x: Fraction, m:int=44)->tuple[Fraction,Fraction]:
    """alternating Taylor, x in [0,1], terms from k=0..m inclusive."""
    assert 0<x<=1
    s=Fraction(0)
    for k in range(m):
        s+=(-1)**k*x**(2*k+1)/(2*k+1)
    term=x**(2*m+1)/(2*m+1)
    return (s,s+term) if m%2 else (s-term,s)

def pi_bounds()->tuple[int,int]:
    # pi = 16 arctan(1/5) - 4 arctan(1/239)
    lo5,hi5=arctan_alt(Fraction(1,5))
    lo239,hi239=arctan_alt(Fraction(1,239))
    lo=16*lo5-4*hi239
    hi=16*hi5-4*lo239
    return lo.numerator*SCALE//lo.denominator,ceildiv(hi.numerator*SCALE,hi.denominator)

def log_ratio_bounds(x: Fraction,terms:int=55)->tuple[Fraction,Fraction]:
    assert 1<=x<=2
    y=(x-1)/(x+1)
    s=Fraction(0)
    for k in range(terms): s+=2*y**(2*k+1)/(2*k+1)
    r=2*y**(2*terms+1)/((2*terms+1)*(1-y*y))
    return s,s+r

@lru_cache(maxsize=None)
def log_int_bounds(n:int)->tuple[int,int]:
    assert n>=1
    if n==1:return 0,0
    k=n.bit_length()-1
    l2,h2=log_ratio_bounds(Fraction(2))
    l1,h1=log_ratio_bounds(Fraction(n,2**k))
    lo=k*l2+l1;hi=k*h2+h1
    return lo.numerator*SCALE//lo.denominator,ceildiv(hi.numerator*SCALE,hi.denominator)

def unit_exp_bounds(zlo:Fraction,zhi:Fraction,m:int=20)->tuple[int,int]:
    """Alternating Taylor for exp(-z), 0<=zlo<=zhi<=1.
    Include both partial sums 2m-1 and 2m (lower and upper)."""
    assert 0<=zlo<=zhi<=1 and m>=1
    def polynomial(z:Fraction,K:int)->Fraction:
        s=Fraction(1);term=Fraction(1)
        for k in range(1,K+1):
            term=term*(-z)/k
            s+=term
        return s
    low=polynomial(zhi,2*m-1)
    high=polynomial(zlo,2*m)
    return max(0,low.numerator*SCALE//low.denominator),min(SCALE,ceildiv(high.numerator*SCALE,high.denominator))

def exp_neg_bounds(xlo:int,xhi:int,m:int=20)->tuple[int,int]:
    """Outward fixed-point bounds for exp(-x) with x measured in SCALE units."""
    assert 0<=xlo<=xhi
    if xhi==0:return SCALE,SCALE
    q=1
    while xhi>q*SCALE:q*=2
    lo,hi=unit_exp_bounds(Fraction(xlo,q*SCALE),Fraction(xhi,q*SCALE),m)
    while q>1:
        lo=lo*lo//SCALE
        hi=ceildiv(hi*hi,SCALE)
        q//=2
    return lo,hi

def positive_product(a:tuple[int,int],b:tuple[int,int])->tuple[int,int]:
    assert 0<=a[0]<=a[1] and 0<=b[0]<=b[1]
    return (a[0]*b[0]//SCALE,ceildiv(a[1]*b[1],SCALE))

def sqrt_scaled(a:tuple[int,int])->tuple[int,int]:
    assert 0<a[0]<=a[1]
    l=math.isqrt(a[0]*SCALE)
    h=math.isqrt(a[1]*SCALE)
    if h*h<a[1]*SCALE:h+=1
    return l,h

def inverse_positive(a:tuple[int,int])->tuple[int,int]:
    assert 0<a[0]<=a[1]
    return SCALE*SCALE//a[1],ceildiv(SCALE*SCALE,a[0])

def spf_sieve(n:int):
    x=list(range(n+1))
    for p in range(2,math.isqrt(n)+1):
        if x[p]==p:
            for k in range(p*p,n+1,p):
                if x[k]==k:x[k]=p
    return x

def power_base(n:int,spf:list[int]):
    p=spf[n]
    while n%p==0:n//=p
    return p if n==1 else None

def fmt(lo:int,hi:int,digits:int=30)->list[str]:
    from decimal import Decimal,localcontext
    with localcontext() as d:
        d.prec=digits+35
        return [format(Decimal(i)/Decimal(SCALE),f'.{digits}f') for i in (lo,hi)]

def prime_gaussian_certificate(N:int=64)->dict:
    """Enclose P_a(N)=(pi*a)^-1/2 sum_{n<=N} Lambda(n)/sqrt(n)*exp(-(ln n)^2/(4a)).

    a=1/50, N>=64. For n>N, ln n>ln N>=4 and monotonic f(x)=ln x/sqrt x*exp(-ln(x)^2/(4a)):
    tail <= 4 sqrt(a/pi) exp(-(ln N)^2/(8a))
         <= e^-100/3 < 1/SCALE (for N>=64).
    All bounds follow from rational Taylor series and integer directed rounding.
    """
    if N<64:raise ValueError('certified global prime tail requires N>=64')
    a=Fraction(1,50)
    pilo,pihi=pi_bounds()
    assert pilo>3*SCALE and pihi<4*SCALE
    pi_a=(pilo//50,ceildiv(pihi,50))
    pref=inverse_positive(sqrt_scaled(pi_a))
    spf=spf_sieve(N)
    total=(0,0)
    details=[]
    for n in range(2,N+1):
        p=power_base(n,spf)
        if p is None:continue
        ll,lu=log_int_bounds(n)
        # x=(ln n)^2/(4a)=(25/2)*(ln n)^2
        xlo=25*ll*ll//(2*SCALE)
        xhi=ceildiv(25*lu*lu,2*SCALE)
        decay=exp_neg_bounds(xlo,xhi)
        logp=log_int_bounds(p)
        w=positive_product(logp,inverse_positive(sqrt_scaled((n*SCALE,n*SCALE))))
        term=positive_product(pref,positive_product(w,decay))
        total=(total[0]+term[0],total[1]+term[1])
        if n<=16:
            details.append({'n':n,'p':p,'P_term_interval':fmt(*term, digits=26)})
    # Rigorous tail at N=64 and beyond: logN>=4, sqrt(a/pi)<1/12,
    # 4 sqrt(a/pi) exp(-100) < exp(-100)/3 < 1/SCALE.
    assert log_int_bounds(N)[0]>=4*SCALE
    assert exp_neg_bounds(100*SCALE,100*SCALE)[1]<=1
    tail_upper_scaled=1
    total_all=(total[0],total[1]+tail_upper_scaled)
    return {'a':'1/50','N':N,'scale':str(SCALE),'pi_interval_scaled':list(map(str,(pilo,pihi))),
            'positive_prime_term_finite_interval':fmt(*total,32),
            'positive_prime_term_all_n_interval':fmt(*total_all,32),
            'prime_term_all_n_exact_scaled':list(map(str,total_all)),
            'prime_tail_bound':'0 <= P(a)-P_N(a) < 10^-40 (rigorous for N>=64)',
            'bound_proof':'Lambda(n)<=ln(n); decreasing integrand for n>=64; integral after y=ln x <= 4*sqrt(a/pi)*exp(-(ln N)^2/(8a)); ln 64>4, sqrt(a/pi)<1/12, exp(-100)<10^-40',
            'first_prime_powers':details,
            'certificate_type':'outward integer/rational arithmetic; analytic prime-tail enclosure ONLY (not archimedean integral or zero sum)'}

def prime_powers_mpmath(N:int,a:mp.mpf)->tuple[mp.mpf,mp.mpf]:
    spf=spf_sieve(N)
    values=[];derivs=[]
    for n in range(2,N+1):
        p=power_base(n,spf)
        if p is None:continue
        u=mp.log(n)
        val=(mp.log(p)/mp.sqrt(n))/mp.sqrt(mp.pi*a)*mp.exp(-u*u/(4*a))
        values.append(val)
        derivs.append(val*(u*u/(4*a*a)-1/(2*a)))
    # values = +P(a); actual prime contribution is -P(a).
    # derivative prime contribution = -P'(a)
    return mp.fsum(values),-mp.fsum(derivs)

def analytic_parts(a:mp.mpf,N:int=128,T:int=128,derivative:bool=False):
    pole=2*mp.exp(a/4)
    dpole=mp.exp(a/4)/2
    logpi=mp.log(mp.pi)
    def density(t):return (mp.re(mp.digamma(mp.mpf(1)/4+1j*t/2))-logpi)/mp.pi
    cuts=[mp.mpf(0),mp.mpf(1),mp.mpf(2),mp.mpf(4),mp.mpf(8),mp.mpf(16),mp.mpf(32),mp.mpf(64),mp.mpf(T)]
    cuts=sorted(set(x for x in cuts if x<=T))
    if cuts[-1]<T:cuts.append(mp.mpf(T))
    arch=mp.quad(lambda t:mp.exp(-a*t*t)*density(t),cuts)
    P,dP=prime_powers_mpmath(N,a)
    ans={'pole':pole,'arch':arch,'prime':-P,'total':pole+arch-P}
    if derivative:
        darch=mp.quad(lambda t:-t*t*mp.exp(-a*t*t)*density(t),cuts)
        ans.update({'d_pole':dpole,'d_arch':darch,'d_prime':dP,'d_total':dpole+darch+dP})
    return ans

def run_exploration(a_values=('0.02','0.04'),J:int=24,N:int=128,dps:int=55,T:int=128):
    if J<5 or N<64 or dps<35 or T<64:raise ValueError('J>=5,N>=64,dps>=35,T>=64')
    with mp.workdps(dps):
        gammas=[mp.im(mp.zetazero(k)) for k in range(1,J+1)]
        cases=[]
        for raw in a_values:
            a=mp.mpf(raw)
            if not 0<a<=mp.mpf('0.1'):raise ValueError('0<a<=0.1')
            part=analytic_parts(a,N,T,derivative=True)
            z=2*mp.fsum(mp.exp(-a*g*g) for g in gammas)
            dz=-2*mp.fsum(g*g*mp.exp(-a*g*g) for g in gammas)
            # Adversarial spectral model ES4: the positions cannot be replaced by matching the smooth count.
            def N0(y):return y/(2*mp.pi)*(mp.log(y/(2*mp.pi))-1)+mp.mpf(7)/8
            model=[]
            for k in range(1,J+1):
                q=mp.mpf(k)-mp.mpf(7)/8
                model.append(2*mp.pi*q/mp.lambertw(q/mp.e))
            model_sum=2*mp.fsum(mp.exp(-a*t*t) for t in model)
            cases.append({'a':raw,'poles':mp.nstr(part['pole'],45),
                          'archimedean_numeric':mp.nstr(part['arch'],45),
                          'prime_numeric':mp.nstr(part['prime'],45),
                          'rhs_numeric':mp.nstr(part['total'],45),
                          'known_critical_zero_sum_numeric':mp.nstr(z,45),
                          'rhs_minus_known_zero_sum_numeric':mp.nstr(part['total']-z,20),
                          'first_J_ES4_model_sum_numeric':mp.nstr(model_sum,30),
                          'ES4_minus_zeta_known_zero_sum_numeric':mp.nstr(model_sum-z,20),
                          'derivative_rhs_numeric':mp.nstr(part['d_total'],40),
                          'derivative_known_zero_sum_numeric':mp.nstr(dz,40),
                          'derivative_discrepancy_numeric':mp.nstr(part['d_total']-dz,20),
                          'zero_side':'only first J known on-critical zeros paired; full formula requires ALL nontrivial zeros without RH assumption',
                          'numeric_status':'mpmath NON-CERTIFIED; arch quadrature and possible off-critical zeros not certified'})
    return {'schema':'NUMERIA.ES6.GAUSSIAN_EXPLICIT_FORMULA.1',
            'scientific_status':'unconditional classical Guinand-Weil identity; finite certified prime term; exploratory full numeric crosscheck; NOT RH proof',
            'source':'Chirre-Goncalves 2019 Guinand-Weil formula and elementary Euler product; full formula unconditional',
            'convention':{'h':'h_a(t)=exp(-a*t^2)',
                          'g':'g_a(u)=(1/(2*pi))*sqrt(pi/a)*exp(-u^2/(4*a))',
                          'zero_variable':'z_rho=(rho-1/2)/i may be complex; no RH assumed',
                          'formula':'sum_rho h(z_rho)=2*exp(a/4)+(1/(2*pi))*int_R h(t)*(Re psi(1/4+it/2)-log(pi))dt - (1/sqrt(pi*a))*sum_{n>=2}Lambda(n)/sqrt(n)*exp(-(ln n)^2/(4*a))',
                          'differential':'derivative in a: d/da sum h(z_rho) = - sum z_rho^2 h(z_rho); derivative of RHS termwise, with justified Gaussian decay'},
            'finite_prime_interval_certificate':prime_gaussian_certificate(64),
            'exploration':{'digits':dps,'zero_count_J':J,'arithmetic_cutoff_N':N,'archimedean_cutoff_T':T,
                           'cases':cases,'limitations':['zero side lists known critical-line zeros only; not all zeros certified',
                          'archimedean quadrature not interval-certified',
                          'ES4 operator has no arithmetic explicit trace identity',
                          'Gaussian test is Schwartz, not compactly supported: cannot be substituted for E7/V24 finite matrices without new bounds',
                          'finite numeric discrepancy is not a certificate of RH or global Weil positivity']}}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--dps',type=int,default=55)
    p.add_argument('--zeros',type=int,default=24)
    p.add_argument('--N',type=int,default=128)
    p.add_argument('--T',type=int,default=128)
    p.add_argument('--out',default='RESULTADOS_ES6.json')
    args=p.parse_args()
    data=run_exploration(J=args.zeros,N=args.N,dps=args.dps,T=args.T)
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(data,sort_keys=True,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('ES6',out,'| certified prime term, numeric explicit formula and parameter differential; NOT RH PROOF')

if __name__=='__main__':main()
