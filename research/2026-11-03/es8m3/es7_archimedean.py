"""NUMERIA ES7: validated enclosure of the Gaussian archimedean integral.

Every reported enclosing endpoint derives from directed rational arithmetic,
proved digamma-series tail bounds, a composite midpoint error bound, and
an explicit Gaussian tail. This does NOT enclose the full Weil zero side,
and it does NOT identify any V24 matrix entry.

A(a) = pi^-1 int_0^infty exp(-a*t*t) [Re psi(1/4+it/2)-log(pi)] dt.
"""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
import argparse, hashlib, json, math, platform
from pathlib import Path

S=10**28

def up(n:int,d:int)->int:
    assert d>0
    return -((-n)//d)

def enc(x:F)->tuple[int,int]:
    return (x.numerator*S//x.denominator,up(x.numerator*S,x.denominator))

def add(x,y):return x[0]+y[0],x[1]+y[1]
def sub(x,y):return x[0]-y[1],x[1]-y[0]
def mul(x,y):
    vals=[a*b for a in x for b in y]
    return min(vals)//S,up(max(vals),S)

def interval_fmt(x, digits=12):
    with localcontext() as c:
        c.prec=digits+40
        step=Decimal(1).scaleb(-digits)
        lo=(Decimal(x[0])/Decimal(S)).quantize(step,rounding=ROUND_FLOOR)
        hi=(Decimal(x[1])/Decimal(S)).quantize(step,rounding=ROUND_CEILING)
        return [format(lo,f'.{digits}f'),format(hi,f'.{digits}f')]

def atanh_log_ratio(x:F, terms:int=46):
    """Log(x) for 1<=x<=2 using a positive series with bounded tail."""
    assert 1<=x<=2
    y=(x-1)/(x+1)
    val=F(0);power=y
    for j in range(terms):
        val+=2*power/(2*j+1)
        power*=y*y
    rem=2*power/((2*terms+1)*(1-y*y))
    return val,val+rem

def pi_interval():
    def atan_series(x,terms=34):
        s=F(0)
        for k in range(terms):s+=(-1)**k*x**(2*k+1)/(2*k+1)
        rem=x**(2*terms+1)/(2*terms+1)
        return (s-rem,s) if terms%2 else (s,s+rem)
    a=atan_series(F(1,5));b=atan_series(F(1,239))
    return enc(F(16)*a[0]-F(4)*b[1])[0], enc(F(16)*a[1]-F(4)*b[0])[1]

def ln2_interval():
    l,u=atanh_log_ratio(F(2))
    return enc(l)[0],enc(u)[1]

def lnpi_interval(pi):
    """pi=2*(pi/2), and pi/2 is enclosed by exact rational endpoints."""
    a=atanh_log_ratio(F(pi[0],2*S));b=atanh_log_ratio(F(pi[1],2*S))
    return add(ln2_interval(),(enc(a[0])[0],enc(b[1])[1]))


def gamma_interval(M:int=4096):
    assert M>=2 and M&(M-1)==0
    H=(sum(S//k for k in range(1,M+1)),sum(up(S,k) for k in range(1,M+1)))
    exponent=M.bit_length()-1
    logM=(ln2_interval()[0]*exponent,ln2_interval()[1]*exponent)
    raw=sub(H,logM)
    # gamma = raw - theta, with 1/(2(M+1)) < theta < 1/(2M)
    return raw[0]-up(S,2*M), raw[1]-(S//(2*(M+1)))

def omega0_interval():
    pi=pi_interval(); gamma=gamma_interval()
    log2=ln2_interval();logpi=lnpi_interval(pi)
    # psi(1/4) = -gamma-pi/2-3*log 2
    b=(pi[0]//2,up(pi[1],2))
    return sub(sub(sub((0,0),gamma),b),add((3*log2[0],3*log2[1]),logpi))


def log1p_interval(u:F, terms=24):
    """log(1+u), 0<=u<=1: atanh with a positive geometric remainder."""
    assert 0<=u<=1
    v=u/(2+u)
    s=F(0);power=v
    for k in range(terms):
        s+=2*power/(2*k+1)
        power*=v*v
    r=2*power/((2*terms+1)*(1-v*v))
    return (enc(s)[0],enc(s+r)[1])


def omega_interval(t:F, terms:int=192, baseline=None):
    assert t>=0 and terms>=2
    base=baseline if baseline is not None else omega0_interval()
    tn=t.numerator;td=t.denominator
    num=16*tn*tn
    partial_lo=partial_hi=0
    for n in range(terms):
        u=4*n+1
        den=u*(u*u*td*td + 4*tn*tn)
        q=num*S//den
        partial_lo+=q
        partial_hi+=up(num*S,den)
    # b_n(t) positive, decreasing in n; integral test for n>=terms.
    # int_M^inf b_x dx = (1/2)*log(1+4t²/(4M+1)^2).
    # int_M^inf <= sum_{n>=M} b_n <= int_(M-1)^inf.
    tail_lo=log1p_interval(F(4)*t*t/(4*terms+1)**2)[0]//2
    tail_hi=up(log1p_interval(F(4)*t*t/(4*terms-3)**2)[1],2)
    return base[0]+partial_lo+tail_lo,base[1]+partial_hi+tail_hi


def exp_negative_interval(x:F, terms:int=17):
    """Outward-rational exp(-x) from alternating Taylor on [0,1], squaring."""
    assert x>=0
    if x==0:return S,S
    q=1
    while x>q:q*=2
    y=x/q
    p=F(1);term=F(1)
    lower=upper=None
    for k in range(1,2*terms+1):
        term*=(-y)/k
        p+=term
        if k==2*terms-1:lower=p
        elif k==2*terms:upper=p
    val=(max(0,enc(lower)[0]),min(S,enc(upper)[1]))
    while q>1:
        val=mul(val,val)
        q//=2
    return val

def gaussian_tail(a:F,T:int,pi_lower:int):
    """|A_tail| <= exp(-aT²)/pi * [7/(2aT)+18(T/(2a)+1/(4a²T))]."""
    exp_hi=exp_negative_interval(a*T*T)[1]
    factor=(F(7)/(2*a*T)+18*(F(T)/(2*a)+1/(4*a*a*T)))
    # exp_hi/S * factor * (S/pi_lower), returned in fixed units
    return up(exp_hi*factor.numerator*S, factor.denominator*pi_lower)


def arch_certificate(a:F=F(1,50),M:int=8192,terms:int=192,T:int=48, derivative:bool=True):
    if a<=0 or a>F(1,20) or M<64 or M>131072 or T<32 or terms<16:raise ValueError('unsafe numerical or analytic configuration')
    pi=pi_interval();base=omega0_interval();h=F(T,M)
    if not (-7*S<base[0]<base[1]<0):raise AssertionError('Omega(0) envelope exceeds |Omega0|<7')
    low_sum=hi_sum=0
    dlow_sum=dhi_sum=0
    Csum=Dsum=0
    prev_left_exp=exp_negative_interval(F(0))
    for j in range(M):
        l=h*j;r=h*(j+1);m=(l+r)/2
        Om=omega_interval(m,terms,base)
        gm=exp_negative_interval(a*m*m)
        fm=mul(Om,gm)
        low_sum+=fm[0];hi_sum+=fm[1]
        if derivative:
            dm=mul(fm,enc(-m*m))
            dlow_sum+=dm[0];dhi_sum+=dm[1]
        # Certified second-derivative bound for f=e^(-a t²) Omega(t).
        # |Omega|<=7+18t², |Omega'|<=36t, |Omega''|<=36.
        # |f''(t)| <= exp(-a*t²)*(36+180*a*t²+72*a²*t⁴+28*a²*t²+14*a).
        y=a*r*r
        C=36+180*y+72*y*y+28*a*y+14*a
        left_hi=prev_left_exp[1]
        C_hi=up(left_hi*C.numerator,C.denominator)
        Csum+=C_hi
        if derivative:
            # |(t²f)''| <= 2|f|+4t|f'|+t²|f''|, with
            # |f'|<=exp(-a t²)*(36t+2at(7+18t²)).
            D=(2*(7+18*r*r)+4*r*(36*r+2*a*r*(7+18*r*r))+r*r*C)
            Dsum+=up(left_hi*D.numerator,D.denominator)
        prev_left_exp=exp_negative_interval(a*r*r)
    # Signed ratio: (h/pi)*sum(f(mid)).
    def scaled_signed_sum(xlo,xhi):
        vv=[(xlo*T*S)//(M*pi[1]),up(xhi*T*S,M*pi[0])]
        # For negative x: lower divides by smaller pi, upper by larger pi.
        if xlo<0: vv[0]=xlo*T*S//(M*pi[0])
        if xhi<0: vv[1]=up(xhi*T*S,M*pi[1])
        return vv
    approx=scaled_signed_sum(low_sum,hi_sum)
    miderr=up(Csum*T*T*T*S,24*M*M*M*pi[0])
    tail=gaussian_tail(a,T,pi[0])
    result=(approx[0]-miderr-tail,approx[1]+miderr+tail)
    out={'a':str(a),'M':M,'T':T,'digamma_terms':terms,'scale':str(S),
         'omega0_exact_scaled':list(map(str,base)),
         'interval_scaled_exact':list(map(str,result)),
         'integral_interval_decimal':interval_fmt(result,12),
         'midpoint_interval_decimal':interval_fmt(approx,12),
         'midpoint_error_bound_decimal':interval_fmt((0,miderr),12)[1],
         'gaussian_tail_bound_decimal':interval_fmt((0,tail),16)[1],
         'status':'CERTIFIED rational enclosure of the archimedean Gaussian term ONLY',
         'analytic_preconditions':['DLMF digamma series and special value psi(1/4)',
             'harmonic-number bound on Euler gamma',
             'integral test for positive decreasing digamma summands',
             'composite midpoint theorem using explicit bound on |f_double_prime|',
             'Gaussian tail bounded by |Omega(t)| <= 7+18*t*t']}
    if derivative:
        approx_d=scaled_signed_sum(dlow_sum,dhi_sum)
        miderr_d=up(Dsum*T*T*T*S,24*M*M*M*pi[0])
        # |dA_tail| <= exp(-aT²)/pi*[7*(T/(2a)+1/(4a²T))+18*(T³/(2a)+3T/(4a²)+3/(8a³T))]
        exp_hi=exp_negative_interval(a*T*T)[1]
        factor=7*(F(T)/(2*a)+1/(4*a*a*T))+18*(F(T**3)/(2*a)+3*F(T)/(4*a*a)+3/(8*a*a*a*T))
        tail_d=up(exp_hi*factor.numerator*S,factor.denominator*pi[0])
        dr=(approx_d[0]-miderr_d-tail_d,approx_d[1]+miderr_d+tail_d)
        out.update({'derivative_interval_scaled_exact':list(map(str,dr)),
          'derivative_interval_decimal':interval_fmt(dr,12),
          'derivative_midpoint_error_bound_decimal':interval_fmt((0,miderr_d),12)[1],
          'derivative_tail_bound_decimal':interval_fmt((0,tail_d),12)[1],
          'derivative_status':'CERTIFIED dA/da for the Gaussian archimedean term ONLY'})
    return out


def run(M=8192,terms=192,T=48):
    values=[arch_certificate(F(1,50),M,terms,T),arch_certificate(F(1,25),M,terms,T)]
    return {'schema':'NUMERIA.ES7.ARCHIMEDEAN_GAUSSIAN_RATIONAL.1',
      'scope':'Weil Gaussian archimedean component and its a derivative; not complete zero-side certification',
      'environment':{'python':platform.python_version()},
      'cases':values,
      'limitations':['Gaussian is not compactly supported; not V24 projected matrices',
          'Prime and pole terms are not reclamped or certified here',
          'No complete zero sum or RH determination',
          'Analytic lemmas used (DLMF series, Euler gamma bound, midpoint theorem) are not formalized in Lean']}


def main():
    p=argparse.ArgumentParser();p.add_argument('--M',type=int,default=8192);p.add_argument('--terms',type=int,default=192);p.add_argument('--T',type=int,default=48);p.add_argument('--out',default='RESULTADOS_ES7.json')
    args=p.parse_args();data=run(args.M,args.terms,args.T)
    output=Path(args.out);output.write_text(json.dumps(data,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')
    print('ES7',output,[(e['a'],e['integral_interval_decimal'],e['derivative_interval_decimal']) for e in data['cases']])
if __name__=='__main__':main()
