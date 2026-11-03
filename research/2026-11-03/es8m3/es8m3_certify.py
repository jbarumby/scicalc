"""NUMERIA ES8-M3: rational enclosures for Gaussian Weil moments k=0..5.

The global explicit formula is used unconditionally. Does not imply universal
Weil positivity, RH, or an identification with E7/V24 projected matrices.

New gates: 4th-order composite Simpson error, Euler-Maclaurin digamma tail,
accurate Euler-gamma bound, 5th order differentiated prime term and tail.
"""
from __future__ import annotations
import argparse, json, math, platform, time
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from functools import lru_cache
import es7_archimedean as es7
import es6_weil_gaussian_baseline as es6

A=Q(1,50); C=400; S=es7.S; S2=es6.SCALE

def interval(lo,hi=None):
    if hi is None:hi=lo
    assert lo<=hi
    return Q(lo),Q(hi)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def sub(x,y):return x[0]-y[1],x[1]-y[0]
def mul(x,y):
    z=[u*v for u in x for v in y]
    return min(z),max(z)
def scl(x,k):return mul(x,interval(Q(k)))
def exact(x):return [[str(v.numerator),str(v.denominator)] for v in x]
def from_exact(x):return tuple(Q(int(n),int(d)) for n,d in x)
def fmt(x,d=15):
    with localcontext() as ctx:
        ctx.prec=d+65;step=Decimal(1).scaleb(-d)
        return [format((Decimal(q.numerator)/Decimal(q.denominator)).quantize(step,rounding=ROUND_FLOOR if i==0 else ROUND_CEILING),f'.{d}f') for i,q in enumerate(x)]
def fixed(x,scale=S):return (x[0]*Q(1,scale),x[1]*Q(1,scale))
def enc(x,scale=S):
    return (x.numerator*scale//x.denominator,es7.up(x.numerator*scale,x.denominator))

@lru_cache(maxsize=1)
def gamma_precise():
    """Euler-Maclaurin bounds: 1/(2n)-1/(12n²) < Hn-log(n)-gamma
    < 1/(2n)-1/(12n²)+1/(120n⁴), n>=1.
    Proof by alternating Bernoulli remainder for harmonic numbers.
    """
    n=4096
    H=(sum(S//j for j in range(1,n+1)),sum(es7.up(S,j) for j in range(1,n+1)))
    lo2,hi2=es7.ln2_interval()
    raw=es7.sub(H,(lo2*12,hi2*12))
    d=Q(1,2*n)-Q(1,12*n*n)
    rlo,rhi=enc(d)[0],enc(d+Q(1,120*n**4))[1]
    return raw[0]-rhi,raw[1]-rlo

@lru_cache(maxsize=1)
def omega0():
    pi=es7.pi_interval()
    log2=es7.ln2_interval()
    logpi=es7.lnpi_interval(pi)
    a=(pi[0]//2,es7.up(pi[1],2))
    return es7.sub(es7.sub(es7.sub((0,0),gamma_precise()),a), (3*log2[0]+logpi[0],3*log2[1]+logpi[1]))

@lru_cache(maxsize=50000)
def omega_interval_em(t,terms=128):
    """Digamma series plus Euler-Maclaurin enclosure for its positive tail.

    f(q)=t²/[q(4q²+t²)], q=n+1/4. f'' and f'''' positive,
    f''' negative for q>0. With p=2 Euler-Maclaurin,
      sum_{n=m}∞ f(n+1/4)
         = I + f(m)/2 - f'(m)/12 + f'''(m)/720 + R,
      |R| <= (1/720) int_m∞ f''''(x) dx = -f'''(m)/720.
    Hence inclusive bounds: I+f/2-f'/12+f'''/360 <= tail <= I+f/2-f'/12.
    All endpoints computed via rational operations and log1p enclosure.
    """
    assert t>=0 and terms>=8
    v=t*t
    partial_lo=partial_hi=0
    tn=t.numerator;td=t.denominator
    num=16*tn*tn
    for n in range(terms):
        u=4*n+1
        den=u*(u*u*td*td+4*tn*tn)
        partial_lo+=num*S//den
        partial_hi+=es7.up(num*S,den)
    q=Q(terms)+Q(1,4)
    fq=v/(q*(4*q*q+v))
    fp=-v*(12*q*q+v)/(q*q*(4*q*q+v)**2)
    fppp=-6*v*(640*q**6+80*q**4*v+16*q*q*v*v+v**3)/(q**4*(4*q*q+v)**4)
    # Integral exactly = (1/2)*log(1+v/(4*q²)).
    if v>4*q*q: raise ValueError("For log1p bound use terms at least t/2")
    il,ih=es7.log1p_interval(v/(4*q*q),24)
    I=(Q(il,2*S),Q(ih,2*S))
    center=fq/2-fp/12
    tail=(I[0]+center+fppp/360, I[1]+center)
    base=omega0()
    return base[0]+partial_lo+enc(tail[0])[0],base[1]+partial_hi+enc(tail[1])[1]

@lru_cache(maxsize=1)
def omega_derivative_bounds():
    # For l>=1: |b_n^(l)(t)| <= l!/(2^l q_n^(l+1)),
    # and sum q_n^(-p) <= 4^p + 4^(p-1)/(p-1), p>1.
    v=[None]
    for l in range(1,5):
        v.append(Q(math.factorial(l),2**l)*(4**(l+1)+Q(4**l,l)))
    return v

def fourth_bound(r,k):
    """Nonnegative rational polynomial bound for d^4/dt^4[(t²/C)^k e^-at² Omega]."""
    omega=[7+18*r*r]+omega_derivative_bounds()[1:]
    g=[Q(1),2*A*r,4*A*A*r*r+2*A,8*A**3*r**3+12*A*A*r,16*A**4*r**4+48*A**3*r*r+12*A*A]
    out=Q(0)
    for i in range(min(4,2*k)+1):
        dpoly=Q(math.factorial(2*k),math.factorial(2*k-i))*r**(2*k-i)/C**k
        for j in range(5-i):
            l=4-i-j
            coeff=Q(math.factorial(4),math.factorial(i)*math.factorial(j)*math.factorial(l))
            out+=coeff*dpoly*g[j]*omega[l]
    return out

def tail_gaussian_coeff(m,a,T):
    """bound int_T∞ t^m exp(-a t²) dt <= exp(-aT²)*C_m, m even."""
    assert m>=0 and m%2==0
    v=Q(1,2*a*T)
    for j in range(2,m+1,2):v=Q(T**(j-1),2*a)+Q(j-1,2*a)*v
    return v

def arch_simpson(M=8192,terms=128,T=48):
    if M<1024 or M%2 or M>131072 or terms<16 or T<40 or T>64:
        raise ValueError('M must be even >=1024; 16<=terms; 40<=T<=64')
    h=Q(T,M)
    pi=es7.pi_interval()
    invpi=(Q(S,pi[1]),Q(S,pi[0]))
    sums=[[0,0] for _ in range(6)]
    for j in range(M+1):
        t=h*j;w=1 if j in (0,M) else 4 if j%2 else 2
        om=omega_interval_em(t,terms)
        e=es7.exp_negative_interval(A*t*t)
        v=es7.mul(om,e)
        factor=Q(1)
        for k in range(6):
            p=es7.mul(v,es7.enc(factor*w))
            sums[k][0]+=p[0];sums[k][1]+=p[1]
            factor*=t*t/C
    arch=[]
    # Per-block Simpson error <= h^5 / 90 * sup |f''''|, block width 2h.
    errors=[0]*6
    expL=es7.exp_negative_interval(Q(0))[1]
    for j in range(M//2):
        l=2*j*h;r=l+2*h
        for k in range(6):
            errors[k]+=es7.up(expL*fourth_bound(r,k).numerator,fourth_bound(r,k).denominator)
        expL=es7.exp_negative_interval(A*r*r)[1]
    expTail=es7.exp_negative_interval(A*T*T)[1]
    for k in range(6):
        approx=mul(mul(fixed(sums[k]),interval(h/3)),invpi)
        err=Q(errors[k],S)*h**5/90*invpi[1]
        # |Omega|<=7+18t² from ES7.
        tail=(Q(expTail,S)*invpi[1]
              *(7*tail_gaussian_coeff(2*k,A,T)+18*tail_gaussian_coeff(2*k+2,A,T))/C**k)
        full=add(approx,interval(-err-tail,err+tail))
        arch.append({'k':k,'interval_rational':exact(full),'display':fmt(full,14),
                     'simpson_error_upper':fmt(interval(0,err),12)[1],
                     'gaussian_tail_upper':fmt(interval(0,tail),12)[1]})
    return {'M':M,'T':T,'terms':terms,'a':'1/50','scale_base':C,
            'method':'interval Simpson with local fourth-derivative bound; Euler-Maclaurin digamma tail; Bernoulli gamma bound',
            'moments':arch}

@lru_cache(maxsize=1)
def p_coeffs():
    """p_k(x) coefficient array ascending; p_(k+1)=(x-k-1/2)p_k-x p_k'."""
    out=[];p=[Q(1)]
    for k in range(6):
        out.append(p)
        nxt=[Q(0)]*(len(p)+1)
        for i,a in enumerate(p):
            nxt[i]+=-(Q(k)+Q(1,2)+i)*a
            nxt[i+1]+=a
        p=nxt
    return out

def eval_poly_interval(coeffs,x):
    z=interval(coeffs[-1])
    for c in reversed(coeffs[:-1]):z=add(mul(z,x),interval(c))
    return z

def prime_derivatives(N=64):
    if N<64:raise ValueError('N >=64 is required for the tail proof')
    # If x>=200, |p_k(x)|/x^k <= 2, proven for each k=0..5.
    pols=p_coeffs()
    assert all(sum(abs(c)*Q(200)**(i-k) for i,c in enumerate(p))<2 for k,p in enumerate(pols))
    # Tail: each term ≤10*625^5 * u^(11) exp(-25u²/2+u/2);
    # since 2k log u<=u² (u>=4,k<=5), its sum ≤
    # 10*625^5 ∫_4∞ u exp(-91 u²/8)du = (40/91)*625^5*exp(-182).
    bound=Q(40*625**5,91)*Q(10,27)**182
    assert bound<Q(1,10**64)
    pi=es6.pi_bounds()
    pi_a=(pi[0]//50,es6.ceildiv(pi[1],50))
    pref=fixed(es6.inverse_positive(es6.sqrt_scaled(pi_a)),S2)
    sieve=es6.spf_sieve(N)
    sums=[interval(0) for _ in range(6)]
    for n in range(2,N+1):
        p=es6.power_base(n,sieve)
        if p is None:continue
        un=fixed(es6.log_int_bounds(n),S2)
        up=fixed(es6.log_int_bounds(p),S2)
        x=scl(mul(un,un),Q(25,2))
        xl=x[0].numerator*S2//x[0].denominator
        xu=es6.ceildiv(x[1].numerator*S2,x[1].denominator)
        expv=fixed(es6.exp_neg_bounds(xl,xu),S2)
        invsqrt=fixed(es6.inverse_positive(es6.sqrt_scaled((n*S2,n*S2))),S2)
        base=mul(pref,mul(up,mul(invsqrt,expv)))
        for k,pco in enumerate(pols):
            sums[k]=add(sums[k],scl(mul(base,eval_poly_interval(pco,x)),50**k))
    sums=[add(z,interval(-Q(1,10**64),Q(1,10**64))) for z in sums]
    return {'N':N,'a':'1/50','tail_bound_per_derivative':'1e-64',
            'derivatives':[{'k':k,'interval_rational':exact(z),'display':fmt(z,18)} for k,z in enumerate(sums)],
            'method':'exact rational interval prime powers with Laguerre recurrence and uniform analytic tail k<=5'}

def moments_and_matrices(arch,prime):
    # e^(a/4) = 1/e^(-1/200) with exact rational enclosure.
    ee=es6.exp_neg_bounds(S2//200,S2//200)
    e=fixed(es6.inverse_positive(ee),S2)
    m=[]
    for k in range(6):
        ak=from_exact(arch['moments'][k]['interval_rational'])
        pk=from_exact(prime['derivatives'][k]['interval_rational'])
        pole=scl(e,Q(2*(-1)**k,(4*C)**k))
        pterm=scl(pk,Q((-1)**(k+1),C**k))
        m.append(add(add(pole,ak),pterm))
    def det2(a,b,c):return sub(mul(a,c),mul(b,b))
    def det3(vals):
        a,b,c,d,e=vals
        return add(sub(mul(a,det2(c,d,e)),mul(b,sub(mul(b,e),mul(c,d)))),mul(c,sub(mul(b,d),mul(c,c))))
    out=[]
    for sh in (0,1):
        mm=m[sh:sh+5]
        d1=mm[0];d2=det2(*mm[:3]);d3=det3(mm)
        status='CERTIFIED_POSITIVE_DEFINITE' if all(z[0]>0 for z in (d1,d2,d3)) else 'CERTIFIED_NOT_PSD' if any(z[1]<0 for z in (d1,d2,d3)) else 'INCONCLUSIVE_INTERVAL'
        out.append({'type':'Hankel_3x3' if sh==0 else 'Stieltjes_shift_3x3',
                    'moments_used':[sh+i for i in range(5)],'leading_principal_minors':[
                        {'size':j+1,'interval_rational':exact(z),'display':fmt(z,18)} for j,z in enumerate((d1,d2,d3))],
                    'status':status,'does_not_establish_RH':True})
    return m,out

def run(M=8192,terms=128,T=48,N=64):
    start=time.monotonic()
    arch=arch_simpson(M,terms,T)
    prime=prime_derivatives(N)
    mom,matrices=moments_and_matrices(arch,prime)
    return {'schema':'NUMERIA.ES8M3.MOMENT_0_TO_5.RATIONAL.1','date':'2026-11-03',
       'status':'finite certified enclosures conditional on classical explicit formula and explicitly derived analytic bounds; no RH proof',
       'parameters':{'a':'1/50','M':M,'digamma_terms':terms,'T':T,'prime_cutoff':N,'normalization':'m_k=mu_k/400^k'},
       'environment':{'python':platform.python_version()},
       'moments':[{'k':k,'interval_rational':exact(z),'display':fmt(z,15)} for k,z in enumerate(mom)],
       'matrix_certificates':matrices,'archimedean':arch,'primes':prime,
       'runtime_seconds_approximate':round(time.monotonic()-start,3),
       'limitations':['No proof of Weil positivity for all test functions or of RH',
         'No direct E7/V24 matrix-entry certification or projected LDL comparison',
         'Finite Gaussian Hankel/Stieltjes moments do not imply a complete positive moment measure of the zeta zeros',
         'Classical analytic tail/derivative lemmas require independent line-by-line review prior to any strong theorem claim',
         'No assumption that all zeta zeros lie on the critical line']}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--M',type=int,default=8192)
    p.add_argument('--terms',type=int,default=128)
    p.add_argument('--T',type=int,default=48)
    p.add_argument('--N',type=int,default=64)
    p.add_argument('--out',default='RESULTADOS_ES8M3.json')
    z=p.parse_args()
    r=run(z.M,z.terms,z.T,z.N)
    Path(z.out).write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('ES8M3',[(i['k'],i['display']) for i in r['moments']])
    print('matrix gates',[(d['type'],d['status'],d['leading_principal_minors'][-1]['display']) for d in r['matrix_certificates']])

if __name__=='__main__':main()
