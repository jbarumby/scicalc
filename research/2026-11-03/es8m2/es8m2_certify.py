"""NUMERIA ES8-M2: rigorous rational enclosures of Gaussian Weil moments 0..3.

Based on classical Guinand-Weil: F(a)=2e^(a/4)+A(a)-P(a),
mu_k=(-1)^k (d/da)^k F(a), where P positive prime-power term.
The theorem is unconditional, but finite moment PSD does not prove RH.

A_k=(1/pi) integral_0^infty t^(2k) exp(-a t*t) Omega(t)dt.
For k=0..3, certified midpoint quadrature + Gaussian tail, carried
from the ES7 kernel; a=1/50 only. Output contains exact rational endpoints.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,platform
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
import es7_archimedean as es7
import es6_weil_gaussian_baseline as es6
import mpmath as mp

A=Q(1,50)
C=400
S=es7.S
S2=es6.SCALE


def interval(lo,hi=None):
    if hi is None:hi=lo
    assert lo<=hi
    return (Q(lo),Q(hi))

def plus(x,y):return (x[0]+y[0],x[1]+y[1])
def minus(x,y):return (x[0]-y[1],x[1]-y[0])
def mul(x,y):
    v=[a*b for a in x for b in y]
    return (min(v),max(v))
def scale(x,c):return mul(x,interval(Q(c)))
def square_pos(x):
    assert x[0]>=0
    return x[0]**2,x[1]**2

def fmt(x,digits=16):
    with localcontext() as ctx:
        ctx.prec=digits+55
        unit=Decimal(1).scaleb(-digits)
        return [format((Decimal(v.numerator)/Decimal(v.denominator)).quantize(unit,rounding=ROUND_FLOOR if j==0 else ROUND_CEILING),f'.{digits}f') for j,v in enumerate(x)]

def exact(x):return [[str(v.numerator),str(v.denominator)] for v in x]
def from_scaled(x,s):return (Q(int(x[0]),s),Q(int(x[1]),s))

def laguerre_derivative_factor(k,x):
    """Exact interval for p_k(x) in P^(k)=P*a^-k*p_k(x), k<=3."""
    if k==0:return interval(1)
    if k==1:return minus(x,interval(Q(1,2)))
    x2=square_pos(x)
    if k==2:return plus(minus(x2,scale(x,3)),interval(Q(3,4)))
    if k==3:
        x3=mul(x2,x)
        return minus(plus(minus(x3,scale(x2,Q(15,2))),scale(x,Q(45,4))),interval(Q(15,8)))
    raise ValueError('only k=0,1,2,3 supported')


def prime_moment_intervals(N=64):
    """Rigorous P^(k) for k=0..3, a=1/50, absolute tail <=1e-70.

    For u>=ln 64>4, x=u²/(4a)=12.5u²>200, |p_k(x)|<=x^k.
    Since Lambda(n)<=ln n, sum tail <= pref*625^k * integral_4^inf
    u^(2k+1)*exp(-12.5u²+u/2) du <= 5*625³*(4/91)*exp(-182)
    < 10^-70. The integral bound uses u/2<=u²/8 and 2k ln u<=u²
    for u>=4 and k<=3. The f(n) majorant is decreasing for n>=64.
    Here e>27/10, so exp(-182)<(10/27)^182, exact rational.
    """
    if N<64:raise ValueError('N>=64 needed for uniform tail proof')
    tail_proof_bound=Q(5*625**3*4,91)*Q(10,27)**182
    assert tail_proof_bound<Q(1,10**70)
    pilo,pihi=es6.pi_bounds()
    pia=(pilo//50,es6.ceildiv(pihi,50))
    pref_s=es6.inverse_positive(es6.sqrt_scaled(pia))
    pref=from_scaled(pref_s,S2)
    spf=es6.spf_sieve(N)
    sums=[interval(Q(0)) for _ in range(4)]
    count=0
    for n in range(2,N+1):
        p=es6.power_base(n,spf)
        if p is None:continue
        count+=1
        un=from_scaled(es6.log_int_bounds(n),S2)
        up=from_scaled(es6.log_int_bounds(p),S2)
        x=scale(square_pos(un),Q(25,2))
        # e^(-x), evaluated with rigorously rounded Taylor at the interval endpoints.
        xlo=x[0].numerator*S2//x[0].denominator
        xhi=es6.ceildiv(x[1].numerator*S2,x[1].denominator)
        expval=from_scaled(es6.exp_neg_bounds(xlo,xhi),S2)
        nsqrt=from_scaled(es6.inverse_positive(es6.sqrt_scaled((n*S2,n*S2))),S2)
        common=mul(pref,mul(up,mul(nsqrt,expval)))
        for k in range(4):
            term=mul(common,scale(laguerre_derivative_factor(k,x),50**k))
            sums[k]=plus(sums[k],term)
    # Enclose signed derivatives of the entire, not truncated, positive P(a)
    sums=[plus(x,interval(-Q(1,10**70),Q(1,10**70))) for x in sums]
    return {'N':N,'prime_powers_count':count,'P_derivatives_rational':[exact(x) for x in sums],
            'P_derivatives_decimals':[fmt(x,18) for x in sums],
            'tail_absolute_upper':'1/10^70 for each derivative order k=0..3',
            'tail_proof_upper_rational':[str(tail_proof_bound.numerator),str(tail_proof_bound.denominator)],
            'method':'exact rational interval evaluation with directed 1e-40 fixed point log, exp, sqrt; analytic tail proof'}


def upper_gaussian_integral_coeff(m,a,T):
    """I_m=integral_T^infty t^m exp(-a t²)dt <= exp(-aT²)*C_m.
    m even, IBP recurrence C_0=1/(2aT), C_m=T^(m-1)/(2a)+(m-1)C_(m-2)/(2a).
    """
    if m<0 or m%2:raise ValueError('m even >=0')
    v=Q(1,2*a*T)
    for j in range(2,m+1,2):v=Q(T**(j-1),2*a)+Q(j-1,2*a)*v
    return v


def arch_moments(M=16384,terms=192,T=48,scale_base=C):
    """Prove rigorous intervals for A_k/C^k, k=0..3 using ES7 Omega bounds."""
    if M<512 or M>131072 or terms<24 or T<40 or T>128:raise ValueError('unsupported certified parameters')
    h=Q(T,M)
    pilo,pihi=es7.pi_interval()
    assert pilo>3*S
    invpi=(S*S//pihi,es7.up(S*S,pilo))
    base=es7.omega0_interval()
    sums=[[0,0] for _ in range(4)]
    errors=[0]*4
    derivative_max_records=[0]*4
    previous_exp=es7.exp_negative_interval(Q(0))
    for j in range(M):
        left=j*h;right=(j+1)*h; mid=(left+right)/2
        om=es7.omega_interval(mid,terms,base)
        e=es7.exp_negative_interval(A*mid*mid)
        f=es7.mul(om,e)
        for k in range(4):
            q=mid**(2*k)/scale_base**k
            wt=es7.enc(q)
            g=es7.mul(f,wt)
            sums[k][0]+=g[0]; sums[k][1]+=g[1]
        # Monotone polynomial bounds for f'' and t^m f derivatives on [left,right].
        y=A*right*right
        B=7+18*right*right
        Cp=36+180*y+72*y*y+28*A*y+14*A
        Dp=36*right+2*A*right*B
        for k in range(4):
            m=2*k
            poly=(right**m*Cp
                  +(2*m*right**(m-1)*Dp if m else 0)
                  +(m*(m-1)*right**(m-2)*B if m>=2 else 0))/scale_base**k
            val=es7.up(previous_exp[1]*poly.numerator,poly.denominator)
            errors[k]+=val
        previous_exp=es7.exp_negative_interval(A*right*right)
    end_exp_hi=es7.exp_negative_interval(A*T*T)[1]
    out=[]
    for k in range(4):
        approx=es7.mul(es7.mul((sums[k][0],sums[k][1]),es7.enc(h)),invpi)
        # midpoint composite error <= h^3/(24*pi) sum sup second derivatives
        m_error=es7.up(errors[k]*h.numerator**3*S,24*h.denominator**3*pilo)
        tail_poly=(7*upper_gaussian_integral_coeff(2*k,A,T)+18*upper_gaussian_integral_coeff(2*k+2,A,T))/scale_base**k
        t_error=es7.up(end_exp_hi*tail_poly.numerator*S,tail_poly.denominator*pilo)
        total=(approx[0]-m_error-t_error,approx[1]+m_error+t_error)
        out.append({'order':k,'rational':exact(from_scaled(total,S)),
                    'display':fmt(from_scaled(total,S),14),
                    'approximate_midpoint':fmt(from_scaled(approx,S),14),
                    'error_midpoint_upper':fmt(interval(0,Q(m_error,S)),12)[1],
                    'error_tail_upper':fmt(interval(0,Q(t_error,S)),12)[1]})
    return {'a':'1/50','M':M,'T':T,'digamma_terms':terms,'scale_base':scale_base,
            'arch_moments_scaled_by_C':out,
            'integration_method':'composite midpoint plus analytic second derivative bound and Gaussian tails; interval operations exact fixed point'}


def extract_rational(raw):return tuple(Q(int(n),int(d)) for n,d in raw)

def combine_moments(arch,prime):
    """mu_k/C^k = [(-1)^k (pole^(k) + archA^(k) - P^(k))] / C^k.
    arch_moment(k) = (-1)^k A^(k) / C^k.
    pole moment = 2*(-1)^k*exp(a/4)/(4C)^k.
    prime contribution = (-1)^(k+1) P^(k)/C^k.
    """
    e_neg=es6.exp_neg_bounds(S2//200,S2//200)
    e_pos=es6.inverse_positive(e_neg)
    E=from_scaled(e_pos,S2)
    M=[]
    for k in range(4):
        Ak=extract_rational(arch['arch_moments_scaled_by_C'][k]['rational'])
        Pk=extract_rational(prime['P_derivatives_rational'][k])
        pole=scale(E,Q(2*(-1)**k,(4*C)**k))
        pp=scale(Pk,Q((-1)**(k+1),C**k))
        mk=plus(plus(pole,Ak),pp)
        M.append(mk)
    return M


def determinant_2x2_intervals(m):
    # Compute det H1=mu0*mu2-mu1² and det Stieltjes H1'=mu1*mu3-mu2².
    determinants=[]
    for shift in (0,1):
        d=minus(mul(m[shift],m[shift+2]),square_pos(m[shift+1]))
        # no false interpretation of tiny/numerically zero eigenvalues.
        diag_pos=m[shift][0]>0 and m[shift+2][0]>0
        state=('CERTIFIED_POSITIVE_DEFINITE' if diag_pos and d[0]>0
               else 'CERTIFIED_INDEFINITE' if d[1]<0
               else 'INCONCLUSIVE_INTERVAL')
        determinants.append({'type':'Hankel' if shift==0 else 'shifted_Stieltjes',
                             'matrix_intervals':[fmt(m[shift],14),fmt(m[shift+1],14),fmt(m[shift+2],14)],
                             'det_rational':exact(d),'det_decimal':fmt(d,14),'status':state,
                             'necessary_condition_only':True})
    return determinants


def classify_refs(folder):
    high_terms=('trace formula','riemann zeros and eigenvalue','li coefficients','spacings between zeros')
    medium_terms=('riemann zeta-function','random matrix','quantum chaos','adelic string')
    rows=[]; seen=set()
    for num in (7,8):
        with (folder/f'References_LeapSpace ({num}).csv').open(encoding='utf-8-sig',newline='') as fp:
            for r in csv.DictReader(fp):
                key=(r['DOI'].strip().lower() or 'title:'+r['Title'].lower().strip())
                if key in seen:continue
                seen.add(key)
                title=r['Title'].strip(); low=title.lower()
                grade='HIGH_CONTEXT' if any(t in low for t in high_terms) else ('SECONDARY' if any(t in low for t in medium_terms) else 'NOT_DIRECT_ES8M2')
                rows.append({'Ref':r['Reference'],'Title':title,'DOI':r['DOI'],'Grade':grade,
                             'Source CSV':str(num),'Caution':'Imported bibliographic metadata only; not independently validated theorem'})
    return rows


def numerical_crosscheck(m, dps=55):
    """High precision diagnostics of the exact full analytic expression, not certificate."""
    with mp.workdps(dps):
        a=mp.mpf(1)/50; cut=400;N=128
        def om(t):return (mp.re(mp.digamma(mp.mpf(1)/4+1j*t/2))-mp.log(mp.pi))/mp.pi
        shifts=[mp.mpf(0),mp.mpf(1),mp.mpf(2),mp.mpf(4),mp.mpf(8),mp.mpf(16),mp.mpf(32),mp.mpf(48),mp.mpf(64)]
        spf=es6.spf_sieve(N)
        vals=[]
        for k in range(4):
            ak=mp.quad(lambda t: (t*t/cut)**k *mp.exp(-a*t*t)*om(t),shifts)
            poles=2*(-1)**k*mp.exp(a/4)/(4*cut)**k
            pterms=[]
            for n in range(2,N+1):
                p=es6.power_base(n,spf)
                if p is None:continue
                u=mp.log(n);x=u*u/(4*a)
                polynomial= (1 if k==0 else x-mp.mpf(1)/2 if k==1 else
                             x*x-3*x+mp.mpf(3)/4 if k==2 else
                             x**3-mp.mpf(15)/2*x*x+mp.mpf(45)/4*x-mp.mpf(15)/8)
                prime=(mp.log(p)/mp.sqrt(n))/mp.sqrt(mp.pi*a)*mp.exp(-x)*(a**(-k))*polynomial
                pterms.append(prime)
            total=ak+poles+((-1)**(k+1))*mp.fsum(pterms)/cut**k
            vals.append(mp.nstr(total,35))
        return vals


def run(M=16384,terms=192,T=48,N=64,include_diagnostic=False):
    archive=arch_moments(M,terms,T)
    prime=prime_moment_intervals(N)
    moments=combine_moments(archive,prime)
    dets=determinant_2x2_intervals(moments)
    base=Path(__file__).resolve().parent
    refs=classify_refs(base/'sources')
    outcome={'schema':'NUMERIA.ES8M2.MOMENT_0_TO_3.RATIONAL.1',
            'status':'CERTIFIED finite rational enclosures conditional on established explicit formula and analytic lemmas; NOT RH proof',
            'a':'1/50','variable':'x=z_rho^2/400','moments_normalized_by_400_power':[
                {'k':k,'interval_rational':exact(v),'display':fmt(v,14),
                 'certified_sign':'POSITIVE' if v[0]>0 else 'NEGATIVE' if v[1]<0 else 'INCONCLUSIVE'}
                for k,v in enumerate(moments)],
            'hankel_and_stieltjes_2x2':dets,
            'archimedean':archive,'prime_derivatives':prime,
            'references':{'n_unique':len(refs), 'n_high_context':sum(r['Grade']=='HIGH_CONTEXT' for r in refs),'source':'two user-provided LeapSpace CSVs; bibliographic triage, no article claims adopted'},
            'limitations':['No proof of universal Weil positivity or RH',
             'Gaussian moments are not entries of E7/V24 projected matrices',
             'Classification of bibliography is metadata-only, and citations may contain errors',
             'No exhaustive verification of off-critical zeros or all dimensions',
             'Only first two 2x2 moment matrices checked; higher Hankel dimensions require moments up to at least mu_4 and mu_5']}
    if include_diagnostic:outcome['floating_point_crosscheck_NOT_CERTIFIED']=numerical_crosscheck(moments)
    return outcome,refs


def main():
    p=argparse.ArgumentParser();p.add_argument('--M',type=int,default=16384);p.add_argument('--terms',type=int,default=192);p.add_argument('--T',type=int,default=48);p.add_argument('--N',type=int,default=64)
    p.add_argument('--out',default='RESULTADOS_ES8M2.json');p.add_argument('--refs-out',default='REFERENCIAS_TRIADAS_ES8M2.csv');p.add_argument('--diagnostic',action='store_true')
    x=p.parse_args();out,refs=run(x.M,x.terms,x.T,x.N,x.diagnostic)
    Path(x.out).write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    with Path(x.refs_out).open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=refs[0].keys());writer.writeheader();writer.writerows(refs)
    print('ES8M2',len(refs),'references',[(m['k'],m['display']) for m in out['moments_normalized_by_400_power']],[(d['type'],d['status']) for d in out['hankel_and_stieltjes_2x2']])

if __name__=='__main__':main()
