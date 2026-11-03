"""NUMERIA ES8-M: exact rational Stieltjes fixtures + limited genuine Gaussian moments.

Validated: finite rational model PSD/indefiniteness, exact rational enclosures
for full mu_0 and mu_1 of Guinand-Weil from prior ES6/ES7 certified pieces.
Exploratory: matrices from a finite list of known critical-line zeros; NOT
certificate of all zeros or of complete Weil positivity.
"""
from __future__ import annotations
import argparse, hashlib, json, platform
from fractions import Fraction as Q
from pathlib import Path
import mpmath as mp

from es6_weil_gaussian_baseline import (
    SCALE, ceildiv, exp_neg_bounds, inverse_positive, positive_product,
    sqrt_scaled, pi_bounds, log_int_bounds, spf_sieve, power_base,
    prime_gaussian_certificate,
)
ROOT=Path(__file__).resolve().parent

def interval_from_scaled(v,scale):
    lo,hi=map(int,v)
    S=int(scale)
    assert lo<=hi
    return (Q(lo,S),Q(hi,S))

def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul_pos(a,q):
    assert q>=0
    return (q*a[0],q*a[1])
def fmt_fraction(x,digits=22):
    from decimal import Decimal,localcontext
    with localcontext() as ctx:
        ctx.prec=digits+50
        return format(Decimal(x.numerator)/Decimal(x.denominator),f'.{digits}f')
def fmt_interval(a,digits=22):
    # Exact rational endpoints also stored separately. Decimal display is NOT relied upon
    # as an outward-rounded enclosure.
    return [fmt_fraction(x,digits) for x in a]
def exact_interval(a):return [[str(x.numerator),str(x.denominator)] for x in a]

def prime_derivative_cert(N=64):
    """Rigorous fixed-point enclosure for P'(a), a=1/50, P positive prime term.
    For n>=2 coefficient c_n = -25 +625 (ln n)^2 >0.
    For N>=64, bound whole omitted tail by 1e-37 (see report).
    """
    if N<64: raise ValueError('N>=64 required for derivative tail proof')
    pi_lo,pi_hi=pi_bounds()
    pref=inverse_positive(sqrt_scaled((pi_lo//50,ceildiv(pi_hi,50))))
    sieve=spf_sieve(N)
    tot=(0,0)
    for n in range(2,N+1):
        p=power_base(n,sieve)
        if p is None:continue
        lo,hi=log_int_bounds(n)
        xlo=25*lo*lo//(2*SCALE)
        xhi=ceildiv(25*hi*hi,2*SCALE)
        decay=exp_neg_bounds(xlo,xhi)
        w=positive_product(log_int_bounds(p),inverse_positive(sqrt_scaled((n*SCALE,n*SCALE))))
        term=positive_product(pref,positive_product(w,decay))
        coeff=(625*lo*lo//SCALE-25*SCALE,ceildiv(625*hi*hi,SCALE)-25*SCALE)
        assert coeff[0]>0
        der=positive_product(term,coeff)
        tot=tot[0]+der[0],tot[1]+der[1]
    assert log_int_bounds(N)[0]>4*SCALE
    # Prime tail: for x=ln t >=4, x^2 exp(-12.5*x^2) <= exp(-6.25*x^2)
    # and exp(x/2)<=exp(x*x/8).  Integrating y exp(-(49/8)y^2)
    # yields exp(-98)/(49/4); 1/sqrt(pi/50)<5. Therefore
    # 0 <= P'_tail <= 625*5*(4/49)*exp(-98) < 10^-37.
    assert exp_neg_bounds(98*SCALE,98*SCALE)[1]<=1
    return {
        'a':'1/50','N':N,'scaled_exact':[str(tot[0]),str(tot[1]+1000)],
        'display':fmt_interval((Q(tot[0],SCALE),Q(tot[1]+1000,SCALE)),24),
        'tail_proof':'0<=tail <= (12500/49) exp(-98) < 10^-37 for N>=64; positivity and monotonicity hold for log n >=4',
        'tail_upper_scaled':'1000','status':'CERTIFIED PRIME DERIVATIVE P_prime(a), not the full moment'
    }

def full_moments_0_1_certified(baseline:dict,N=64):
    row=next(c for c in baseline['cases'] if c['a']=='1/50')
    s=row['scale']
    a0=interval_from_scaled(row['interval_scaled_exact'],s)
    a1=interval_from_scaled(row['derivative_interval_scaled_exact'],s)
    pc=prime_gaussian_certificate(N)
    p0=interval_from_scaled(pc['prime_term_all_n_exact_scaled'],pc['scale'])
    pdc=prime_derivative_cert(N)
    p1=interval_from_scaled(pdc['scaled_exact'],SCALE)
    # exp(a/4)=exp(1/200)=1/exp(-1/200), with rational directed rounding.
    e_neg=exp_neg_bounds(SCALE//200,SCALE//200)
    e_pos=inverse_positive(e_neg)
    e=(Q(e_pos[0],SCALE),Q(e_pos[1],SCALE))
    # F(a)=2e^(a/4)+A(a)-P(a)
    # mu_0=F(a); mu_1= -F'(a)=P'(a)-A'(a)-(1/2)e^(a/4)
    m0=sub(add(mul_pos(e,2),a0),p0)
    m1=sub(sub(p1,a1),mul_pos(e,Q(1,2)))
    assert m0[0]>0 and m1[0]>0
    return {
        'a':'1/50','mu0':{'rational_interval':exact_interval(m0),'display_not_outward':fmt_interval(m0,20),'certified_positive':True},
        'mu1':{'rational_interval':exact_interval(m1),'display_not_outward':fmt_interval(m1,20),'certified_positive':True},
        'prime_derivative':pdc,
        'precision_scope':'CERTIFIED ONLY mu0 and mu1 of full explicit formula; n=0 scalar Hankel/Stieltjes tests, NOT n>=1 matrices',
        'components':{'arch_integral':row['integral_interval_decimal'],'arch_derivative':row['derivative_interval_decimal'],
                      'prime_P':pc['positive_prime_term_all_n_interval'],'pole_exp':fmt_interval(e,20)}
    }

def moments_of_finite_rational_positive_atoms(nodes,weights,kmax):
    assert len(nodes)==len(weights)>0 and all(x>=0 for x in nodes) and all(w>0 for w in weights)
    return [sum((w*x**k for x,w in zip(nodes,weights)),Q(0)) for k in range(kmax+1)]

def hankel(m,degree,shift=0):
    return [[m[i+j+shift] for j in range(degree+1)] for i in range(degree+1)]

def exact_ldl(A):
    """Exact symmetric LDL. Zero pivot returns inconclusive; negative pivot disproves PSD."""
    n=len(A)
    assert all(len(row)==n for row in A)
    assert all(A[i][j]==A[j][i] for i in range(n) for j in range(n))
    L=[[Q(int(i==j)) for j in range(n)] for i in range(n)]
    ds=[]
    for k in range(n):
        d=A[k][k]-sum((L[k][j]*L[k][j]*ds[j] for j in range(k)),Q(0))
        ds.append(d)
        if d==0:
            return {'state':'INCONCLUSIVE_ZERO_PIVOT','pivots':[str(x) for x in ds]}
        if d<0:
            return {'state':'CERTIFIED_NOT_PSD','pivots':[str(x) for x in ds]}
        for i in range(k+1,n):
            num=A[i][k]-sum((L[i][j]*L[k][j]*ds[j] for j in range(k)),Q(0))
            L[i][k]=num/d
    return {'state':'CERTIFIED_PD_RATIONAL_FIXTURE','pivots':[str(x) for x in ds]}

def toy_benchmarks():
    nodes=list(map(Q, [1,2,3,4,5]))
    weights=list(map(Q, [1,2,3,5,8]))
    moments=moments_of_finite_rational_positive_atoms(nodes,weights,7)
    positive={str(n):{'H':exact_ldl(hankel(moments,n,0)),
                      'H_shift':exact_ldl(hankel(moments,n,1))} for n in (0,1,2,3)}
    assert all(v[k]['state']=='CERTIFIED_PD_RATIONAL_FIXTURE' for v in positive.values() for k in ('H','H_shift'))
    # Adversarial conjugate non-real atoms 1±i, each weight 1:
    # m0=2,m1=2,m2=0,m3=-4, H1=[[2,2],[2,0]] indefinite.
    anti=[Q(2),Q(2),Q(0),Q(-4)]
    negative={'H':exact_ldl(hankel(anti,1,0)),'H_shift':exact_ldl(hankel(anti,1,1))}
    assert negative['H']['state']=='CERTIFIED_NOT_PSD'
    return {'positive_atoms':{'nodes':list(map(str,nodes)),'weights':list(map(str,weights)),
                               'moment_first_8':list(map(str,moments)),'matrices':positive},
            'offaxis_fixture':{'description':'conjugate complex atom locations 1 +/- i with equal positive weights; NOT zeta zeros',
                               'first_four_moments':list(map(str,anti)),'H1':negative['H'],'H1_shift':negative['H_shift']},
            'status':'exact rational arithmetic, applies ONLY to these artificial finite atomic fixtures'}

def laguerre_minus_half(k,x):
    if k==0:return mp.mpf(1)
    old=mp.mpf(1);curr=mp.mpf('0.5')-x
    for n in range(1,k):
        nxt=((2*n+mp.mpf('0.5')-x)*curr-(n-mp.mpf('0.5'))*old)/(n+1)
        old,curr=curr,nxt
    return curr

def exploratory_known_critical_zeros_and_full_formula(a='0.02',J=32,nmax=5,dps=45):
    """Non-interval comparison; full formula side is analytic, zeros are incomplete."""
    if J<10 or nmax>8 or dps<30: raise ValueError('J>=10,nmax<=8,dps>=30')
    with mp.workdps(dps):
        a=mp.mpf(a)
        gammas=[mp.im(mp.zetazero(i)) for i in range(1,J+1)]
        scale=mp.mpf(400)
        moments_zeros=[2*mp.fsum((g*g/scale)**k *mp.exp(-a*g*g) for g in gammas) for k in range(nmax+1)]
        # Derivatives of known Guinand-Weil analytic side, via integral against t^(2k).
        N=128
        spf=spf_sieve(N)
        powers=[(n,power_base(n,spf)) for n in range(2,N+1)]
        powers=[(n,p) for n,p in powers if p is not None]
        def kernel(t):return (mp.re(mp.digamma(mp.mpf('0.25')+1j*t/2))-mp.log(mp.pi))/mp.pi
        cuts=[0,1,2,4,8,16,32,64,96]
        full=[]
        for k in range(nmax+1):
            arch=mp.quad(lambda t:t**(2*k)*mp.exp(-a*t*t)*kernel(t),cuts)
            pole=2*((-mp.mpf(1)/4)**k)*mp.exp(a/4)
            pp=[]
            for n,p in powers:
                u=mp.log(n);x=u*u/(4*a)
                P=(mp.log(p)/mp.sqrt(n))*mp.exp(-x)/mp.sqrt(mp.pi*a)
                pp.append(P*a**(-k)*mp.factorial(k)*laguerre_minus_half(k,x))
            mu=(pole+arch-mp.fsum(pp))/(scale**k)
            full.append(mu)
        diffs=[abs(p-q) for p,q in zip(full,moments_zeros)]
        matrices=[]
        for n in (0,1,2):
            if 2*n+1>nmax:continue
            h0=mp.matrix(hankel(moments_zeros,n,0))
            h1=mp.matrix(hankel(moments_zeros,n,1))
            ev0=mp.eigsy(h0,eigvals_only=True);ev1=mp.eigsy(h1,eigvals_only=True)
            matrices.append({'n':n,'matrix_size':n+1,
                             'min_eigen_known_zeros_H':mp.nstr(min(ev0),18),
                             'min_eigen_known_zeros_shifted':mp.nstr(min(ev1),18),
                             'status':'NONCERTIFIED finite critical-zero sample only'})
        return {'a':mp.nstr(a,10),'zero_count':J,'work_dps':dps,'monomial_scale':'y=gamma^2/400',
                'moments_full_explicit_rhs_numeric':[mp.nstr(x,28) for x in full],
                'moments_known_critical_zeros_numeric':[mp.nstr(x,28) for x in moments_zeros],
                'absolute_difference':[mp.nstr(x,15) for x in diffs],
                'Hankel_finite_critical_zero_eigenvalues':matrices,
                'note':'ALL entries and eigenvalues in this section are NONCERTIFIED; zero list not completeness-certified; no info about possible off-line zeros'}

def run(zeros=32,degree=5,dps=45):
    baseline=json.loads((ROOT/'ES7_RESULTADOS_BASELINE.json').read_text(encoding='utf-8'))
    r={'schema':'NUMERIA.ES8M.MOMENTS_STIELTJES.1',
       'scientific_status':'certified full mu0,mu1 only; rational fixtures; all higher actual moments exploratory; NOT RH proof, not full Weil PSD',
       'environment':{'python':platform.python_version(),'mpmath':mp.__version__},
       'certified_full_gaussian_low_moments':full_moments_0_1_certified(baseline),
       'exact_rational_stieltjes_control':toy_benchmarks(),
       'exploratory_higher_moments':exploratory_known_critical_zeros_and_full_formula(J=zeros,nmax=degree,dps=dps),
       'gates':{'H_n_actual_zeta_for_n_ge_1':'NOT_CERTIFIED',
                'Stieltjes_shift_actual_zeta_for_n_ge_1':'NOT_CERTIFIED',
                'Carleman_equivalence_to_RH':'NOT_ESTABLISHED',
                'Weil_universal_positivity':'NOT_PROVED',
                'E7_V24_matrix_identification':'NOT_ESTABLISHED'},
       'original_sources':{'DeepSeek_tex':'sources/deepseek_latex_20261103_74b66a.tex',
                           'Qwen_txt':'sources/Qwen_latex_20261103_pag4a82de.txt'}}
    return r

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',default='RESULTADOS_ES8M.json')
    p.add_argument('--zeros',type=int,default=32)
    p.add_argument('--degree',type=int,default=5)
    p.add_argument('--dps',type=int,default=45)
    a=p.parse_args()
    r=run(a.zeros,a.degree,a.dps)
    Path(a.out).write_text(json.dumps(r,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print('ES8M certified mu0 and mu1 positive; exact rational surrogate controls pass; higher Hankel NOT CERTIFIED; RH UNPROVED')

if __name__=='__main__':main()
