"""NUMERIA ES8-P3: rational matrix enclosures for V17 archimedean term.

No mpmath, numpy or zeta-zero input is used for the certificates.  The
Guinand--Weil formula and V17 Fourier normalization are EXTERNAL premises.
The exponential-kernel identity and explicit sine-kernel integrals reduce
real half-line digamma/Fourier quadrature to a positive-term series; both
its finite terms and its Euler--Maclaurin tails are enclosed by Fractions.

Tests of the tapered matrix use the CONDITIONAL V24 absolute error envelope.
No identification with the historical, originally shipped E7/V16 matrix is
asserted or automatically validated; no universal positivity or RH claim.
"""
from __future__ import annotations
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial, isqrt, floor
from pathlib import Path
import argparse, json, platform, hashlib

from es8p2_prime_cert import Iv, I, L, SCALE, constants, projected_prime, projection_column, format_bound

BERN={2:Q(1,6),4:Q(-1,30),6:Q(1,42),8:Q(-1,30),10:Q(5,66),
      12:Q(-691,2730),14:Q(7,6),16:Q(-3617,510),18:Q(43867,798)}

def rise(s,k):
    z=1
    for j in range(k):z*=s+j
    return z

@lru_cache(maxsize=None)
def ln_pi() -> Iv:
    """log(pi) = log 2 + 2atanh((pi/2-1)/(pi/2+1)), monotone endpoints."""
    pi=constants()['pi']
    def log_x(x:Q):
        assert 1<x<2
        v=(x-1)/(x+1)
        s=sum((2*v**(2*j+1)/Q(2*j+1) for j in range(55)),Q())
        rem=2*v**111/(111*(1-v*v))
        return s,s+rem
    alo,_=log_x(pi.lo/2)
    _,bhi=log_x(pi.hi/2)
    return constants()['ln2']+Iv.rounded(alo,bhi)

@lru_cache(maxsize=1)
def gamma_rational()->Iv:
    """Euler constant by H_256-ln(256)-1/(512)+Bernoulli EM corrections.

    Standard real Euler--Maclaurin alternating remainder for harmonic sums
    is bounded by the magnitude of the next Bernoulli term (at N>0).
    This is an analytic, non-mechanized lemma explicitly recorded in README.
    """
    N=256
    H=sum((Q(1,k) for k in range(1,N+1)),Q(0))
    S=H-Q(1,2*N)
    for k in range(1,9):S+=BERN[2*k]/Q(2*k*N**(2*k))
    # Conservative periodic-Bernoulli EM remainder, using LAST retained
    # correction size (rather than assuming next-term alternating sharpness).
    rem=abs(BERN[16])/Q(16*N**16)
    return Iv.rounded(S-rem,S+rem)-8*constants()['ln2']

@lru_cache(maxsize=1)
def omega_zero()->Iv:
    pi=constants()['pi'];lg2=constants()['ln2']
    return -gamma_rational()-pi/2-3*lg2-ln_pi()

@lru_cache(maxsize=None)
def zeta_power_enclosure(s:int,q:Q,corrections:int=3)->Iv:
    """Hurwitz zeta(s,q) for integer s>1 and rational q>=1.

    Euler--Maclaurin via periodic Bernoulli remainder. Apply B2,..,B6,
    with the conservative periodic-Bernoulli remainder
    |R| <= |B6|/6! * (s)_5 * q^(-s-5). All exact rational.
    """
    assert s>=2 and q>=1 and corrections==3
    val=1/(Q(s-1)*q**(s-1))+1/(2*q**s)
    for k in range(1,corrections+1):
        val+=BERN[2*k]*rise(s,2*k-1)/Q(factorial(2*k)*q**(s+2*k-1))
    err=abs(BERN[6])*rise(s,5)/Q(factorial(6)*q**(s+5))
    return Iv.rounded(val-err,val+err)

@lru_cache(maxsize=None)
def exp_minus_interval(x:Q)->Iv:
    """Elementary rational exp(-x) bounds: alternating Taylor for x/2^k<=1,
    then repeated squaring. For x>=40, 0<=e^-x<=2^-floor(x) suffices.
    """
    assert x>=0
    if x==0:return I(1)
    if x>=40:return Iv.rounded(0,Q(1,2**(x.numerator//x.denominator)))
    k=1
    while x>k:k*=2
    y=x/k
    v=Q(0);term=Q(1)
    lo=hi=None
    for j in range(35):
        if j:term*=(-y)/j
        v+=term
        if j==33:lo=v
        if j==34:hi=v
    a=Iv.rounded(max(Q(0),lo),min(Q(1),hi))
    while k>1:
        a*=a;k//=2
    return a

@lru_cache(maxsize=None)
def mode(k:int):
    if not 1<=k<=10:raise ValueError('finite modes 1..10 only')
    w=I(k)*constants()['pi']/(2*L)
    return w,w*w/4

def ipow(v, k):
    """Nonnegative integer power for fixed-point interval."""
    if k<0:raise ValueError
    result=I(1)
    for _ in range(k):result=result*v
    return result

def power_h(i,j,k):
    """Complete homogeneous symmetric polynomial h_k(x_i,x_j)."""
    x,y=mode(i)[1],mode(j)[1]
    return sum((ipow(x,v)*ipow(y,k-v) for v in range(k+1)),I(0))

@lru_cache(maxsize=None)
def arch_tail(i:int,j:int,M:int,P:int=4)->Iv:
    """Rigorous tail from index n>=M. Series expansion in 1/q^2 with
    alternating geometric remainder; zeta-power enclosed via EM. The
    separate exponential factor tail is <= 2^-floor(4L(M+1/4)).
    """
    if M<32 or P!=4:raise ValueError('M>=32, P=4 required')
    if (i-j)%2:return I(0)
    q=Q(M)+Q(1,4)
    assert (mode(i)[1]+mode(j)[1]).hi < q*q, 'alternating tail requires (xi+xj)/q²<1'
    zcache={s:zeta_power_enclosure(s,q) for s in range(3,16)}
    out=I(0)
    if i==j:
        x=mode(i)[1]
        dsum=sum(((-1)**k*ipow(x,k+1)*zcache[3+2*k] for k in range(P+1)),I(0))
        rem=(ipow(x,P+2)*zcache[5+2*P]).hi
        out+=dsum+Iv.rounded(-rem,rem)
    wi,_=mode(i);wj,_=mode(j)
    coeff=wi*wj/(8*L)
    ksum=sum(((-1)**k*power_h(i,j,k)*zcache[4+2*k] for k in range(P+1)),I(0))
    rem=(power_h(i,j,P+1)*zcache[6+2*P]).hi
    exp_upper=Q(1,2**((4*L*q).numerator//(4*L*q).denominator))
    exp_err=(coeff*zcache[4]).hi*exp_upper
    ksum+=Iv.rounded(-rem-exp_err,rem+exp_err)
    return out-coeff*ksum

@lru_cache(maxsize=None)
def arch_raw(i:int,j:int,M:int=128)->Iv:
    if not 1<=i<=10 or not 1<=j<=10:raise ValueError('modes 1..10 only')
    if i>j:return arch_raw(j,i,M)
    if (i-j)%2:return I(0)
    if M<32 or M>2048:raise ValueError('32<=M<=2048 only')
    wi,xi=mode(i);wj,xj=mode(j)
    v=omega_zero() if i==j else I(0)
    coeff=wi*wj/(8*L)
    for n in range(M):
        q=Q(4*n+1,4)
        if i==j:v+=xi/(q*(q*q+xi))
        exp=exp_minus_interval(4*L*q)
        numer=I(1)-(-1)**i*exp
        v-=coeff*numer/((q*q+xi)*(q*q+xj))
    return v+arch_tail(i,j,M)

@lru_cache(maxsize=None)
def arch_projected(i:int,j:int,N:int,M:int=128)->Iv:
    if not 3<=i<=N<=10 or not 3<=j<=N:raise ValueError('bad modes')
    if i>j:return arch_projected(j,i,N,M)
    if (i-j)%2:return I(0)
    ai=projection_column(i,N);aj=projection_column(j,N)
    return sum((coef*arch_raw(p,q,M)*beta for p,coef in ai.items() for q,beta in aj.items()),I(0))

@lru_cache(maxsize=None)
def full_projected(i:int,j:int,N:int,M:int=128)->Iv:
    if (i-j)%2:return I(0)
    return arch_projected(i,j,N,M)+projected_prime(i,j,N)

def ldl_certify(matrix:list[list[Iv]],shift:Q=Q(0))->dict:
    """Validated unpivoted interval LDL. All lower pivots >0 proves SPD.
    Applying it to (A-shift I) gives a *certified spectral lower bound*.
    """
    n=len(matrix)
    if n<1 or any(len(row)!=n for row in matrix):raise ValueError('square matrix required')
    D=[];U=[[I(0) for j in range(n)] for i in range(n)]
    for i in range(n):
        d=matrix[i][i]-shift-sum((U[i][k]*U[i][k]*D[k] for k in range(i)),I(0))
        D.append(d)
        if d.lo<=0:
            return {'status':'INCONCLUSIVE_OR_NONPOSITIVE','failed_pivot':i+1,
                    'pivots':[v.as_json()['decimal_outward'] for v in D]}
        for j in range(i+1,n):
            t=matrix[j][i]-sum((U[j][k]*U[i][k]*D[k] for k in range(i)),I(0))
            U[j][i]=t/d
    return {'status':'CERTIFIED_POSITIVE_DEFINITE','shift':str(shift),
            'pivots':[v.as_json()['decimal_outward'] for v in D],
            'min_lower_pivot':format_bound(min(d.lo for d in D),True,14)}

@lru_cache(maxsize=1)
def sqrt_L()->Iv:
    v=(L.numerator*SCALE*SCALE)//L.denominator
    k=isqrt(v)
    return Iv(Q(k,SCALE),Q(k+1,SCALE))

@lru_cache(maxsize=None)
def perturb_bound_raw(i:int,j:int,delta:Q)->Q:
    """V24 *conditional* absolute bound E_ij: full arch + poles + 2-prime.
    All upper constants rounded outward from rational intervals.
    """
    if not 1<=i<=10 or not 1<=j<=10 or not 0<delta<L:raise ValueError
    pi=constants()['pi'];log2=constants()['ln2']
    s=sqrt_L();slo,shi=s.lo,s.hi
    def params(k):
        w=Q(k)*pi.hi/(2*L)
        lip=w/slo
        curv=w*w/slo
        return 2*delta*delta*lip, 2*L*curv+288*lip
    ei,di=params(i);ej,dj=params(j)
    B=2*shi; poleq=1/(1-L/2)
    arch=(12*B*(ei+ej)+16*(ei*dj+ej*di))/pi.lo
    pole=2*B*poleq*poleq*(ei+ej)
    prime=2*log2.hi/slo*(ei+ej)
    return arch+pole+prime

@lru_cache(maxsize=None)
def taper_projected_row_bound(N:int,delta:Q)->Q:
    """Elementwise |B^t E B|, then infinity norm. Projection uses the
    *same* B as the interval matrix. No direct use of E7/V16 pivots.
    """
    def absupper(v:Iv):return max(abs(v.lo),abs(v.hi))
    projected=[]
    for i in range(3,N+1):
        row=[]
        c=projection_column(i,N)
        for j in range(3,N+1):
            d=projection_column(j,N)
            row.append(sum((absupper(ci)*absupper(dj)*perturb_bound_raw(p,q,delta)
                            for p,ci in c.items() for q,dj in d.items()),Q(0)))
        projected.append(row)
    return max((sum(row,Q(0)) for row in projected))

def run(M:int=128,delta:Q=Q(1,1000000)):
    if M<32 or M>2048:raise ValueError('M out of bounds')
    outputs={};shift=Q(1,20)
    for N in (6,8,10):
        vals=[[full_projected(i,j,N,M) for j in range(3,N+1)] for i in range(3,N+1)]
        arch=[[arch_projected(i,j,N,M) for j in range(3,N+1)] for i in range(3,N+1)]
        assert all(vals[i][j]==vals[j][i] for i in range(N-2) for j in range(N-2))
        ldl=ldl_certify(vals,shift)
        perturb=taper_projected_row_bound(N,delta)
        outputs[str(N)]={
            'dimension':N-2,'block_sizes':[(N-2+1)//2,(N-2)//2],
            'arch_projected_intervals':[[v.as_json()['decimal_outward'] for v in row] for row in arch],
            'full_projected_intervals':[[v.as_json()['decimal_outward'] for v in row] for row in vals],
            'full_projected_scaled_integer_bounds':[[[str(int(v.lo*SCALE)),str(int(v.hi*SCALE))] for v in row] for row in vals],
            'certified_shifted_LDL':ldl,
            'taper_projected_row_error_upper':format_bound(perturb,False,18),
            'taper_error_smaller_than_shift':perturb<shift,
            'conditional_tapered_positive':ldl['status']=='CERTIFIED_POSITIVE_DEFINITE' and perturb<shift,
            'maximum_arch_entry_width':format_bound(max(v.width() for row in arch for v in row),False,14),
            'cross_parity_exact_zero':all(vals[i][j]==I(0) for i in range(N-2) for j in range(N-2) if (i-j)%2)}
    return {'schema':'NUMERIA.ES8P3.ARCH_PROJECTED_INTERVAL_LDL.1',
        'status':'CONDITIONAL_FINITE_PROJECTED_WEIL_POSITIVITY_CERTIFICATE_NOT_RH',
        'L':str(L),'M':M,'delta':str(delta),'shift_for_LDL':str(shift),
        'omega0_enclosure':omega_zero().as_json(),
        'gamma_enclosure':gamma_rational().as_json(),
        'methods':{'arch':'exact exponential kernel Fourier identity + sine overlap + rational finite terms + Hurwitz zeta Euler-Maclaurin remainder',
                   'projection':'V17 real parity and two-pole-annihilating columns B',
                   'prime':'ES8-P2 prime= -sqrt(2)log(2) symmetric convolution',
                   'taper':'V24 conditional absolute bounds, propagated through the same B',
                   'spectral':'rational outward interval LDL of M-(1/20)I'},
        'environment':{'python':platform.python_version(),'exact_fraction_arithmetic':True,'scale':str(SCALE)},
        'finite_matrices':outputs,
        'premises_not_automatically_verified':['Guinand-Weil explicit formula specialization to the real projected V17 basis',
           'V24 analytic derivative and Fourier tail bounds for the C-infinity taper',
           'Agreement of this freshly reconstructed matrix with original archived E7/V16 source code and files'],
        'not_proven':['RH, universal Weil positivity or every N',
           'Exact identity with historical E7/V16 matrix outputs',
           'Full Lean/Coq mechanization of Euler-Maclaurin and cutoff lemmas',
           'Certified spectral identification with an adelic Hilbert space operator']}

def main():
    p=argparse.ArgumentParser();p.add_argument('--M',type=int,default=128)
    p.add_argument('--delta',default='1/1000000');p.add_argument('--out',default='RESULTADOS_ES8P3.json')
    a=p.parse_args()
    if '/' in a.delta:
        x,y=a.delta.split('/');delta=Q(int(x),int(y))
    else:delta=Q(a.delta)
    result=run(a.M,delta)
    dest=Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    for N,r in result['finite_matrices'].items():
        print(f'ES8-P3 N={N} LDL={r["certified_shifted_LDL"]["status"]} gap>1/20; taper error <= {r["taper_projected_row_error_upper"]}; positive={r["conditional_tapered_positive"]}')
    print('STATUS: CONDITIONAL FINITE ONLY; historical E7/V16 comparison still open; NOT RH')

if __name__=='__main__':main()
