#!/usr/bin/env python3
"""P245 Gate B: exact-integer directed fixed-point interval replay of selected GL60 entries.

No Arb, MPFR, SciPy or mpmath arithmetic is used for certification.  mpmath
is used solely to propose Legendre-root brackets; EVERY bracket is proved by
integer-rational polynomial sign tests.  All special-function bounds are
analytic and all rounding is outward, by integer floor/ceiling.

Scope: 48 GL60 sample-entry inclusions across 8 index pairs, two parities,
panels 0,503,1007. High-degree enclosures are coarse Bessel majorants.  This does not redo all 142506 matrix entries or certify
identification with the Weil/Zhu quadratic form.
"""
from __future__ import annotations
import argparse, hashlib, json, math, os, sys
from fractions import Fraction as Q
from pathlib import Path
from functools import lru_cache

PREC = int(os.environ.get('NUMERIA_P245_FIXED_PRECISION', '256'))
if PREC not in (256, 512):
    raise ValueError('supported fixed-point precisions are 256 or 512 bits')
S = 1 << PREC
ROOT_PREC = 224
ROOT_SCALE = 1 << ROOT_PREC
TAYLOR_N = 100
LOG_N = 115 if PREC == 256 else 205
BERNOULLI_K = 30 if PREC == 256 else 45
SHIFT = 512
PAIRS = ((0, 0), (0, 1), (1, 1), (0, 16), (16, 16),
         (0, 376), (1, 376), (376, 376))
PANELS = (0, 503, 1007)
PILOT_HASHES = {
    'P236_even_0.json': '1b75938e12fca8d117443d347d714ee98b4b2c2dde540f50d5a18940f90cf5c2',
    'P236_even_503_900.json': 'b028262eaa3ca7255210573a50ccef9f39326889533e647bb4635e7d4d4436d7',
    'P236_even_1007_1152.json': '2d78d723aca037f50a4a8d845016ac43e179264b16261765423239eb3c05655a',
    'P236_odd_0.json': '59d3258a3d4130dd122972aec6b27fd404ce9007a036eb8375b1a9902f67e0c4',
    'P236_odd_503_900.json': 'fc897a58700b3f70ce59988f38efa97e41e5d592c0743bc5def483644700061b',
    'P236_odd_1007_1152.json': 'e06c138750b236958752c20b6fbed1d24ef3a91ea552e446801a27f987f5ba96',
}


def ceildiv(a: int, b: int) -> int:
    assert b > 0
    return -((-a)//b)


class I:
    __slots__ = ('l', 'u')
    def __init__(self, lo: int, hi: int | None = None):
        self.l = int(lo)
        self.u = self.l if hi is None else int(hi)
        if self.l > self.u:
            raise ValueError('reversed interval')

    @classmethod
    def rational(cls, q):
        q = Q(q)
        return cls((q.numerator * S)//q.denominator,
                   ceildiv(q.numerator * S, q.denominator))

    @classmethod
    def endpoints(cls, lo: Q, hi: Q):
        return cls(lo.numerator*S//lo.denominator,
                   ceildiv(hi.numerator*S,hi.denominator))

    def __add__(self, other):
        other=iv(other)
        return I(self.l+other.l, self.u+other.u)
    __radd__ = __add__

    def __neg__(self): return I(-self.u, -self.l)
    def __sub__(self,other): return self + (-iv(other))
    def __rsub__(self,other): return iv(other) - self

    def __mul__(self,other):
        other=iv(other)
        ps=(self.l*other.l,self.l*other.u,self.u*other.l,self.u*other.u)
        return I(min(ps)//S,ceildiv(max(ps),S))
    __rmul__=__mul__

    def rec(self):
        if self.l <= 0 <= self.u:
            raise ZeroDivisionError('interval contains zero')
        return I((S*S)//self.u,ceildiv(S*S,self.l))
    def __truediv__(self,other): return self * iv(other).rec()
    def __rtruediv__(self,other): return iv(other) / self

    def __pow__(self,n):
        if n<0: return (self**(-n)).rec()
        if n==0:return iv(1)
        acc=iv(1); b=self
        while n:
            if n&1:acc=acc*b
            n>>=1
            if n:b=b*b
        return acc

    def sqrt(self):
        if self.l<0: raise ValueError('sqrt negative')
        floor=math.isqrt(self.l*S)
        maxn=self.u*S; ceil=math.isqrt(maxn)
        if ceil*ceil<maxn:ceil+=1
        return I(floor,ceil)

    def pad(self,ticks=1):return I(self.l-ticks,self.u+ticks)
    def abs_upper(self):return max(abs(self.l),abs(self.u))
    def contain_fraction(self,q:Q):
        q=Q(q)
        return self.l*q.denominator<=q.numerator*S<=self.u*q.denominator
    def contain_interval(self,other):
        other=iv(other)
        return self.l<=other.l and other.u<=self.u
    def desc(self):return {'lower':str(self.l),'upper':str(self.u),'scale_power_of_two':PREC}
    def approx(self):return f'{(self.l+self.u)/(2*S):.17g}'
    def __repr__(self):return f'I({self.approx()}, rad<{(self.u-self.l)/(2*S):.3g})'


def iv(x):return x if isinstance(x,I) else I.rational(x)


def pi_interval():
    def atan_bracket(m,N=100):
        s=Q(0)
        for k in range(N):
            t=Q(1,(2*k+1)*pow(m,2*k+1))
            s += t if k%2==0 else -t
        next_t=Q(1,(2*N+1)*pow(m,2*N+1))
        return (s, s+next_t) if N%2==0 else (s-next_t,s)
    a,b=atan_bracket(5);c,d=atan_bracket(239)
    return I.endpoints(16*a-4*d,16*b-4*c)

PI=pi_interval()


def _series_odd_log(z):
    if z.abs_upper()*8>3*S:raise ArithmeticError('log series |z| >3/8')
    z2=z*z;term=z;val=iv(0)
    for k in range(LOG_N):
        val=val+term/Q(2*k+1)
        term=term*z2
    # A priori remainder for |z|<=1/3
    rem=Q(2,1)*Q(3,8)**(2*LOG_N+1)/Q(2*LOG_N+1)/(1-Q(9,64))
    if rem>=Q(1,S):raise ArithmeticError('log remainder exceeds one unit')
    return (2*val).pad(1)

LN2=_series_odd_log(iv(Q(1,3)))


def log_interval(x):
    x=iv(x)
    if x.l<=0:raise ArithmeticError('log nonpositive')
    # Exact binary range reduction without relying on floating point
    k=x.l.bit_length()-1-PREC
    v=x/(iv(2)**k) if k>=0 else x*(iv(2)**(-k))
    # v between 1 and 2 for non-point endpoints except on a binary boundary
    z=(v-1)/(v+1)
    if z.abs_upper()*8>3*S:
        # use another k to place v in [1/2,2]
        if z.l<0: k-=1;v=v*2
        else:k+=1;v=v/2
        z=(v-1)/(v+1)
    return _series_odd_log(z)+k*LN2


def _sin_cos_small(t):
    t=iv(t)
    if t.abs_upper()>4*S:raise ArithmeticError('Taylor input not reduced to [-4,4]')
    y2=t*t
    cosv=iv(1);sinv=t
    c=iv(1);s=t
    for k in range(1,TAYLOR_N+1):
        c=-c*y2/Q((2*k-1)*(2*k))
        s=-s*y2/Q((2*k)*(2*k+1))
        cosv=cosv+c
        sinv=sinv+s
    # Taylor with Lagrange derivative bound <=1, for all real arguments.
    rcos=Q(4)**(2*TAYLOR_N+1)/math.factorial(2*TAYLOR_N+1)
    rsin=Q(4)**(2*TAYLOR_N+2)/math.factorial(2*TAYLOR_N+2)
    if max(rcos,rsin)>=Q(1,S):raise ArithmeticError('trigonometric Taylor remainder too large')
    return sinv.pad(1),cosv.pad(1)


def sin_cos(t):
    t=iv(t)
    period=2*PI
    # A deterministic integer k is a range-reduction choice, not a proof decision
    # Integer midpoint quotient avoids float and the remainder is checked.
    k=(((t.l+t.u)//2)+((period.l+period.u)//4))//((period.l+period.u)//2)
    small=t - k*period
    if small.abs_upper()>4*S:
        raise ArithmeticError('trig reduction failed')
    return _sin_cos_small(small)


@lru_cache(1)
def bernoulli_even():
    max_order=2*BERNOULLI_K
    nums=[Q(1)]
    for n in range(1,max_order+1):
        bn=-sum(Q(math.comb(n+1,k))*nums[k] for k in range(n))/Q(n+1)
        nums.append(bn)
    return [nums[2*k] for k in range(1,BERNOULLI_K+1)]


def re_digamma(y):
    """Rigorous real part of psi(1/4+i*y), shifted to 512+1/4+i*y.

    Binet / Euler-Maclaurin terms through 2*K; the omitted Binet integral is
    bounded by 4*(2K+1)!/(6^(2K+2)*a^(2K)*(a^2-ymax^2)), since pi>3 and
    zeta(s)<2 for s>=2.  Here 0<=y<=252.
    """
    y=iv(y); a=iv(Q(4*SHIFT+1,4));a2=a*a; den=a2+y*y
    if y.abs_upper()>252*S:raise ValueError('digamma y >252')
    inv_re=a/den; inv_im=-y/den
    invsq_re=inv_re*inv_re-inv_im*inv_im
    invsq_im=2*inv_re*inv_im
    pw_re,pw_im=invsq_re,invsq_im
    v=log_interval(den)/2-inv_re/2
    for k,bk in enumerate(bernoulli_even(),start=1):
        v=v+(-bk/Q(2*k))*pw_re
        if k!=BERNOULLI_K:
            pw_re,pw_im=pw_re*invsq_re-pw_im*invsq_im,pw_re*invsq_im+pw_im*invsq_re
    k=BERNOULLI_K
    remainder=Q(4*math.factorial(2*k+1),6**(2*k+2))*Q(1,SHIFT**(2*k))*Q(1,Q(4*SHIFT+1,4)**2-Q(252)**2)
    if remainder>=Q(1,S):raise AssertionError('Binet remainder > 1 ulp')
    v=v.pad(1)
    # psi(z+n) = psi(z) + sum_{j=0}^{n-1} 1/(z+j)
    for j in range(SHIFT):
        aj=iv(Q(4*j+1,4))
        v=v-aj/(aj*aj+y*y)
    return v


def spherical_j_small(x):
    x=iv(x)
    if x.l<=0:raise ValueError('positive x required')
    s,c=sin_cos(x)
    j0=s/x
    j1=s/(x*x)-c/x
    j2=3*j1/x-j0
    j3=5*j2/x-j1
    return (j0,j1,j2,j3)


def double_factorial_odd(n):
    return math.prod(range(1,2*n+2,2))


def bessel_small_arg_series(x,n):
    """j_n on 0<x<1 via a signed power series with proven tail <1 ulp."""
    x=iv(x)
    if x.l<=0 or x.u>=S: raise ValueError('small-argument series only for 0<x<1')
    K=32 if PREC==256 else 64
    x2=x*x
    term=iv(1); sm=term
    for j in range(1,K):
        term=-term*x2/Q(2*j*(2*n+2*j+1))
        sm=sm+term
    upper=Q(1)
    for j in range(1,K+1):upper/=2*j*(2*n+2*j+1)
    # x<1 so the first omitted term is <= this exact rational.
    if upper>=Q(1,S):raise ArithmeticError('Bessel power-series tail not small')
    return (x**n)*sm.pad(1)/double_factorial_odd(n)


def bessel_forward(x,n,previous):
    """Forward 3-term recurrence for large x and moderate n."""
    arr=list(previous[:2])
    for k in range(1,n):arr.append((2*k+1)*arr[k]/x-arr[k-1])
    return arr[n]


def legendre_at(x,n):
    x=Q(x)
    p0,p1=Q(1),x
    if n==0:return p0
    if n==1:return p1
    for j in range(2,n+1):
        p0,p1=p1,((2*j-1)*x*p1-(j-1)*p0)/j
    return p1


def legendre_iv(x,n):
    p0,p1=iv(1),iv(x)
    if n==0:return p0
    if n==1:return p1
    for j in range(2,n+1):
        p0,p1=p1,((2*j-1)*x*p1-(j-1)*p0)/j
    return p1


def certified_nodes(root_file:Path|None=None):
    if root_file is None:
        # mpmath is a ROOT GUESSER ONLY, and its values are NOT trusted as bounds.
        import mpmath as mp
        mp.mp.dps=110
        guesses,_=mp.gauss_quadrature(60,'legendre')
        integers=[int(mp.floor(g*ROOT_SCALE)) for g in guesses]
    else:
        obj=json.loads(root_file.read_text())
        if (obj.get('schema')!='NUMERIA.ES8P245B.RATIONAL_GL60_ROOT_BRACKETS.1'
            or obj.get('root_degree')!=60 or obj.get('root_denominator_power_of_two')!=ROOT_PREC
            or len(obj.get('roots',[]))!=60):
            raise ValueError('invalid root-bracket file header')
        integers=[]
        for idx,item in enumerate(obj['roots']):
            if (item.get('index')!=idx
                or int(item['upper_numerator'])!=int(item['lower_numerator'])+1):
                raise ValueError('invalid root-bracket entry')
            integers.append(int(item['lower_numerator']))
    roots=[]
    for idx,k in enumerate(integers):
        lower=Q(k,ROOT_SCALE);upper=Q(k+1,ROOT_SCALE)
        sl=legendre_at(lower,60);su=legendre_at(upper,60)
        if not sl or not su or (sl>0)==(su>0):
            raise ArithmeticError(f'not isolated root {idx}')
        x=I.endpoints(lower,upper)
        p59=legendre_iv(x,59)
        # Exact at the root: w=2(1-x^2)/(60^2 P_59(x)^2)
        w=2*(1-x*x)/(3600*p59*p59)
        if w.l<=0 or x.l<=-S or x.u>=S:
            raise ArithmeticError(f'bad weight {idx}')
        roots.append((x,w,{'index':idx,'a':str(k),'b':str(k+1),'denominator_power_of_two':ROOT_PREC,
                            'P60_a_sign':1 if sl>0 else -1,'P60_b_sign':1 if su>0 else -1,
                            'weight_interval':w.desc()}))
    # 60 isolated roots in (-1,1) are automatically all roots, by degree.
    if len(roots)!=60 or any(roots[i][0].u>=roots[i+1][0].l for i in range(59)):
        raise ArithmeticError('Legendre root coverage/order invalid')
    if not sum((w for x,w,meta in roots),iv(0)).contain_fraction(Q(2)):
        raise ArithmeticError('GL60 weight-sum enclosure does not contain 2')
    return roots


def source_ball_to_fraction(b):
    (m,e),(r,k)=b
    def d(a,b):return Q(int(a)*(2**int(b))) if int(b)>=0 else Q(int(a),2**(-int(b)))
    return d(m,e)-d(r,k),d(m,e)+d(r,k)


def verify_against_pilot(panel,parity,value_map,source_dir:Path):
    sector='odd' if parity else 'even'
    src=source_dir/f'P236_{sector}_{panel}{"" if panel==0 else "_900" if panel==503 else "_1152"}.json'
    b=src.read_bytes()
    if hashlib.sha256(b).hexdigest()!=PILOT_HASHES[src.name]:
        raise ArithmeticError('pilot-source SHA-256 not pinned: '+src.name)
    data=json.loads(b)
    assert data['sector']==sector and data['panel_index']==panel
    samples={(i,j):ball for i,j,ball in data['gauss_interval_matrix_sample']}
    records=[]
    for i,j in PAIRS:
        lo,hi=source_ball_to_fraction(samples[i,j]); test=value_map[i,j]
        if not (test.contain_fraction(lo) and test.contain_fraction(hi)):

            raise ArithmeticError(f'P245B external pilot not contained: {sector} {panel} {i},{j}; {test}')
        records.append({'pair':[i,j],'degree_pair':[2*i+parity,2*j+parity],
                        'certified_interval':test.desc(),'encloses_arb_sample_interval':True})
    return {'source_file':src.name,'source_sha256':hashlib.sha256(b).hexdigest(),
            'sector':sector,'panel':panel,'comparisons':records}


def compute_all(source_dir:Path,root_file:Path|None=None):
    roots=certified_nodes(root_file)
    log_pi=log_interval(PI)
    logs=[LN2,log_interval(3),2*LN2,log_interval(5)]
    coeff=[2*LN2/iv(2).sqrt(),2*logs[1]/iv(3).sqrt(),LN2,2*logs[3]/iv(5).sqrt()]
    beta=log_interval(iv(504)/(2*PI))-iv(Q(1,504))-sum(coeff,iv(0))
    L=iv(Q(161,200))
    norms=[(2*L*(2*n+1)).sqrt() for n in (0,1,2,3,32,33)]
    report={'schema':'NUMERIA.ES8P245B.EXACT_FIXEDPOINT_GL60_SPOTCERT.1',
            'proof_arithmetic':f'Python arbitrary-size integers directed {PREC}-bit fixed-point; exact Fraction sign/root isolation',
            'root_input':'pinned integer root brackets (no mpmath at audit time)' if root_file is not None else 'mpmath guesses, rigorously isolated',
            'root_brackets_sha256':hashlib.sha256(root_file.read_bytes()).hexdigest() if root_file else None,
            'parameters':{'L':'161/200','T':504,'GL_nodes':60,'precision_bits':PREC,
                          'root_bracket_bits':ROOT_PREC,'digamma_shift':SHIFT,
                          'digamma_Binet_terms':BERNOULLI_K},
            'beta_interval':beta.desc(),
            'nodes':[tup[2] for tup in roots],
            'pilots':[], 'certified_selected_entries':0,
            'all_GL60_panels_certified':False,'Weil_Zhu_identity_certified':False,'RH_proved':False}
    for panel in PANELS:
        vals={p:{pair:iv(0) for pair in PAIRS} for p in (0,1)}
        x_max=Q(161*(panel+1),400)
        node_abs_sum=Q(0)
        high_bounds={n:Q(50)*x_max**n/double_factorial_odd(n) for n in (752,753)}
        for idx,(root,weight,_) in enumerate(roots):
            t=(2*panel+1+root)/4
            x=L*t
            dig=re_digamma(t/2)
            ps=dig-log_pi-sum((coef*sin_cos(t*ln)[1] for coef,ln in zip(coeff,logs)),iv(0))
            z=weight*(ps-beta)/(4*PI)
            node_abs_sum+=Q(z.abs_upper(),S)
            jb=spherical_j_small(x)
            for par in (0,1):
                n16=32+par
                j16=bessel_small_arg_series(x,n16) if panel==0 else bessel_forward(x,n16,jb)
                n376=752+par
                hu=high_bounds[n376]
                u={0:norms[par]*jb[par],
                   1:-norms[2+par]*jb[2+par],
                   16:norms[4+par]*j16,
                   376:I.endpoints(-hu,hu)}
                for i,j in PAIRS:
                    vals[par][i,j]=vals[par][i,j]+z*u[i]*u[j]
        for par in (0,1):
            rec=verify_against_pilot(panel,par,vals[par],source_dir)
            # Exact-rational high-degree magnitude bounds are FAR sharper
            # than the common zero-centered intervals for tiny modes.
            src_map=json.loads((source_dir/rec['source_file']).read_text())
            sampled={(i,j):b for i,j,b in src_map['gauss_interval_matrix_sample']}
            rec['high_degree_magnitude_bounds']=[]
            h=high_bounds[752+par]
            for i,j in PAIRS:
                if 376 not in (i,j):continue
                # |U_0|,|U_1|<4; high-mode |U_376|<50*xmax^n/(2n+1)!!.
                ui=h if i==376 else Q(4)
                uj=h if j==376 else Q(4)
                bound=node_abs_sum*ui*uj
                low,up=source_ball_to_fraction(sampled[i,j])
                if not (-bound<low<=up<bound):
                    raise ArithmeticError('high-mode source not inside exact magnitude bound')
                k=max(0,((bound.denominator.bit_length()-bound.numerator.bit_length())*30102999566)//100000000000-3)
                while bound<Q(1,10**(k+1)):k+=1
                if not bound<Q(1,10**k):raise ArithmeticError('bad high magnitude exponent')
                rec['high_degree_magnitude_bounds'].append({
                    'pair':[i,j], 'absolute_value_strictly_below_10_to_negative_power':k,
                    'verified_exact_rational':True})
            report['pilots'].append(rec)
            report['certified_selected_entries']+=len(PAIRS)
        print('P245-B rigorous GL60 panel',panel,'verified both parities',flush=True)
    if report['certified_selected_entries']!=48:raise AssertionError('incomplete selected matrix entries')
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--pilot-dir',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--roots',type=Path,default=None,help='verified rational root brackets; do not run a numerical root finder')
    a=p.parse_args()
    result=compute_all(a.pilot_dir,a.roots)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print('CERTIFIED',result['certified_selected_entries'],'GL60 entry inclusions independent of Arb',flush=True)

if __name__=='__main__':main()
