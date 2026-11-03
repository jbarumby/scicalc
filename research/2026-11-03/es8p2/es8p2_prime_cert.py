"""NUMERIA ES8-P2: rigorous rational intervals for the V17 finite prime-2 matrix.

This module DOES NOT compute the full Weil matrix and DOES NOT prove RH.
Certificate targets: exact real half-line symmetrization, exact parity zeros,
finite prime overlap kernel (only n=2 for L=9/20), its pole-annihilating
projection for N=6,8,10. Rational evaluations of pi/log2/sqrt2/trig are
outward, including explicit Taylor remainder and argument-Lipschitz bounds.

Mathematical normalization ASSUMED from the already-audited V17 model:
 Prime_ij = -2*log(2)/sqrt(2) * (Cij(log2)+Cji(log2))/2.
This is not an independent rederivation of Guinand-Weil from first principles.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as Q
from math import factorial, isqrt
from functools import lru_cache
import argparse, json, platform
from pathlib import Path

L=Q(9,20)
SCALE=10**40
M=80

def down(x:Q)->Q:
    x=Q(x);return Q(x.numerator*SCALE//x.denominator,SCALE)
def up(x:Q)->Q:
    x=Q(x);return Q(-((-x.numerator*SCALE)//x.denominator),SCALE)

@dataclass(frozen=True)
class Iv:
    lo: Q
    hi: Q
    def __post_init__(self):
        if self.lo>self.hi:raise ValueError('invalid interval')
    @staticmethod
    def point(v):
        v=Q(v);return Iv(v,v)
    @staticmethod
    def rounded(lo,hi):return Iv(down(lo),up(hi))
    def __add__(self,other):
        other=I(other);return Iv.rounded(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __neg__(self):return Iv(-self.hi,-self.lo)
    def __sub__(self,other):return self+-I(other)
    def __rsub__(self,other):return I(other)+-self
    def __mul__(self,other):
        other=I(other); a=[x*y for x in (self.lo,self.hi) for y in (other.lo,other.hi)]
        return Iv.rounded(min(a),max(a))
    __rmul__=__mul__
    def __truediv__(self,other):
        other=I(other)
        if other.lo<=0<=other.hi:raise ZeroDivisionError('denominator interval crosses zero')
        return self*Iv.rounded(1/other.hi,1/other.lo)
    def __rtruediv__(self,other):return I(other)/self
    def contains(self,x):return self.lo<=Q(x)<=self.hi
    def width(self):return self.hi-self.lo
    def as_json(self):return {'lo':str(self.lo),'hi':str(self.hi),
        'decimal_outward':[format_bound(self.lo,True),format_bound(self.hi,False)]}

def I(v):return v if isinstance(v,Iv) else Iv.point(v)

def format_bound(x:Q,lower:bool,digits:int=18)->str:
    s=10**digits
    q=(x.numerator*s//x.denominator) if lower else -((-x.numerator*s)//x.denominator)
    sign='-' if q<0 else ''
    q=abs(q);return f'{sign}{q//s}.{q%s:0{digits}d}'

def atan_bounds(q:Q,terms=45)->Iv:
    assert 0<q<1 and terms>1
    S=sum(((-1)**j*q**(2*j+1)/Q(2*j+1) for j in range(terms)),Q())
    R=q**(2*terms+1)/Q(2*terms+1)
    return Iv.rounded(S,S+R) if terms%2==0 else Iv.rounded(S-R,S)

def atanh_bounds(q:Q,terms=90)->Iv:
    assert 0<q<1 and terms>1
    S=sum((2*q**(2*j+1)/Q(2*j+1) for j in range(terms)),Q())
    R=2*q**(2*terms+1)/(Q(2*terms+1)*(1-q*q))
    return Iv.rounded(S,S+R)

@lru_cache(maxsize=1)
def constants():
    pi=16*atan_bounds(Q(1,5))-4*atan_bounds(Q(1,239))
    lg2=atanh_bounds(Q(1,3));lg3=atanh_bounds(Q(1,2))
    assert lg2.hi<2*L<lg3.lo
    scale_sqrt=isqrt(2*SCALE*SCALE)
    rt2=Iv(Q(scale_sqrt,SCALE),Q(scale_sqrt+1,SCALE))
    assert rt2.lo**2<=2<=rt2.hi**2
    assert pi.lo>3 and pi.hi<4
    return {'pi':pi,'ln2':lg2,'ln3':lg3,'sqrt2':rt2}


def sincos_rational_at(t:Q):
    """Taylor at exact rational t; Lagrange remainder with |f^(n)|<=1."""
    t=Q(t)
    sinpoly=Q(0);cospoly=Q(0)
    termsin=t;termcos=Q(1)
    for j in range(M+1):
        sinpoly+=termsin;cospoly+=termcos
        termsin=-termsin*t*t/Q((2*j+2)*(2*j+3))
        termcos=-termcos*t*t/Q((2*j+1)*(2*j+2))
    sinrem=abs(t)**(2*M+3)/Q(factorial(2*M+3))
    cosrem=abs(t)**(2*M+2)/Q(factorial(2*M+2))
    return (Iv.rounded(sinpoly-sinrem,sinpoly+sinrem),
            Iv.rounded(cospoly-cosrem,cospoly+cosrem))

@lru_cache(maxsize=10)
def trig_mode(k:int):
    if not 1<=k<=10:raise ValueError('ES8-P2 allows modes k=1..10')
    c=constants()
    arg=I(k)*c['pi']*c['ln2']/I(2*L)
    # Use a short fixed-denominator rational centre; enclosure widened via
    # sin/cos global 1-Lipschitz bound, so no assumption about endpoints.
    denom=10**24
    n=(((arg.lo+arg.hi)/2).numerator*denom//((arg.lo+arg.hi)/2).denominator)
    t=Q(n,denom)
    dist=max(abs(arg.lo-t),abs(arg.hi-t))
    s,co=sincos_rational_at(t)
    return Iv.rounded(s.lo-dist,s.hi+dist),Iv.rounded(co.lo-dist,co.hi+dist)

@lru_cache(maxsize=None)
def overlap_sym(i:int,j:int)->Iv:
    """C_sym_ij(log 2), exact parity zeros and analytic convolution formulas.

    If i==j: (1-u/(2L))*cos(w_i*u) + sin(w_i*u)/(i*pi).
    If i!=j same parity: 2/(pi*(i²-j²))*(i*sin(w_j*u)-j*sin(w_i*u)).
    If opposite parity: 0 exactly by reflection.
    """
    if not 1<=i<=10 or not 1<=j<=10:raise ValueError('modes 1..10')
    if i>j:return overlap_sym(j,i)
    if (i-j)%2:return Iv.point(0)
    c=constants()
    si,ci=trig_mode(i)
    if i==j:return (I(1)-c['ln2']/I(2*L))*ci+si/(I(i)*c['pi'])
    sj,_=trig_mode(j)
    return (2/(I(i*i-j*j)*c['pi']))*(I(i)*sj-I(j)*si)

@lru_cache(maxsize=None)
def raw_prime(i:int,j:int)->Iv:
    if i>j:return raw_prime(j,i)
    c=constants()
    return -c['sqrt2']*c['ln2']*overlap_sym(i,j)

def projection_column(k:int,N:int)->dict[int,Iv]:
    if not 3<=k<=N<=10:raise ValueError('projection modes 3..N<=10')
    c=constants()
    p=1 if k%2 else 2
    x=c['pi']*c['pi']
    coeff=-I(Q(k,p))*(I(p*p)*x+I(L*L))/(I(k*k)*x+I(L*L))
    return {p:coeff,k:Iv.point(1)}

def projected_prime(i:int,j:int,N:int)->Iv:
    """Entry of B^T Prime B in the poles-annihilating real subspace."""
    if not 3<=i<=N or not 3<=j<=N:raise ValueError('bad projected modes')
    if i>j:return projected_prime(j,i,N)
    if (i-j)%2:return Iv.point(0)
    a=projection_column(i,N);b=projection_column(j,N)
    val=Iv.point(0)
    for r,ai in a.items():
        for s,bj in b.items():val+=ai*raw_prime(r,s)*bj
    return val

def audit(N_values=(6,8,10)):
    c=constants()
    outputs={}
    for N in N_values:
        mat=[[projected_prime(i,j,N) for j in range(3,N+1)] for i in range(3,N+1)]
        assert all(mat[i][j]==mat[j][i] for i in range(N-2) for j in range(N-2))
        assert all(mat[i][j].lo==0==mat[i][j].hi for i in range(N-2) for j in range(N-2) if (i-j)%2)
        outputs[str(N)]={'dimension':N-2,
          'projected_prime_intervals':[[v.as_json()['decimal_outward'] for v in row] for row in mat],
          'projected_prime_intervals_exact_scaled':[[[str(int(v.lo*SCALE)),str(int(v.hi*SCALE))] for v in row] for row in mat],
          'sample_entries':{f'{i},{j}':mat[i-3][j-3].as_json() for i,j in [(3,3),(3,4),(3,5),(4,4)] if i<=N and j<=N},
          'projected_cross_parity_exact_zeros':sum((i-j)%2==1 for i in range(3,N+1) for j in range(3,N+1)),
          'max_width':format_bound(max(v.width() for row in mat for v in row),False,24)}
    return {'schema':'NUMERIA.ES8P2.PRIME_PROJECTION_INTERVAL.1',
      'scientific_status':'rigorous rational finite prime term enclosures for V17 convention; analytic real-part parity gate; NOT full Weil or RH proof',
      'constants':{k:v.as_json() for k,v in c.items()},
      'source_convention':'V17 f_k, FT forward +irx; arch=(1/pi)int_0^inf Omega Re(F_i conjugate(F_j)); prime=-2ln2/sqrt2 * symmetric convolution',
      'mathematical_identity':'real integral (-inf,+inf) of even Omega(r) F_i(r)conj(F_j(r)) equals twice halfline Omega(r) Re(F_i(r)conj(F_j(r))) by F_k(-r)=conj(F_k(r)); corrected real halfline integrand',
      'condition':'finite support [-L,L], 2L=9/10, only prime power 2; pole projection coefficients are dependent on pi^2',
      'finite_prime_matrices':outputs,
      'limitations':['not global Guinand-Weil re-derivation','not interval-certified archimedean entries','not independent validation of original E7/V16 full G+P source','not a positivity test for full matrix','no RH proof'],
      'environment':{'python':platform.python_version(),'stdlib_exact':'fractions.Fraction for rational intervals','Taylor_M':M,'scale':str(SCALE)}}

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='RESULTADOS_ES8P2.json')
    a=p.parse_args();res=audit()
    dest=Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(res,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('ES8-P2 PASS rational prime projections N=6,8,10, cross parity zeros; FULL_WEIL=NOT_ESTABLISHED')
    for N in [6,8,10]:print('N',N,'max width',res['finite_prime_matrices'][str(N)]['max_width'])

if __name__=='__main__':main()
