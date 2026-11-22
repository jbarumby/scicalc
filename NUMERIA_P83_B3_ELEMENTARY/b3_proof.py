"""NUMERIA P83: elementary independent analytic proof of the SHORT-RANGE b3=0.605.

The conclusion is valid for EVERY REAL X>=10**9:
  sum_{X<n<5X} (sum_{d|n,d<=X} mu(d))^2/n^2 <= 0.605/X.

No use of Ramaré Lemma 3.1's sharp squarefree error bound, nor
Ramaré Lemma 5.6 itself; only elementary Möbius inversion, Euler product
zeta(2)=pi²/6, Archimedes/Machin pi lower bound and partial summation.
The exhaustive finite checks are regression tests, NOT proof for X>=10**9.
Ramaré's separate Hilbert (Montgomery–Vaughan) mean-value inequality
still remains an imported external analytic theorem.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path

SCHEMA="NUMERIA.P83.INDEPENDENT_B3.SHORT_RANGE.1"
X0=10**9
B3=Q(121,200)
PI_LOWER=Q(157,50)
SQFREE_ERROR_CONSTANT=Q(201,100)
SQRT_X_LOWER=30000
S1_ERROR_FACTOR=Q(65,24)
SOURCE_R="https://ramare-olivier.github.io/Maths/DensityEstimate-34.pdf"
SOURCE_KLN="https://www.cs.uleth.ca/~kadiri/articles/Explicit-zero-density-for-the-Riemann-zeta-function-JMAA-May2018.pdf"

def mu_sieve(limit):
    if not isinstance(limit,int) or limit<1:
        raise ValueError("limit must be positive integer")
    lp=[0]*(limit+1)
    primes=[]
    mob=[0]*(limit+1)
    mob[1]=1
    for i in range(2,limit+1):
        if lp[i]==0:
            lp[i]=i;primes.append(i);mob[i]=-1
        for p in primes:
            if p>lp[i] or p*i>limit:break
            lp[p*i]=p
            if p==lp[i]:
                mob[p*i]=0
                break
            mob[p*i]=-mob[i]
    return mob

def coefficient_sieve(X):
    """For integer X, a_X(n)=sum_{d|n,d<=X}mu(d); all n<5X."""
    if not isinstance(X,int) or X<1:raise ValueError("positive integer X")
    M=5*X
    mob=mu_sieve(X)
    a=[0]*(M+1)
    for d in range(1,X+1):
        if mob[d]:
            for n in range(d,M+1,d):
                a[n]+=mob[d]
    return a

def exact_finite_short_sum(X):
    """Diagnostic only; rational values, X small."""
    a=coefficient_sieve(X)
    s1=sum((Q(a[n]*a[n],n*n) for n in range(X+1,2*X+1)),Q())
    s2=sum((Q(a[n]*a[n],n*n) for n in range(2*X+1,5*X)),Q())
    if any(abs(a[n])>1 for n in range(2*X+1,5*X)):
        raise ArithmeticError("false coefficient absolute bound in (2X,5X)")
    return {
        "X":X,
        "S1":str(s1),
        "S2":str(s2),
        "total":str(s1+s2),
        "normalized_total":str(X*(s1+s2)),
        "at_most_one_for_2X_to_5X":True,
    }

def proof():
    """Compute every numerical constant in the general proof as exact rational.

    Euler product and the partial summation integral are proved in README.
    pi bound: Machin identity pi=16 atan(1/5)-4 atan(1/239);
    atan(x) alternating bounds give pi > 16(1/5-1/(3*5^3))-4/239>157/50.
    """
    pi_from_machin_lower=16*(Q(1,5)-Q(1,375))-Q(4,239)
    if not pi_from_machin_lower>PI_LOWER:
        raise ArithmeticError("pi lower bound not verified")

    # sqrt(X0)>30000; from m=floor(sqrt x) >=sqrt x-1:
    # |Q(x)-6x/pi²| <= m+x/m <= sqrt x*(2+1/(sqrt x-1)).
    if SQRT_X_LOWER*SQRT_X_LOWER>=X0:
        raise ArithmeticError("incorrect sqrt X0 lower")
    if not Q(2)+Q(1,SQRT_X_LOWER-1) < SQFREE_ERROR_CONSTANT:
        raise ArithmeticError("squarefree elementary discrepancy bound invalid")

    # Partial summation: X*S1 <= 3/pi²+
    # C/sqrt(X)*[1+sqrt(2)/4+(4/3)*(1-2^-3/2)]
    # sqrt 2 <3/2 and (1-2^-3/2)<1.
    if not (Q(3,2)**2 > 2 and
            Q(1)+Q(3,8)+Q(4,3)==S1_ERROR_FACTOR):
        raise ArithmeticError("S1 integral error factor incorrect")

    principal=Q(3)/PI_LOWER**2
    s1_error=SQFREE_ERROR_CONSTANT*S1_ERROR_FACTOR/SQRT_X_LOWER
    s2_principal=Q(3,10)
    s2_edge=Q(1,4*X0)
    total_upper=principal+s1_error+s2_principal+s2_edge
    slack=B3-total_upper
    if not slack>0:
        raise ArithmeticError("independent b3 bound failed")
    return {
        "schema":SCHEMA,
        "proof_scope":"ALL REAL X>=1e9",
        "theorem":"sum_(X<n<5X) a_X(n)^2/n^2 <= (121/200)/X",
        "a_X_n":"sum_(d|n,d<=X) mu(d); equals lambda_X(n) for n>1",
        "elementary_inputs":[
            "mu(n)^2=sum_(d^2|n)mu(d)",
            "sum_(n>=1)mu(n)/n^2=1/zeta(2)=6/pi^2",
            "|Q(x)-(6/pi^2)x|<=floor(sqrt(x))+x/floor(sqrt(x))",
            "partial summation for S1 on X<n<=2X",
            "|a_X(n)|<=1 for 2X<n<5X by mu-divisor case distinctions",
            "S2<= 3/(10X)+1/(4X^2) for decreasing n^-2",
            "pi>157/50 via Machin identity and alternating arctan",
        ],
        "X0":X0,
        "pi_machin_rational_lower":str(pi_from_machin_lower),
        "pi_conservative_lower":str(PI_LOWER),
        "squarefree_error_constant":str(SQFREE_ERROR_CONSTANT),
        "squarefree_error_domain":"x>=1e9",
        "sqrt_X0_lower":SQRT_X_LOWER,
        "S1_error_factor_rational_upper":str(S1_ERROR_FACTOR),
        "normalized_S1_majorant":str(principal+s1_error),
        "normalized_S2_majorant":str(s2_principal+s2_edge),
        "normalized_total_majorant":str(total_upper),
        "b3":str(B3),
        "positive_rational_slack":str(slack),
        "b3_inequality_proved_from_elementary_standard_identities":True,
        "ramare_short_range_lemma_required_as_external":False,
        "mv_hilbert_inequality_independently_proved":False,
        "RZ_S_epsilon_bounds_independently_proved":False,
        "RH_globally_proved":False,
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--diagnose-X",type=int,default=0)
    args=p.parse_args()
    doc=proof()
    if args.diagnose_X:
        if not 1<=args.diagnose_X<=200:
            raise ValueError("diagnose-X must be 1..200; finite only")
        doc["finite_diagnostic"]=exact_finite_short_sum(args.diagnose_X)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(doc,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({k:doc[k] for k in
        ("proof_scope","normalized_total_majorant","b3","positive_rational_slack",
         "ramare_short_range_lemma_required_as_external")},indent=2))

if __name__=="__main__":
    main()
