"""NUMERIA P82: independently check arithmetic stages of CS v2 Lemma 6.

Distinguishes exact proofs (Dirichlet coefficient convolution, minus-one
identity, splitting at 5X, exponent/log algebra) from imported analytic
lemmas (KLN mean-value theorem; Ramaré short-range b3=0.605; RZ S_eps).
Pure rational tests here do NOT prove those external theorems.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from math import gcd
from pathlib import Path

from NUMERIA_P81_TRANSFER.transfer import audit as audit_p81

SCHEMA="NUMERIA.P82.FX.MEAN_SQUARE_AUDIT.1"
HERE=Path(__file__).resolve().parent
ORIGINAL_ARTICLE="https://arxiv.org/html/2507.15184v2"
KLN_ARTICLE="https://www.cs.uleth.ca/~kadiri/articles/Explicit-zero-density-for-the-Riemann-zeta-function-JMAA-May2018.pdf"
SOURCE_SHA256={
 "passagem_completa_3_pdf":"a52ae9cac7777130c5383b4c70267663bbf0594c05798d6b3b001e32c06d90ab",
 "proposicao_auxiliar_b1_pdf":"acd532d7f48c828577f3d1080f19f6741ab23af6594412b34ac5dbe7ed246841",
 "proposicao_auxiliar_capa_html":"a81879bc525fb4dd6ea4ed2d930dfbc2d2b3d30e011a07e0f6a9d2ce9a096994",
}

def mu(n):
    if not isinstance(n,int) or n<1:
        raise ValueError("mu: positive integer required")
    if n==1:return 1
    m=n
    out=1
    p=2
    while p*p<=m:
        if m%p==0:
            m//=p
            if m%p==0:return 0
            out=-out
        while m%p==0:m//=p
        p+=1
    if m>1:out=-out
    return out

def divisor_raw(n,X):
    """a_X(n) = sum_{d|n,d<=X} mu(d), including a_X(1)=1."""
    return sum(mu(d) for d in range(1,min(n,X)+1) if n%d==0)

def lam(n,X):
    """Dirichlet coefficient of zeta(s)M_X(s)-1."""
    return divisor_raw(n,X)-(1 if n==1 else 0)

def truncated_diagonal(X,N,tau):
    """Finite coefficient convolution, without zeta evaluation.

    sum_{n<=N} a_X(n)^2 /n^tau =
      sum_{d,e<=X} mu(d)mu(e)/lcm(d,e)^tau *
             sum_{k<=N/lcm(d,e)} 1/k^tau.
    For N>=X, lhs = 1 + sum_{n<=N} lambda_X(n)^2/n^tau.
    All expressions evaluated as exact fractions for integer tau>=2.
    """
    if type(X)!=int or type(N)!=int or not 1<=X<=N:
        raise ValueError("require 1<=X<=N")
    if not isinstance(tau,int) or tau<2:
        raise ValueError("tau must be integer >=2")
    direct=sum((Q(divisor_raw(n,X)**2,n**tau) for n in range(1,N+1)),Q(0))
    mol=sum((Q(lam(n,X)**2,n**tau) for n in range(1,N+1)),Q(0))
    mob=[mu(d) for d in range(1,X+1)]
    harmonic=[Q(0)]*(N+1)
    for k in range(1,N+1):
        harmonic[k]=harmonic[k-1]+Q(1,k**tau)
    conv=Q(0)
    for d in range(1,X+1):
        if not mob[d-1]:
            continue
        for e in range(1,X+1):
            if not mob[e-1]:
                continue
            l=d*e//gcd(d,e)
            if l>N: continue
            conv+=Q(mob[d-1]*mob[e-1],l**tau)*harmonic[N//l]
    if direct!=conv or direct!=mol+1:
        raise ArithmeticError("Dirichlet convolution or -1 identity failed")
    if any(lam(n,X)!=0 for n in range(1,X+1)):
        raise ArithmeticError("lambda_X nonzero below X")
    return {"X":X,"N":N,"tau":tau,
            "raw_diagonal":str(direct),
            "lambda_square_sum":str(mol),
            "difference_exact":str(direct-mol),
            "double_sum_truncation":str(conv),
            "all_equalities_exact":True}

def rational_mv_weight_identity(X,N,q,T,P):
    """Imported KLN mean value bound has sum lambda² n^-q (T+P(n+1)).
    This *algebraic re-grouping* is exact; KLN's inequality is external.
    """
    T,P=Q(T),Q(P)
    if min(T,P)<0 or q<3:
        raise ValueError("invalid weight")
    a=sum((Q(lam(n,X)**2,n**(q-1)) for n in range(1,N+1)),Q(0))
    b=sum((Q(lam(n,X)**2,n**q) for n in range(1,N+1)),Q(0))
    direct=sum((Q(lam(n,X)**2,n**q)*(T+P*(n+1))
                for n in range(1,N+1)),Q(0))
    rebuilt=P*a+(T+P)*b
    if direct!=rebuilt:
        raise ArithmeticError("KLN weight splitting incorrect")
    return {"X":X,"N":N,"q":q,"direct":str(direct),
            "term_first_sum":str(P*a),
            "term_second_sum":str((T+P)*b),
            "reassembled":str(rebuilt)}

def finite_splitting(X,N,a):
    """Spot-test basic power monotonicities, exact for integer a>=1.
    No false extrapolation of small-X tests to X>=3e12 is allowed.
    """
    if type(a)!=int or a<1 or not 1<=X<=N:
        raise ValueError("invalid exponent or range")
    first=Q(0);second=Q(0);rhs_first=Q(0);rhs_second=Q(0)
    for n in range(X+1,N+1):
        val=lam(n,X)**2
        if n<=5*X:
            first+=Q(val,n**(2+2*a))
            rhs_first+=Q(val,X**(2*a)*n**2)
        else:
            second+=Q(val,n**(2+2*a))
            rhs_second+=Q(val,(5*X)**(1+a)*n**(1+a))
    if first>rhs_first or second>rhs_second:
        raise ArithmeticError("KLN n>X or n>5X inequality wrong")
    return {"X":X,"N":N,"a":a,
            "first_piece_satisfied":first<=rhs_first,
            "second_piece_satisfied":second<=rhs_second,
            "combined_satisfied":first+second<=rhs_first+rhs_second}

def main_audit():
    p81=audit_p81()
    if not p81["all_checks_pass"] or p81["total_branch_zone_checks"]!=96:
        raise ArithmeticError("P81 public S_epsilon transfer not reproduced")
    identities=[]
    for X,N in ((2,45),(3,48),(5,75),(8,96),(12,144)):
        for tau in (2,3):
            identities.append(truncated_diagonal(X,N,tau))
    mv=[]
    for X,N in ((2,45),(5,75),(12,144)):
        for q in (3,4):
            mv.append(rational_mv_weight_identity(X,N,q,Q(13,7),Q(22,9)))
    splits=[]
    for X,N in ((2,45),(5,75),(12,144)):
        for a in (1,2):
            splits.append(finite_splitting(X,N,a))
    if not(len(identities)==10 and len(mv)==6 and len(splits)==6):
        raise ArithmeticError("incomplete P82 exact tests")
    return {
       "schema":SCHEMA,
       "source": {"CS_v2":ORIGINAL_ARTICLE,"KLN_2018":KLN_ARTICLE,
                  "RZ_v3":"https://arxiv.org/html/2603.25961v3"},
       "uploaded_source_sha256":SOURCE_SHA256,
       "coefficients":{
          "f_X":"zeta(s)*sum_{d<=X}mu(d)/d^s-1",
          "lambda_X_n":"sum_{d|n,d<=X}mu(d) minus indicator(n=1)",
          "lambda_vanishes_n_le_X":True,
          "exact_infinite_sum_identity":"sum lambda_X(n)^2/n^tau = zeta(tau)*S_(tau-1)(X) - 1 (tau>1)",
          "CS_v2_uses_safe_inequality":"sum lambda_X^2/n^tau <= zeta(tau)*S_(tau-1)(X)"
       },
       "analytic_integration":{
          "F_X":"integral_0^T |f_X(sigma+it)|^2 dt",
          "KLN_3_14_imported":True,
          "KLN_m0_squared":"1+(2/3)*sqrt(6/5)",
          "inequality":"F_X <= pi*m0*sum lambda^2/n^(2sigma-1) + (T+pi*m0)*sum lambda^2/n^(2sigma)",
          "proof":"substitute u_n=lambda_X(n)*n^-sigma in KLN 3.14, regroup n+1 exactly",
          "source_conditions":"real coefficients; sigma>1; T>=0; absolute weighted square summability",
       },
       "tail_splitting":{
          "short_range":"X<n<=5X: n^(-2-2a)<=X^(-2a)*n^(-2)",
          "long_range":"n>5X: n^(-2-2a)<=(5X)^(-1-a)*n^(-1-a)",
          "b3_short_range_imported":"sum_(X<n<=5X) lambda^2/n^2 <= 0.605/X (KLN lemma 3.4 eq 3.8)",
          "RZ_source_imported":"sum lambda^2/n^(1+a) <= b4/a from conditional P81, X>=3e12",
          "a":"delta/log X, with 0<delta<=0.303; epsilon_2=2a and epsilon_1=a",
          "combined_exact_form":"sum lambda^2/n^(2+2a) <= (0.605*exp(-2delta)+b4*logX*exp(-delta)*5^(-a)/(5delta))/X",
       },
       "final_constants":{
          "C3":"pi*m0*b4*exp(2delta*gamma/logX0)/(2delta)",
          "C4":"b4/(5delta*exp(delta))+0.605*exp(-2delta)/logX0",
          "F_X_bound":"F_X(1+delta/logX,0,T,1/2) <= logX*(C3+C4*(T+pi*m0)/X)",
          "scope":"X>=X0=3e12; T>=0; 0<delta<=0.303; external KLN b3 and MV + RZ S_eps",
          "for_CS_corollary_1":"X=h*T, h=1, T>=3e12 (CS proof of Corollary 1)",
          "not_full_CS_lemma6_scope":"original X0>=3e10; P82 public RZ branch only X0>=3e12",
       },
       "finite_rational_mobius_examples":identities,
       "finite_rational_mean_weight_examples":mv,
       "finite_rational_tail_examples":splits,
       "status":{
          "dirichlet_convolution_proved_algebraically":True,
          "CS_FX_inequality_derived_conditionally":True,
          "KLN_mean_value_theorem_independently_proved":False,
          "KLN_b3_short_range_independently_proved":False,
          "RZ_S_epsilon_theorems_independently_proved":False,
          "Riemann_hypothesis_proved":False,
       }
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    data=main_audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({
        "exact_diagonal_tests":len(data["finite_rational_mobius_examples"]),
        "exact_mv_weight_tests":len(data["finite_rational_mean_weight_examples"]),
        "exact_tail_tests":len(data["finite_rational_tail_examples"]),
        "conclusion":"conditional F_X verified, three external analytic gates remain open"
    },indent=2))

if __name__=="__main__":
    main()
