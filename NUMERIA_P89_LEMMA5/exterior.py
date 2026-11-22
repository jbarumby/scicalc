"""NUMERIA P89: an independent all-spacing theorem for EXTERIOR Lemma6 tails.

For sorted distinct reals x and a<b, h=x_b-x_a>0, delta_i nearest gaps,
 L = sum_{k<a} delta_k / ((x_a-x_k)^2 (x_b-x_k)^2)
 R = sum_{k>b} delta_k / ((x_k-x_a)^2 (x_k-x_b)^2).
We prove for ALL finite REAL configurations (README):
 L <= delta_a sum_{j>=1} 1/[(j delta_a)^2 (h+j delta_a)^2],
 R <= delta_b sum_{j>=1} 1/[(j delta_b)^2 (h+j delta_b)^2].
The proof uses the already-reconstructed Shan convex one-sided theorem in
P87, with f(t)=1/(t²(t+h)²) and monotonicity of g*f(j*g).
No Preissmann Lemma5 or Lemma6 is used as an input for THESE exterior tails.

The central interval a<k<b is still open at the SHARP Lemma6 bound.
The exact rational K-term upper bounds below are diagnostics, not a substitute
for the analytical all-real-node Shan argument.
"""
import argparse
import json
import hashlib
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
from NUMERIA_P86_PREISSMANN.spectral_rational import PAIRS,rational_nodes,gaps
from NUMERIA_P88_LEMMA6.compensation import machin_pi_lower

SCHEMA="NUMERIA.P89.PREISSMANN.LEMMA5.EXTERIOR.1"
K=48

def kernel(t,h):
    t,h=F(t),F(h)
    if t<=0 or h<=0: raise ValueError("positive t and h required")
    return 1/(t*t*(t+h)**2)

def majorant(delta,h,K=K):
    """Fully rational finite majorant for delta*sum_{j>=1} kernel(j delta,h).

    Sum j=1..K is exact; for j>K, rectangle/integral comparison
       delta sum_{j>K} f(j delta) <= integral_{K delta}^infty f(t)dt
                                     <= 1/[3 (K delta)^3].
    """
    delta,h=F(delta),F(h)
    if delta<=0 or h<=0 or type(K)!=int or K<1:
        raise ValueError("positive delta,h and integer K required")
    s=delta*sum((kernel(j*delta,h) for j in range(1,K+1)),F())
    tail=1/(3*(K*delta)**3)
    return s,tail,s+tail

def exact_parts(values,a,b):
    x=tuple(sorted(rational_nodes(values)))
    if type(a)!=int or type(b)!=int or not 0<=a<b<len(x):
        raise ValueError("valid ordered endpoint indices required")
    delta=gaps(x)
    h=x[b]-x[a]
    left=sum((delta[k]*kernel(x[a]-x[k],h) for k in range(a)),F())
    right=sum((delta[k]*kernel(x[k]-x[b],h) for k in range(b+1,len(x))),F())
    middle=sum((delta[k]/((x[k]-x[a])**2*(x[k]-x[b])**2)
                for k in range(a+1,b)),F())
    direct=sum((delta[k]/((x[k]-x[a])**2*(x[k]-x[b])**2)
                for k in range(len(x)) if k not in (a,b)),F())
    if left+right+middle!=direct:
        raise ArithmeticError("three-region decomposition is wrong")
    _, lt, lu=majorant(delta[a],h)
    _, rt, ru=majorant(delta[b],h)
    if not left<=lu or not right<=ru:
        raise ArithmeticError("exterior Shan envelope regression failure")
    return {
       "a_index":a+1,"b_index":b+1,"h":str(h),
       "delta_a":str(delta[a]),"delta_b":str(delta[b]),
       "left_exact":str(left),"middle_exact":str(middle),
       "right_exact":str(right),"total_exact":str(direct),
       "left_bound_rational":str(lu),"right_bound_rational":str(ru),
       "left_tail_K_rational":str(lt),"right_tail_K_rational":str(rt),
       "left_gap_rational":str(lu-left),
       "right_gap_rational":str(ru-right),
       "bounds_verified_exact":True,
    }

def source_errata():
    """Exact identities that correct the source P87 note and its B5 gate."""
    # For t>0:
    # psi_1(t) = integral s exp(-t s)/(1-exp(-s)) ds.
    # u=ts => t psi_1(t) = integral u exp(-u)/
    #                            [t(1-exp(-u/t))] du.
    # For fixed u, derivative of denominator is
    # d/dt[t(1-exp(-u/t))] = 1-(1+u/t)exp(-u/t)>0.
    #
    # The original note erroneously omits 1/t AND reverses the inequality
    # for 1-exp(-u/t). Both errors are substantive.
    #
    # Its Python B5 uses the UPPER bound of pi for Q<pi²/3:
    # this is directionally invalid. A LOWER pi bound is mandatory.
    pi_lower=machin_pi_lower()
    lower=pi_lower*pi_lower/3
    # Exact rational demonstration of upper-bound logical fallacy:
    # Q=22/7^2/3? just pick Q satisfying actual bound< Q<upper bound;
    # no external theorem is claimed by this simple diagnostic.
    return {
      "schema":"NUMERIA.P89.P87.NOTE.ERRATA.1",
      "uploaded_sha256":{
        "p87_teorema_preissmann.pdf":"616de178ba794135deafe3dccd041c50f60f6553a370c6b0bd93e99362096b05",
        "p87_teorema_preissmann.tex":"d16771dab750ef4f00e00bf841ddcf5f179e412ef6392914ea819a0e7b9d706d",
        "p87_verificacao.py":"5089f7016827dffb5a6202cd783ca70ab099ba5a87db37e954e2192a1b8ff14a",
        "p87_ci.yml":"d6fccbdfaee9e8bf1bd5116f69d38ebba9e31049a87a00aac52b51360ed2be2a",
        "p87_capa.html":"a1d420bdf90d56009d9fdb2574b1344b50262ddd1a6041ed8d5dee3153558135"
      },
      "laplace_original_missing_factor_one_over_t":True,
      "laplace_original_reversed_exponential_inequality":True,
      "correct_laplace":"t*psi1(t)=integral_0^inf u*exp(-u)/[t*(1-exp(-u/t))]du",
      "correct_denominator_derivative":"1-(1+u/t)*exp(-u/t)>0",
      "B5_original_used_upper_pi_bound_not_sufficient":True,
      "B5_correct_use_lower_pi_bound":str(pi_lower),
      "B5_correct_Q_bound_lower":str(lower),
      "B3_B4_B6_B7_floating_grid_are_diagnostic_only":True,
      "full_24_gates_are_not_24_exact_proofs":True,
      "energy_lemma_proof_repaired_in_errata_md":True,
      "independent_Preissmann_Lemma6_proved":False,
      "RH_proved":False,
    }

def audit():
    groups={
       "uniform_seventeen":[str(i) for i in range(17)],
       "uniform_thirtyone":[str(i) for i in range(31)],
       "clusters":["0","1/10000","1/5000","1/4000","3/1000","1","7","8","80"],
       "skewed":["0","1/4","1/2","3","40","50","101","110","700"],
       "geometric":["0","1","3","7","15","31","63","127","255"]
    }
    outcomes=[]
    for name,points in groups.items():
        n=len(points)
        selected=sorted({
            (0,1),(0,n-1),(n//3,2*n//3),(n//2,n-1),
            (1,n-2),(n//3,n//3+1),
        })
        records=[exact_parts(points,a,b) for a,b in selected]
        outcomes.append({"name":name,"N":n,"pairs_checked":len(records),"records":records})
    return {
       "schema":SCHEMA,
       "source":{"P87_Shans_theorem":"NUMERIA_P87_PREISSMANN/README.md",
                 "Preissmann_Lemma5":"https://doi.org/10.5169/seals-53823",
                 "Preissmann_Lemma6":"NUMERIA_P88_LEMMA6/README.md"},
       "theorem_exterior_all_real_spacings":
          "L<=delta_a*sum_{j>=1}[(j delta_a)^2(h+j delta_a)^2]^-1, similarly R with delta_b",
       "general_proof":"Shan one-sided positive convex decreasing f(t)=1/[t²(t+h)²]; monotone g*f(jg) for g>0",
       "three_regions_split_exact":True,
       "exterior_bound_proved_for_all_finite_real_nodes":True,
       "interior_middle_sharp_Lemma6_proved":False,
       "full_Lemma5_variational_proof_reproduced":False,
       "full_Lemma6_preissmann_reproved":False,
       "RZ_double_sum_source_extended_below_3e12":False,
       "RH_proved":False,
       "diagnostic_K":K,
       "finite_sets":outcomes,
       "finite_checked_pairs":sum(q["pairs_checked"] for q in outcomes),
       "p87_errata":source_errata(),
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    x=p.parse_args()
    d=audit()
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"finite_checked_pairs":d["finite_checked_pairs"],
        "exterior_all_real_proven":True,
        "full_Lemma6_proven":False},indent=2))

if __name__=="__main__":main()
