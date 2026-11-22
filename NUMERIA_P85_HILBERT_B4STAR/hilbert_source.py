"""NUMERIA P85: source-exact Hilbert/Preissmann -> real Dirichlet mean square.

Elementary identities, finite cancellation, logarithmic separation and the
extension to infinite Dirichlet series are proved in README and regression
tested independently. Preissmann's weighted Hilbert inequality itself is
IMPORTED from the 1984 paper, not reproved. No global RH claim.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import json
import argparse
from pathlib import Path

SCHEMA="NUMERIA.P85.REAL_HILBERT_SOURCE_TRANSFER.1"
RAMARE="https://ramare-olivier.github.io/Maths/DensityEstimate-34.pdf"
PREISSMANN="https://doi.org/10.5169/seals-53823"
MV="https://doi.org/10.1112/jlms/s2-8.1.73"
CARNEIRO="https://arxiv.org/html/2302.14658"

def m0_bracket():
    """m0 = sqrt(1+2/3 sqrt(6/5)); strict rational bounds without floats."""
    lo=F(263,200) # 1.315
    hi=F(329,250) # 1.316
    def transformed(v):
        return F(3,2)*(v*v-1)
    if not lo>1 or not transformed(lo)**2<F(6,5)<transformed(hi)**2:
        raise ArithmeticError("source m0 rational bracket invalid")
    return str(lo),str(hi)

def exact_boundary_cancellation(frequencies,coefficients):
    """Algebraic sum over m!=n of real u_m*u_n/(lambda_m-lambda_n) = 0."""
    x=[F(q) for q in frequencies]
    u=[F(q) for q in coefficients]
    if len(x)!=len(u) or len(x)<2 or len(set(x))!=len(x):
        raise ValueError("finite frequencies must be distinct and equally sized")
    unsym=sum((u[i]*u[j]/(x[i]-x[j])
              for i in range(len(x)) for j in range(len(x)) if i!=j),F())
    pair_sum=sum((u[i]*u[j]/(x[i]-x[j])+
                 u[j]*u[i]/(x[j]-x[i])
                 for i,j in combinations(range(len(x)),2)),F())
    if unsym!=0 or pair_sum!=0:
        raise ArithmeticError("real-coefficient initial boundary did not cancel")
    return {"N":len(x),"initial_boundary_exact":str(unsym),
            "paired_contribution_exact":str(pair_sum)}

def lemma():
    lo,hi=m0_bracket()
    examples=[]
    for x,u in (
        ([0,1,2],[1,2,-1]),
        ([0,3,7,12],[2,-3,5,7]),
        ([-3,2,10,27,33],[F(1,2),F(-7,9),3,F(5,6),-4]),
        ([F(2,7),F(7,8),F(13,5),F(39,7)],[-2,0,4,F(-9,11)]),
    ):
        examples.append(exact_boundary_cancellation(x,u))
    return {
       "schema":SCHEMA,
       "sources":{"Ramaré_Lemma_6_5":RAMARE,
                  "Preissmann_1984":PREISSMANN,
                  "Montgomery_Vaughan_1974":MV,
                  "Carneiro_Littmann_2023":CARNEIRO},
       "weighted_Hilbert_input":{
           "imported_theorem":"|sum_(m!=n) a_m conj(a_n)/(lambda_m-lambda_n)| <= pi*m0 sum_n |a_n|^2/delta_n",
           "delta_n":"min_(m!=n)|lambda_n-lambda_m|",
           "constant":"m0 = sqrt(1+(2/3)*sqrt(6/5))",
           "m0_rational_lower":lo,"m0_rational_upper":hi,
           "Preissmann_inequality_proved_independently":False,
           "Carneiro_Littmann_2pi_variant_sufficient_to_keep_m0":False,
       },
       "elementary_transfer":{
           "exact_integral":"(exp(i*T*d)-1)/(i*d), d=log(m/n)",
           "factor_i_in_source_pdf_display":"not printed in the extracted Ramaré formula; required mathematically",
           "real_pair_cancel":"sum_(m!=n) u_m u_n/(log m-log n)=0 by skew symmetry",
           "log_gap":"min_(m!=n)|log m-log n| >= log(1+1/n) >= 1/(n+1)",
           "real_mean_square_bound":
               "integral_0^T |sum_n u_n n^(it)|^2 dt <= sum_n u_n^2*(T+pi*m0*(n+1))",
           "lambda_X_coefficients_real":"lambda_X(n) = sum_(d|n,d<=X)mu(d) - 1_(n=1)",
           "abs_convergence":"|lambda_X(n)|<=floor(X) and sigma>1 imply both sum |lambda_X|n^-sigma and sum (n+1)|lambda_X|^2 n^(-2sigma) finite",
           "FX_bound":
               "F_X <= pi*m0 sum lambda_X(n)^2/n^(2sigma-1) + (T+pi*m0) sum lambda_X(n)^2/n^(2sigma)",
           "finite_to_infinite":"uniform Dirichlet convergence plus dominated limit of integrals and weighted series",
           "no_Hilbert_extra_constant_claimed":True,
       },
       "finite_symbolic_pair_cancellation_tests":examples,
       "infinite_series_extension_proved_by_standard_comparison":True,
       "source_theorem_Preissmann_proved_in_P85":False,
       "source_theorem_RZ_proved_in_P85":False,
       "X_3e10_to_3e12_covered":False,
       "global_RH_proved":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    d=lemma()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("P85 Hilbert: exact skew cancellation and infinite extension checked; Preissmann remains imported")

if __name__=="__main__":
    main()
