"""NUMERIA P88 — sharpen the exact open gate in Preissmann Lemma 6.

New UNIVERSAL results, proved in README:
 (1) A weaker all-real-nodes mixed kernel bound using Shan Lemma1;
 (2) an exact equivalence of sharp Lemma6 to a *mixed compensation*
     inequality E_ab>=0, exposing the obstruction to a diagonal-only proof;
 (3) an exact infinite uniform lattice equality, including the coefficient -3,
     with a non-asymptotic finite-truncation tail <= 2/(3*M^3).

Finite Fraction checks with a rational lower bound on pi are diagnostics.
The full arbitrary-spacing compensation inequality E_ab>=0 is NOT
independently proved in P88; Preissmann Lemma6 remains imported.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations

from NUMERIA_P86_PREISSMANN.spectral_rational import PAIRS, rational_nodes, gaps
from NUMERIA_P87_PREISSMANN.lemma6_bridge import machin_pi_lower

SCHEMA = "NUMERIA.P88.PREISSMANN.LEMMA6.MIXED_COMPENSATION.1"
SOURCE_YANGJIT = "https://arxiv.org/pdf/2203.14950"
SOURCE_PREISSMANN = "https://www.e-periodica.ch/cntmng?pid=ens-001:1984:30::31"

def reciprocal_kernel_identity(u, v):
    """Exactly check 1/(u²v²) = [u^-2+v^-2-2/(uv)]/(u-v)²."""
    u,v = F(u),F(v)
    if not u or not v or u==v:
        raise ValueError("nonzero distinct real offsets required")
    lhs = 1/(u*u*v*v)
    rhs = (1/(u*u)+1/(v*v)-2/(u*v))/(u-v)**2
    if lhs!=rhs:
        raise ArithmeticError("partial fractions not identical")
    return {"u":str(u),"v":str(v),"lhs":str(lhs),"rhs":str(rhs)}

def harmonics(m):
    if not isinstance(m,int) or m<0:
        raise ValueError("nonnegative integer")
    return sum((F(1,k) for k in range(1,m+1)),F())

def uniform_cross_parts(m):
    """Compute each infinite cross-series component symbolically via H_m.

    left tail = H_m/m, right tail = H_m/m,
    interior = -2H_(m-1)/m; total = 2/m².
    Each underlying series is absolutely convergent.
    """
    if not isinstance(m,int) or m<1:
        raise ValueError("m>=1 required")
    left=harmonics(m)/m
    right=left
    interior=-2*harmonics(m-1)/m
    total=left+right+interior
    if total!=F(2,m*m):
        raise ArithmeticError("uniform cross sum not 2/m²")
    return {"m":m,"left_tail":str(left),"right_tail":str(right),
            "between_endpoints":str(interior),
            "cross_sum_total":str(total)}

def lattice_finite(m,M,spacing=1):
    """Exact finite truncation of an infinite equispaced lattice.

    Nodes k*d for k in [-M,m+M], omit k=0,m in the cross kernel.
    Difference from infinite equals missing positive tails, bounded by
    2/(3 M³*d³). Never use finite lattice alone as uniform proof.
    """
    if type(m)!=int or m<1 or type(M)!=int or M<1:
        raise ValueError("positive integer m and M required")
    d=F(spacing)
    if d<=0:
        raise ValueError("positive spacing needed")
    s=sum((F(1,((k*d)**2)*((k-m)*d)**2)*d
           for k in range(-M,m+M+1) if k not in (0,m)),F())
    cross=sum((F(1,k*(k-m)) for k in range(-M,m+M+1)
               if k not in (0,m)),F())
    if not s>0:
        raise ArithmeticError("lattice sum must be positive")
    d3=d**3
    tail=F(2,3*M**3)/d3
    pi_lo=machin_pi_lower()
    formal_lower=2*pi_lo*pi_lo/(3*m*m*d3)-F(6,m**4)/d3
    if not formal_lower>s:
        raise ArithmeticError("finite lattice inequality not verified against rigorous lower pi")
    # The analytic upper bound on remainder can be checked against
    # rational upper bound on pi: using 12-term alternating arctangents,
    # pi_upper = 16*(sum_{k=0..10} arctan(1/5) terms)
    #             -4*(sum_{k=0..11} arctan(1/239) terms).
    # Lower and upper are strictly one-sided by the alternating series.
    from NUMERIA_P87_PREISSMANN.lemma6_bridge import atan_partial
    pi_up=16*atan_partial(F(1,5),10)-4*atan_partial(F(1,239),11)
    if pi_up<=pi_lo:
        raise ArithmeticError("incorrect Machin bounds")
    formal_upper=2*pi_up*pi_up/(3*m*m*d3)-F(6,m**4)/d3
    if not 0<formal_lower-s<formal_upper-s<=tail:
        raise ArithmeticError("exact uniform-lattice finite tail bound failed")
    return {
       "m":m, "radius_M":M, "spacing":str(d),
       "finite_left":str(s),
       "cross_sum_finite":str(cross),
       "limiting_formula":"d^-3 * (4*zeta(2)/m^2 - 6/m^4)",
       "pi_lower_rational":str(pi_lo),"pi_upper_rational":str(pi_up),
       "limiting_minus_finite_lower":str(formal_lower-s),
       "limiting_minus_finite_upper":str(formal_upper-s),
       "positive_tail_upper_exact":str(tail),
       "exact_tail_bound_verified":True,
    }

def mixed_balance(values,i,j):
    """The sharp Lemma6 is EXACTLY equivalent to E_ij>=0.

    E = Δ_i + Δ_j + 2C - 2(δ_i+δ_j)/h²,
    Δ_i=pi²/(3δ_i)-Σ_{k!=i}δ_k/(x_k-x_i)²,
    C=Σ_{k!=i,j}δ_k/[(x_k-x_i)(x_k-x_j)].
    h²*(RHS_Lemma6 - LHS_Lemma6) = E identically.
    """
    x=tuple(sorted(rational_nodes(values)))
    d=gaps(x)
    if type(i)!=int or type(j)!=int or not (0<=i<j<len(x)):
        raise ValueError("require 0 <= i < j < len(nodes)")
    h=x[j]-x[i]
    na=sum((d[k]/(x[k]-x[i])**2 for k in range(len(x)) if k!=i),F())
    nb=sum((d[k]/(x[k]-x[j])**2 for k in range(len(x)) if k!=j),F())
    cross=sum((d[k]/((x[k]-x[i])*(x[k]-x[j]))
               for k in range(len(x)) if k!=i and k!=j),F())
    mixed=sum((d[k]/((x[k]-x[i])**2*(x[k]-x[j])**2)
               for k in range(len(x)) if k!=i and k!=j),F())
    pi_lo=machin_pi_lower()
    positive_factor=F(1,d[i])+F(1,d[j])
    # exact rational lower bound on the symbolic compensation
    E_lo=(pi_lo*pi_lo/3)*positive_factor -na-nb+2*cross-F(2)*(d[i]+d[j])/h**2
    R_lo=(pi_lo*pi_lo/3)*positive_factor/h**2 -F(3)*(d[i]+d[j])/h**4
    if h*h*(R_lo-mixed)!=E_lo:
        raise ArithmeticError("sharp compensation equivalence failed")
    if not E_lo>0:
        raise ArithmeticError(f"finite Lemma6 mixed compensation negative: indices {i},{j}")
    # The weaker, proven all-N theorem is LHS <= 2/h² * (pi²/3
    #*(1/δ_i+1/δ_j) - (δ_i+δ_j)/h²).
    relaxed=2/h**2 * ((pi_lo*pi_lo/3)*positive_factor-(d[i]+d[j])/h**2)
    if not relaxed>=mixed:
        raise ArithmeticError("universal relaxed bound finite diagnostic failed")
    return {"i":i+1,"j":j+1,
            "endpoint_distance_exact":str(h),
            "diagonal_i_exact":str(na),"diagonal_j_exact":str(nb),
            "mixed_cross_exact":str(cross),
            "sharp_original_LHS_exact":str(mixed),
            "sharp_RHS_with_rational_pi_lower":str(R_lo),
            "compensation_with_pi_lower_exact":str(E_lo),
            "E_equals_h_squared_times_sharp_slack":True,
            "relaxed_RHS_with_pi_lower":str(relaxed),
            "all_finite_signs_strict":True}

def audit():
    kernel=[
        reciprocal_kernel_identity("1/7","-13/11"),
        reciprocal_kernel_identity("-2","4/9"),
        reciprocal_kernel_identity("9/13","7/3"),
    ]
    lattice=[lattice_finite(m,M,d) for m,M,d in
       ((1,8,1),(2,9,1),(3,16,F(3,2)),(5,20,F(1,3)),
        (8,18,F(11,7)),(13,24,F(2,5)),(20,32,1))]
    configs={**PAIRS,
        "nine_exponential":[str(F(2**k)) for k in range(9)],
        "strongly_clustered":["0","1/10000","1/5000","1/2500","1","5","5.001","100"],
        "uniform_31":[str(k) for k in range(31)]
    }
    triples=[]
    for name,values in configs.items():
        records=[mixed_balance(values,i,j)
                 for i,j in combinations(range(len(values)),2)]
        triples.append({
            "configuration":name, "N":len(values),
            "pairs":len(records),
            "minimum_positive_compensation_with_pi_lower":str(
                min(F(r["compensation_with_pi_lower_exact"]) for r in records)),
            "negative_cross_count":sum(F(r["mixed_cross_exact"])<0 for r in records),
            "records":records,
        })
    # A decisive counterexample to assuming mixed cross >=0:
    simple=mixed_balance(["0","1","2","3"],0,3)
    if F(simple["mixed_cross_exact"])!=-1:
        raise ArithmeticError("cross-term obstruction not demonstrated")
    return {
        "schema":SCHEMA,
        "sources":{
            "Yangjit_Lemma3_2_eq3_6":"https://arxiv.org/pdf/2203.14950",
            "Preissmann_original_Lemme6":SOURCE_PREISSMANN,
            "P87_general_diagonal":"NUMERIA_P87_PREISSMANN/README.md",
        },
        "theorem_weak_all_real_nodes":
           "For all finite distinct real nodes and distinct indices i,j, mixed S <= 2/h²*[pi²/3*(1/delta_i+1/delta_j)-(delta_i+delta_j)/h²]",
        "proof_weak":"(u-v)²<=2(u²+v²) and P87 diagonal bound, subtract the other endpoint exactly",
        "sharp_lemma6_equivalent_gate":
           "E_ij=(pi²/(3δ_i)-Σ_(k!=i)δ_k/(x_k-x_i)²)+(pi²/(3δ_j)-Σ_(k!=j)δ_k/(x_k-x_j)²)+2Σ_(k!=i,j)δ_k/((x_k-x_i)(x_k-x_j))-2(δ_i+δ_j)/h²>=0",
        "uniform_lattice_exact":"For x_k=k*d, endpoints 0 and m*d, infinite sum = d^-3*(4*zeta(2)/m²-6/m⁴), equality in Lemma6",
        "uniform_cross_series_exact":"Σ_(k in Z, k!=0,m)1/[k(k-m)] = 2/m² from harmonic telescoping",
        "uniform_sharp_coefficient_negative_three":True,
        "infinite_uniform_identity_proved_analytically_in_README":True,
        "universal_relaxed_mixed_inequality_proved":True,
        "arbitrary_spacing_sharp_compensation_proved":False,
        "independent_full_Preissmann_Lemma6_proved":False,
        "RZ_S_epsilon_reproved":False,
        "RH_proved":False,
        "partial_fraction_examples":kernel,
        "uniform_cross_harmonic":[uniform_cross_parts(m) for m in (1,2,3,5,8,13,20)],
        "finite_lattice_tail_examples":lattice,
        "cross_sign_counterexample":simple,
        "finite_mixed_configurations":triples,
        "finite_mixed_pair_count":sum(r["pairs"] for r in triples)
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    data=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({
      "status":"UNIFORM_LATTICE_SHARPNESS_PLUS_UNIVERSAL_RELAXED_BOUND",
      "finite_checked_pairs":data["finite_mixed_pair_count"],
      "sharp_all_spacings_Lemma6_proved":False,
      "uniform_sharp_coefficient_negative_three":True
    },indent=2))

if __name__=="__main__":
    main()
