"""NUMERIA P90: sharp Preissmann Lemma 6 for nearby pairs, all real spacings.

MATHEMATICAL THEOREMS (full universal proofs in README):
  - The original SHARP Lemma 6 holds for all finite real configurations
    whenever b-a=1 (adjacent endpoints), or b-a=2 (one interior point).
    Exterior nodes can be arbitrarily numerous and irregular.
  - Equivalent sufficient budget for ANY b>a, h=x_b-x_a:
       M_(a,b) <= D(delta_a/h)+D(delta_b/h) after setting h=1,
       D(s)=(2*zeta(2)-S(s))/s-3*s,
       S(s)=sum_(j>=1) 1/[j²(1+j*s)²].
    This suffices for sharp Lemma6 by P89 exterior tails.
    It is proven for b-a<=2, NOT for general b-a.
  - D(s) strictly decreases: exact rational polynomial
       1+3z+6z²+2z³ > 0 controls each summand's derivative.
  - Equality in the budget for normalized UNIFORM m-gap interior
    for each integer m>=1, via finite harmonic sums.

RATIONAL TESTS are DIAGNOSTIC; the universal proofs in README are not
deduced by finite testing. No proof of RH or full all-pair Lemma6.
"""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path

from NUMERIA_P86_PREISSMANN.spectral_rational import rational_nodes,gaps
from NUMERIA_P88_LEMMA6.compensation import machin_pi_lower
from NUMERIA_P89_LEMMA5.exterior import majorant

SCHEMA="NUMERIA.P90.PREISSMANN.INTERIOR.NEIGHBOR_SHARP.1"
CHORD_SCHEMA="NUMERIA.P90.CHORD.KNOTS.EXACT.1"

def zeta2_lower():
    p=machin_pi_lower()
    if not p>F(3):raise ArithmeticError("pi lower must exceed 3")
    return p*p/6

def H(n,p=1):
    if type(n)!=int or n<0 or type(p)!=int or p<1:
        raise ValueError("valid harmonic indices required")
    return sum((F(1,j**p) for j in range(1,n+1)),F())

def uniform_budget_identity(m):
    """Exact equality M_uniform(m)=2*D(1/m), no transcendental floats."""
    if type(m)!=int or m<1:raise ValueError("m>=1")
    D=F(m)*H(m,2)+2*H(m,1)-F(3,m)
    M=sum((F(m**3,k*k*(m-k)**2) for k in range(1,m)),F())
    if M!=2*D:
        raise ArithmeticError(f"uniform harmonic identity failed for m={m}")
    return {"m":m,"D_exact":str(D),"interior_exact":str(M),
            "two_sided_budget_exact":str(2*D),
            "identity_true":True}

def chord_knot_values():
    """Knot values S(0),S(1/2),S(1) in the Q+Q*zeta(2) module."""
    # For S(s)=sum_(j>=1) 1/[j²(1+js)²]:
    # S(0)=Z, S(1/2)=2Z-11/4, S(1)=2Z-3.
    # Derived via partial fractions, not numerical approximations.
    # Exact m=2 finite-harmonic computation:
    # S(1/m)=2Z-H_m^(2)-(2/m)*H_m.
    b2=-H(2,2)-F(2,2)*H(2)
    b1=-H(1,2)-F(2)*H(1)
    if b2!=-F(11,4) or b1!=-F(3):
        raise ArithmeticError("S knot constants incorrect")
    return {
        "schema":CHORD_SCHEMA,
        "S_at_zero":"Z",
        "S_at_half":"2*Z-11/4",
        "S_at_one":"2*Z-3",
        "Z":"zeta(2)=pi^2/6 > 3/2",
        "convex_chord_0_to_1":"S(s) <= (1-s)*Z+s*(2*Z-3), 0<=s<=1",
        "convex_chord_0_to_half":"S(t) <= (1-2t)*Z+2t*(2*Z-11/4), 0<=t<=1/2",
        "convex_chord_half_to_one":"S(1-t) <= 2t*(2*Z-11/4)+(1-2t)*(2*Z-3), 0<=t<=1/2",
        "D_monotonic_rational_polynomial":"1+3*z+6*z^2+2*z^3 > 0 for z>0",
        "D_for_adjacent":"D(s) >= (1-s)*(Z/s+3) >= 0",
        "one_interior_terminal_polynomial":
            "(3/2)*(1-2t)*(1-t)^2 + t*(1-t)*(11/2-3t) - 1 = (1-2t)*(1+t)/2 >=0",
    }

def one_interior_poly(t):
    t=F(t)
    if t<=0 or t>F(1,2):
        raise ValueError("t in (0,1/2] required")
    lhs=F(3,2)*(1-2*t)*(1-t)**2+t*(1-t)*(F(11,2)-3*t)-1
    rhs=(1-2*t)*(1+t)/2
    if lhs!=rhs or rhs<0:
        raise ArithmeticError("two-gap convex chord terminal gate failed")
    return {"t":str(t),"positive_polynomial_exact":str(rhs),
            "factorization_exact":True}

def mid_scaled(gaps_interior):
    """Normalized interior kernel for a fixed pair of endpoints 0 and 1."""
    g=[F(v) for v in gaps_interior]
    if not g or any(v<=0 for v in g):
        raise ValueError("at least one positive gap required")
    total=sum(g,F())
    a=[v/total for v in g]
    t=F(0); mid=F(0)
    for j in range(1,len(a)):
        t+=a[j-1]
        mid+=min(a[j-1],a[j])/(t*t*(1-t)**2)
    return a,mid

def D_lower(s,K=96):
    """PROVEN rational lower enclosure of D(s), not a floating guess.

    Z >= pi_lower²/6.
    S(s) <= sum_(j=1..K)1/[j²(1+js)²] + 1/(3*s²*K³).
    Also S(s)<=Z implies D(s)>=Z/s-3s. Max of two valid
    lower bounds mitigates very small s.
    """
    s=F(s)
    if not F(0)<s<=1 or type(K)!=int or K<1:
        raise ValueError("0<s<=1 and K>=1 required")
    Z=zeta2_lower()
    prefix=sum((F(1,j*j*(1+j*s)**2) for j in range(1,K+1)),F())
    upper=prefix+F(1,3*s*s*K**3)
    first=(2*Z-upper)/s-3*s
    second=Z/s-3*s
    return max(first,second)

def budget_instance(g):
    a,mid=mid_scaled(g)
    da,db=a[0],a[-1]
    # δa,δb may be smaller due to OUTSIDE nodes; D decreasing, so
    # using the endpoint interior gaps is the minimal available budget.
    bound=D_lower(da)+D_lower(db)
    if len(a)>2 and bound<=mid:
        raise ArithmeticError("selected nonuniform interior budget not certified")
    return {"gaps_exact":[str(x) for x in g],
            "normalized_gaps":[str(x) for x in a],
            "interior_exact":str(mid),
            "D_endpoint_lower_sum":str(bound),
            "certified_positive_margin":str(bound-mid),
            "strict_for_selected_example":bound>mid}

def sharp_pair_direct(nodes,a,b):
    """Independent Fraction scan using a *lower* exact Machin pi bound."""
    x=tuple(sorted(rational_nodes(nodes)))
    if type(a)!=int or type(b)!=int or not 0<=a<b<len(x):
        raise ValueError("bad indices")
    d=gaps(x);h=x[b]-x[a]
    term=sum((d[k]/((x[k]-x[a])**2*(x[k]-x[b])**2)
              for k in range(len(x)) if k!=a and k!=b),F())
    B=(machin_pi_lower()**2/F(3))*(
        F(1,d[a])+F(1,d[b]))/(h*h)-3*(d[a]+d[b])/(h**4)
    if b-a<=2 and not B>term:
        raise ArithmeticError("nearby-pair sharp Lemma6 finite regression failed")
    return {"a":a+1,"b":b+1,"interior_count":b-a-1,
            "delta_a":str(d[a]),"delta_b":str(d[b]),"h":str(h),
            "sharp_LHS_exact":str(term),
            "sharp_RHS_pi_lower_exact":str(B),
            "strict_pi_lower_margin":str(B-term),
            "all_real_theorem_covered":b-a<=2}

def proof_ledger():
    examples=[]
    for name,values in {
       "uniform_7":[str(j) for j in range(7)],
       "uniform_19":[str(j) for j in range(19)],
       "alternating":["0","1","4","5","8","10","21","22","25"],
       "very_clustered":["0","1/2500","1/1250","1/400","1","5","100","101","201"],
       "irregular":["0","1/3","2","10/3","16","67/3","100"],
       "wide_spacings":["0","1/20000","1/10000","3/10000","1","1000"],
    }.items():
        certs=[]
        for a in range(len(values)-1):
            for step in (1,2):
                b=a+step
                if b<len(values):
                    certs.append(sharp_pair_direct(values,a,b))
        examples.append({"name":name,"node_count":len(values),
                         "nearby_pairs":certs,"count":len(certs)})
    gap_sets=[
       ["1","2","1","4"],
       ["2","3","1","5"],
       ["1","1","2","2","3"],
       ["2","2","3","4","5"],
       ["1","4","1","1","3","2"],
       ["1","2","3"],
       ["2","1","4"],
       ["1","2","1","2","1"],
    ]
    interiors=[budget_instance(v) for v in gap_sets]
    if not all(x["strict_for_selected_example"] for x in interiors):
        raise ArithmeticError("strict nonuniform interior diagnostics failed")
    return {
       "schema":SCHEMA,
       "theorem":"Sharp Preissmann Lemma6 for ANY finite real configuration if b-a<=2",
       "analytic_source":{
          "P89_uniform_exterior_tails":"NUMERIA_P89_LEMMA5/README.md",
          "Preissmann_original_Lemma5":"https://www.e-periodica.ch/cntmng?pid=ens-001:1984:30::29",
          "Yangjit_Lemma3_2":"https://arxiv.org/pdf/2203.14950"
       },
       "universal_adjacent_pairs_proved_in_README":True,
       "universal_one_interior_pairs_proved_in_README":True,
       "D_strictly_decreasing_proved_by_polynomial":True,
       "arbitrary_interior_count_Lemma6_proved":False,
       "full_Preissmann_Lemma5_variations_reconstructed":False,
       "full_Preissmann_Lemma6_independently_reproved":False,
       "m0_without_external_Lemma6_proved":False,
       "RZ_S_epsilon_reproved":False,
       "RH_proved":False,
       "uploaded_P87":{
          "ledger_SHA256":"f593dc80743a5510c064aaa6ad3d4b2149f9950389c2294e16e425121b728da7",
          "log_SHA256":"1b594dda439681ab78a286dca2aa810c3dab4b5f4568da47ecce380959d11da3",
          "tex_SHA256":"d16771dab750ef4f00e00bf841ddcf5f179e412ef6392914ea819a0e7b9d706d",
          "verificacao_SHA256":"5089f7016827dffb5a6202cd783ca70ab099ba5a87db37e954e2192a1b8ff14a",
          "ledger_gates_reported":24,"ledger_gates_failed_reported":0,
          "ledger_B5_wrong_upper_pi_and_float_margins":True,
          "cross_control_imported_in_upload":True
       },
       "chord_knots":chord_knot_values(),
       "one_interior_symbolic_regressions":[one_interior_poly(F(k,100)) for k in (1,7,17,25,33,49,50)],
       "uniform_all_m_closed_forms":[uniform_budget_identity(m) for m in range(1,41)],
       "exact_nearby_pair_cases":examples,
       "exact_nearby_pair_count":sum(x["count"] for x in examples),
       "interior_sufficient_budget_finite_checks":interiors,
       "original_mixed_compensation_all_spacing_gate_remains":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    d=proof_ledger()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("P90 closed nearby-pair theorem: zero/one interior node, any real spacings")
    print("Finite independent checks:",d["exact_nearby_pair_count"])
    print("Full all-pair L6 still external:",not d["full_Preissmann_Lemma6_independently_reproved"])

if __name__=="__main__":main()
