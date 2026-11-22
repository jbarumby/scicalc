"""P87: ledger for Preissmann Lemma 6 and exact spectral combination.

The all-N two-point/two-denominator LEMMA 6 is NOT established here.
We verify exact rational finite test instances of its statement using a
strict rational lower bound on pi, and prove the purely ALGEBRAIC
implication:
 Lemma1 (diagonal, now independently proved in P87) + Lemma6
 => U/V <= (pi²/3)(1+sqrt(6/5))
 => ||B||² <= pi²(1+(2/3)sqrt(6/5)).
No finite tests should be promoted to all-N theorem.
"""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
from NUMERIA_P86_PREISSMANN.spectral_rational import PAIRS,rational_nodes,gaps

SCHEMA="NUMERIA.P87.PREISSMANN.LEMMA6.CONDITIONAL.BRIDGE.1"
LEMMA6_SOURCE="https://arxiv.org/pdf/2203.14950"
ORIGINAL="https://www.e-periodica.ch/cntmng?pid=ens-001:1984:30::29"

def atan_partial(x,n):
    if x<=0 or x>=1 or n<0:raise ValueError("x in (0,1), n>=0 required")
    return sum(((-1 if k%2 else 1)*x**(2*k+1)/F(2*k+1)
                for k in range(n+1)),F(0))

def machin_pi_lower():
    """Rational strictly smaller than pi, with proven alternating signs."""
    # atan(1/5) lower with last term odd k=11
    # atan(1/239) upper with last term even k=10
    p=16*atan_partial(F(1,5),11)-4*atan_partial(F(1,239),10)
    if not F(157,50)<p:
        raise ArithmeticError("pi bound too coarse")
    return p

def lemma6_finite(nodes,selected=None):
    """Finite exact instances of Lemma6 (3.6) with rigorous pi lower.

    For all finite samples, we prove a STRONGER inequality than the
    original RHS using pi_lower^2; this cannot replace the universal
    all-configuration theorem in Preissmann 1984.
    """
    x=tuple(sorted(rational_nodes(nodes)))
    d=gaps(x)
    pairs=list(combinations(range(len(x)),2))
    if selected is not None:
        pairs=[pair for pair in pairs if pair in selected]
    if not pairs:raise ValueError("no pairs")
    p=machin_pi_lower()
    slack=None
    rows=[]
    for l,m in pairs:
        h=x[m]-x[l]
        term=sum((d[k]/((x[k]-x[l])**2*(x[k]-x[m])**2)
                  for k in range(len(x)) if k!=l and k!=m),F(0))
        rhs=p*p*(d[l]+d[m])/(3*d[l]*d[m]*h*h)-3*(d[l]+d[m])/(h**4)
        diff=rhs-term
        if diff<=0:
            raise ArithmeticError(f"finite Lemma6 fails with pi_lower: ({l},{m})")
        if slack is None or diff<slack:slack=diff
        rows.append({"i":l+1,"j":m+1,"left_exact":str(term),
                     "rhs_with_pi_lower_exact":str(rhs),
                     "slack_exact":str(diff),"finite_strict":True})
    return {"size":len(x),"x":[str(q) for q in x],
            "pi_lower":str(p),"number_pairs":len(rows),
            "minimum_exact_slack":str(slack),"pairs":rows,
            "all_finite_gates_pass":True}

def quadratic_surd_chain():
    """Exact coefficient algebra in Q[sqrt(6/5)] (no float pi/zeta).

    Let c=(1+sqrt(6/5))/3. Then
    c²-(2/3)c=1/45. So u²<=(1/45)+(2/3)u implies
    u <= c (u>=0). The Preissmann triple-sum Lemma6 supplies this
    inequality after dividing by pi^4 and V^2.
    Also 1/3+2*c = 1+(2/3)sqrt(6/5).
    """
    t2=F(6,5)
    def add(x,y):return (x[0]+y[0],x[1]+y[1])
    def mul(x,y):return (x[0]*y[0]+x[1]*y[1]*t2,
                        x[0]*y[1]+x[1]*y[0])
    c=(F(1,3),F(1,3))
    test=add(mul(c,c),mul((F(-2,3),F(0)),c))
    target=(F(1,45),F(0))
    if test!=target:
        raise ArithmeticError("Preissmann quadratic coefficient identity failed")
    norm=add((F(1,3),F(0)),mul((F(2),F(0)),c))
    if norm!=(F(1),F(2,3)):
        raise ArithmeticError("m0 square coefficient chain failed")
    return {
        "c":"(1+sqrt(6/5))/3",
        "coefficient_check":"c^2-(2/3)c=1/45",
        "exact_quadratic_components_c":[str(z) for z in c],
        "exact_quadratic_identity_lhs":[str(z) for z in test],
        "exact_m0_squared":[str(z) for z in norm],
        "conditional_U_bound":"U^2 <= V*(pi^4/45*V + (2*pi^2/3)*U)",
        "conditional_spectral":"||B||_2^2 <= pi^2/3+2U/V <= pi^2*(1+(2/3)*sqrt(6/5))",
        "algebra_closed_given_Lemma6":True,
    }

def audit():
    samples={}
    for name,values in PAIRS.items():
        samples[name]=lemma6_finite(values)
    # 41 nodes generates 820*39 rational factors; inspect
    # representative distant and adjacent pairs without pretending exhaustive
    extra=lemma6_finite(range(41),selected={(0,1),(10,20),(19,20),(20,21),(0,40)})
    samples["uniform41_representative_pairs"]=extra
    total=sum(x["number_pairs"] for x in samples.values())
    return {
        "schema":SCHEMA,
        "source":{"Yangjit_Theorem1_2_Lemma3_2":LEMMA6_SOURCE,
                  "Preissmann_1984_Lemme6":ORIGINAL},
        "precise_uniform_Lemma6":
            "Σ_(k!=l,m) δ_k/((x_k-x_l)^2(x_k-x_m)^2) <= pi^2*(δ_l+δ_m)/(3δ_lδ_m(x_l-x_m)^2) - 3*(δ_l+δ_m)/(x_l-x_m)^4",
        "pi_lower_rational":str(machin_pi_lower()),
        "finite_samples":samples,"finite_exact_pair_count":total,
        "spectral_bridge":quadratic_surd_chain(),
        "Preissmann_Lemma1_all_N_reconstructed_elsewhere_in_P87":True,
        "Preissmann_Lemma6_all_N_proved_here":False,
        "Preissmann_constant_m0_unconditionally_independent":False,
        "finite_samples_are_not_a_universal_proof":True,
        "C3_C4_B1_recertification_repeated":False,
        "global_RH_proved":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    cert=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(cert,indent=2,ensure_ascii=False)+"\n")
    print("P87 finite Lemma6 pair checks",cert["finite_exact_pair_count"])
    print("P87 source state: Lemma1 complete via Shan; Lemma6 all N not proved")

if __name__=="__main__":
    main()
