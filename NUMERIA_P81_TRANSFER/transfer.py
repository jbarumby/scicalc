"""NUMERIA P81: conditional exact-rational transfer S_epsilon -> zeta*S -> Lemma 6.

Analytic input (not proved here): public Ramaré–Zúñiga bounds S_epsilon(X)
in P80's three height zones, including v1 Sec. 6 and v3 Thm. 1.2.
We prove algebraically that the selected P80 b4 controls BOTH
epsilon_2=2delta/log(X) and epsilon_1=delta/log(X), uniformly X>=3e12,
using only zeta(1+epsilon) <= 1 + 1/epsilon (integral test).

The C3/C4 transfer is conditional on the KLN/Chourasiya–Simonič
preliminary F_X inequality, specifically CS v2 Lemma 6 proof line before
the final substitution. No independent proof of the external lemma or RH.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from NUMERIA_P80_SOURCE_GATE.v3_gate import (
    audit as audit_p80, best_zone, fraction_repr, ZONES,
)

HERE=Path(__file__).resolve().parent
SCHEMA="NUMERIA.P81.S_EPS_TO_LEMMA6.CONDITIONAL.1"
SOURCE_PDF_SHA256="69baebd7d0d79a00ed76109f3a5aacdc1be2e79dd66489900c82d813566a2b43"
SOURCE_JSON_SHA256="0e821f2f285b3d15fe127d2728836cbc3af5b044ca916bb9209741becd73ce83"

def check_positive_bound(epsilon_max, source_cap, b4):
    """Conditional theorem:
       S_eps(X)<=source_cap for 0<eps<=epsilon_max
       zeta(1+eps)<=1+1/eps
       implies eps*zeta(1+eps)*S_eps(X) <= (1+eps)*source_cap.
       A b4 greater than this EXACT rational expression closes the gate.
    """
    e=F(epsilon_max); cap=F(source_cap); b=F(b4)
    if not (0<e<=F(1,25) and 0<cap and b>0):
        raise ArithmeticError("parameters outside imported source theorem scope")
    bound=cap*(1+e)
    if not b>bound:
        raise ArithmeticError("published b4 not above public zeta*S upper bound")
    return b-bound

def audit():
    baseline=audit_p80()
    rows=[]
    smallest=None
    for row in baseline["rows"]:
        i=row["row"]
        delta=F(row["delta"])
        b4old=F(row["b4_P78"])
        b4=F(row["b4_P80_suggested"])
        if not F(0)<delta<=F("0.303") or not F(0)<b4<=b4old:
            raise ArithmeticError("bad row delta or proposed b4")
        zones=[]
        for zi,(zone_label,log_floor,opts) in enumerate(ZONES):
            source=best_zone(delta,zi)
            cap=F(source["S_cap"])
            valid_eps=F(source["epsilon_cap"])
            zone=[]
            for branch,k in (("epsilon_2",2),("epsilon_1",1)):
                eps_max=k*delta/F(log_floor)
                if eps_max>valid_eps:
                    raise ArithmeticError("public S_eps bound outside its epsilon scope")
                slack=check_positive_bound(eps_max,cap,b4)
                # e2 branch must be as tightly controlled as the P80 source already is.
                if branch=="epsilon_2" and F(source["product_upper"]) != cap*(1+eps_max):
                    raise ArithmeticError("P80 upper normalization mismatch")
                if smallest is None or slack<smallest:
                    smallest=slack
                zone.append({
                    "branch":branch,
                    "epsilon_zone_max":fraction_repr(eps_max),
                    "source_epsilon_cap":fraction_repr(valid_eps),
                    "source_S_cap":fraction_repr(cap),
                    "zeta_product_bound":fraction_repr(cap*(1+eps_max)),
                    "slack_under_b4":fraction_repr(slack),
                    "imported_from":source["selected_version"],
                    "source_location":source["citation_location"],
                })
            zones.append({
                "zone":zone_label,
                "log_X_lower":log_floor,
                "branches":zone,
                "both_branches_strict":True,
            })
        rows.append({
            "row":i,
            "delta":fraction_repr(delta),
            "b4_from_P80":fraction_repr(b4),
            "b4_previous_P78":fraction_repr(b4old),
            "uniform_lambda_bound":"epsilon * zeta(1+epsilon) * S_epsilon(X) < b4",
            "zones":zones,
            "both_epsilon_branches_closed":True,
        })
    if len(rows)!=16 or smallest is None or not smallest>0:
        raise ArithmeticError("missing source branch or nonpositive slack")
    # For C3, with sigma=1+delta/logX, eps2=2delta/logX:
    # pi*m0*Lambda(1+eps2) <= pi*m0*b4*logX/(2delta).
    # Original C3 includes the harmless exp(2gamma delta/logX0) >=1.
    # For C4a, eps1=delta/logX, (5X)^(-eps1) = 5^(-eps1)*e^(-delta):
    # Lambda(1+eps1)/(5(5X)^eps1)
    # <= [b4*logX/(5delta e^delta)]*5^(-eps1)
    # <= b4*logX/(5delta e^delta), since 5>1 and eps1>0.
    # The independent b3*e^(-2delta)/logX0 term is unchanged.
    return {
        "schema":SCHEMA,
        "public_source_proof_status":"conditional on Ramaré–Zúñiga v1/v3 S_epsilon estimates",
        "user_note_pdf_sha256":SOURCE_PDF_SHA256,
        "user_note_json_sha256":SOURCE_JSON_SHA256,
        "source_P78_sha256":baseline["P78_source_sha256"],
        "region":"X>=3e12; 16 rows, 3 height zones, 2 epsilon branches",
        "zeta_upper_used":"zeta(1+eps) <= 1 + 1/eps (integral comparison)",
        "unneeded_unsupported_claims":[
            "G_prime<0 from an unbounded O(epsilon^2) remainder",
            "dropping exp(gamma*epsilon_1) from an upper bound without 5^(-epsilon_1)",
            "fitting L numerically proves it equals log(T0) exactly",
        ],
        "rows":rows,
        "checked_analytic_implications":{
            "eps2_C3":"pi*m0*Lambda(1+eps2) <= pi*m0*b4*logX/(2*delta) <= C3*logX",
            "eps1_C4":"Lambda(1+eps1)/(5*(5X)^eps1) <= b4*logX*e^(-delta)*5^(-eps1)/(5*delta) <= C4a*logX",
            "C4_residual":"b3*X^(-2delta/logX)=b3*e^(-2delta) <= b3*e^(-2delta)*logX/logX0",
            "final_F_X":"F_X bound follows if the *independent* preliminary CS/KLN mean-square inequality is accepted",
        },
        "smallest_rational_b4_slack":fraction_repr(smallest),
        "total_branch_zone_checks":len(rows)*3*2,
        "all_checks_pass":True,
        "external_S_epsilon_theorems_proved_here":False,
        "independent_KLN_F_X_lemma_proved_here":False,
        "original_CS_Lemma6_full_X_range_proved_here":False,
        "Suzuki_new_window_proved_here":False,
        "global_RH_proved":False,
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    r=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({
        "checks":r["total_branch_zone_checks"],
        "smallest_rational_b4_slack":r["smallest_rational_b4_slack"],
        "all_checks_pass":r["all_checks_pass"],
        "scope":r["region"],
    },indent=2))

if __name__=="__main__":
    main()
