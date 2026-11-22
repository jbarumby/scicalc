"""NUMERIA P80: exact-rational v1 vs v3 Möbius-double-sum source comparison.

NO analytic theorem is proved here: all upper bounds on S_epsilon(X)
are imported from Ramaré–Zúñiga Alterman, arXiv:2603.25961v1/v3.
The only unconditional deductions performed by this module are rational
interval/parameter bookkeeping and monotonicity of imported majorants.
No global RH statement, no replacement of existing P78/P79 certificates.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from math import factorial
from pathlib import Path

HERE=Path(__file__).resolve().parent
P78=HERE.parent/"NUMERIA_P78_RAMARE_PUBLIC"/"certificate.json"
SCHEMA="NUMERIA.P80.RAMARE_V1_V3.PIECEWISE.1"
SOURCES={
    "v1":"https://arxiv.org/pdf/2603.25961v1",
    "v3":"https://arxiv.org/pdf/2603.25961v3",
    "cs_v2":"https://arxiv.org/html/2507.15184v2",
}
# Entries: (source version, epsilon cap, S_epsilon(X) upper cap, location).
# The v1 <= 0.526 for [1e12,1e16), epsilon <= 1/100,
# comes from SECTION SIX calculations, printed PDF page 21,
# and is NOT a displayed case in the v1 Theorem 1.2 statement!
ZONES=(
    ("3e12<=X<1e16",28,(
        ("v1","1/100","0.526","section 6, printed page 21"),
        ("v1","1/50","0.5503","theorem 1.2"),
        ("v1","1/25","0.593","theorem 1.2"),
        ("v3","1/50","0.5501","theorem 1.2, printed page 2"),
        ("v3","1/25","0.5936","theorem 1.2, printed page 2"),
    )),
    ("1e16<=X<1e33",36,(
        ("v1","1/100","0.4901","theorem 1.2"),
        ("v1","1/50","0.5104","theorem 1.2"),
        ("v3","1/100","0.4901","theorem 1.2"),
        ("v3","1/50","0.5103","theorem 1.2"),
    )),
    ("X>=1e33",75,(
        ("v1","1/100","0.4794","theorem 1.2"),
        ("v1","1/50","0.4905","theorem 1.2"),
        ("v1","1/25","0.5099","theorem 1.2"),
        ("v3","1/100","0.4794","theorem 1.2"),
        ("v3","1/50","0.4905","theorem 1.2"),
        ("v3","1/25","0.51","theorem 1.2"),
    )),
)
DELTA_LIMIT=F("0.303")
ORIGINAL_B4_FLOOR=F("0.54")

def fraction_repr(f):
    f=F(f)
    return f"{f.numerator}/{f.denominator}"

def verify_log_lower():
    """(11/4)>e with an exact rational series bound; powers upper-limit e^n."""
    e_upper=sum((F(1,factorial(k)) for k in range(7)), F(0))
    e_upper += F(8,7*factorial(7))
    if not (e_upper<F(11,4) and
            F(11,4)**28<F(3*10**12) and
            F(11,4)**36<F(10**16) and
            F(11,4)**75<F(10**33)):
        raise ArithmeticError("elementary e/log height proof invalid")
    return fraction_repr(e_upper)

def ceil_strict_1000(value):
    v=F(value)*1000
    return F(v.numerator//v.denominator+1,1000)

def best_zone(delta,zone_index,versions=("v1","v3")):
    label,log_lower,opts=ZONES[zone_index]
    eta=2*F(delta)/log_lower
    if not 0<eta<=F(1,25):
        raise ArithmeticError("eta outside preprint parameter domain")
    eligible=[]
    for version,epsilon_cap,S_bound,location in opts:
        if version not in versions or eta>F(epsilon_cap):
            continue
        # zeta(1+eta) <= 1+1/eta, and exp(-gamma eta) <= 1.
        # Hence eta exp(-gamma eta) zeta(1+eta) S_eta
        # <= (1+eta) * S_bound.
        upper=F(S_bound)*(1+eta)
        eligible.append((upper,version,F(epsilon_cap),F(S_bound),location))
    if not eligible:
        raise ArithmeticError(f"no authorized literature bound: {label} delta={delta}")
    best=min(eligible,key=lambda x:x[0])
    upper,version,epsilon_cap,S_bound,location=best
    return {
        "zone":label,"log_X_strict_lower":log_lower,
        "eta_max":fraction_repr(eta),"selected_version":version,
        "S_cap":fraction_repr(S_bound),
        "epsilon_cap":fraction_repr(epsilon_cap),
        "citation_location":location,
        "product_upper":fraction_repr(upper)
    }

def audit():
    verify_log_lower()
    source_bytes=P78.read_bytes()
    p78=json.loads(source_bytes)
    if p78.get("schema")!="NUMERIA.P78.RESTRICTED_B4.1" or len(p78.get("rows",[]))!=16:
        raise ArithmeticError("unexpected source P78 certificate")
    rows=[]
    for i,src in enumerate(p78["rows"],1):
        if src["row"]!=i:
            raise ArithmeticError("source row mismatch")
        delta=F(src["delta"])
        if not 0<delta<=DELTA_LIMIT:
            raise ArithmeticError("invalid delta in P78")
        old=[best_zone(delta,j,("v1",)) for j in range(3)]
        new=[best_zone(delta,j) for j in range(3)]
        old_worst=max(F(x["product_upper"]) for x in old)
        new_worst=max(F(x["product_upper"]) for x in new)
        old_b4=F(src["b4"])
        if old_worst!=F(src["worst_upper"]):
            raise ArithmeticError(f"P78 baseline not reproduced: row {i}")
        if not old_b4>old_worst or not new_worst<=old_worst:
            raise ArithmeticError("non-conservative source comparison")
        suggested=max(ORIGINAL_B4_FLOOR,ceil_strict_1000(new_worst))
        if not new_worst<suggested<=old_b4:
            raise ArithmeticError("suggested b4 not an upper bound within old cap")
        rows.append({
            "row":i,"delta":fraction_repr(delta),
            "b4_P78":fraction_repr(old_b4),
            "b4_P80_suggested":fraction_repr(suggested),
            "worst_v1_upper":fraction_repr(old_worst),
            "worst_min_v1_v3_upper":fraction_repr(new_worst),
            "old_zone_sources":old,"optimized_zone_sources":new,
            "v3_strictly_improves_a_zone":any(
                F(n["product_upper"])<F(o["product_upper"]) for n,o in zip(new,old)),
            "P79_existing_b4_still_valid":bool(old_b4>new_worst),
            "candidate_does_not_worsen_P79":bool(suggested<=old_b4)
        })
    changed=[r["row"] for r in rows if F(r["b4_P80_suggested"])<F(r["b4_P78"])]
    if changed != [10]:
        raise ArithmeticError(f"unexpected suggested b4 changes: {changed}")
    return {
        "schema":SCHEMA,"source":SOURCES,
        "P78_source_sha256":hashlib.sha256(source_bytes).hexdigest(),
        "methods":"exact fractions; bounds on S_epsilon imported from theorem 1.2 and v1 section 6",
        "scope":"X>=3e12; eta=2delta/log(X); 0<delta<=0.303; rows 1..16",
        "rows":rows,"changed_b4_rows":changed,
        "preserve_P79_B1_table":True,
        "reason":"new b4 <= old b4, and P79 B1(b4) increasing; old certified B1 still upper-bounds",
        "numerical_Suzuki_replay_not_performed_here":True,
        "external_lemmas_fully_verified":False,
        "RH_globally_proved":False,
    }

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--verify",type=Path)
    args=parser.parse_args()
    data=audit()
    if args.verify and json.loads(args.verify.read_text())!=data:
        raise SystemExit("P80 source certificate mismatch")
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
      "rows":len(data["rows"]),"changed_b4_rows":data["changed_b4_rows"],
      "row10_old_b4":data["rows"][9]["b4_P78"],
      "row10_suggested_b4":data["rows"][9]["b4_P80_suggested"],
      "all_candidate_non_worsening":all(r["candidate_does_not_worsen_P79"] for r in data["rows"]),
    },indent=2))

if __name__=="__main__":
    main()
