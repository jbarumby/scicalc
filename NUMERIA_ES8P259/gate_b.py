"""NUMERIA ES8-P259 Gate B: staged certified odd 377x377 banded-floor test.

Data: P253 authenticated GL60 aggregate, P258 authenticated artifact.
Mathematics: R_259 = R_beta + (1/4-beta)(I-G_504) + (1/4)(I-G_850).

Proof is **conditional on** the proven P253 GL60 error budget, and the
rigorous endpoint Bessel formulas in P258. The code computes actual interval
matrices using python-flint/Arb. It never infers SPD from numerical float.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import file_digest
import json
from pathlib import Path

from flint import arb, arb_mat, ctx

from NUMERIA_ES8P251.gate_b.finite_head_certify import (
    assemble_finite_head, certify_center_weyl,
)
from NUMERIA_ES8P251.gate_b.finite_head import _dyadic_fraction
from NUMERIA_ES8P251.gate_b.finite_head_quadrature import fixed_beta
from NUMERIA_ES8P258.audit import (
    OLD_AGGREGATE, OLD_SHA, HALF_WIDTH, MODES, CUTOFF, FLOOR,
    _old_matrices, update_head, _trial_after_failed_pivot, _schur,
)
from NUMERIA_ES8P258.concentration import analytic_concentration, matrix_rows
from NUMERIA_ES8P258.independent_rayleigh import replay_rayleigh

SOURCE_SHA = "1aa765357cfe1fa46f864d5b57fc407157375e75b3fa8b23b77f0db5c2e4985c"
CUTOFF2 = 850
FLOOR2 = F(1,2)
TARGET_PREFIXES = (29, 36, 48, 64, 96, 128, 192, 259, 377)


def require(test, message):
    if not test:
        raise ArithmeticError(message)


def validate_parameters():
    require(CUTOFF==504 and MODES==377 and FLOOR==F(1,4),
            "P258 model parameters changed")
    require(HALF_WIDTH==F("0.972955074527656652552676371721590"),
            "P258 rational halfwidth changed")
    # For x>1, log(x) > 2(x-1)/(x+1); pi enters neither test.
    bound=F(2*(CUTOFF2-CUTOFF), CUTOFF2+CUTOFF)
    require(bound>F(1,2),"floor 1/2 not certified from Binet lower bound")
    return {"cutoffs":[CUTOFF,CUTOFF2],"floors":[str(FLOOR),str(FLOOR2)],
            "log_850_504_lower":str(bound)}


def check_p258_archive(path):
    p=Path(path)
    with p.open("rb") as handle:
        digest=file_digest(handle,"sha256").hexdigest()
    require(digest==SOURCE_SHA,"P258 artifact SHA mismatch")
    j=json.loads(p.read_text())
    require(j["schema"]=="NUMERIA.ES8P258.RAISED_AUXILIARY.1", "P258 schema mismatch")
    require(F(j["L_exact"])==HALF_WIDTH and j["cutoff"]==CUTOFF
            and j["modes_per_sector"]==MODES,"P258 geometry mismatch")
    require(j["sectors"]["even"]["validation"]["certified_positive"],
            "P258 even margin not established")
    require(j["sectors"]["odd"]["negative_witness"] is not None,
            "P258 odd negative witness missing")
    return j


def prefix(M, size):
    if not 1<=size<=M.nrows() or M.nrows()!=M.ncols():
        raise ValueError("incorrect prefix size")
    return arb_mat(size,size,[M[i,j] for i in range(size) for j in range(size)])


def rational_negative_witness(matrix, k, errors):
    """Verify sign with independent exact rational matrix-box replay.

    Does NOT claim transfer to any other half-width; 850 cutoff changes
    the model and invalidates blindly reusing P258's Lipschitz factor.
    """
    if k is None:
        return None
    try:
        trial=_trial_after_failed_pivot(matrix,k)
    except (ZeroDivisionError, ValueError):
        return None
    if not 0<len(trial)<=matrix.nrows():
        return None
    dyadic=[_dyadic_fraction(x.man_exp()) for x in trial]
    m=len(dyadic)
    out=replay_rayleigh(dyadic,{"P259_banded":(matrix_rows(prefix(matrix,m)),1)},
                        operator_errors=errors)
    if out["witness_strict_sign"]!="negative":
        return None
    return {"coordinate_count":m,"coordinates_dyadic":
            [[str(q.numerator),str(q.denominator)] for q in dyadic],
            "exact_rational_replay":out,"scope":"negative at exact P258 rational halfwidth only"}


def sector_prefix(whole_old, size, bits):
    with ctx.workprec(bits):
        quad=prefix(whole_old,size)
        old,assembly=assemble_finite_head(quad,"odd",cutoff=CUTOFF,
                                           half_width=HALF_WIDTH)
        gram504=analytic_concentration(cutoff=CUTOFF,modes=size,
                                          half_width=HALF_WIDTH,bits=bits)[1]
        gram850=analytic_concentration(cutoff=CUTOFF2,modes=size,
                                          half_width=HALF_WIDTH,bits=bits)[1]
        delta=arb(1)/4-fixed_beta()
        require(0<delta<arb(1)/4,"invalid base floor ordering")
        iden=arb_mat(size,size,[arb(i==j) for i in range(size) for j in range(size)])
        assembled=update_head(old,gram504,delta)+(iden-gram850)/4
        errors={"P253_weighted_GL60":assembly["gl60_error_budget"]["row_sum_bound"]}
        validation=certify_center_weyl(assembled,sum(errors.values()))
        result={"dimension":size,"status":validation["status"],
                "certified_positive":validation["certified_positive"],
                "positive_pivot_count":validation["positive_pivot_count"],
                "first_inconclusive_pivot_index":validation.get("first_inconclusive_pivot_index"),
                "margin_lower_dyadic":validation.get("margin_lower_dyadic"),
                "analytic_GL60_error_exact":str(errors["P253_weighted_GL60"]),
                "negative_witness":None}
        if not result["certified_positive"]:
            result["negative_witness"]=rational_negative_witness(
                assembled,validation.get("first_inconclusive_pivot_index"),errors)
        return result


def audit(artifact, max_dimension, bits):
    params=validate_parameters()
    archive=check_p258_archive(artifact)
    require(bits>=1024,"precision too low for tiny spectral margins")
    require(29<=max_dimension<=MODES,"invalid maximal prefix")
    old_odd=_old_matrices(OLD_AGGREGATE)[1]
    records=[]
    detected=False
    for n in TARGET_PREFIXES:
        if n>max_dimension:
            break
        record=sector_prefix(old_odd,n,bits)
        records.append(record)
        if record["negative_witness"] is not None:
            detected=True
            break
        # 'inconclusive' is neither positive nor negative; diagnose separately.
        if not record["certified_positive"]:
            break
    overall="CERTIFIED_NEGATIVE_DIRECTION" if detected else (
        "CERTIFIED_ODD_HEAD_POSITIVE" if records and records[-1]["dimension"]==377
        and records[-1]["certified_positive"] else "PARTIAL_OR_INCONCLUSIVE")
    schur=None
    if overall=="CERTIFIED_ODD_HEAD_POSITIVE":
        # Even block in P259 >= P258 even block due to positive high-pass update.
        schur=_schur({"even":archive["sectors"]["even"],
                     "odd":{"validation":{"certified_positive":True,
                      "margin_lower_dyadic":records[-1]["margin_lower_dyadic"]}}})
    return {"schema":"NUMERIA.ES8P259.GATE_B.1","source_P258_sha256":SOURCE_SHA,
            "model":"four prime, floor 1/4 from 504 and 1/2 from 850",
            "proven_floor":params,"bits":bits,"checked_prefixes":records,
            "status":overall,"modal_Schur_and_local_Weil":schur,
            "finite_auxiliary_global_coercivity_proved":bool(schur),
            "Riemann_hypothesis_proved":False,
            "notes":"If a prefix is inconclusive, increase Arb bits or improve preconditioner. "
                    "A failed interval LDL pivot is not proof of negativity. "
                    "No P258 dilation bound reused after shifting the 850 cutoff."}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--p258",type=Path,required=True)
    p.add_argument("--max-dimension",type=int,default=96)
    p.add_argument("--bits",type=int,default=1536)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    report=audit(args.p258,args.max_dimension,args.bits)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print("P259 GATE B",report["status"])
    for a in report["checked_prefixes"]:
        print(a["dimension"],a["status"],"negative?",a["negative_witness"] is not None)

if __name__=="__main__":
    main()
