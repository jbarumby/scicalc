#!/usr/bin/env python3
"""Exact rational Sylvester audit from two separately calculated P92 Arb entries.

Does not infer positivity from a green CI job, printed decimals, or P89/P95.
"""
import json,sys
from fractions import Fraction
from pathlib import Path
def exact_pair(x):
    a,b=Fraction(x["lo"]),Fraction(x["hi"])
    if a>b:raise AssertionError("inverted interval")
    if not (Fraction(x["im_lo"])<=0<=Fraction(x["im_hi"])):
        raise AssertionError("imaginary interval does not contain zero")
    return a,b
def assemble(paths):
    if len(paths)!=2:raise ValueError("expected W11 and W12 results")
    found={}
    for path in paths:
        d=json.loads(Path(path).read_text(encoding="utf-8"))
        if d.get("status")!="FINITE_INTERVAL_ENTRY_CONDITIONAL_ON_FORMULA":
            raise RuntimeError("incomplete/nonfinite Arb entry "+str(path))
        n=d.get("entry_name")
        if n not in ("W11","W12") or n in found:
            raise ValueError("missing or duplicate pair")
        if d.get("epsilon_denominator")!=2048 or d.get("tolerance")!="1e-8":
            raise ValueError("incompatible experimental settings")
        if d.get("normalization")!="P92":
            raise ValueError("incompatible normalization")
        found[n]=d
    if set(found)!={"W11","W12"}:raise ValueError("missing entry")
    if found["W11"]["precision_bits"]!=found["W12"]["precision_bits"]:
        raise ValueError("different Arb precision")
    diag=exact_pair(found["W11"]["entry"]["W"])
    cross=exact_pair(found["W12"]["entry"]["W"])
    lower_diag=diag[0]
    abs_cross=max(abs(cross[0]),abs(cross[1]))
    # W11=W22 by translation invariance at the fixed test family.
    det_lower=lower_diag**2 - abs_cross**2 if lower_diag>0 else None
    eig_lower=lower_diag-abs_cross
    positive=bool(lower_diag>0 and det_lower is not None and det_lower>0 and eig_lower>0)
    return {
        "status": "INTERVAL_POSITIVE_DEFINITE_CONDITIONAL" if positive
          else "INCONCLUSIVE_FINITE_MATRIX",
        "precision_bits": found["W11"]["precision_bits"],
        "normalization":"P92",
        "symmetric_matrix":"[[W11,W12],[W12,W11]]",
        "diagonal_translation_invariance":True,
        "W11": {"lo":str(diag[0]),"hi":str(diag[1])},
        "W22": {"lo":str(diag[0]),"hi":str(diag[1])},
        "W12": {"lo":str(cross[0]),"hi":str(cross[1])},
        "lambda_lower":str(eig_lower),
        "determinant_lower":str(det_lower) if det_lower is not None else None,
        "determinant_positive":positive,
        "warning":"Only a finite-dimensional conditional result; no universal Weil positivity/RH."
    }
if __name__=="__main__":
    d=assemble(sys.argv[1:3])
    print(json.dumps(d,indent=2))
