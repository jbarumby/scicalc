#!/usr/bin/env python3
"""Pure-rational, fail-closed audit for one ES8-P99 direct Arb 2x2 output."""
import json,sys
from fractions import Fraction as F

def evaluate(path,epsilon,bits):
    j=json.load(open(path,encoding="utf-8"))
    status=j.get("status")
    if status not in ("INCONCLUSIVE_DETERMINANT_NOT_POSITIVE","INTERVAL_POSITIVE_DEFINITE_CONDITIONAL"):
        raise RuntimeError(f"No usable finite Arb matrix: {status}: {j.get('reason')}")
    if j.get("epsilon_denominator") != epsilon or j.get("precision_bits") != bits or j.get("centers") != [1,2]:
        raise AssertionError("Metadata mismatch")
    entries={}
    for k in ("W11","W22","W12"):
        raw=j["entries"][k]["W"]
        lo,hi=F(raw["lo"]),F(raw["hi"])
        if not lo<=hi:raise AssertionError(f"Inverted {k}")
        if not (F(raw["im_lo"]) <= 0 <= F(raw["im_hi"])):
            raise AssertionError(f"Nonreal {k}")
        entries[k]={"lo":str(lo),"hi":str(hi),"width":str(hi-lo),
                    "lo_float":float(lo),"hi_float":float(hi),"width_float":float(hi-lo)}
    l11,l22=F(entries["W11"]["lo"]),F(entries["W22"]["lo"])
    a12=max(abs(F(entries["W12"]["lo"])),abs(F(entries["W12"]["hi"])))
    diag_positive=l11>0 and l22>0
    det_lower=l11*l22-a12*a12 if diag_positive else None
    positive=bool(diag_positive and det_lower>0)
    if positive != (status=="INTERVAL_POSITIVE_DEFINITE_CONDITIONAL"):
        raise AssertionError("Wrong decision")
    return {"epsilon_denominator":epsilon,"precision_bits":bits,
            "tolerance":j["tolerance"],"status":status,
            "entries":entries,"lower_diagonal_both_positive":diag_positive,
            "determinant_lower":str(det_lower) if det_lower is not None else None,
            "determinant_lower_positive":positive,
            "elapsed_seconds":j.get("elapsed_seconds")}
if __name__=="__main__":
    fn,epsilon,bits=sys.argv[1],int(sys.argv[2]),int(sys.argv[3])
    print(json.dumps(evaluate(fn,epsilon,bits),indent=2,ensure_ascii=False))
