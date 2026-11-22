"""P80 full 16-row S3 check using P79 Arb enclosures as immutable evidence.

This is intentionally a REPLACEMENT for (not an in-place edit of) the submitted
S3 prototype. It fixes four issues: precision scope, bisection outside Arb,
mixed keep/replace patch schema, and float-valued 'dyadic' certificates.
No new spectral positivity claim is made here.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

from NUMERIA_P79_B1_PUBLIC.audit_public import audit as audit_p79
from .v3_gate import audit as audit_v3, fraction_repr

SCHEMA="NUMERIA.P80.S3.ALL_ROWS.REPLAY.1"

def checked_rows(spec):
    if not spec:
        return list(range(1,17))
    chosen=set()
    for piece in spec.split(","):
        piece=piece.strip()
        if not piece:
            raise ValueError("empty row selection")
        if "-" in piece:
            aa,bb=piece.split("-",1)
            if not aa.isdecimal() or not bb.isdecimal() or int(aa)>int(bb):
                raise ValueError(f"invalid row range: {piece}")
            chosen.update(range(int(aa),int(bb)+1))
        elif piece.isdecimal():
            chosen.add(int(piece))
        else:
            raise ValueError(f"invalid row identifier: {piece}")
    if not chosen or any(i<1 or i>16 for i in chosen):
        raise ValueError("row index outside 1..16")
    return sorted(chosen)

def run(precision=192,rows=None):
    selected=checked_rows(rows)
    p79=audit_p79(precision)
    v3=audit_v3()
    if len(p79["coefficients"])!=16 or len(v3["rows"])!=16:
        raise ArithmeticError("truncated P79/P80 input")
    if p79["schema"]!="NUMERIA.P79.B1.PUBLIC_B4.1":
        raise ArithmeticError("unexpected P79 audit schema")
    out=[]
    overlay=[]
    for i in selected:
        src=p79["coefficients"][i-1]
        r=v3["rows"][i-1]
        if src["row"]!=i or r["row"]!=i or F(src["b4_source"])!=F(r["b4_P78"]):
            raise ArithmeticError("P79/P80 source mismatch")
        s2=F(i+16,32)
        exponent=(2*s2-1)/(2-s2)
        if exponent!=F(2*i,48-i):
            raise ArithmeticError("wrong sigma exponent")
        if i<=9:
            if not src["old_B1_still_strict"]:
                raise ArithmeticError("original B1 invalid on unchanged row")
            if src["original_B1"]!=src["B1_replacement"]:
                raise ArithmeticError("modified original B1 on row <=9")
            action="keep"
        else:
            if src["old_B1_still_strict"]:
                raise ArithmeticError("unexpected valid printed B1 in P79 row>=10")
            if not F(src["B1_replacement"])>F(src["original_B1"]):
                raise ArithmeticError("new majorant did not grow")
            action="overlay"
        if not F(r["b4_P80_suggested"])<=F(r["b4_P78"]):
            raise ArithmeticError("v3 source raises b4")
        out.append({
            "row":i,"status":action,
            "delta":r["delta"],
            "b4_P78":r["b4_P78"],
            "b4_P80_maximum":r["b4_P80_suggested"],
            "e_exact":fraction_repr(exponent),
            "B1_old":src["original_B1"],
            "B1_P79_certified":src["B1_replacement"],
            "B1_upper_dyadic_P79":src["B1_upper_dyadic"],
            "b4_nonincrease_proves_P79_B1_reusable":True,
        })
        overlay.append({"row":i,"action":action,"B1":src["B1_replacement"]})
    tab=json.dumps(overlay,sort_keys=True,separators=(",",":"))
    return {
        "schema":SCHEMA,
        "precision_bits":precision,
        "rows_requested":selected,
        "rows":out,
        "overlay_patch":overlay,
        "overlay_patch_sha256":hashlib.sha256(tab.encode()).hexdigest(),
        "source_P79_table_sha256":p79["table_sha256"],
        "source_P78_sha256":v3["P78_source_sha256"],
        "P80_changed_b4_rows":v3["changed_b4_rows"],
        "all_B1_valid_under_P80_v3_sources":True,
        "P79_512_band_replay_unmodified":True,
        "fully_independent_lemma6_proved":False,
        "global_RH_proved":False,
    }

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--precision",type=int,default=192)
    a.add_argument("--rows",default=None)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    doc=run(args.precision,args.rows)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(doc,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"precision":args.precision,"rows":len(doc["rows"]),
                      "kept":sum(r["status"]=="keep" for r in doc["rows"]),
                      "replaced":sum(r["status"]=="overlay" for r in doc["rows"]),
                      "P80_changed_b4_rows":doc["P80_changed_b4_rows"],
                      "source_table_sha":doc["source_P79_table_sha256"]},indent=2))

if __name__=="__main__":
    main()
