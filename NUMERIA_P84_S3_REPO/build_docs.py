"""Build fresh, source-provenanced Doc1/Doc2 reports from the P84 Arb runs.

No user-uploaded float envelope, fitted alpha/beta/D or G grid is promoted
to a rigorous interval. The outputs explicitly mark the portions not proved.
"""
import argparse
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json

from .audit_s3 import exact, DOCS, DOC_SHA, ORIG
from NUMERIA_P81_TRANSFER.transfer import audit as audit_p81
from NUMERIA_P83_B3_ELEMENTARY.b3_proof import proof as audit_b3

BITS=("192","256","384")

def build(folder):
    runs={}
    for b in BITS:
        p=folder/f"s3-interval-{b}.json"
        d=json.loads(p.read_text(encoding="utf-8"))
        if d["precision_bits"]!=int(b) or d["status"]!="REPO_INTERVAL_CERTIFIED_CONDITIONALLY" or len(d["rows"])!=16:
            raise ArithmeticError(f"failed or incomplete Arb run: {b}")
        if d["source_sha256"]["DOC1_upload_unverified"]!=DOC_SHA[DOCS[0]] or d["source_sha256"]["DOC2_upload_unverified"]!=DOC_SHA[DOCS[1]]:
            raise ArithmeticError("document input hash mismatch")
        runs[b]=d
    if len({json.dumps(x["source_sha256"],sort_keys=True) for x in runs.values()})!=1:
        raise ArithmeticError("source changes across precisions")
    results1=[]
    results2=[]
    for i in range(16):
        rb={b:runs[b]["rows"][i] for b in BITS}
        if [r["row"] for r in rb.values()]!=[i+1]*3:
            raise ArithmeticError("mismatched row ordering")
        if len({r["B1_S2_adopted"] for r in rb.values()})!=1:
            raise ArithmeticError("adopted table changed")
        d1={"row":i+1,"b4":rb["256"]["b4_public_exact"],
            "B1_robust":rb["256"]["B1_S2_adopted"],"by_precision":{}}
        d2={"row":i+1,"by_precision":{}}
        for bits,r in rb.items():
            margin=F(r["adopted_margin_exact_rational"])
            if margin<=0 or (i>=9 and margin<F(1,1000)):
                raise ArithmeticError("weak B1 margin")
            for field in ("K1_integral_lower_dyadic","K1_integral_upper_dyadic",
                          "K1_tail_upper_dyadic","K2_upper_dyadic",
                          "K3_lower_dyadic","B1_upper_dyadic",
                          "dK3_db4_lower_dyadic","dlog_B1_db4_lower_dyadic"):
                if exact(r[field])<=0:
                    raise ArithmeticError(f"negative or missing interval endpoint: {field}")
            d1["by_precision"][bits]={
                "B1_upper_dyadic":r["B1_upper_dyadic"],
                "margin_exact_rational":r["adopted_margin_exact_rational"],
                "dK3_db4_lower_dyadic":r["dK3_db4_lower_dyadic"],
                "dlog_B1_db4_lower_dyadic":r["dlog_B1_db4_lower_dyadic"],
            }
            d2["by_precision"][bits]={
                "K1_integral_lower_dyadic":r["K1_integral_lower_dyadic"],
                "K1_integral_upper_dyadic":r["K1_integral_upper_dyadic"],
                "K1_tail_upper_dyadic":r["K1_tail_upper_dyadic"],
                "K2_upper_dyadic":r["K2_upper_dyadic"],
                "K3_lower_dyadic":r["K3_lower_dyadic"],
                "B1_upper_dyadic":r["B1_upper_dyadic"],
            }
        results1.append(d1)
        results2.append(d2)
    p81=audit_p81()
    if not p81["all_checks_pass"] or p81["total_branch_zone_checks"]!=96:
        raise ArithmeticError("P81 source transfer not confirmed")
    b3=audit_b3()
    if not b3["b3_inequality_proved_from_elementary_standard_identities"] or F(b3["positive_rational_slack"])<=0:
        raise ArithmeticError("P83 elementary b3 proof not confirmed")
    provenance=runs["256"]["source_sha256"]
    doc1={
        "schema":"NUMERIA.P84.S3.DOC1.REPO_INTERVAL.1",
        "status":"REPO_INTERVAL_CERTIFIED_CONDITIONALLY",
        "source":"user DOC1 is unverified and is NOT itself this certificate",
        "precisions_bits":[192,256,384],
        "rows":results1,
        "monotonicity_proved_algebraically_where_K3_gt1":True,
        "b4_star_root_from_user_report_recertified":False,
        "G_monotonicity_from_user_report_recertified":False,
        "first_nine_rows_P78_source_present_and_B1_checked":True,
        "full_original_Lemma6_domain_including_3e10_recertified":False,
        "provenance_sha256":provenance,
        "global_RH_proved":False
    }
    doc2={
        "schema":"NUMERIA.P84.S3.DOC2.REPO_INTERVAL.1",
        "status":"REPO_INTERVAL_CERTIFIED_CONDITIONALLY",
        "source":"user DOC2 float-precision envelopes are never used as proof endpoints",
        "precisions_bits":[192,256,384],
        "rows":results2,
        "P81_public_S_epsilon_transfer_exact_rational_gates":96,
        "P83_independent_b3_rational_slack":b3["positive_rational_slack"],
        "C3_C4_K3_B1_computed_by_P79_interval_formulas":True,
        "user_DOC2_G_monotonicity_zeta_28_checks_recertified":False,
        "external_mean_value_preissmann_theorem_reproved":False,
        "external_RZ_S_epsilon_theorems_reproved":False,
        "full_X_from_3e10_proved":False,
        "provenance_sha256":provenance,
        "global_RH_proved":False
    }
    return doc1,doc2

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input-dir",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    x=p.parse_args()
    d1,d2=build(x.input_dir)
    x.output_dir.mkdir(parents=True,exist_ok=True)
    for name,v in (("DOC1_repo_interval.json",d1),("DOC2_repo_interval.json",d2)):
        (x.output_dir/name).write_text(json.dumps(v,ensure_ascii=False,indent=2)+"\n")
        print(name,v["status"],"rows",len(v["rows"]))

if __name__=="__main__":
    main()
