"""NUMERIA P84: independently replay S2 table with P79 Arb quadrature.

Explicit computational certification only, conditional on original P77
quadrature method and P78/P80/P81 external analytic inputs.  The uploaded
S3 documents are UNVERIFIED claims, even if self-labelled REPO_CERTIFIED.
No binary float or ad-hoc relative enclosure is accepted in proof gates.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

from flint import arb, ctx
from NUMERIA_P76_COEFFICIENTS.table import (
    TABLE, TABLE_CS_V2_PRINTED, TABLE_S2_PUBLIC_ROBUST,
    S2_ROBUST_B1, validate_table_s2_robust,
)
from NUMERIA_P75_CS2.cs_v2 import TABLE_V2_S2_PUBLIC_ROBUST
from NUMERIA_P79_B1_PUBLIC import audit_public as p79
from NUMERIA_P81_TRANSFER.transfer import audit as p81_audit
from NUMERIA_P83_B3_ELEMENTARY.b3_proof import proof as b3_proof

ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
PATCH=ROOT/"NUMERIA_P79_B1_PUBLIC"/"prototype_input"/"p78_s2_table_patch.json"
ORIG=HERE/"input_unverified"
DOCS=("p78_s3_doc1_certificado.json","p78_s3_doc2_certificado.json")
DOC_SHA={
 "p78_s3_doc1_certificado.json":"093dd70b4fc6ddf0385094cd112699034f2fad2c36a3c51758c7f7c044652137",
 "p78_s3_doc2_certificado.json":"ac356d2997e93945e8ce0751dd9e2f00054b45b9284a800d3ed690b57d8c6846",
}
SCHEMA="NUMERIA.P84.S3.ADOPTED.TABLE.INTERVAL.1"
ALLOWED_BITS=(192,256,384)

def exact(pair):
    if not isinstance(pair,(list,tuple)) or len(pair)!=2:
        raise ArithmeticError("missing or invalid dyadic endpoint")
    m,e=pair
    if not isinstance(m,str) or not isinstance(e,int):
        raise ArithmeticError("dyadic endpoint must have string mantissa and integer exponent")
    if abs(e)>100000:
        raise ArithmeticError("dyadic exponent outside guard")
    return F(int(m))*(F(2)**e)

def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def strict_source_checks():
    if not validate_table_s2_robust():
        raise ArithmeticError("versioned S2 TABLE invalid")
    if TABLE!=TABLE_CS_V2_PRINTED:
        raise ArithmeticError("historical TABLE overwritten")
    p78,p78_sha = p79.load_source()
    patch=json.loads(PATCH.read_text(encoding="utf-8"))
    if patch.get("schema")!="NUMERIA.P78.S2.TABLE.PATCH":
        raise ArithmeticError("invalid patch schema")
    if [x.get("row") for x in patch.get("rows",[])]!=list(range(10,17)):
        raise ArithmeticError("patch rows not exactly 10..16")
    for row in patch["rows"]:
        i=row["row"]
        o,a=TABLE[i-1],TABLE_S2_PUBLIC_ROBUST[i-1]
        r=p78[i-1]
        if F(row["delta"])!=F(o.delta) or F(row["b4_public"])!=F(r["b4"]):
            raise ArithmeticError(f"P78 b4/delta mismatch row {i}")
        if row["b1_old"]!=o.b1 or row["b1_new_robust"]!=a.b1:
            raise ArithmeticError(f"non-matching S2 patch row {i}")
        if a.b1!=S2_ROBUST_B1[i]:
            raise ArithmeticError(f"S2 map mismatch row {i}")
    for name in DOCS:
        if not (ORIG/name).is_file() or file_hash(ORIG/name)!=DOC_SHA[name]:
            raise ArithmeticError(f"original user-document integrity mismatch: {name}")
    if tuple((r.b1,r.b2,r.b3) for r in TABLE_S2_PUBLIC_ROBUST)!=TABLE_V2_S2_PUBLIC_ROBUST:
        raise ArithmeticError("P75 consumer overlay differs from adopted TABLE")
    return p78,p78_sha,patch

def reconstruct(bits):
    if bits not in ALLOWED_BITS:
        raise ValueError(f"supported exact workprec values: {ALLOWED_BITS}")
    p78,p78_sha,patch=strict_source_checks()
    rows=[]
    with ctx.workprec(bits):
        cs=p79.constants()
        L,pi,m0=cs["logh"],cs["pi"],cs["m0"]
        cs["c1"]=6/pi**2+p79.q("1.048")/L
        cs["c2"]=pi*m0*p79.q(".62")/L+6*m0/(pi*p79.T0)+pi*m0*p79.q("1.048")/(arb(p79.T0)*L)
        for i,(old,row) in enumerate(zip(TABLE,TABLE_S2_PUBLIC_ROBUST),1):
            b4=F(p78[i-1]["b4"])
            calc=p79.row_audit(row,cs,p79.exact_ball(b4))
            upper=exact(calc["B1_upper_dyadic"])
            margin=F(row.b1)-upper
            if margin<=0 or not calc["strict"]:
                raise ArithmeticError(f"FAIL adopted TABLE at row {i}, bits {bits}")
            if i<=9 and old.b1!=row.b1:
                raise ArithmeticError("cannot alter first nine printed B1")
            if i>=10:
                if not F(old.b1)<upper:
                    raise ArithmeticError(f"original printed B1 unexpectedly passes row {i}")
                if margin<F(1,1000):
                    raise ArithmeticError(f"robust S2 margin < 1e-3: row {i}")
            delta=p79.q(row.delta)
            kap=p79.q(row.kappa)
            A0=(arb(10).log()*p79.q(row.log10_a0)).exp()
            c=3/(4*kap)
            # K2 is an analytic expression, NOT a second quadrature.
            two13=p79.pw(arb(2),p79.rat(1,3))
            k2=two13*c.sqrt()*pi.sqrt()/2 + two13*(4*kap/3)*p79.pw(A0/arb(p79.T0),arb(3))
            # Algebraic derivative of K3 with respect to b4:
            coupling=pi.sqrt()/(2*(2*kap).sqrt())+pi*m0/p79.T0
            alpha=pi*m0/(2*delta)*(2*delta*arb.const_euler()/L).exp() + coupling/(5*delta*delta.exp())
            beta=coupling*p79.q(".605")*(-2*delta).exp()/L
            k3=alpha*p79.exact_ball(b4)+beta
            s2=p79.rat(row.sigma_start+1,32)
            exponent=(2*s2-1)/(2-s2)
            derivative=exponent*alpha/k3
            if not (k2>0 and alpha>0 and beta>0 and k3>1 and exponent>0 and derivative>0):
                raise ArithmeticError(f"monotonicity structural gate failed: row {i}")
            rows.append({
                "row":i, "sigma":calc["sigma"],"delta":row.delta,
                "b4_public_exact":str(b4),
                "B1_CS_original":old.b1,
                "B1_S2_adopted":row.b1,
                "B1_upper_dyadic":calc["B1_upper_dyadic"],
                "adopted_margin_exact_rational":str(margin),
                "adopted_margin_positive":True,
                "old_B1_fails":bool(i>=10),
                "margin_lower_dyadic":calc["margin_lower_dyadic"],
                "K1_integral_lower_dyadic":calc["k1_integral_lower_dyadic"],
                "K1_integral_upper_dyadic":calc["k1_integral_upper_dyadic"],
                "K1_tail_upper_dyadic":calc["k1_tail_upper_dyadic"],
                "K2_upper_dyadic":p79.dyadic(k2.upper()),
                "K3_lower_dyadic":p79.dyadic(k3.lower()),
                "dK3_db4_lower_dyadic":p79.dyadic(alpha.lower()),
                "dlog_B1_db4_lower_dyadic":p79.dyadic(derivative.lower()),
                "symbolic_monotonicity_domain":"k3>1, e=(2sigma2-1)/(2-sigma2)>0; remaining K3 negative-power maxima saturate at 1",
            })
    if len(rows)!=16:
        raise ArithmeticError("incomplete 16-row audit")
    return {
        "schema":SCHEMA,
        "precision_bits":bits,
        "table_variant":"S2_ROBUST_EXPLICIT",
        "scope":"X>=3e12, external Ramaré–Zúñiga v1/v3 and KLN/MV mean-value inputs imported",
        "rows":rows,
        "source_sha256":{
            "P78_certificate":p78_sha,
            "S2_patch":file_hash(PATCH),
            "DOC1_upload_unverified":file_hash(ORIG/DOCS[0]),
            "DOC2_upload_unverified":file_hash(ORIG/DOCS[1]),
            "P77_method":file_hash(ROOT/"NUMERIA_P77_B1"/"audit_b1.py"),
            "P79_method":file_hash(ROOT/"NUMERIA_P79_B1_PUBLIC"/"audit_public.py"),
            "TABLE_module":file_hash(ROOT/"NUMERIA_P76_COEFFICIENTS"/"table.py"),
            "P75_cs_v2_consumer":file_hash(ROOT/"NUMERIA_P75_CS2"/"cs_v2.py"),
        },
        "doc1_reconstructed":"B1 exact Arb endpoints, structural positive derivative of K3 and B1",
        "doc1_b4star_root_not_recertified":True,
        "doc2_reconstructed":"P77 formula K1 interval, K2 formula, C3/C4 K3, B1; external source transfer P81",
        "doc2_original_G_monotonicity_and_28_zeta_claims_not_recertified":True,
        "analytic_external_S_epsilon_reproved":False,
        "analytic_external_MV_reproved":False,
        "X_3e10_to_3e12_covered":False,
        "S2_first_nine_B1_checked_by_Arb":True,
        "status":"REPO_INTERVAL_CERTIFIED_CONDITIONALLY",
        "global_RH_proved":False,
    }

def replay_suzuki():
    from NUMERIA_P75_CS2 import cs_v2
    from NUMERIA_P75_CS2.certify import run
    original=cs_v2.TABLE_V2
    try:
        cs_v2.TABLE_V2=TABLE_V2_S2_PUBLIC_ROBUST
        result=run(start="18",end="79",count=500,prec=256,segments=512,mesh="1/64")
    finally:
        cs_v2.TABLE_V2=original
    result["source_table"]="TABLE_V2_S2_PUBLIC_ROBUST (not P79 six-decimal table)"
    result["analytic_status"]="conditional_on_imported_public_bounds_and_mean_value_theorem"
    result["global_RH_proved"]=False
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits",type=int,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--replay-output",type=Path)
    a=parser.parse_args()
    res=reconstruct(a.bits)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(res,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"bits":a.bits,"rows":len(res["rows"]),"status":res["status"],
        "minimum_margin":str(min(F(x["adopted_margin_exact_rational"]) for x in res["rows"][9:]))},indent=2))
    if a.replay_output:
        if a.bits!=256:
            raise SystemExit("--replay-output only supported for 256-bit run")
        d=replay_suzuki()
        a.replay_output.parent.mkdir(parents=True,exist_ok=True)
        a.replay_output.write_text(json.dumps(d,ensure_ascii=False,indent=2)+"\n")
        print(json.dumps({"replay_status":d["status"],"end":d["certified_through_exact"],"cells":d["cell_count"]},indent=2))

if __name__=="__main__":
    main()
