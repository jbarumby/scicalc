#!/usr/bin/env python3
"""V22 arithmetic audit of E7 midpoint bounds. All decisions are exact Fraction.
Assumes source M dominates the true second derivative; does NOT certify that fact.
"""
import argparse,hashlib,json,zipfile
from fractions import Fraction as F
from pathlib import Path
E7_SHA="4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7"
TICK=F(1,10**45)
SOURCES={"G11":"corrigidos/RESULTADO_B3_G11_E2_CORRIGIDO.json","G22":"corrigidos/RESULTADO_B3_G22_E2_CORRIGIDO.json"}
def number(d):return F(int(d["numerator"]),int(d["denominator"]))
def round_up(x):
    n=x/TICK
    return F(-(-n.numerator//n.denominator))*TICK
def component(j,which):
    bands=j["mesh"];prev=F(0);old=F(0);new=F(0);rows=[]
    for b in bands:
        a,z=F(b["from"]),F(b["to"]);n=b["cells"]
        m=number(b["f2_absolute_bound"] if which=="G11" else b["f2_bound"])
        old_error=number(b["midpoint_error_bound"] if which=="G11" else b["error_bound"])
        expected=m*(z-a)**3/(24*n*n)
        assert a==prev and a<z and n>0 and m>=0 and abs(old_error-expected)<=TICK
        fixed=max(old_error,round_up(expected))
        assert fixed>=expected and fixed>=old_error
        rows.append({"from":str(a),"to":str(z),"cells":n,"M_declared":str(m),
            "midpoint_formula_exact":str(expected),"historical":str(old_error),
            "corrected":str(fixed),"historical_enclosed":old_error>=expected,
            "added_slack":str(fixed-old_error)})
        old+=old_error;new+=fixed;prev=z
    assert prev==160
    key="quadrature_absolute_error_upper" if which=="G11" else "G22_quadrature_error_upper"
    assert old==number(j[key])
    return {"component":which,"bands":rows,"count":len(rows),
        "points":sum(x["cells"] for x in rows),"old_total":str(old),
        "correct_total":str(new),"correction":str(new-old),
        "historical_underbounds":sum(not x["historical_enclosed"] for x in rows)}
def audit(path):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==E7_SHA,"E7 SHA mismatch"
    with zipfile.ZipFile(path) as z:
        c={name:component(json.loads(z.read(source)),name) for name,source in SOURCES.items()}
        b5=json.loads(z.read("b5/RESULTADO_B5.json"))
    assert sum(x["points"] for x in b5["interval_results"])==15754
    assert c["G11"]["points"]==9524 and c["G22"]["points"]==6230
    assert c["G11"]["count"]==6 and c["G22"]["count"]==5
    assert c["G11"]["historical_underbounds"]==0 and c["G22"]["historical_underbounds"]==1
    assert c["G22"]["correction"]==str(TICK)
    assert [(a["from"],a["to"]) for a in c["G22"]["bands"] if not a["historical_enclosed"]]==[("96","160")]
    return {"schema":"NUMERIA.V22.MIDPOINT_ARITHMETIC.1",
      "status":"G22_ONE_TICK_OUTWARD_ROUNDING_CORRECTED_CONDITIONAL",
      "e7_sha256":E7_SHA,"error_formula":"M*(b-a)^3/(24*n^2)",
      "unit_tick":"1/10^45","b5_records":15754,"components":c,
      "not_certified":["M >= true sup|f''| independently for each band",
                       "independent verification of B5 nodal transcendental enclosures",
                       "global Weil positivity or the Riemann Hypothesis"]}
def main():
    p=argparse.ArgumentParser();p.add_argument("--e7",required=True,type=Path);p.add_argument("--out",required=True,type=Path);a=p.parse_args()
    try:r=audit(a.e7)
    except Exception as e:print("FAIL",type(e).__name__,e);return 2
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("PASS: 11 B3 bands, 15754 B5 records, G22 last band corrected by one upward tick")
    return 0
if __name__=="__main__":raise SystemExit(main())
