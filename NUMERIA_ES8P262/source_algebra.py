"""Exact rational audit of the KPS AIF (4.8),(4.17),(4.18) consistency.

Reads an independently verified Arb-moment JSON, then works entirely in
fractions. Detects a discrepancy *as printed*: it does NOT infer RH status,
and does NOT claim the intended source formula has been recovered.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from math import factorial, comb
from hashlib import sha256
from pathlib import Path
import json

SCHEMA = "NUMERIA.ES8P262.KPS_THETA_BERNSTEIN_SIGN_GATE.1"

def need(ok,msg):
    if not ok: raise ArithmeticError(msg)

def decode_dyadic(pair):
    """Convert signed binary [mantissa, exponent] to exact Fraction."""
    m,e=pair
    return Q(int(m))*Q(2)**int(e)

def interval(lo, hi):
    l,h=Q(lo),Q(hi)
    need(l<=h,"reversed enclosure")
    return (l,h)

def add(a,b):
    return (a[0]+b[0],a[1]+b[1])

def subtract(a,b):
    return (a[0]-b[1],a[1]-b[0])

def product(a,b):
    p=[x*y for x in a for y in b]
    return min(p),max(p)

def quotient(a,b):
    need(not b[0]<=0<=b[1],"division by an enclosure containing zero")
    inv=(Q(1,b[1]),Q(1,b[0]))
    return product(a,inv)

def scale(a,value):
    return product(a,(Q(value),Q(value)))

def from_dyadic_box(raw):
    return interval(*[decode_dyadic(x) for x in raw])

def eq417_rhs(n, moments):
    """Literal RHS of printed (4.17), n>=1:
    n!/(2n)! * [32*binom(2n,2)*F(2n-2)-F(2n)]/2^(2n-1).
    """
    need(n>=1,"n>=1 only; F(-2) undefined for n=0")
    fn=scale(subtract(scale(moments[2*n-2],32*comb(2*n,2)),
                      moments[2*n]),Q(factorial(n),factorial(2*n)*2**(2*n-1)))
    return fn

def verify(path,expected_sha=None):
    raw=Path(path).read_bytes()
    digest=sha256(raw).hexdigest()
    if expected_sha is not None:
        need(digest==expected_sha,"P262 moment SHA mismatch")
    data=json.loads(raw)
    need(data.get("schema")==SCHEMA,"wrong moment source schema")
    need(data.get("precision_bits",0)>=128,"insufficient source precision")
    moments={int(k):from_dyadic_box(v["arb_enclosure_dyadic"])
             for k,v in data["moments"].items()}
    need(set((0,2,4))<=set(moments),"moment sources F0,F2,F4 missing")
    need(all(lo>0 for lo,hi in moments.values()),"moments not positive")
    g2=subtract(scale(quotient(moments[0],moments[2]),64), (Q(1),Q(1)))
    g4=subtract(scale(quotient(moments[2],moments[4]),384), (Q(1),Q(1)))
    need(g2[0]>0 and g4[0]>0,"G(2),G(4) not strictly positive")
    phi_printed=quotient(scale(g2,Q(-1,16)),g4)
    need(phi_printed[1]<0, "printed Eq4.18 not demonstrably negative")
    a1=eq417_rhs(1,moments)
    a2=eq417_rhs(2,moments)
    need(a1[0]>0 and a2[0]>0,"printed coefficients not positive")
    phi_recurrence=quotient(a1,a2)
    need(phi_recurrence[0]>0, "the W recurrence did not yield positive ratio")
    need(phi_printed[1]<phi_recurrence[0],
         "printed formulas might overlap; mismatch not certified")
    return {
        "schema":"NUMERIA.ES8P262.PRINTED_SOURCE_CONSISTENCY.1",
        "moment_source_sha256":digest,
        "paper_pdf":"https://www.numdam.org/item/10.5802/aif.3600.pdf",
        "formula_pages":{"J_and_W_definition":395,"theta_coefficients_and_phi":397},
        "exact_arithmetic":"fractions.Fraction on outward binary Arb moment endpoints",
        "G2_interval":[str(v) for v in g2],
        "G4_interval":[str(v) for v in g4],
        "paper_eq418_phi2_interval":[str(v) for v in phi_printed],
        "paper_eq417_a1_interval":[str(v) for v in a1],
        "paper_eq417_a2_interval":[str(v) for v in a2],
        "W_recurrence_phi2_interval":[str(v) for v in phi_recurrence],
        "phi2_intervals_disjoint":True,
        "printed_eq417_and_eq418_cannot_both_follow_from_W_recurrence_for_these_moments":True,
        "source_ambiguity":"check signed Taylor coefficients, xi normalization and any errata before attributing an error",
        "intended_Bernstein_candidate_ruled_out":False,
        "Riemann_hypothesis_proved_or_disproved":False,
        "Weil_positivity_proved_or_disproved":False,
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--moments",type=Path,required=True)
    parser.add_argument("--expected-sha")
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    result=verify(args.moments,args.expected_sha)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("P262 source algebra: PRINTED phi2 DISJOINT from W recurrence")
if __name__=="__main__":
    main()
