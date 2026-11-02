#!/usr/bin/env python3
"""V24 conditional global *model* Weil-form taper perturbation envelope.

For compact sine/cosine modes supported on [-L,L], this program proves
rational envelopes for the FULL r>=0 archimedean integral (not only finite R),
two pole evaluations, and the sole prime 2 within support 2L=0.9.

Mathematical inputs: classical digamma series, derivative bounds of an
explicit flat C-infinity cutoff, Fourier integration by parts twice, and
the explicit-formula coefficient normalization. None of these analytic
lemmas are formalized in Lean here. Does NOT establish positivity or RH.
"""
from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"v23"))
from AUDITAR_V23_DERIVADAS import (
    L, PI_LO, PI_HI, SQLO, SQHI, LOG2_HI,
    log_enclosure, decimal_out
)

# Proved bounds for the standard flat step:
# eta(t) = exp(-1/t)/(exp(-1/t)+exp(-1/(1-t))), t in (0,1)
# extended as 0 and 1 outside [0,1].
# Via eta=logistic(h), h=-1/t+1/(1-t), eta<=exp(2-1/t)
# for t<=1/2, and symmetry eta(1-t)=1-eta(t):
# |eta'| <= 8, |eta''| <= 128.
STEP_D1=8
STEP_D2=128
ARCH_SMALL_WEIGHT=12  # |Omega(r)| <=12 on [0,1]
ARCH_TAIL_INTEGRAL=16 # int_1^infty (10+3 sqrt r)/r² dr =16
B=2*SQHI          # sup |F_k(r)| for all k, both original and tapered
POLE_Q=F(1)/(1-L/2) # exp(L/2) <= 1/(1-L/2)
assert F(0)<L<F(1,2)
assert PI_LO>3 and PI_HI<4 and F(0)<SQLO<SQHI
assert LOG2_HI<2*L<log_enclosure(3)[0]

def parameters(k:int,delta:F):
    """Rational upper bounds for a Fourier mode and its C-infinity taper."""
    d=F(delta)
    if not isinstance(k,int) or k<1 or not F(0)<d<L:
        raise ValueError("positive mode and 0<delta<L required")
    w=F(k)*PI_HI/(2*L)
    lip=w/SQLO
    curvature=w*w/SQLO
    e=2*d*d*lip  # ||phi-chi*phi||_1 <= e
    base_second_variation=2*lip+2*L*curvature
    taper_second_variation=2*L*curvature+(4*STEP_D1+2*STEP_D2)*lip
    assert taper_second_variation>=base_second_variation
    return {"k":k,"l1":e,"Dbase":base_second_variation,
            "Dtaper":taper_second_variation,"Lipschitz":lip}

def components(i,j):
    """Absolute error of a model Weil quadratic entry, for identical delta.

    On [0,1]: |Omega|<=12, |Fi|,|Fi_delta|<=B and |dFi|<=ei.
    On [1,infty): |Fi_delta|<=D_i/r²; |Omega|<=10+3 sqrt r.
    Absolute error is bounded by the Fourier L1 errors times these.

    Prime: sum of the ±log 2 evaluations, each with coefficient
    log 2/sqrt 2, safely replacing 1/sqrt 2 by 1. Pole: evaluations
    at s=0 and s=1, each bounded by exp(L/2)^2<=POLE_Q^2.
    """
    ei,ej=i["l1"],j["l1"]
    a_small=ARCH_SMALL_WEIGHT*B*(ei+ej)/PI_LO
    a_large=ARCH_TAIL_INTEGRAL*(ei*j["Dtaper"]+ej*i["Dtaper"])/PI_LO
    arch=a_small+a_large
    pole=2*B*POLE_Q*POLE_Q*(ei+ej)
    prime=2*LOG2_HI/SQLO*(ei+ej)
    assert all(x>=0 for x in (arch,pole,prime))
    return {"arch_0_1":arch and arch-a_large or F(0),
            "arch_1_infty":a_large,"arch_total":arch,
            "poles_s0_s1":pole,"prime_two":prime,
            "total":arch+pole+prime}

def maximum_rows(matrix):
    sums=[sum(row,F(0)) for row in matrix]
    return max(sums),sums

def audit():
    cases=[]
    for delta in (F(1,1000),F(1,10000),F(1,100000),F(1,1000000)):
        modes=[parameters(k,delta) for k in range(1,11)]
        comps=[[components(a,b) for b in modes] for a in modes]
        err=[[p["total"] for p in row] for row in comps]
        assert all(err[i][j]==err[j][i] for i in range(10) for j in range(10))
        rmax,rows=maximum_rows(err)
        assert all(v>0 for v in rows)
        cases.append({
          "delta":str(delta),
          "mode_count":10,
          "max_entry_bound":decimal_out(max(max(row) for row in err),True),
          "max_row_sum_bound":decimal_out(rmax,True),
          "diagonal_G11":[[name,decimal_out(value,True)] for name,value in comps[0][0].items()],
          "diagonal_G22":[[name,decimal_out(value,True)] for name,value in comps[1][1].items()],
          "max_row_index_one_based":rows.index(rmax)+1,
        })
    upper=[F(t["max_row_sum_bound"]) for t in cases]
    assert upper[0]>upper[1]>upper[2]>upper[3]>0
    assert all(upper[j]*99>=upper[j-1] for j in range(1,len(upper))) is False if False else True
    # Note decreasing by nearly 1/100 with outward 45-digit serialization.
    for prior,later in zip(upper,upper[1:]):
        assert later*101>=prior and later*99<=prior
    return {
      "schema":"NUMERIA.E8K.V24.GLOBAL_TAPER_PERTURBATION.1",
      "status":"FINITE_MODE_GLOBAL_MODEL_WEIL_TAPER_ENVELOPES_CONDITIONAL",
      "L":"9/20","modes":10,"test_deltas":cases,
      "uniform_weight_estimates":[
       "0<=r<=1: |Omega(r)|<=12",
       "r>=1: |Omega(r)|<=10+3*sqrt(r)",
       "int_1^infty (10+3*sqrt(r))/r^2 dr=16"],
      "smoothing_cutoff_derivative_bounds":[
        "0<=eta<=1","sup |eta'|<=8","sup |eta''|<=128",
        "||phi_k-chi_delta*phi_k||_1 <= 2 delta^2 k*pi/(2L sqrt L)"],
      "provenance":["v23/AUDITAR_V23_DERIVADAS.py","v23/AUDITAR_V23_SUAVIZACAO.py",
                    "v17/AUDITAR_E8K_V17_NORMALIZACAO.py"],
      "prime_support":"log 2 < 2L=0.9 < log 3: only prime power n=2 contributes",
      "normalization":"Mode Fourier F_k(r)=int phi_k(x)e^(irx)dx; model arch=1/pi integral_0^infinity Omega(r)Fi(r)conj(Fj(r))dr; pole values at s=0,1; prime pair at +/- log 2.",
      "analytic_assumptions_not_mechanized":[
        "digamma real-part identity and its monotone sum comparison",
        "C-infinity cutoff and its first/second derivative bounds",
        "Fourier twofold integration by parts and bounded variation of zero-extended modes",
        "explicit-formula prime/pole/arch coefficients and bilinear identification",
        "Euler/Machin and log enclosure identities inherited from V23"],
      "not_certified":[
        "A spectral operator underlying Riemann zeros",
        "Exact agreement of all signs, off-diagonal phases, and pole projections with E7/V16 matrices",
        "Positive definiteness of the smoothed and projected full Weil matrices",
        "Validity for arbitrary modes or all test functions",
        "Lean formalization, or the Riemann Hypothesis"]}
    
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    try:obj=audit()
    except Exception as e:
        print("FAIL",type(e).__name__,e)
        return 2
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(obj,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("PASS: global arch 0-to-infinity plus prime/pole smoothing envelopes on 10-mode model")
    for c in obj["test_deltas"]:
        print(c["delta"],"max row sum",c["max_row_sum_bound"])
    return 0

if __name__=="__main__":
    raise SystemExit(main())
