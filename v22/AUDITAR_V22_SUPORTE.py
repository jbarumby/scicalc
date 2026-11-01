#!/usr/bin/env python3
"""Independent exact support/smoothness warning for Burnol vs NUMERIA."""
from fractions import Fraction as F
import argparse,json
from pathlib import Path
L=F(9,20)
def enclosure(x,N=80):
 t=F(x-1,x+1)
 s=sum((2*t**(2*k+1)/F(2*k+1) for k in range(N)),F(0))
 return s,s+2*t**(2*N+1)/(F(2*N+1)*(1-t*t))
def audit():
 l2,h2=enclosure(2);l3,h3=enclosure(3)
 assert L<l2<h2<2*L<l3<h3
 return {"schema":"NUMERIA.V22.BURNOL_SUPPORT.1",
  "status":"PRIME_2_PRESENT_AND_SMOOTHNESS_OBLIGATION_OPEN",
  "L":"9/20","double_L":"9/10",
  "log2_enclosure":[str(l2),str(h2)],"log3_enclosure":[str(l3),str(h3)],
  "prime_powers_with_log_n_le_2L":[2],
  "sharp_sine_C_infinity":False,
  "explanation":"Truncated sin modes have nonzero one-sided derivatives at endpoints; Burnol's compact smooth-test hypothesis not automatically satisfied.",
  "not_certified":["Burnol's unspecified sufficiently-small support includes L=0.45",
                   "Extension of Burnol positivity beyond its support hypothesis",
                   "Approximation by smooth tests with an explicit convergence/error bound"]}
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
 j=audit();a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(j,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
 print("PASS: log 2<2L<log 3, prime 2 present, sine cutoff is not C-infinity");return 0
if __name__=="__main__":raise SystemExit(main())
