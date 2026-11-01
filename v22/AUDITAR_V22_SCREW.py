#!/usr/bin/env python3
"""Exact finite PSD test of a toy screw-function kernel, not the RH-specific one."""
import argparse,json
from itertools import combinations
from fractions import Fraction as F
from pathlib import Path
def phase(n):return ((1,0),(0,1),(-1,0),(0,-1))[n%4]
def det(m):
 n=len(m)
 if n==0:return F(1)
 if n==1:return F(m[0][0])
 return sum((F(-1 if j%2 else 1)*m[0][j]*det([r[:j]+r[j+1:] for r in m[1:]]) for j in range(n)),F(0))
def kernel(times,freq,weights):
 g=[]
 for t in times:
  row=[]
  for u in times:
   val=F(0)
   for w,a in zip(freq,weights):
    ct,st=phase(w*t);cu,su=phase(w*u)
    val+=a*((ct-1)*(cu-1)+st*su)
   row.append(val)
  g.append(row)
 return g
def audit():
 t=[0,1,2,3,4];g=kernel(t,[1,2,3],[F(2),F(3),F(5)])
 minors=[]
 for m in range(1,len(t)+1):
  for ids in combinations(range(len(t)),m):
   d=det([[g[i][j] for j in ids] for i in ids])
   assert d>=0;minors.append(d)
 assert len(minors)==31 and kernel([0,1],[1],[F(-1)])[1][1]==-2
 return {"schema":"NUMERIA.V22.SCREW_FINITE_TOY.1",
 "status":"FINITE_TOY_SCREW_GRAM_PSD_EXACT_NOT_ZETA",
 "grid":"times multiples of π/2","frequencies":[1,2,3],"weights":[2,3,5],
 "principal_minors":31,"negative_weight_detected":True,
 "not_certified":["The specific zeta screw function of Matsumoto-Suzuki (2026)",
                  "Positivity for arbitrary real times and arbitrary dimension",
                  "Riemann Hypothesis"]}
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
 j=audit();a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(j,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
 print("PASS: 31 exact finite PSD principal minors and negative-weight adversary");return 0
if __name__=="__main__":raise SystemExit(main())
