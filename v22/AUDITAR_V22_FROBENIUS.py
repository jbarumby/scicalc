#!/usr/bin/env python3
"""V22: exact arithmetic point counts for E:y²=x³-x over F5 and F25."""
from pathlib import Path
import json,argparse
P=5
def mul(x,y):
 a,b=x;c,d=y
 return ((a*c+2*b*d)%P,(a*d+b*c)%P)
def add(x,y):return ((x[0]+y[0])%P,(x[1]+y[1])%P)
def power(x,n):
 v=(1,0)
 for _ in range(n):v=mul(v,x)
 return v
def count(ext):
 elems=[(a,b) for a in range(P) for b in range(P)] if ext else [(a,0) for a in range(P)]
 total=1
 for x in elems:
  target=add(power(x,3),((-x[0])%P,(-x[1])%P))
  total+=sum(power(y,2)==target for y in elems)
 return total
def audit():
 assert all((a*a)%5!=2 for a in range(5))
 n1=count(False);n2=count(True)
 a1=5+1-n1;a2=a1*a1-10
 assert (n1,n2,a1,a2)==(8,32,-2,-6)
 assert a1*a1<=4*P and (-1)**2+2**2==P
 return {"schema":"NUMERIA.V22.FROBENIUS_F5_F25.1",
 "status":"EXACT_FINITE_FIELD_FROBENIUS_MODEL_CHECKED",
 "curve":"y²=x³-x","points":{"F5":n1,"F25":n2},
 "frobenius_polynomial":"X²+2X+5","frobenius_eigenvalues":["-1+2i","-1-2i"],
 "modulus_squared":5,"note":"Finite-field model, not the classical zeta-function RH."}
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
 j=audit();a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(j,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
 print("PASS F5=8 F25=32, |Frobenius eigenvalues|²=5");return 0
if __name__=="__main__":raise SystemExit(main())
