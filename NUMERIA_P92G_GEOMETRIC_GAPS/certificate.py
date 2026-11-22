"""P92-G: exact rational regressions for the geometric-gap Delta theorem.

The *universal* statement for arbitrary real q>1 and any dimension is proved
in README.md, not by this finite computational regression.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
from math import isqrt
import argparse
import hashlib
import json

SCHEMA = "NUMERIA.P92G.DELTA.GEOMETRIC_GAPS.SCALE_FREE.1"
TARGET = F(15, 16)


def points_and_deltas(r, n):
    """r=sqrt(q)>1 rational; gaps 1,q,q²,... ."""
    r = F(r)
    if r <= 1 or not isinstance(n, int) or n < 2:
        raise ValueError("need rational sqrt(q)>1 and n>=2")
    q = r*r
    gaps = [q**j for j in range(n-1)]
    x = [F(0)]
    for g in gaps:
        x.append(x[-1]+g)
    delta = [gaps[0]]+[min(gaps[j-1],gaps[j]) for j in range(1,len(gaps))]+[gaps[-1]]
    if delta != [F(1)] + [q**(j-1) for j in range(1,n)]:
        raise ArithmeticError("nearest-neighbor distances do not match theory")
    return x, delta


def geometric_kernel(r, k):
    r = F(r)
    if r<=1 or not isinstance(k,int) or k<1:
        raise ValueError("require r>1, k>=1")
    q=r*r
    return (q-1)**2/(2*q*q)*r**k/(1+r**k)**2


def kernel_one_side_ceiling(r):
    """Exact rational majorant S<=t(1-t)(1+t)^2/2, t=1/r."""
    r=F(r)
    if r<=1:
        raise ValueError("r>1 required")
    t=1/r
    return t*(1-t)*(1+t)**2/2


def entries(r,n):
    r=F(r)
    x,d=points_and_deltas(r,n)
    roots=[F(1)]+[r**(j-1) for j in range(1,n)]
    return [[F(0) if i==j else
             roots[i]*roots[j]*(roots[i]-roots[j])**2/(2*(x[i]-x[j])**2)
             for j in range(n)] for i in range(n)]


def exact_ldlt_positive(threshold, delta_matrix):
    """Exact LDL test for threshold*I-Delta, independent of Schur code."""
    threshold=F(threshold)
    n=len(delta_matrix)
    a=[[threshold*(i==j)-delta_matrix[i][j] for j in range(n)] for i in range(n)]
    least=None
    for k in range(n):
        piv=a[k][k]
        if piv<=0:
            raise ArithmeticError("non-positive exact LDL pivot")
        least=piv if least is None else min(least,piv)
        for i in range(k+1,n):
            for j in range(i,n):
                a[j][i]-=a[i][k]*a[j][k]/piv
                a[i][j]=a[j][i]
    return least


def audit_case(r,n,ldlt=False):
    r=F(r)
    x,d=points_and_deltas(r,n)
    matrix=entries(r,n)
    s=kernel_one_side_ceiling(r)
    # The rational inequality 3S<15/16 is proved analytically for ALL real r>1.
    if not (0<=3*s<TARGET):
        raise ArithmeticError("invalid upper envelope")
    if n>=3:
        for i in range(1,n):
            for j in range(i+1,n):
                if matrix[i][j]!=geometric_kernel(r,j-i):
                    raise ArithmeticError("Toeplitz kernel exact identity failed")
        for j in range(2,n):
            if not (0<=matrix[0][j]<=matrix[1][j]):
                raise ArithmeticError("anchor dominated by first Toeplitz row failed")
    if matrix[0][1]!=0:
        raise ArithmeticError("first anchor interaction must vanish")
    row_sums=[sum(row,F(0)) for row in matrix]
    if any(row>3*s for row in row_sums):
        raise ArithmeticError("Schur row sum violates scale-free theorem")
    pivot=exact_ldlt_positive(TARGET,matrix) if ldlt else None
    return {"sqrt_gap_ratio":str(r),"gap_ratio":str(r*r),
            "points":n, "distinct_nearest_neighbor_gaps":len(set(d)),
            "actual_max_row_sum":str(max(row_sums)),
            "universal_rational_ceiling":str(3*s),
            "margin_to_15_over_16":str(TARGET-max(row_sums)),
            "exact_ldlt_min_pivot":None if pivot is None else str(pivot)}


def certify():
    roots=(F(101,100),F(11,10),F(6,5),F(3,2),F(2),F(3),F(5))
    sizes=(2,3,4,6,9,13,17,25)
    rows=[]
    for r,n in product(roots,sizes):
        rows.append(audit_case(r,n,ldlt=n<=13))
    large=audit_case(F(11,10),60,ldlt=False)
    rows.append(large)
    if len(rows)!=57 or sum(x["exact_ldlt_min_pivot"] is not None for x in rows)!=42:
        raise ArithmeticError("wrong P92G regression census")
    if max(x["distinct_nearest_neighbor_gaps"] for x in rows)!=59:
        raise ArithmeticError("unbounded-scale direction not tested")
    return {"schema":SCHEMA,
            "theorem":"For every real q>1 and every finite geometric-gap configuration, ||Delta||_2 < 15/16",
            "proof":"Toeplitz interaction, anchor dominated by first row, Schur, exact max of quartic on (0,1)",
            "all_real_geometric_ratios_proved_analytically":True,
            "number_of_distinct_scales_unbounded_in_theorem":True,
            "arbitrary_point_configurations_proved":False,
            "general_Preissmann_Lemma6_proved":False,
            "Riemann_Hypothesis_proved":False,
            "rational_regressions":len(rows),
            "independent_rational_ldlt":sum(x["exact_ldlt_min_pivot"] is not None for x in rows),
            "cases":rows}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    obj=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("P92G PASS:",obj["rational_regressions"],"exact cases,",
          obj["independent_rational_ldlt"],"exact LDL; q>1 theorem is analytic")
    print("Full Preissmann, arbitrary configurations, and RH remain open")
    print("SHA256:",hashlib.sha256(args.out.read_bytes()).hexdigest())

if __name__=="__main__":
    main()
