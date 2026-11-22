"""P92-K: all-Fraction regression for quadratically separated multiclusters.

The real-parameter, arbitrary-dimensional theorem is proved in README.md.
Finite Fraction checks below certify specific arrays, not the universal theorem.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
from NUMERIA_P92J_MULTIANCHOR.certificate import ldlt_positive

SCHEMA = "NUMERIA.P92K.DELTA.MULTICLUSTER.QUADRATIC_SEPARATION.1"
ZERO = F(0)
BLOCK_FACTOR = F(5, 14)
TARGET = F(3, 2)

def construct(roots, sizes):
    roots=tuple(F(v) for v in roots)
    if not roots or len(roots)!=len(sizes) or any(v<=0 for v in roots):
        raise ValueError("positive roots with matching lengths required")
    if any(not isinstance(n,int) or n<2 for n in sizes):
        raise ValueError("each cluster must contain >=2 points")
    steps=tuple(s*s for s in roots)
    H=max(steps)
    x=[]
    groups=[]
    for r,(h,n) in enumerate(zip(steps,sizes)):
        first=ZERO if r==0 else x[-1]+H*r*r
        group=[first+j*h for j in range(n)]
        groups.append(tuple(range(len(x),len(x)+n)))
        x.extend(group)
    return tuple(x),tuple(groups),steps,roots

def nearest(x):
    return tuple(min(abs(x[i]-x[j]) for j in range(len(x)) if j!=i)
                 for i in range(len(x)))

def matrix(x,roots):
    n=len(x)
    return tuple(tuple(
        ZERO if i==j else roots[i]*roots[j]*(roots[i]-roots[j])**2
             /(2*(x[i]-x[j])**2)
        for j in range(n)) for i in range(n))

def audit_case(roots,sizes, full_ldlt=False):
    x,groups,h,s=construct(roots,sizes)
    obs=nearest(x)
    if any(obs[i]!=h[r] for r,g in enumerate(groups) for i in g):
        raise ArithmeticError("nearest neighbor gap mismatch")
    nodal=tuple(s[r] for r,g in enumerate(groups) for _ in g)
    if any(nodal[i]**2!=obs[i] for i in range(len(x))):
        raise ArithmeticError("inconsistent square roots")
    delta=matrix(x,nodal)
    K=len(groups)
    A=[[ZERO]*K for _ in range(K)]
    pairs=[]
    for r in range(K):
        if any(delta[i][j]!=0 for i in groups[r] for j in groups[r]):
            raise ArithmeticError("intrablock Delta not zero")
        for t in range(r+1,K):
            D=x[groups[t][0]]-x[groups[r][-1]]
            rho=D/max(h[r],h[t])
            if rho<(t-r)**2:
                raise ArithmeticError("quadratic separation failed")
            hs=sum((delta[i][j]**2 for i in groups[r] for j in groups[t]),ZERO)
            cap=F(1,8)/rho**2
            if not 0<=hs<=cap:
                raise ArithmeticError("Hilbert-Schmidt inequality failed")
            A[r][t]=A[t][r]=BLOCK_FACTOR/rho
            pairs.append(dict(r=r,s=t,rho=str(rho),hs_square=str(hs),
                              hs_square_ceiling=str(cap)))
    rows=[sum(row,ZERO) for row in A]
    if not all(v<F(10,7) for v in rows):
        raise ArithmeticError("comparison row sum too large")
    reduced=[[TARGET*(i==j)-A[i][j] for j in range(K)] for i in range(K)]
    ldlt_positive(reduced)
    if full_ldlt:
        n=len(x)
        ldlt_positive([[TARGET*(i==j)-delta[i][j] for j in range(n)]
                       for i in range(n)])
    return dict(group_count=K,point_count=len(x),pair_count=len(pairs),
                max_to_min_gap=str(max(obs)/min(obs)),
                full_matrix_fraction_ldlt=bool(full_ldlt),
                comparison_fraction_ldlt=True,
                max_rational_comparison_row_sum=str(max(rows)),
                pairs=pairs)

def certify():
    roots=(F(1),F(2),F(1,2),F(3),F(4,3),F(5,2),F(1,3),F(4),
           F(2,3),F(7,2),F(5,4),F(6,5),F(3,4),F(7),
           F(1,5),F(9,2),F(3,5),F(11,3),F(1,7),
           F(13,4),F(2,5),F(8,3),F(1,8),F(5))
    plans=((2,2,2),(3,2,3),(4,2,3),(5,3,2),(6,2,3),
           (8,2,2),(10,2,3),(12,2,2),(16,2,2),
           (20,2,2),(24,2,2))
    records=[]
    for k,a,b in plans:
        sizes=tuple(a if r%2==0 else b for r in range(k))
        records.append(audit_case(roots[:k],sizes,full_ldlt=k<=5))
    total=sum(x["pair_count"] for x in records)
    full=sum(x["full_matrix_fraction_ldlt"] for x in records)
    if (len(records),total,full)!=(11,760,4):
        raise ArithmeticError("unexpected certificate census")
    return dict(
        schema=SCHEMA, fraction_configurations=11,
        pairwise_interactions=760, exact_full_LDL=4,
        exact_comparison_LDL=11,
        theorem="For arbitrary real parameters satisfying D_rs >= |r-s|^2 max(h_r,h_s), ||Delta|| <= zeta(2)/sqrt(2) < 10/7",
        all_real_parameters_with_quadratic_separation_proved=True,
        arbitrary_configurations_proved=False,
        general_Preissmann_Lemma6_proved=False,
        Riemann_Hypothesis_proved=False,
        cases=records)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    d=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("P92K PASS: 11 exact arrays, 760 exact block HS bounds, 4 full LDL, 11 block LDL")
    print("Global unrestricted Delta bound, full Preissmann, RH not proved")
    print("SHA256",hashlib.sha256(args.out.read_bytes()).hexdigest())

if __name__=="__main__":
    main()
