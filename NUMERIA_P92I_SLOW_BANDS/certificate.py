"""P92-I exact Fraction regression: bounded-nearest-gap bands and irregular block gluing.

Universal real-parameter proofs: README.md. Finite Fraction samples only check
implementation; they do NOT prove general Preissmann or Riemann hypothesis.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from math import isqrt
from pathlib import Path

SCHEMA = "NUMERIA.P92I.DELTA.SLOW_BANDS.BLOCK_GLUE.1"
ZETA2_UP = F(5,3)
RESIDUAL_LOWER = F(25,8)**2/F(3)*F(109,100)  # 2725/768
RATIO_ROOT_EXTRA = F(2037,1000)

def exact_root(q):
    q=F(q)
    if q<0: raise ValueError("negative radicand")
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    if a*a!=q.numerator or b*b!=q.denominator:
        raise ArithmeticError("expected rational square")
    return F(a,b)

def root_upper(q,digits=20):
    q=F(q)
    if q<0 or digits<1: raise ValueError("invalid square root request")
    if q==0: return F(0)
    scale=10**digits
    lower=isqrt(q.numerator*scale*scale//q.denominator)
    return F(lower if F(lower*lower,scale*scale)==q else lower+1,scale)

def coefficient(s):
    s=F(s)
    if not 1<=s<=3: raise ValueError("ratio outside proven range")
    return s*(s-1)**2

def sample_geometry(root_gaps,bridges=()):
    root_gaps=[list(map(F,a)) for a in root_gaps]
    bridges=tuple(map(F,bridges))
    if not root_gaps or len(bridges)!=len(root_gaps)-1 or any(
            not a or any(v<=0 for v in a) for a in root_gaps) or any(v<=0 for v in bridges):
        raise ValueError("invalid blocks or bridge count")
    x=[F(0)]
    blocks=[]
    for r,gaps in enumerate(root_gaps):
        if r: x.append(x[-1]+bridges[r-1])
        start=len(x)-1
        for v in gaps: x.append(x[-1]+v*v)
        blocks.append(tuple(range(start,len(x))))
    delta=[min(abs(x[i]-x[j]) for j in range(len(x)) if i!=j) for i in range(len(x))]
    roots=[exact_root(d) for d in delta]
    return x,delta,roots,blocks

def delta_matrix(x,roots):
    n=len(x)
    matrix=[[F(0)]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            v=roots[i]*roots[j]*(roots[i]-roots[j])**2/(2*(x[j]-x[i])**2)
            matrix[i][j]=matrix[j][i]=v
    return matrix

def integral_bound(S,h,g):
    S,h,g=map(F,(S,h,g))
    if min(S,h,g)<=0:raise ValueError("positive geometry required")
    return S**-4+(1/h+1/g)/(3*S**3)+1/(6*h*g*S*S)

def exact_ldlt(A):
    A=[list(map(F,row)) for row in A]
    n=len(A); smallest=None
    for k in range(n):
        p=A[k][k]
        if p<=0:raise ArithmeticError("nonpositive exact LDL pivot")
        smallest=p if smallest is None else min(smallest,p)
        for i in range(k+1,n):
            for j in range(i,n):
                A[j][i]-=A[i][k]*A[j][k]/p
                A[i][j]=A[j][i]
    return smallest

def verify(root_gaps,bridges=(),ldlt=False):
    x,d,roots,blocks=sample_geometry(root_gaps,bridges)
    matrix=delta_matrix(x,roots)
    local=[]; comparator_rows=[]
    for block in blocks:
        rr=[roots[i] for i in block]
        h,H=min(v*v for v in rr),max(v*v for v in rr)
        sr=max(rr)/min(rr)
        if sr>2:raise ArithmeticError("band spread exceeds four")
        for p in range(len(block)-1):
            if x[block[p+1]]-x[block[p]]<h:
                raise ArithmeticError("internal gap below nearest-neighbor minimum")
        c=coefficient(sr)
        for i in block:
            for j in block:
                if i!=j and matrix[i][j]>c/(2*(i-j)**2):
                    raise ArithmeticError("band entry envelope invalid")
        comparator_rows.append(ZETA2_UP*c)
        local.append({"members":len(block),"h":str(h),"H":str(H),
                      "sqrt_spread":str(sr),"local_coefficient":str(c),
                      "local_zeta_upper":str(ZETA2_UP*c),
                      "first_index":block[0],"last_index":block[-1],
                      "max_root":str(max(rr))})
    cross=[]
    for r in range(len(blocks)):
        for s in range(r+1,len(blocks)):
            S=x[blocks[s][0]]-x[blocks[r][-1]]
            h,g=F(local[r]["h"]),F(local[s]["h"])
            A,B=F(local[r]["max_root"]),F(local[s]["max_root"])
            kappa=A*B*(A+B)**2/2
            hs_upper=kappa*kappa*integral_bound(S,h,g)
            actual=sum((matrix[i][j]**2 for i in blocks[r] for j in blocks[s]),F(0))
            norm_upper=root_upper(hs_upper)
            if actual>hs_upper or norm_upper**2<hs_upper:
                raise ArithmeticError("cross-block HS envelope invalid")
            comparator_rows[r]+=norm_upper
            comparator_rows[s]+=norm_upper
            cross.append({"blocks":[r,s],"separation":str(S),
                          "kappa_upper":str(kappa),"actual_HS_squared":str(actual),
                          "integral_HS_squared_upper":str(hs_upper),
                          "block_norm_upper":str(norm_upper)})
    upper=max(comparator_rows)
    pivots=None
    if ldlt:
        pivots=[]
        for sign in (1,-1):
            n=len(x)
            A=[[RESIDUAL_LOWER*(i==j)-sign*matrix[i][j]
               for j in range(n)] for i in range(n)]
            pivots.append(str(exact_ldlt(A)))
    return {"point_count":len(x),"blocks":len(blocks),
            "distinct_nearest_gaps":len(set(d)),
            "global_gap_spread":str(max(d)/min(d)),
            "actual_max_row_sum":str(max(sum(row,F(0)) for row in matrix)),
            "rational_block_norm_upper":str(upper),
            "strict_residual_lower":str(RESIDUAL_LOWER),
            "block_test_passes":upper<RESIDUAL_LOWER,
            "block_data":local,"cross_block_data":cross,
            "two_sided_exact_ldlt_min_pivots":pivots}

def certify():
    s=RATIO_ROOT_EXTRA;c=coefficient(s)
    if not 2<s<3 or not c*c<F(24,5) or RESIDUAL_LOWER!=F(2725,768):
        raise ArithmeticError("refined spread or rational PI lower invalid")
    A=[F(1),F(6,5),F(7,5),F(5,4),F(3,2)]
    B=[F(20),F(22),F(25),F(24),F(26)]
    inputs=[
      ("flat_7_irregular",[[F(1),F(6,5),F(7,5),F(3,2),F(4,3),F(5,4)]],(),True),
      ("flat_11_many_values",[[F(1)+F(k,12) for k in range(10)]],(),True),
      ("two_bands_contrast_676",[A,B],(F(1600),),True),
      ("two_bands_huge_separation",[A,B],(F(100000),),False),
      ("two_bands_low_contrast",[[F(1),F(4,3),F(3,2)],[F(3),F(10,3),F(7,2)]],(F(300),),True),
      ("three_bands_contrast_4225",[A,[F(5),F(11,2),F(6),F(13,2)],
                                      [F(50),F(55),F(60),F(65)]],(F(1000),F(100000)),False),
    ]
    rows=[]
    for name,gaps,bridges,ldlt in inputs:
        row=verify(gaps,bridges,ldlt)
        row["name"]=name
        if not row["block_test_passes"]:raise ArithmeticError("block gate inconclusive")
        rows.append(row)
    if len(rows)!=6 or sum(r["two_sided_exact_ldlt_min_pivots"] is not None for r in rows)!=4:
        raise ArithmeticError("incorrect P92-I sample census")
    return {"schema":SCHEMA,
            "all_real_bounded_spread_proved_analytically":True,
            "bounded_spread_factor_four_proved":True,
            "refined_ratio_sqrt":str(s),"refined_ratio":str(s*s),
            "refined_polynomial_margin":str(F(24,5)-c*c),
            "irregular_block_comparator_proved":True,
            "arbitrary_configurations_global_residual_proved":False,
            "general_Preissmann_Lemma6_proved":False,
            "Riemann_Hypothesis_proved":False,
            "rational_regressions":len(rows),"two_sided_LDL_cases":4,"cases":rows}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    d=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("P92I PASSED: 6 Fraction families; 4 exact two-sided LDL; all global statements delimited")
    print("SHA256:",hashlib.sha256(args.out.read_bytes()).hexdigest())

if __name__=="__main__":
    main()
