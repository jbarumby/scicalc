"""NUMERIA P92-J: exact rational multisource Delta bipartite certificates.

The all-real, all-dimensions theorem is derived in README.md; these
Fraction computations are independent finite regressions and LDL checks.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from itertools import product
from pathlib import Path

SCHEMA = "NUMERIA.P92J.DELTA.MULTIANCHOR.SCHUR.1"
BOUND_SQUARE = F(5, 48)
EASY_BOUND_SQUARE = F(1, 9)
P88_RESIDUAL_LOWER = F(2725, 768)
ZERO = F(0)


def coordinates(s, anchors, cluster):
    """h=1, L=s*s: r equally spaced anchors followed by M dense points."""
    s = F(s)
    if s < 1 or not isinstance(anchors, int) or not isinstance(cluster, int) or anchors < 1 or cluster < 2:
        raise ValueError("s>=1, anchors>=1, cluster>=2 required")
    lam = s*s
    x = tuple(lam*i for i in range(anchors)) + tuple(lam*anchors + j for j in range(cluster))
    d = tuple(lam for _ in range(anchors)) + tuple(F(1) for _ in range(cluster))
    for i in range(len(x)):
        observed = min(abs(x[i]-x[j]) for j in range(len(x)) if i != j)
        if observed != d[i]:
            raise ArithmeticError("nearest gap mismatch")
    return x, d


def rectangular_weights(s, anchors, cluster):
    s=F(s)
    if s<1 or anchors<1 or cluster<2:
        raise ValueError("invalid parameters")
    lam=s*s
    coeff=s*(s-1)**2/2
    return tuple(tuple(coeff/(m*lam+j)**2 for j in range(cluster)) for m in range(1,anchors+1))


def sqrt_upper_envelope(s):
    """Rational bound ||Delta||_HS^2 <= 5/48 (1-u)^4 (1+3u^2)."""
    s=F(s)
    if s<1:
        raise ValueError("s>=1 required")
    u=1/s
    return BOUND_SQUARE*(1-u)**4*(1+3*u*u)


def ldlt_positive(a):
    """Fraction-only LDL; return pivot bit lengths, not vast decimal strings."""
    a=[list(map(F,row)) for row in a]
    n=len(a)
    if any(len(row)!=n for row in a):
        raise ValueError("must be square")
    for i in range(n):
        for j in range(i):
            if a[i][j]!=a[j][i]:
                raise ArithmeticError("non-symmetric input")
    lens=[]
    for k in range(n):
        pivot=a[k][k]
        if pivot<=0:
            raise ArithmeticError(f"nonpositive exact LDL pivot {k}")
        lens.append(pivot.numerator.bit_length())
        for i in range(k+1,n):
            for j in range(i,n):
                a[j][i]-=a[i][k]*a[j][k]/pivot
                a[i][j]=a[j][i]
    return lens


def delta_matrix_independent(s, anchors, cluster):
    """Compute original Delta formula from all coordinates and nearest gaps."""
    x,d=coordinates(s,anchors,cluster)
    roots=[F(s)]*anchors+[F(1)]*cluster
    n=len(x)
    return [[ZERO if i==j else roots[i]*roots[j]*(roots[i]-roots[j])**2/(2*(x[i]-x[j])**2)
             for j in range(n)] for i in range(n)]


def schur_core_reduce(core, cross, t):
    """For [[core,cross],[cross^T,0]], tI-Delta SPD iff tI-core-WW^T/t SPD."""
    t=F(t)
    if t<=0:
        raise ValueError("positive threshold required")
    r=len(core)
    if not r or len(cross)!=r:
        raise ValueError("bad core/cross")
    if any(len(row)!=r for row in core):
        raise ValueError("core must be square")
    m=len(cross[0])
    if any(len(row)!=m for row in cross):
        raise ValueError("cross rows must have same length")
    reduced=[[t*(i==j)-F(core[i][j])-sum((F(cross[i][k])*F(cross[j][k]) for k in range(m)),ZERO)/t
              for j in range(r)] for i in range(r)]
    return ldlt_positive(reduced)


def audit_case(s, anchors, cluster, full_ldlt=False):
    s=F(s)
    weights=rectangular_weights(s,anchors,cluster)
    # reversed rows: weights enumerate distance-to-cluster 1..r; irrelevant to norms.
    square_norm_hs=sum((v*v for row in weights for v in row), ZERO)
    bound=sqrt_upper_envelope(s)
    if not (ZERO<=square_norm_hs<=bound<=BOUND_SQUARE<EASY_BOUND_SQUARE):
        raise ArithmeticError("Hilbert-Schmidt analytic envelope regression failed")
    # Core of equally scaled anchors vanishes identically; exact Schur complement.
    pivots=schur_core_reduce([[ZERO]*anchors for _ in range(anchors)],weights,F(1,3))
    full_pivots=None
    if full_ldlt:
        m=delta_matrix_independent(s,anchors,cluster)
        # In the exact original matrix, the anchor and cluster subblocks are zero.
        for i in range(anchors):
            for j in range(anchors):
                if m[i][j]!=0:
                    raise ArithmeticError("anchor block not zero")
            for j in range(cluster):
                if m[i][anchors+j]!=weights[anchors-1-i][j]:
                    raise ArithmeticError("cross-block original formula mismatch")
        if any(m[i][j]!=0 for i in range(anchors,len(m)) for j in range(anchors,len(m))):
            raise ArithmeticError("dense-cluster block not zero")
        full_pivots=[]
        for sign in (1,-1):
            a=[[F(1,3)*(i==j)-sign*m[i][j] for j in range(len(m))] for i in range(len(m))]
            full_pivots.append(ldlt_positive(a))
    last_anchor_row=sum(weights[0],ZERO)
    return {
        "sqrt_scale_ratio":str(s),"gap_ratio":str(s*s),
        "anchors":anchors,"cluster_points":cluster,"total_points":anchors+cluster,
        "hs_square_exact":str(square_norm_hs),
        "analytic_rational_envelope":str(bound),
        "nearest_anchor_row_sum_exact":str(last_anchor_row),
        "exact_reduced_schur_ldlt_positive":True,
        "reduced_schur_pivot_numerator_bits":pivots,
        "full_two_sided_ldlt_positive":full_pivots is not None,
        "full_two_sided_ldlt_pivot_numerator_bits":full_pivots,
    }


def certify():
    roots=(F(1),F(6,5),F(3,2),F(2),F(3),F(8),F(16))
    shapes=((1,2),(1,9),(2,4),(3,7),(5,10),(15,25))
    records=[audit_case(s,r,m,full_ldlt=(r+m<=10))
             for s,(r,m) in product(roots,shapes)]
    stress=audit_case(F(16),4,280,full_ldlt=False)
    # The row-sum failure is established WITHOUT relying on summing 280 entries:
    # decreasing integral bound sum_{k=0}^{M-1} c/(lambda+k)^2 > integral_0^M.
    s=F(16);lam=s*s;M=280
    coefficient=s*(s-1)**2/2
    exact_integral_lower=coefficient*(1/lam-1/(lam+M))
    if exact_integral_lower!=F(7875,2144) or not (
            stress["nearest_anchor_row_sum_exact"] and
            F(stress["nearest_anchor_row_sum_exact"])>exact_integral_lower>P88_RESIDUAL_LOWER):
        raise ArithmeticError("stress test did not defeat row-sum bound")
    stress["row_sum_lower_by_exact_integral"]=str(exact_integral_lower)
    stress["fails_original_row_sum_sufficient_test"]=True
    records.append(stress)
    n_full=sum(int(x["full_two_sided_ldlt_positive"]) for x in records)
    if len(records)!=43 or n_full!=28:
        raise ArithmeticError("wrong exact certificate census")
    return {
        "schema":SCHEMA,
        "analytical_family_theorem": "For all reals L>=h>0 and all integers r>=1,M>=2: ||Delta||_2^2 < 5/48, hence ||Delta||_2 < 1/3",
        "general_core_leaf_Schur_reduction_proved":True,
        "multi_anchor_family_all_real_parameters_proved":True,
        "arbitrary_point_configurations_proved":False,
        "all_block_partitions_pass_Schur":False,
        "full_Preissmann_Lemma6_proved":False,
        "Riemann_Hypothesis_proved":False,
        "fraction_regressions":len(records),
        "reduced_exact_ldlt_cases":len(records),
        "independent_full_two_sided_ldlt_cases":n_full,
        "row_sum_stress_example":{"anchors":4,"cluster_points":280,
            "lower_bound":str(exact_integral_lower),
            "P88_rational_residual_lower":str(P88_RESIDUAL_LOWER),
            "row_sum_fails_but_spectral_passes":True},
        "records":records,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    cert=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(cert,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("P92J PASS:",cert["fraction_regressions"],"fraction cases;",
          cert["reduced_exact_ldlt_cases"],"reduced LDL;",
          cert["independent_full_two_sided_ldlt_cases"],"independent full two-sided LDL")
    print("Matrix row-sum bound fails for 4+280 points while exact spectral bound <1/3")
    print("General Delta estimate, full Preissmann and RH remain OPEN")
    print("SHA256:",hashlib.sha256(args.out.read_bytes()).hexdigest())

if __name__=="__main__":
    main()
