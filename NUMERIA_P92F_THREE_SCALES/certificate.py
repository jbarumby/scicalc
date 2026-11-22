"""P92-F exact rational regressions for analytic scale-count Schur theorem.

The universal theorem is proven with real-variable inequalities in README.md;
this script tests rational-square scale examples and exact LDL pivots.
No float, CAS or unsupported global-Delta status.
"""
from fractions import Fraction as F
from itertools import product
from math import isqrt
from pathlib import Path
import argparse
import hashlib
import json

SCHEMA='NUMERIA.P92F.DELTA.SCALE_COUNT.SCHUR.1'
ZETA2_RATIONAL_UPPER=F(5,3)
PI_RATIONAL_LOWER=F(25,8)
SQRT_6_OVER_5_RATIONAL_LOWER=F(109,100)

def square_root_exact(q):
    q=F(q)
    if q<0:
        raise ValueError('negative square')
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    if a*a!=q.numerator or b*b!=q.denominator:
        raise ValueError('rational gap not a rational square; regression only')
    return F(a,b)

def points_deltas_from_gaps(gaps):
    gaps=tuple(map(F,gaps))
    if not gaps or any(t<=0 for t in gaps):
        raise ValueError('positive gaps required')
    x=[F(0)]
    for g in gaps:
        x.append(x[-1]+g)
    delta=[gaps[0]] + [min(gaps[i-1],gaps[i]) for i in range(1,len(gaps))] + [gaps[-1]]
    return x,delta

def delta_matrix(x,d,roots):
    n=len(x)
    return [[F(0) if i==j else
             roots[i]*roots[j]*(roots[i]-roots[j])**2/(2*(x[i]-x[j])**2)
             for j in range(n)] for i in range(n)]

def weighted_schur_rows(x,d,roots):
    n=len(x)
    return [sum((d[j]*(roots[i]-roots[j])**2/(2*(x[i]-x[j])**2)
                 for j in range(n) if j!=i), F(0)) for i in range(n)]

def check_class_shells(x,d,roots):
    for i,a in enumerate(d):
        for b in sorted(set(d)):
            if a==b:
                continue
            term=sum((d[j]*(roots[i]-roots[j])**2/(2*(x[i]-x[j])**2)
                      for j in range(len(d)) if j!=i and d[j]==b),F(0))
            if a>=b:
                u=square_root_exact(b/a)
                ceiling=(1-u)**2*(1+u*u)
            else:
                u=square_root_exact(a/b)
                ceiling=ZETA2_RATIONAL_UPPER*(1-u)**2
            if not (F(0)<=term<=ceiling):
                raise ArithmeticError('shell comparison failed')

def exact_positive_definite_ldlt(diag, matrix):
    """Independent exact Sylvester/LDLT check on diag*I-Delta."""
    n=len(matrix)
    a=[[F(diag) * (i==j)-matrix[i][j] for j in range(n)] for i in range(n)]
    pivots=[]
    for k in range(n):
        pivot=a[k][k]
        if pivot<=0:
            raise ArithmeticError('nonpositive exact LDL pivot')
        pivots.append(pivot)
        for i in range(k+1,n):
            for j in range(i,n):
                a[j][i]-=a[i][k]*a[j][k]/pivot
                a[i][j]=a[j][i]
    return min(pivots)

def audit_gap_sequence(gaps, independent_ldlt=False):
    x,d=points_deltas_from_gaps(gaps)
    roots=tuple(square_root_exact(v) for v in d)
    vals=sorted(set(d))
    rows=weighted_schur_rows(x,d,roots)
    matrix=delta_matrix(x,d,roots)
    check_class_shells(x,d,roots)
    for i in range(len(d)):
        if rows[i]!=sum((matrix[i][j]*roots[j]/roots[i]
                         for j in range(len(d))),F(0)):
            raise ArithmeticError('weighted Schur row does not match Delta')
        if rows[i]>(len(vals)-1)*ZETA2_RATIONAL_UPPER:
            raise ArithmeticError('unexpected Schur bound violation')
    if len(vals)<=3:
        if not (2*ZETA2_RATIONAL_UPPER <
                PI_RATIONAL_LOWER**2/3*SQRT_6_OVER_5_RATIONAL_LOWER):
            raise ArithmeticError('incorrect rational threshold comparison')
    pivot=None
    if independent_ldlt:
        pivot=exact_positive_definite_ldlt((len(vals)-1)*ZETA2_RATIONAL_UPPER+F(1,100),matrix)
    return {'gaps':[str(t) for t in gaps], 'node_count':len(x),
            'distinct_nearest_gap_values':len(vals),
            'gap_values':[str(t) for t in vals],
            'max_weighted_schur_row':str(max(rows)),
            'rational_zeta_upper_bound':str((len(vals)-1)*ZETA2_RATIONAL_UPPER),
            'ldlt_min_pivot':str(pivot) if pivot is not None else None}

def certify():
    patterns=[(1,4,4,9,9), (9,9,4,4,1), (1,4,9,9,4,1,4,9),
              (9,4,1,4,9,4,1,4,9), (1,)*25,
              (1,4,9)*7, (9,9,4,4,1,1)*5,
              (1,1,1,4,4,4,9,9,9)*3]
    rows=[]
    for n in range(2,6):
        for gap in product((1,4,9),repeat=n):
            rows.append(audit_gap_sequence(gap,independent_ldlt=n<=3))
    base_count=len(rows)
    for gap in patterns:
        rows.append(audit_gap_sequence(gap,independent_ldlt=len(gap)<=10))
    if len(rows)!=368 or base_count!=360:
        raise ArithmeticError('unexpected exact rational coverage')
    if max(row['distinct_nearest_gap_values'] for row in rows)!=3:
        raise ArithmeticError('missing three-scale regression')
    return {'schema':SCHEMA,'universal_analytic_statement':
            'for any finite distinct real nodes with at most m distinct nearest-neighbor distances, ||Delta||_2 <= (m-1)*zeta(2)',
            'three_scale_global_class_bound': '||Delta||_2 <= pi^2/3, provided at most three distinct nearest-neighbor gaps',
            'scale_class_unrestricted_number_of_nodes':True,
            'three_scales_imply_P88_residual_bound':True,
            'arbitrarily_many_scales_proved':False,
            'full_Preissmann_Lemma5_or_6_proved':False,
            'Riemann_Hypothesis_proved':False,
            'regression_count':len(rows),
            'ldlt_certified_count':sum(row['ldlt_min_pivot'] is not None for row in rows),
            'regressions':rows}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    o=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(o,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print('P92F PASS:',o['regression_count'],'rational exact gap-sequence regressions;',
          o['ldlt_certified_count'],'exact positive definite LDL checks')
    print('Scope: at most THREE distinct nearest-neighbor gap values, NOT arbitrary scales or RH')
    print('SHA256',hashlib.sha256(args.out.read_bytes()).hexdigest())

if __name__=='__main__':
    main()
