"""P92-H exact rational regressions for variable-ratio lacunary nearest-neighbor gaps.

The real-parameter theorem is proved in README.md, not by finite regressions.
No floats, CAS, numerical eigensolvers or tolerance tests in the verifier.
"""
from fractions import Fraction as Q
from itertools import cycle
from pathlib import Path
import argparse
import hashlib
import json

SCHEMA = "NUMERIA.P92H.DELTA.LACUNARY.GAP_RATIOS.1"
THRESHOLD_GAP_RATIO = Q(64, 25)
THRESHOLD_SQRT_RATIO = Q(8, 5)
T = Q(5, 8)
PHI_MAX = Q(4, 27)
ROW_BOUND = Q(4, 9) + T**4/(1-T) + Q(2, 27)
TARGET = Q(15, 16)


def phi(s):
    s=Q(s)
    if not 0 <= s <= 1:
        raise ValueError('phi domain requires 0<=s<=1')
    return s*(1-s)**2


def point_data(root_multipliers):
    """Gap roots h_k=sqrt(g_k), h_0=1, h_{k+1}/h_k>=8/5."""
    mults=tuple(map(Q,root_multipliers))
    if any(r<THRESHOLD_SQRT_RATIO for r in mults):
        raise ValueError('each consecutive gap ratio must be at least 64/25')
    gap_roots=[Q(1)]
    for ratio in mults:
        gap_roots.append(gap_roots[-1]*ratio)
    gaps=[r*r for r in gap_roots]
    xs=[Q(0)]
    for gap in gaps:
        xs.append(xs[-1]+gap)
    # 2+(len(gaps)-1)=len(xs) nearest-neighbor gap roots
    delta_roots=[gap_roots[0]]+gap_roots
    if len(delta_roots)!=len(xs):
        raise ArithmeticError('bad nearest-neighbor root count')
    actual=[min(abs(xs[i]-xs[j]) for j in range(len(xs)) if j!=i)
            for i in range(len(xs))]
    if actual != [r*r for r in delta_roots]:
        raise ArithmeticError('nearest gap identities fail')
    return xs,delta_roots


def delta_matrix(xs, delta_roots):
    xs=list(map(Q,xs))
    root=list(map(Q,delta_roots))
    if len(xs)!=len(root) or len(xs)<2 or any(xs[i+1]<=xs[i] for i in range(len(xs)-1)):
        raise ValueError('invalid points and roots')
    n=len(xs)
    M=[[Q(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            M[i][j]=M[j][i]=root[i]*root[j]*(root[i]-root[j])**2/(2*(xs[j]-xs[i])**2)
    return M


def index_envelope(k):
    if not isinstance(k,int) or k<1:
        raise ValueError('positive integer gap in index required')
    return min(PHI_MAX,T**k)


def exact_ldlt_min_pivot(M, signed=1):
    """Exact LDL positive definiteness of ROW_BOUND I - signed*Delta."""
    n=len(M)
    a=[[ROW_BOUND*(i==j)-signed*M[i][j] for j in range(n)] for i in range(n)]
    least=None
    for k in range(n):
        pivot=a[k][k]
        if pivot<=0:
            raise ArithmeticError('nonpositive exact LDL pivot')
        least=pivot if least is None else min(least,pivot)
        for i in range(k+1,n):
            aik=a[i][k]
            for j in range(i,n):
                a[j][i] -= aik*a[j][k]/pivot
                a[i][j] = a[j][i]
    return least


def verify_example(root_multipliers,ldlt=False):
    xs, roots=point_data(root_multipliers)
    M=delta_matrix(xs,roots)
    n=len(xs)
    for i in range(1,n):
        for j in range(i+1,n):
            s=roots[i]/roots[j]
            k=j-i
            if s>T**k or M[i][j]>phi(s)/2 or M[i][j]>index_envelope(k)/2:
                raise ArithmeticError('pair envelope incorrect')
    if M[0][1]!=0 or any(M[0][j]>M[1][j] for j in range(2,n)):
        raise ArithmeticError('anchor-row domination violated')
    sums=[sum(row,Q(0)) for row in M]
    if max(sums)>ROW_BOUND:
        raise ArithmeticError('Schur bound violated')
    smallest_pivots=None
    if ldlt:
        smallest_pivots=[str(exact_ldlt_min_pivot(M,+1)),str(exact_ldlt_min_pivot(M,-1))]
    return {'point_count':n,'distinct_nearest_gaps':len(set(roots)),
            'root_multipliers':[str(v) for v in map(Q,root_multipliers)],
            'max_exact_row_sum':str(max(sums)),
            'margin_to_rational_bound':str(ROW_BOUND-max(sums)),
            'two_sided_exact_ldlt_min_pivots':smallest_pivots}


def certify():
    if not (T==Q(5,8) and ROW_BOUND==Q(12793,13824) and ROW_BOUND<TARGET):
        raise ArithmeticError('wrong universal rational constant')
    if not (PHI_MAX==Q(4,27) and
            ROW_BOUND==3*PHI_MAX+T**4/(1-T)+PHI_MAX/2):
        raise ArithmeticError('wrong series majorant')
    multipliers=(Q(8,5),Q(13,8),Q(5,3),Q(7,4),Q(2),Q(9,4),Q(3),Q(5),Q(10))
    cases=[]
    for size in range(2,15):
        for variant in range(5):
            seq=tuple(multipliers[(variant+size+3*k*k+k)%len(multipliers)]
                      for k in range(size-2))
            cases.append(verify_example(seq,ldlt=size<=8))
    # Large contrast, near-threshold, and alternating-ratio cases.
    for n in (16,20,28):
        cases.append(verify_example((Q(8,5),)*(n-2),ldlt=False))
        cases.append(verify_example(tuple(Q(8,5) if i%2==0 else Q(9,4)
                                          for i in range(n-2)),ldlt=False))
    if len(cases)!=71:
        raise ArithmeticError('wrong regression count')
    ldlt_count=sum(c['two_sided_exact_ldlt_min_pivots'] is not None for c in cases)
    if ldlt_count!=35:
        raise ArithmeticError(f'unexpected LDL case count {ldlt_count}')
    return {'schema':SCHEMA,
            'theorem':'If consecutive gaps g[j+1]/g[j] >=64/25, then ||Delta||_2 <=12793/13824 <15/16 for all N>=2.',
            'real_parameter_proof':'Index envelope phi(s)<=min(4/27,s); anchor-row domination; symmetric nonnegative Schur norm',
            'minimum_gap_ratio':str(THRESHOLD_GAP_RATIO),
            'rational_norm_upper':str(ROW_BOUND),
            'strict_target':str(TARGET),
            'exact_margin_to_target':str(TARGET-ROW_BOUND),
            'all_lacunary_variable_ratios_proved_analytically':True,
            'geometric_exact_ratio_required':False,
            'all_gap_ratios_above_one_proved_to_meet_15_over_16':False,
            'arbitrary_point_configurations_proved':False,
            'full_Preissmann_Lemma6_proved':False,
            'Riemann_Hypothesis_proved':False,
            'exact_rational_regression_count':len(cases),
            'two_sided_exact_ldlt_cases':ldlt_count,
            'cases':cases}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',required=True,type=Path)
    args=p.parse_args()
    cert=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(cert,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('P92H PASS:',cert['exact_rational_regression_count'],'Fraction examples,',
          cert['two_sided_exact_ldlt_cases'],'two-sided LDL checks')
    print('Universal RATIO-LACUNARY theorem proven analytically. General Preissmann and RH open.')
    print('SHA256:',hashlib.sha256(args.out.read_bytes()).hexdigest())


if __name__=='__main__':
    main()
