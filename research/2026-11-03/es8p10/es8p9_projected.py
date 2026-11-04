#!/usr/bin/env python3
"""Exact rational, finite projected matrix checks using NEW ES8-P9 integrals.

Prime-2 convolution and pole-killing B basis preserve the V16 / V17
conventions, but are *not* a new proof of the global Weil criterion.
"""
from __future__ import annotations
from fractions import Fraction as Q
from decimal import Decimal
from pathlib import Path
import argparse,json
from AUDITAR_FASE3_PROJECAO import I, PI, LN2, SQ2, L, overlap, exact_ldl_pivots
from es8p6_independent_audit import sharp_tail_upper

BASE=Path(__file__).resolve().parent

def d(x):return Q(Decimal(x))

def projected_matrix(N:int,records:dict):
    if N not in (6,8,10):raise ValueError('only N 6,8,10 audited')
    G=[[I(0) for _ in range(N)] for _ in range(N)]
    P=[[I(0) for _ in range(N)] for _ in range(N)]
    for i in range(1,N+1):
        for j in range(i,N+1):
            if (i-j)%2:continue
            record=records[f'{i},{j}']
            lo,hi=map(d,record['independent_finite_with_cauchy'])
            t=sharp_tail_upper(i,j)
            G[i-1][j-1]=G[j-1][i-1]=I(lo-t,hi+t)
            P[i-1][j-1]=P[j-1][i-1]=(-2*LN2/SQ2)*overlap(LN2,i,j)
    B=[[I(0) for _ in range(N-2)] for _ in range(N)]
    pi2=PI*PI
    for col,k in enumerate(range(3,N+1)):
        p=1 if k%2 else 2
        B[p-1][col]=-(I(Q(k,p))*(p*p*pi2+L*L)/(k*k*pi2+L*L))
        B[k-1][col]=I(1)
    W=[[G[i][j]+P[i][j] for j in range(N)] for i in range(N)]
    dim=N-2
    WB=[[sum((W[i][k]*B[k][j] for k in range(N)),I(0)) for j in range(dim)] for i in range(N)]
    A=[[sum((B[k][i]*WB[k][j] for k in range(N)),I(0)) for j in range(dim)] for i in range(dim)]
    for i in range(dim):
        for j in range(i+1,dim):
            q=I(min(A[i][j].min(),A[j][i].min()),max(A[i][j].max(),A[j][i].max()))
            A[i][j]=A[j][i]=q
    mid=[[x.mid() for x in row]for row in A]
    delta=max(sum((x.width()/2 for x in row),Q()) for row in A)
    margin=Q(1,20)
    shifted=[[mid[i][j]-(delta+margin if i==j else 0) for j in range(dim)] for i in range(dim)]
    piv,fail=exact_ldl_pivots(shifted)
    return {'N':N,'dimension':dim,'certified_M_minus_0p05I_positive':fail is None,
            'min_shifted_pivot_approx':float(min(piv)),
            'perturbation_row_sum_upper_approx':float(delta),
            'perturbation_row_sum_upper_exact':str(delta),
            'matrix_max_width_approx':float(max(x.width() for row in A for x in row)),
            'first_projected_diagonal_approx':float(mid[0][0]),
            'failed_pivot_index':fail,
            'subspace_basis_not_orthonormal':True,
            'prime_and_projection':'V16/V17 Fourier convention; internal interval arithmetic (not independently re-proved)' }

def run(filepath):
    obj=json.loads(Path(filepath).read_text(encoding='utf-8'))
    if obj['status']!='INTEGRATED_30_ENTRIES_CONDITIONALLY_CERTIFIED':raise ValueError('require complete 30-entry integration')
    out={'schema':'NUMERIA.ES8P9.FINITE_PROJECTION.1','input_path':Path(filepath).name,
         'conditional_on':'special-function, Cauchy interpolation and tail lemmas; V16 prime normalization',
         'matrices':[projected_matrix(N,obj['entries']) for N in (6,8,10)]}
    out['status']='FINITE_PROJECTED_POSITIVITY_CONDITIONAL' if all(x['certified_M_minus_0p05I_positive'] for x in out['matrices']) else 'INCONCLUSIVE'
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,default=BASE/'RESULTADOS_ES8P9_INTEGRAIS.json');ap.add_argument('--out',type=Path,default=BASE/'RESULTADOS_ES8P9_PROJECAO.json')
    a=ap.parse_args();v=run(a.input);a.out.write_text(json.dumps(v,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print(v['status']);[print(x['N'],x['perturbation_row_sum_upper_approx'],x['min_shifted_pivot_approx']) for x in v['matrices']]
    return 0 if v['status']=='FINITE_PROJECTED_POSITIVITY_CONDITIONAL' else 2
if __name__=='__main__':raise SystemExit(main())
