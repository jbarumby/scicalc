"""NUMERIA ES8-P19: finite-support forms, block-projection correction and countermodels.

Only exact Fraction arithmetic. This is an ALGEBRAIC audit, not an analytic
certification of the global Weil form or of an unbounded operator.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import product
import json
import hashlib
from pathlib import Path


def matrix_rank(a):
    """Exact Gaussian-elimination rank over Q."""
    if not a:
        return 0
    a = [[F(v) for v in row] for row in a]
    n, m = len(a), len(a[0])
    assert all(len(row) == m for row in a)
    r = 0
    for col in range(m):
        k = next((j for j in range(r, n) if a[j][col]), None)
        if k is None: continue
        a[r], a[k] = a[k], a[r]
        pivot = a[r][col]
        a[r] = [v / pivot for v in a[r]]
        for j in range(n):
            if j != r and a[j][col]:
                factor = a[j][col]
                a[j] = [a[j][i]-factor*a[r][i] for i in range(m)]
        r += 1
        if r == n: break
    return r


def projected_kernel(A00, A0h, Ah0, Ahh, C):
    """B^T A B, B=(-C,I), where C is 2 x n, exact rational matrices."""
    n = len(Ahh)
    assert all(len(row) == n for row in Ahh)
    assert len(A00) == 2 and all(len(row) == 2 for row in A00)
    assert len(A0h) == 2 and all(len(row) == n for row in A0h)
    assert len(Ah0) == n and all(len(row) == 2 for row in Ah0)
    assert len(C) == 2 and all(len(row) == n for row in C)
    return [[F(Ahh[i][j]) - sum((F(C[a][i])*F(A0h[a][j]) for a in range(2)),F())
             - sum((F(Ah0[i][a])*F(C[a][j]) for a in range(2)),F())
             + sum((F(C[a][i])*F(A00[a][b])*F(C[b][j]) for a,b in product(range(2), repeat=2)),F())
             for j in range(n)] for i in range(n)]


def correction_factors(A00, A0h, Ah0, C):
    """4 outer products U V^T give B^T A B - A_HH exactly.
    u_a(i)=C[a][i], v_a(j)=sum_b A00[a][b] C[b][j]-A0h[a][j]
    u_(2+a)(i)=-Ah0[i][a], v_(2+a)(j)=C[a][j]
    No selfadjointness of A needed for this algebraic identity.
    """
    n=len(C[0]); assert len(C)==2
    U=[[F(C[a][i]) for a in range(2)] + [-F(Ah0[i][a]) for a in range(2)] for i in range(n)]
    V=[[sum((F(A00[a][b])*F(C[b][j]) for b in range(2)),F())-F(A0h[a][j])
        for a in range(2)]+[F(C[a][j]) for a in range(2)] for j in range(n)]
    return U,V


def multiply_uvt(U,V):
    return [[sum((U[i][a]*V[j][a] for a in range(4)),F()) for j in range(len(V))] for i in range(len(U))]


def qform_on_support(A,x):
    """A infinite kernel represented as callback A(i,j), x finite sparse dict."""
    return sum((F(x[i])*F(A(i,j))*F(x[j]) for i,j in product(x,x)),F())


def fixed_models(n):
    """Deterministic, nontrivial symmetric rational block, all sizes."""
    A00=[[F(3,2),F(1,7)],[F(1,7),F(5,4)]]
    A0h=[[F((-1)**(a+j)*(a+1)*(j+1), j+3+a) for j in range(n)] for a in range(2)]
    Ah0=[[A0h[a][i] for a in range(2)] for i in range(n)]
    Ahh=[[F(3+i,2) if i==j else F((-1)**(i+j), (i+j+5)*3) for j in range(n)] for i in range(n)]
    C=[[F((-1)**(a+j), (j+2)*(a+2) + a*j) for j in range(n)] for a in range(2)]
    return A00,A0h,Ah0,Ahh,C


def run():
    out=[]
    for n in [1,2,3,4,8,10,16,32]:
        A00,A0h,Ah0,Ahh,C=fixed_models(n)
        P=projected_kernel(A00,A0h,Ah0,Ahh,C)
        U,V=correction_factors(A00,A0h,Ah0,C)
        E=multiply_uvt(U,V)
        assert all(P[i][j]-Ahh[i][j]==E[i][j] for i,j in product(range(n),repeat=2))
        rank=matrix_rank(E)
        assert rank<=min(4,n)
        out.append({'n':n,'rank_of_pole_projection_correction':rank,
                    'rank_upper_bound':4, 'exact_outer_product_equality':True,
                    'symmetric_projected_matrix':all(P[i][j]==P[j][i] for i,j in product(range(n),repeat=2))})
    # All finite sections of 1/(k+1) diag are strictly positive,
    # but its spectral margin is 0 in the infinite-dimensional limit.
    diag_margins={str(n):str(F(1,n)) for n in [1,4,8,16,100,1000]}
    # The tiny-norm infinite rank operator epsilon*I has every truncation rank n.
    eps=F(1,1000)
    eps_rank={str(n):n for n in [1,4,8,16,100]}
    # Finite 8x8 positivity does not imply positivity of an arbitrary extension.
    bad_extend=lambda i,j: (F(-1) if i==j==9 else F(1) if i==j else F())
    assert all(qform_on_support(bad_extend,{i:F(1)})==1 for i in range(8))
    assert qform_on_support(bad_extend,{9:F(1)})==-1
    # correct: an off-diagonal finite-prime kernel can produce nonzero same parity entries.
    data={'schema':'NUMERIA.ES8P19.BLOCK_FINITE_SUPPORT.1',
          'scientific_status':'EXACT_FINITE_ALGEBRA_AND_COUNTERMODELS_ONLY',
          'source_scope':'Qwen finite-rank manuscript and ES8-P14/P18 operator structure',
          'block_correction':{'two_pole_rows':2,'max_rank':4,'cases':out},
          'counterexamples':{
             'small_bounded_perturbation_can_have_infinite_rank':{'operator':'(1/1000) I on ell2','norm':str(eps),'truncation_ranks':eps_rank},
             'all_finite_sections_positive_no_uniform_gap':{'matrix':'diag(1/(k+1))','first_n_minima':diag_margins,'infimum_gap':'0'},
             'finite_prefix_does_not_fix_infinite_extension':{'first_eight_diagonal_entries':'1','ninth_index_diagonal_entry':'-1'},
             'bounded_error_not_finite_rank':True},
          'unproved': ['the actual infinite Weil kernel equals the constructed algebraic form',
                       'analytic boundedness/semiboundedness/closedness of global Weil form',
                       'the full prime matrix is diagonal',
                       'spectral margin 1/20 in all dimensions',
                       'RH']}
    return data


def write(path):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    data=run();p.write_text(json.dumps(data,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    return data

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--out',default='RESULTADOS_ES8P19.json');arg=p.parse_args()
    v=write(arg.out)
    print('ES8-P19 PASSED: rank<=4 all dimensions tested, finite-support countermodels, no global Weil claim')
    print('RANKS',[(x['n'],x['rank_of_pole_projection_correction']) for x in v['block_correction']['cases']])
