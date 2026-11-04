"""Exact Fraction verification of the ES8-P13 LDL/SOS certificates.
Does not replace Lean compilation; computes original P10 interval inputs again.
"""
import json
from fractions import Fraction as Q
from pathlib import Path
import pytest
from generate_es8p13 import make_lean, exact_ldl, MAX_ERROR, SHIFT, S
from es8p10_lemmas import new_intervals
from scratch import matrix

HERE=Path(__file__).resolve().parent
CERT=json.loads((HERE/'RESULTADOS_ES8P13.json').read_text())

def data(n):
    p=CERT[f'N{n+2}']
    C=[[Q(t) for t in row]for row in p['matrix_coarse']]
    D=[Q(t) for t in p['ldl_pivots_exact']]
    U=[[Q(t) for t in row]for row in p['upper_unit_factor']]
    return p,C,D,U

@pytest.mark.parametrize('dim',[4,6,8])
def test_real_exact_ldl(dim):
    p,C,D,U=data(dim)
    assert all(v>0 for v in D)
    assert D==exact_ldl(C)[0]
    assert U==exact_ldl(C)[1]
    assert all(U[i][j]==(1 if i==j else 0) for i in range(dim) for j in range(i+1))
    for i in range(dim):
        for j in range(dim):
            expected=C[i][j]-(SHIFT if i==j else 0)
            actual=sum((D[k]*U[k][i]*U[k][j] for k in range(dim)),Q())
            assert actual==expected

@pytest.mark.parametrize('dim',[4,6,8])
def test_raw_source_intervals_to_coarse_error(dim):
    p,C,_,_=data(dim)
    assert p['dimension']==dim
    intervals=matrix(dim+2,new_intervals())
    maxrow=Q(0)
    for i in range(dim):
        r=Q(0)
        for j in range(dim):
            A=intervals[i][j]
            assert A.lo==p['source_interval_raw_integer_lower'][i][j]
            assert A.hi==p['source_interval_raw_integer_upper'][i][j]
            assert C[i][j].denominator<=S
            r+=max(abs(A.min()-C[i][j]),abs(A.max()-C[i][j]))
        maxrow=max(maxrow,r)
    assert maxrow==Q(p['max_abs_row_error_exact'])
    assert maxrow<MAX_ERROR

def test_parent_prefix_restrictions():
    p8,C8,_,_=data(8)
    for n in (4,6):
        p,C,_,_=data(n)
        assert all(C[i][j]==C8[i][j] for i in range(n) for j in range(n))
        assert all(p['source_interval_raw_integer_lower'][i][j]==p8['source_interval_raw_integer_lower'][i][j]
                   for i in range(n)for j in range(n))

def test_negative_mutation_detected():
    _,C,d,U=data(8)
    C[0][2]=-C[0][2]
    assert sum((d[k]*U[k][0]*U[k][2]for k in range(8)),Q())!=C[0][2]

def test_negative_pivot_mutation_detected():
    _,C,D,U=data(8)
    D[0]=-D[0]
    assert D[0]<0
    assert sum((D[k]*U[k][0]**2 for k in range(8)),Q()) != C[0][0]-SHIFT

def test_statement_scope_limits():
    src=HERE/'ES8P13.lean' if (HERE/'ES8P13.lean').exists() else HERE.parents[2]/'NUMERIA_E8K_FORMAL/ES8P13.lean'
    lean=src.read_text()
    assert 'theorem concrete_positive4' in lean
    assert 'theorem concrete_positive6' in lean
    assert 'theorem concrete_positive8' in lean
    assert 'NO Lean proof of the complete numerical-to-analytic chain' in lean
    assert not any(('sorry' in line or 'admit' in line)for line in lean.splitlines() if not line.strip().startswith('--'))

if __name__=='__main__':
    raise SystemExit(pytest.main(['-q',str(__file__)]))
