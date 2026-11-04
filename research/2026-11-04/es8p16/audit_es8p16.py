"""Independent exact-rational audit of the ES8-P14 concrete source.
This checks finite nestedness, NOT existence or positivity of an infinite Weil matrix.
No use of floating point for any mathematical decision.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
DEFAULT_SOURCE = ROOT / 'sources' / 'ES8P14Concrete.lean'
FIELDS = ('center', 'lower', 'upper', 'error')


def parse_matrix(src: str, field: str, n: int) -> list[list[Fraction]]:
    pat = rf'def {field}{n} \(i j : Fin {n}\) : ℚ :=\s*match i\.val, j\.val with\s*'
    m = re.search(pat, src)
    if not m:
        raise ValueError(f'Missing matrix {field}{n}')
    body = src[m.end():].split('  | _, _ => 0', 1)[0]
    out = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    hits = set()
    pair = re.compile(r'^\s*\|\s*(\d+),\s*(\d+)\s*=>\s*\((-?\d+)\s*:\s*ℚ\)\s*/\s*(\d+)\s*$', re.M)
    for x in pair.finditer(body):
        i, j, a, b = map(int, x.groups())
        if not (0 <= i < n and 0 <= j < n) or (i, j) in hits:
            raise ValueError(f'Bad or duplicate entry in {field}{n}: {(i,j)}')
        hits.add((i,j))
        out[i][j] = Fraction(a, b)
    if len(hits) != len(re.findall(r'^\s*\|\s*\d+,\s*\d+\s*=>', body, flags=re.M)):
        raise ValueError(f'Unparsed entries in {field}{n}')
    return out


def verify(src_path: Path) -> dict:
    text = src_path.read_text(encoding='utf-8')
    tables = {(field, n): parse_matrix(text, field, n) for field in FIELDS for n in (4, 6, 8)}
    nesting = {}
    for field in FIELDS:
        for lo, hi in ((4,6), (6,8)):
            A, B = tables[field, lo], tables[field, hi]
            assert all(A[i][j] == B[i][j] for i in range(lo) for j in range(lo)), (field, lo, hi)
            nesting[f'{field}_{lo}_to_{hi}'] = lo*lo
    rows = {}
    hs = {}
    for n in (4,6,8):
        e = tables['error', n]
        assert all(e[i][j]>=0 and e[i][j]==e[j][i] for i in range(n) for j in range(n))
        rowvals = [sum(row) for row in e]
        maxrow = max(rowvals)
        hs2 = sum(q*q for row in e for q in row)
        rows[str(n)] = {'max_row_fraction':str(maxrow), 'is_below_1_over_100':maxrow < Fraction(1,100)}
        hs[str(n)] = {'frobenius_norm_squared':str(hs2), 'is_below_1_over_100_squared':hs2 < Fraction(1,10000)}
    # Explicit LDL counterexample with both pivots = 1:
    # [[1,t],[t,t*t+1]]=U^T U, U=[[1,t],[0,1]].
    t=10
    x=(-t,1)
    q = (x[0] + t*x[1])**2 + x[1]**2
    norm2=x[0]**2+x[1]**2
    assert q==1 and Fraction(q,norm2)<Fraction(1,20)
    # The LDL pivots of C are NOT obtained by adding the shift to the pivots of C-shift*I.
    shifted=tables['center',4]
    shift=Fraction(3,50)
    pivot2_shifted=(shifted[2][2]-shift)-shifted[2][0]**2/(shifted[0][0]-shift)
    pivot2_center=shifted[2][2]-shifted[2][0]**2/shifted[0][0]
    assert pivot2_center-pivot2_shifted != shift
    # Same prefix, two distinct infinite continuations: differ only at (8,8).
    C8=tables['center',8]
    def A(i,j): return C8[i][j] if i<8 and j<8 else Fraction(0)
    def B(i,j): return A(i,j)+(Fraction(1) if (i,j)==(8,8) else Fraction(0))
    assert all(A(i,j)==B(i,j) for i in range(8) for j in range(8))
    assert A(8,8)!=B(8,8)
    return {
       'schema':'NUMERIA.ES8P16.STRUCTURAL_AUDIT.1',
       'source_sha256':hashlib.sha256(src_path.read_bytes()).hexdigest(),
       'observed_dimensions':[4,6,8],
       'nested_exact_rational_comparisons':nesting,
       'error_row_budgets':rows,
       'frobenius_exploratory_finite':hs,
       'ldl_pivots_not_spectral_gap':{'t':t,'pivots':['1','1'], 'rayleigh_fraction':str(Fraction(q,norm2)), 'fails_margin_1_over_20':True},
       'shifted_ldl_pivot_nonadditivity':{
           'shift':str(shift), 'pivot_at_index2_of_shifted':str(pivot2_shifted),
           'pivot_at_index2_of_center':str(pivot2_center),
           'shift_can_be_added_to_pivots':False},
       'two_distinct_infinite_extensions_agree_on_8x8':True,
       'projection_factor_counterexample':{'B':'2', 'raw_error':'1/100', 'projected_error':'1/25'},
       'scientific_status':'Exact checks of finite certificates; infinite family, Schur, projection norm, uniform spectral gap and analytic Weil enclosure UNPROVED.'
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,default=DEFAULT_SOURCE)
    p.add_argument('--out',type=Path,default=ROOT/'RESULTADOS_ES8P16.json')
    a=p.parse_args()
    r=verify(a.source)
    a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('Exact rational structural checks OK; wrote',a.out)
    print('Frobenius bound 1/100 achieved at N=4,6,8:',[r['frobenius_exploratory_finite'][str(n)]['is_below_1_over_100_squared'] for n in (4,6,8)])

if __name__=='__main__': main()
