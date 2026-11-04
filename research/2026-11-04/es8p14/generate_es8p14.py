#!/usr/bin/env python3
"""Generate Lean source with EXACT 90-digit original intervals and rational error envelopes.
Reads the ES8-P13 JSON only. Does not infer analytic validity of its endpoints.
"""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from math import ceil
import json
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'RESULTADOS_ES8P13.json'
SCALE=10**90
ROUND=10**6


def rat(v):
    v=F(v)
    if v.denominator==1:return f'({v.numerator} : ℚ)'
    return f'({v.numerator} : ℚ) / {v.denominator}'


def matrix_def(name, vals, n):
    out=[f'def {name} (i j : Fin {n}) : ℚ :=', '  match i.val, j.val with']
    for i in range(n):
        for j in range(n):
            v=F(vals[i][j])
            if not v:continue
            out.append(f'  | {i}, {j} => {rat(v)}')
    out.append('  | _, _ => 0')
    return '\n'.join(out)


def main():
    d=json.loads(SOURCE.read_text())
    out=['import Mathlib.Tactic','import ES8P13','import ES8P14','',
    '/-! ES8-P14: exact finite rational interval endpoint witnesses.',
    'Numerical analytical inclusion of the Weil integrals remains an external',
    'premise. Explicit matrices are rational; no claim about all real matrices or RH. -/',
    'namespace Numeria.ES8P14Concrete']
    results={'schema':'NUMERIA.ES8P14.CONCRETE_BRIDGE.1',
             'source_ES8P13_sha256':sha256(SOURCE.read_bytes()).hexdigest(),
             'raw_scale':str(SCALE),'envelope_scale':str(ROUND),'margin':'1/20',
             'shift':'3/50','budget':'1/100','dimensions':{}}
    for N in (6,8,10):
        n=N-2
        a=d[f'N{N}']
        c=[[F(v) for v in row] for row in a['matrix_coarse']]
        lo=[[F(v,SCALE) for v in row] for row in a['source_interval_raw_integer_lower']]
        hi=[[F(v,SCALE) for v in row] for row in a['source_interval_raw_integer_upper']]
        exact=[[max(abs(lo[i][j]-c[i][j]),abs(hi[i][j]-c[i][j])) for j in range(n)] for i in range(n)]
        eps=[[F((v*ROUND).__ceil__(),ROUND) for v in row] for row in exact]
        assert all(eps[i][j]==eps[j][i] for i in range(n) for j in range(n)),N
        assert all(lo[i][j]<=hi[i][j] for i in range(n) for j in range(n)),N
        assert all(sum(eps[i],F(0)) < F(1,100) for i in range(n)),N
        for i in range(n):
            for j in range(n):
                assert lo[i][j]-c[i][j]>=-eps[i][j],(N,i,j)
                assert hi[i][j]-c[i][j]<=eps[i][j],(N,i,j)
        out.extend(['',f'/-! ## N={N}, dimension={n}. -/',
                   matrix_def(f'center{n}',c,n),
                   matrix_def(f'lower{n}',lo,n),
                   matrix_def(f'upper{n}',hi,n),
                   matrix_def(f'error{n}',eps,n)])
        out.append(f'''
set_option maxHeartbeats 3000000 in
theorem center_id{n} (i j : Fin {n}) :
    center{n} i j = Numeria.ES8P13.matrix{n} i j +
        (if i = j then (3 : ℚ) / 50 else 0) := by
  fin_cases i <;> fin_cases j <;> norm_num [center{n}, Numeria.ES8P13.matrix{n}]

set_option maxHeartbeats 4000000 in
theorem raw_box{n} : ∀ i j : Fin {n},
    -error{n} i j ≤ lower{n} i j - center{n} i j ∧
      upper{n} i j - center{n} i j ≤ error{n} i j := by
  intro i j
  fin_cases i <;> fin_cases j <;> norm_num [lower{n}, upper{n}, center{n}, error{n}]

theorem error_symmetric{n} : ∀ i j : Fin {n}, error{n} i j = error{n} j i := by
  intro i j
  fin_cases i <;> fin_cases j <;> norm_num [error{n}]

set_option maxHeartbeats 3000000 in
theorem error_rows{n} : ∀ i : Fin {n},
    (∑ j : Fin {n}, error{n} i j) ≤ (1 : ℚ) / 100 := by
  intro i
  fin_cases i <;> norm_num [error{n}, Fin.sum_univ_succ]

/-- Rational matrix with entries in the original interval boxes is positive
    with margin 1/20, CONDITIONAL on the source intervals being analytic enclosures. -/
theorem box_positive{n} (M : Fin {n} → Fin {n} → ℚ)
    (hbox : ∀ i j : Fin {n}, lower{n} i j ≤ M i j ∧ M i j ≤ upper{n} i j)
    (x : Fin {n} → ℚ) (hx : ∃ i, x i ≠ 0) :
    (1 : ℚ) / 20 * (∑ i : Fin {n}, x i ^ 2) < Numeria.ES8P14.qForm M x := by
  have hentry : ∀ i j : Fin {n},
      -error{n} i j ≤ M i j - center{n} i j ∧
       M i j - center{n} i j ≤ error{n} i j := by
    intro i j
    have hr := raw_box{n} i j
    have hb := hbox i j
    constructor <;> linarith
  have hcoarse :
      0 < Numeria.ES8P14.qForm
         (fun i j => center{n} i j - if i = j then (3 : ℚ) / 50 else 0) x := by
    have heq : (fun i j => center{n} i j - if i = j then (3 : ℚ) / 50 else 0) =
        Numeria.ES8P13.matrix{n} := by
      funext i j
      rw [center_id{n} i j]
      split_ifs <;> ring
    rw [heq]
    exact Numeria.ES8P13.concrete_positive{n} x hx
  exact Numeria.ES8P14.positive_from_LDL_and_box M center{n} error{n} x
    hcoarse error_symmetric{n} hentry error_rows{n}
''')
        results['dimensions'][str(N)]={'projected_dimension':n,
          'max_exact_row_error':str(max((sum(row,F(0)) for row in exact))),
          'max_rounded_row_error':str(max((sum(row,F(0)) for row in eps))),
          'max_rounded_row_error_decimal':str(float(max((sum(row,F(0)) for row in eps)))),
          'max_entry_rounded':str(max(max(row) for row in eps)),
          'number_of_raw_endpoints':2*n*n,
          'source_interval_sha256':sha256(json.dumps({'lower':a['source_interval_raw_integer_lower'],'upper':a['source_interval_raw_integer_upper']},sort_keys=True,separators=(',',':')).encode()).hexdigest()}
    out.append('\nend Numeria.ES8P14Concrete\n')
    (ROOT/'ES8P14Concrete.lean').write_text('\n'.join(out),encoding='utf-8')
    (ROOT/'RESULTADOS_ES8P14.json').write_text(json.dumps(results,sort_keys=True,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('Lean concrete chars',len('\n'.join(out)))
    for N,r in results['dimensions'].items():print('N',N,'row error',r['max_rounded_row_error'])

if __name__=='__main__':main()
