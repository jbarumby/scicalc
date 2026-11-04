#!/usr/bin/env python3
"""Generate a CHECKABLE concrete Lean LDL/SOS witness from ES8-P10.
Inputs are original P9 fixedpoint intervals, P10 analytic constants and code.
No approximate floating point is used for certificate construction or validation.
The scientific result is conditional on the analytic inclusion lemmas and prime normalization.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
import json,hashlib,sys
ROOT=Path(__file__).resolve().parent
UPSTREAM=ROOT if (ROOT/'es8p10_lemmas.py').is_file() else ROOT.parents[2]/'research/2026-11-03/es8p10'
sys.path.insert(0,str(UPSTREAM))
sys.path.insert(0,str(ROOT))
from scratch import matrix
from es8p10_lemmas import new_intervals
S=10000
SHIFT=F(3,50)
MAX_ERROR=F(1,100)

def rational(x):
 x=F(x)
 if x.denominator==1:return str(x.numerator)
 return f'({x.numerator} : ℚ) / {x.denominator}'

def exact_ldl(C):
 n=len(C);L=[[F(int(i==j)) for j in range(n)]for i in range(n)];d=[]
 for i in range(n):
  v=C[i][i]-(SHIFT)-sum((L[i][k]*L[i][k]*d[k] for k in range(i)),F(0))
  assert v>0
  d.append(v)
  for j in range(i+1,n):
   L[j][i]=(C[j][i]-sum((L[j][k]*L[i][k]*d[k]for k in range(i)),F(0)))/v
 U=[[L[j][i] for j in range(n)]for i in range(n)]
 for i in range(n):
  for j in range(n):
   assert sum((d[k]*U[k][i]*U[k][j] for k in range(n)),F(0))==C[i][j]-(SHIFT if i==j else 0)
 return d,U

def lean_matrix_def(name,A,n):
 lines=[f'def {name} (i j : Fin {n}) : ℚ :=','  match i.val, j.val with']
 for i in range(n):
  for j in range(n):
   if A[i][j]!=0:lines.append(f'  | {i}, {j} => {rational(A[i][j])}')
 lines.append('  | _, _ => 0')
 return '\n'.join(lines)

def lean_diag_def(name,d,n):
 lines=[f'def {name} (i : Fin {n}) : ℚ :=','  match i.val with']
 for i,v in enumerate(d):lines.append(f'  | {i} => {rational(v)}')
 lines.append('  | _ => 0')
 return '\n'.join(lines)

def make_lean():
 records=new_intervals()
 inputs={}
 src_files=['RESULTADOS_ES8P9_INTEGRAIS.json','RESULTADOS_ES8P10.json']
 for fn in src_files: inputs[fn]=hashlib.sha256((UPSTREAM/fn).read_bytes()).hexdigest()
 A10=matrix(10,records)
 for n in (6,8):
  sub=matrix(n,records)
  for i in range(n-2):
   for j in range(n-2):
    assert (A10[i][j].lo,A10[i][j].hi)==(sub[i][j].lo,sub[i][j].hi),f'{n},{i},{j} do not share raw intervals'
 C=[[F(round(A10[i][j].mid()*S),S)for j in range(8)]for i in range(8)]
 assert all(C[i][j]==C[j][i]for i in range(8)for j in range(8))
 data={'schema':'NUMERIA.ES8P13.CONCRETE_LDL.1','base_from':'ES8-P10 finite interval matrix + ES8-P9 original nodal sums','source_sha256':inputs,'rounding_scale':S,'shift_exact':str(SHIFT),'margin_original_exact':'1/20','residual_budget_exact':str(MAX_ERROR),'matrix_sizes':[4,6,8],'numerical_approximations_not_used_for_proofs':True}
 text=[]
 text.append('import Mathlib.Tactic\n\n/-!\nNUMERIA–RIEMANN ES8-P13: concrete exact rational LDL certificates.\n\nThis proves only the POSITIVITY OF THE EXPLICIT RATIONAL COARSE MATRICES Cn - 3/50 * I.\nThe exact 90-digit original arithmetic endpoints are checked by a separate Python\nreference replay: max absolute row sum of A_n - Cn <= 1/100.\nThe bridge from the special-function integral to those endpoints, and the\ngeneral perturbation lemma, are not formalized by this file.\nThus there is NO Lean proof of the complete numerical-to-analytic chain or RH.\n-/\n\nnamespace Numeria.ES8P13\n')
 text.append('''
/-- For the actual *finite rational* matrix represented by its LDL witness. -/
def qForm {n : ℕ} (C : Fin n → Fin n → ℚ) (x : Fin n → ℚ) : ℚ :=
  ∑ i : Fin n, ∑ j : Fin n, x i * C i j * x j

def factoredQ {n : ℕ} (d : Fin n → ℚ) (U : Fin n → Fin n → ℚ)
    (x : Fin n → ℚ) : ℚ :=
  ∑ k : Fin n, d k * (∑ i : Fin n, U k i * x i) ^ 2

/-- Pure algebra: a coefficientwise LDL identity gives equality of quadratic forms. -/
theorem qForm_eq_factored {n : ℕ} (C : Fin n → Fin n → ℚ)
    (d : Fin n → ℚ) (U : Fin n → Fin n → ℚ)
    (h : ∀ i j : Fin n, C i j = ∑ k : Fin n, d k * U k i * U k j)
    (x : Fin n → ℚ) : qForm C x = factoredQ d U x := by
  unfold qForm factoredQ
  simp_rw [h]
  calc
    (∑ i : Fin n, ∑ j : Fin n,
        x i * (∑ k : Fin n, d k * U k i * U k j) * x j) =
      ∑ i : Fin n, ∑ j : Fin n, ∑ k : Fin n,
        d k * U k i * x i * U k j * x j := by
          apply Finset.sum_congr rfl
          intro i _
          apply Finset.sum_congr rfl
          intro j _
          rw [Finset.mul_sum, Finset.sum_mul]
          exact Finset.sum_congr rfl (by intro k _; ring)
    _ = ∑ k : Fin n, ∑ i : Fin n, ∑ j : Fin n,
          d k * U k i * x i * U k j * x j := by
            calc
              _ = ∑ i : Fin n, ∑ k : Fin n, ∑ j : Fin n,
                  d k * U k i * x i * U k j * x j := by
                    apply Finset.sum_congr rfl
                    intro i _
                    rw [Finset.sum_comm]
              _ = _ := by rw [Finset.sum_comm]
    _ = ∑ k : Fin n, d k * (∑ i : Fin n, U k i * x i) ^ 2 := by
        apply Finset.sum_congr rfl
        intro k _
        simp only [pow_two, Finset.mul_sum, Finset.sum_mul]
        apply Finset.sum_congr rfl
        intro i _
        apply Finset.sum_congr rfl
        intro j _
        ring

/-- The standard finite SOS strict positivity criterion. -/
theorem factoredQ_pos {n : ℕ} (d : Fin n → ℚ)
    (U : Fin n → Fin n → ℚ) (hd : ∀ i, 0 < d i)
    (hinj : ∀ x : Fin n → ℚ,
       (∀ i, (∑ j : Fin n, U i j * x j) = 0) → ∀ i, x i = 0)
    (x : Fin n → ℚ) (hx : ∃ i, x i ≠ 0) : 0 < factoredQ d U x := by
  have hy : ∃ i : Fin n, (∑ j : Fin n, U i j * x j) ≠ 0 := by
    by_contra h
    push Not at h
    obtain ⟨j, hj⟩ := hx
    exact hj (hinj x h j)
  obtain ⟨i, hi⟩ := hy
  unfold factoredQ
  apply Finset.sum_pos'
  · intro k _
    exact mul_nonneg (le_of_lt (hd k)) (sq_nonneg _)
  · refine ⟨i, Finset.mem_univ _, ?_⟩
    exact mul_pos (hd i) (sq_pos_of_ne_zero hi)
''')
 for n in (4,6,8):
  d,U=exact_ldl([r[:n]for r in C[:n]])
  mid=[r[:n]for r in C[:n]]
  shifted=[[mid[i][j]-(SHIFT if i==j else 0) for j in range(n)]for i in range(n)]
  errors=[[max(abs(A10[i][j].min()-mid[i][j]),abs(A10[i][j].max()-mid[i][j]))for j in range(n)]for i in range(n)]
  row_errors=[sum(row,F(0))for row in errors]
  assert max(row_errors) < MAX_ERROR
  data[f'N{n+2}']={'dimension':n, 'matrix_coarse':[[str(v)for v in row]for row in mid],
         'ldl_pivots_exact':[str(v)for v in d], 'upper_unit_factor':[[str(v)for v in row]for row in U],
         'max_abs_row_error_exact':str(max(row_errors)),'max_abs_row_error_float_for_reading_only':float(max(row_errors)),
         'minimum_shifted_pivot_float_for_reading_only':float(min(d)),
         'source_interval_raw_integer_lower':[[A10[i][j].lo for j in range(n)] for i in range(n)],
         'source_interval_raw_integer_upper':[[A10[i][j].hi for j in range(n)] for i in range(n)]}
  # Define precise matrix with shifted diagonal; avoid if i=j where Lean needs simp by cases.
  text.append(f'\n/-! ## Concrete projected matrix N={n+2}, projected dimension {n} -/\n')
  text.append(lean_matrix_def(f'matrix{n}',shifted,n))
  text.append(lean_diag_def(f'pivots{n}',d,n))
  text.append(lean_matrix_def(f'upper{n}',U,n))
  text.append(f'''\nset_option maxHeartbeats 6000000 in
theorem factor_entries{n} :
    ∀ i j : Fin {n}, matrix{n} i j =
      ∑ k : Fin {n}, pivots{n} k * upper{n} k i * upper{n} k j := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    norm_num [matrix{n}, pivots{n}, upper{n}, Fin.sum_univ_succ]

theorem pivots_pos{n} : ∀ i : Fin {n}, 0 < pivots{n} i := by
  intro i
  fin_cases i <;> norm_num [pivots{n}]
''')
  text.append(f'\ntheorem upper_inj{n} (x : Fin {n} → ℚ)\n    (h : ∀ i : Fin {n}, (∑ j : Fin {n}, upper{n} i j * x j) = 0) :\n    ∀ i : Fin {n}, x i = 0 := by\n')
  for i in reversed(range(n)):
   simpqs=[f'h{x}' for x in reversed(range(i+1,n))]
   text.append(f'  have h{i} : x ⟨{i}, by decide⟩ = 0 := by\n')
   text.append(f'    have hi := h ⟨{i}, by decide⟩\n')
   text.append(f'    simp [upper{n}, Fin.sum_univ_succ'+(', '+', '.join(simpqs) if simpqs else '')+'] at hi\n')
   text.append('    linarith\n')
  text.append('  intro i\n  fin_cases i\n')
  text.append('  all_goals assumption\n')
  text.append(f'''\ntheorem concrete_positive{n} (x : Fin {n} → ℚ)
    (hx : ∃ i, x i ≠ 0) : 0 < qForm matrix{n} x := by
  rw [qForm_eq_factored matrix{n} pivots{n} upper{n} factor_entries{n}]
  exact factoredQ_pos pivots{n} upper{n} pivots_pos{n} upper_inj{n} x hx
''')
 # prove value w.r.t. legacy center using explicit 4/6/8 rational identity? already done
 text.append('\nend Numeria.ES8P13\n')
 out_lean=ROOT/'ES8P13.lean' if (ROOT/'ES8P13.lean').exists() else ROOT.parents[2]/'NUMERIA_E8K_FORMAL/ES8P13.lean'
 out_lean.write_text('\n'.join(text),encoding='utf-8')
 (ROOT/'RESULTADOS_ES8P13.json').write_text(json.dumps(data,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
 print('generated Lean chars:',out_lean.stat().st_size, 'JSON size:',(ROOT/'RESULTADOS_ES8P13.json').stat().st_size)
 for n in (4,6,8):
  r=data[f'N{n+2}'];print('N',n+2,'max row error',r['max_abs_row_error_float_for_reading_only'],'pivot min',r['minimum_shifted_pivot_float_for_reading_only'])

if __name__=='__main__':make_lean()
