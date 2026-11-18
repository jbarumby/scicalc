#!/usr/bin/env python3
"""Emit the Lean P231 finite-prefix bridge from P231_PREFIX4_DATA.json."""
from __future__ import annotations
import argparse, hashlib, json
from fractions import Fraction as F
from pathlib import Path


def q(s: str) -> F:
    return F(s)


def qlean(x: F) -> str:
    x=F(x)
    if x.denominator == 1:
        return f'({x.numerator} : ℚ)'
    return f'({x.numerator} : ℚ) / {x.denominator}'


def match2(name: str, rows, type_sig: str) -> str:
    lines=[f'def {name} {type_sig} : ℚ :=', '  match i.val, j.val with']
    for i,row in enumerate(rows):
        for j,value in enumerate(row):
            lines.append(f'  | {i}, {j} => {qlean(q(value))}')
    lines.append('  | _, _ => 0')
    return '\n'.join(lines)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--data',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    d=json.loads(a.data.read_text())
    if d['sector']!='even' or d['matrix_dimension']!=4 or not all(d['checks'].values()):
        raise ValueError('expected validated even-sector prefix4 data')
    sha=hashlib.sha256(a.data.read_bytes()).hexdigest()
    p193_sha=d['P193_center_gzip_sha256']
    p230_lower=d['P230_200_mode_target_lambda_lower_bound_exact']
    center=d['rounded_center_matrix_exact']
    L=d['shifted_LDL_unit_lower_exact']
    piv=d['shifted_LDL_diagonal_exact']
    upper=[[L[j][i] for j in range(4)] for i in range(4)]
    cdef=match2('center4',center,'(i j : Fin 4)')
    udef=match2('upper4',upper,'(i j : Fin 4)')
    pdef='def pivots4 (i : Fin 4) : ℚ :=\n  match i.val with\n'+''.join(f'  | {i} => {qlean(q(v))}\n' for i,v in enumerate(piv))+'  | _ => 0'

    source=f'''import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
NUMERIA ES8-P231: exact finite robustness test for the first four even modes.
Generated from P231_PREFIX4_DATA.json (SHA-256 {sha}).
P193 center gzip SHA-256: {p193_sha}.
P230 full 200-mode finite lower bound: {p230_lower}.
The rational boxes are validated by the accompanying Python extractor against
P227/P229 and the P230 source hash. This file proves a finite conditional
statement; it does not prove the analytic inclusion of Weil entries, any tail
or global positivity, or RH.
-/
namespace Numeria.ES8P231

{cdef}

{pdef}

{udef}

def shift : ℚ := {qlean(q(d['center_shift_exact']))}
def epsilon : ℚ := {qlean(q(d['entrywise_epsilon_exact']))}
def rho : ℚ := {qlean(q(d['quadratic_error_rho_exact']))}

def shifted4 (i j : Fin 4) : ℚ :=
  center4 i j - if i = j then shift else 0

def err4 (_i _j : Fin 4) : ℚ := epsilon

def qForm {{n : ℕ}} (M : Fin n → Fin n → ℚ) (x : Fin n → ℚ) : ℚ :=
  ∑ i : Fin n, ∑ j : Fin n, x i * M i j * x j

def factoredQ {{n : ℕ}} (d : Fin n → ℚ) (U : Fin n → Fin n → ℚ)
    (x : Fin n → ℚ) : ℚ :=
  ∑ k : Fin n, d k * (∑ j : Fin n, U k j * x j) ^ 2

set_option maxHeartbeats 6000000 in
theorem factor_entries4 : ∀ i j : Fin 4,
    shifted4 i j = ∑ k : Fin 4, pivots4 k * upper4 k i * upper4 k j := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    norm_num [shifted4, center4, shift, pivots4, upper4, Fin.sum_univ_succ]

theorem pivots_pos4 : ∀ i : Fin 4, 0 < pivots4 i := by
  intro i
  fin_cases i <;> norm_num [pivots4]

theorem upper_inj4 (x : Fin 4 → ℚ)
    (h : ∀ i : Fin 4, (∑ j : Fin 4, upper4 i j * x j) = 0) :
    ∀ i : Fin 4, x i = 0 := by
  have h3 : x (3 : Fin 4) = 0 := by
    have hi := h (3 : Fin 4)
    simp [upper4, Fin.sum_univ_succ] at hi
    linarith
  have h2 : x (2 : Fin 4) = 0 := by
    have hi := h (2 : Fin 4)
    simp [upper4, Fin.sum_univ_succ, h3] at hi
    linarith
  have h1 : x (1 : Fin 4) = 0 := by
    have hi := h (1 : Fin 4)
    simp [upper4, Fin.sum_univ_succ, h3, h2] at hi
    linarith
  have h0 : x (0 : Fin 4) = 0 := by
    have hi := h (0 : Fin 4)
    simp [upper4, Fin.sum_univ_succ, h3, h2, h1] at hi
    linarith
  intro i
  fin_cases i <;> assumption

theorem qForm_eq_factored {{n : ℕ}} (C : Fin n → Fin n → ℚ)
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

theorem factoredQ_nonneg4 (x : Fin 4 → ℚ) : 0 ≤ factoredQ pivots4 upper4 x := by
  unfold factoredQ
  apply Finset.sum_nonneg
  intro k hk
  exact mul_nonneg (le_of_lt (pivots_pos4 k)) (sq_nonneg _)

theorem qForm_eq_dotProduct (M : Matrix (Fin 4) (Fin 4) ℚ) (x : Fin 4 → ℚ) :
    qForm M x = dotProduct x (M.mulVec x) := by
  classical
  change (∑ i : Fin 4, ∑ j : Fin 4, x i * M i j * x j) =
    ∑ i : Fin 4, x i * (∑ j : Fin 4, M i j * x j)
  calc
    _ = ∑ i : Fin 4, ∑ j : Fin 4, x i * (M i j * x j) := by
      apply Finset.sum_congr rfl
      intro i _
      apply Finset.sum_congr rfl
      intro j _
      ring
    _ = _ := by
      apply Finset.sum_congr rfl
      intro i _
      rw [Finset.mul_sum]

theorem diagonal_shift4 (C : Fin 4 → Fin 4 → ℚ) (x : Fin 4 → ℚ) (s : ℚ) :
    qForm (fun i j => C i j + if i = j then s else 0) x =
      qForm C x + s * (∑ i : Fin 4, x i ^ 2) := by
  unfold qForm
  simp only [mul_add, add_mul, Finset.sum_add_distrib]
  congr 1
  calc
    (∑ i : Fin 4, ∑ j : Fin 4, x i * (if i = j then s else 0) * x j) =
      ∑ i : Fin 4, s * x i ^ 2 := by
        apply Finset.sum_congr rfl
        intro i _
        simp [pow_two, mul_left_comm, mul_comm]
    _ = _ := by rw [Finset.mul_sum]

theorem center_shift_id4 (x : Fin 4 → ℚ) :
    qForm center4 x = qForm shifted4 x + shift * (∑ i : Fin 4, x i ^ 2) := by
  have hid : ∀ i j : Fin 4,
      center4 i j = shifted4 i j + if i = j then shift else 0 := by
    intro i j
    simp [shifted4]
  calc
    qForm center4 x = qForm (fun i j => shifted4 i j + if i = j then shift else 0) x := by
      congr 1
      funext i
      funext j
      exact hid i j
    _ = _ := diagonal_shift4 shifted4 x shift

theorem center_lower4 (x : Fin 4 → ℚ) :
    shift * (∑ i : Fin 4, x i ^ 2) ≤ qForm center4 x := by
  have hf := qForm_eq_factored shifted4 pivots4 upper4 factor_entries4 x
  rw [center_shift_id4, hf]
  have hn := factoredQ_nonneg4 x
  linarith

/-- A bound for one entry error, converted into a quadratic-form bound. -/
theorem error_pair_lower (e t u v : ℚ) (hlo : -e ≤ t) (hhi : t ≤ e) :
    -(e * (u ^ 2 + v ^ 2)) / 2 ≤ u * t * v := by
  have hp : 0 ≤ (e + t) * (u + v) ^ 2 :=
    mul_nonneg (by linarith) (sq_nonneg _)
  have hm : 0 ≤ (e - t) * (u - v) ^ 2 :=
    mul_nonneg (by linarith) (sq_nonneg _)
  nlinarith [hp, hm]

private theorem sum_div_two (f : Fin 4 → ℚ) :
    (∑ i : Fin 4, f i) / 2 = ∑ i : Fin 4, f i / 2 := by
  rw [div_eq_mul_inv, Finset.sum_mul]
  simp only [div_eq_mul_inv]

private theorem double_sum_div_two (f : Fin 4 → Fin 4 → ℚ) :
    (∑ i : Fin 4, ∑ j : Fin 4, f i j) / 2 =
      ∑ i : Fin 4, ∑ j : Fin 4, f i j / 2 := by
  rw [sum_div_two]
  apply Finset.sum_congr rfl
  intro i _
  exact sum_div_two (f i)

private theorem symmetric_envelope (eps : Fin 4 → Fin 4 → ℚ) (x : Fin 4 → ℚ)
    (hsym : ∀ i j, eps i j = eps j i) :
    (∑ i : Fin 4, ∑ j : Fin 4, eps i j * (x i ^ 2 + x j ^ 2) / 2) =
    ∑ i : Fin 4, (∑ j : Fin 4, eps i j) * x i ^ 2 := by
  have hs :
      (∑ i : Fin 4, ∑ j : Fin 4, eps i j * x j ^ 2) =
      ∑ i : Fin 4, ∑ j : Fin 4, eps i j * x i ^ 2 := by
    calc
      _ = ∑ j : Fin 4, ∑ i : Fin 4, eps i j * x j ^ 2 := by rw [Finset.sum_comm]
      _ = ∑ j : Fin 4, ∑ i : Fin 4, eps j i * x j ^ 2 := by
        apply Finset.sum_congr rfl
        intro j _
        apply Finset.sum_congr rfl
        intro i _
        rw [hsym]
      _ = _ := rfl
  calc
    (∑ i : Fin 4, ∑ j : Fin 4, eps i j * (x i ^ 2 + x j ^ 2) / 2) =
      ∑ i : Fin 4, ∑ j : Fin 4,
        (eps i j * x i ^ 2 / 2 + eps i j * x j ^ 2 / 2) := by
          apply Finset.sum_congr rfl
          intro i _
          apply Finset.sum_congr rfl
          intro j _
          ring
    _ = (∑ i : Fin 4, ∑ j : Fin 4, eps i j * x i ^ 2 / 2) +
        (∑ i : Fin 4, ∑ j : Fin 4, eps i j * x j ^ 2 / 2) := by
          simp only [Finset.sum_add_distrib]
    _ = ((∑ i : Fin 4, ∑ j : Fin 4, eps i j * x i ^ 2) +
        (∑ i : Fin 4, ∑ j : Fin 4, eps i j * x j ^ 2)) / 2 := by
          rw [← double_sum_div_two, ← double_sum_div_two]
          ring
    _ = ∑ i : Fin 4, ∑ j : Fin 4, eps i j * x i ^ 2 := by rw [hs]; ring
    _ = _ := by simp only [Finset.sum_mul]

theorem qForm_row_error_lower (M C eps : Fin 4 → Fin 4 → ℚ)
    (x : Fin 4 → ℚ) (delta : ℚ)
    (hsym : ∀ i j, eps i j = eps j i)
    (hentry : ∀ i j, -eps i j ≤ M i j - C i j ∧ M i j - C i j ≤ eps i j)
    (hrow : ∀ i, (∑ j : Fin 4, eps i j) ≤ delta) :
    qForm C x - delta * (∑ i : Fin 4, x i ^ 2) ≤ qForm M x := by
  have hb :
      -(∑ i : Fin 4, ∑ j : Fin 4, eps i j * (x i ^ 2 + x j ^ 2) / 2) ≤
      ∑ i : Fin 4, ∑ j : Fin 4, x i * (M i j - C i j) * x j := by
    simp only [← Finset.sum_neg_distrib]
    apply Finset.sum_le_sum
    intro i _
    apply Finset.sum_le_sum
    intro j _
    simpa only [neg_div] using
      (error_pair_lower (eps i j) (M i j - C i j) (x i) (x j)
        (hentry i j).1 (hentry i j).2)
  have hr :
      (∑ i : Fin 4, (∑ j : Fin 4, eps i j) * x i ^ 2) ≤
      delta * (∑ i : Fin 4, x i ^ 2) := by
    calc
      _ ≤ ∑ i : Fin 4, delta * x i ^ 2 := by
        apply Finset.sum_le_sum
        intro i _
        exact mul_le_mul_of_nonneg_right (hrow i) (sq_nonneg _)
      _ = _ := by rw [Finset.mul_sum]
  have hd :
      (∑ i : Fin 4, ∑ j : Fin 4, x i * (M i j - C i j) * x j) =
      qForm M x - qForm C x := by
    simp only [qForm, mul_sub, sub_mul, Finset.sum_sub_distrib]
  rw [symmetric_envelope eps x hsym] at hb
  rw [hd] at hb
  linarith

/-- This is P120's finite conditional bridge: hypotheses must be proved elsewhere. -/
theorem matrix_bridge (W C E : Matrix (Fin 4) (Fin 4) ℚ)
    (hW_hermitian : W.IsHermitian)
    (h_decomp : W = C + E)
    (h_margin : rho < shift)
    (h_center : ∀ x : Fin 4 → ℚ, x ≠ 0 →
      shift * (∑ i, x i ^ 2) ≤ dotProduct x (C.mulVec x))
    (h_error : ∀ x : Fin 4 → ℚ, x ≠ 0 →
      -(rho * (∑ i, x i ^ 2)) ≤ dotProduct x (E.mulVec x)) : W.PosDef := by
  classical
  apply Matrix.PosDef.of_dotProduct_mulVec_pos hW_hermitian
  intro x hx
  have hq : 0 < ∑ i, x i ^ 2 := by
    by_contra h
    have hz : ∀ i, x i = 0 := by
      intro i
      by_contra hi
      have hp : 0 < x i ^ 2 := sq_pos_of_ne_zero hi
      have hle : x i ^ 2 ≤ ∑ j, x j ^ 2 :=
        Finset.single_le_sum (fun j _ => sq_nonneg (x j)) (Finset.mem_univ i)
      exact h (lt_of_lt_of_le hp hle)
    exact hx (funext hz)
  have hsplit : dotProduct x (W.mulVec x) =
      dotProduct x (C.mulVec x) + dotProduct x (E.mulVec x) := by
    rw [h_decomp, Matrix.add_mulVec, dotProduct_add]
  have hpos : 0 < dotProduct x (W.mulVec x) := by
    rw [hsplit]
    have hm := h_center x hx
    have he := h_error x hx
    have hmargin : 0 < (shift - rho) * (∑ i, x i ^ 2) :=
      mul_pos (sub_pos.mpr h_margin) hq
    have hsum : (shift - rho) * (∑ i, x i ^ 2) ≤
        dotProduct x (C.mulVec x) + dotProduct x (E.mulVec x) := by
      nlinarith [hm, he]
    exact lt_of_lt_of_le hmargin hsum
  simpa using hpos

/-- Every rational matrix in the P227 prefix box is positive definite,
provided it is Hermitian. P227/P229-to-Weil analytic inclusion is external. -/
theorem P231_box_positive (W : Matrix (Fin 4) (Fin 4) ℚ)
    (hW_hermitian : W.IsHermitian)
    (hbox : ∀ i j : Fin 4, |W i j - center4 i j| ≤ epsilon) : W.PosDef := by
  let C : Matrix (Fin 4) (Fin 4) ℚ := center4
  let E : Matrix (Fin 4) (Fin 4) ℚ := W - C
  have hdecomp : W = C + E := by
    ext i j
    change W i j = center4 i j + (W i j - center4 i j)
    ring
  have hentry : ∀ i j : Fin 4,
      -err4 i j ≤ W i j - center4 i j ∧ W i j - center4 i j ≤ err4 i j := by
    intro i j
    have h := abs_le.mp (hbox i j)
    simpa [err4, epsilon] using h
  have hrows : ∀ i : Fin 4, (∑ j : Fin 4, err4 i j) ≤ rho := by
    intro i
    norm_num [err4, epsilon, rho, Fin.sum_univ_succ]
  have hcenter : ∀ x : Fin 4 → ℚ, x ≠ 0 →
      shift * (∑ i, x i ^ 2) ≤ dotProduct x (C.mulVec x) := by
    intro x hx
    rw [← qForm_eq_dotProduct C x]
    simpa [C] using center_lower4 x
  have herror : ∀ x : Fin 4 → ℚ, x ≠ 0 →
      -(rho * (∑ i, x i ^ 2)) ≤ dotProduct x (E.mulVec x) := by
    intro x hx
    have hpert := qForm_row_error_lower W center4 err4 x rho
      (by intro i j; rfl) hentry hrows
    have hdiff : qForm E x = qForm W x - qForm center4 x := by
      change qForm (fun i j => W i j - center4 i j) x = _
      simp only [qForm, mul_sub, sub_mul, Finset.sum_sub_distrib]
    rw [← qForm_eq_dotProduct E x, hdiff]
    nlinarith [hpert]
  have hm : rho < shift := by norm_num [rho, shift]
  exact matrix_bridge W C E hW_hermitian hdecomp hm hcenter herror

end Numeria.ES8P231
'''
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(source)
    print(f'generated {a.out} from {a.data}; sha256={sha}; bytes={len(source.encode())}')

if __name__=='__main__': main()
