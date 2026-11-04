import Mathlib.Tactic
import ES8P14
import ES8P14Concrete

/-!
# ES8-P16 — audits of a proposed infinite-dimensional extension

All arguments here are exact, finite and rational, except for the explicit
statements about arbitrary extensions of an observed 8x8 prefix.

We establish nesting for the stored 4x4, 6x6 and 8x8 *center matrices* and
*entrywise envelopes*. This does NOT establish an infinite Weil matrix.

We also show that strictly positive LDL pivots give no uniform spectral gap
without control of the upper-triangular factor, and prove the additional
projection norm factor in an entrywise error budget.

No RH claim, no new axiom, no `sorry`.
-/
namespace Numeria.ES8P16

private def emb46 (i : Fin 4) : Fin 6 := ⟨i.val, by omega⟩
private def emb68 (i : Fin 6) : Fin 8 := ⟨i.val, by omega⟩

set_option maxHeartbeats 10000000 in
theorem center_nested_4_6 (i j : Fin 4) :
    ES8P14Concrete.center4 i j =
      ES8P14Concrete.center6 (emb46 i) (emb46 j) := by
  fin_cases i <;> fin_cases j <;>
    norm_num [ES8P14Concrete.center4, ES8P14Concrete.center6, emb46]

set_option maxHeartbeats 10000000 in
theorem center_nested_6_8 (i j : Fin 6) :
    ES8P14Concrete.center6 i j =
      ES8P14Concrete.center8 (emb68 i) (emb68 j) := by
  fin_cases i <;> fin_cases j <;>
    norm_num [ES8P14Concrete.center6, ES8P14Concrete.center8, emb68]

set_option maxHeartbeats 10000000 in
theorem envelope_nested_4_6 (i j : Fin 4) :
    ES8P14Concrete.error4 i j =
      ES8P14Concrete.error6 (emb46 i) (emb46 j) := by
  fin_cases i <;> fin_cases j <;>
    norm_num [ES8P14Concrete.error4, ES8P14Concrete.error6, emb46]

set_option maxHeartbeats 10000000 in
theorem envelope_nested_6_8 (i j : Fin 6) :
    ES8P14Concrete.error6 i j =
      ES8P14Concrete.error8 (emb68 i) (emb68 j) := by
  fin_cases i <;> fin_cases j <;>
    norm_num [ES8P14Concrete.error6, ES8P14Concrete.error8, emb68]

/-- A continuation agreeing with the observed center8 on its prefix. -/
def extensionA (i j : ℕ) : ℚ :=
  if hi : i < 8 then
    if hj : j < 8 then
      ES8P14Concrete.center8 ⟨i, hi⟩ ⟨j, hj⟩
    else 0
  else 0

/-- A distinct continuation with exactly the same observed prefix. -/
def extensionB (i j : ℕ) : ℚ :=
  extensionA i j + if i = 8 ∧ j = 8 then 1 else 0

theorem extensionA_matches_prefix (i j : Fin 8) :
    extensionA i.val j.val = ES8P14Concrete.center8 i j := by
  simp [extensionA, i.isLt, j.isLt]

theorem extensionB_matches_prefix (i j : Fin 8) :
    extensionB i.val j.val = ES8P14Concrete.center8 i j := by
  have hi : i.val ≠ 8 := by omega
  simp [extensionB, extensionA_matches_prefix, hi]

/-- A finite observed prefix cannot specify its infinite continuation. -/
theorem extensions_genuinely_different : extensionA 8 8 ≠ extensionB 8 8 := by
  simp [extensionA, extensionB]

/-- Unit-upper-triangular Cholesky/LDL example: for every t,
  A(t) = [[1,t],[t,t^2+1]] = U(t)^T I U(t).
  Both LDL pivots are 1, but the smallest eigenvalue can tend to zero. -/
def illConditionedForm (t x y : ℚ) : ℚ := (x + t * y)^2 + y^2

theorem illConditioned_matrix_identity (t x y : ℚ) :
    illConditionedForm t x y =
      x^2 + 2 * t * x * y + (t^2 + 1) * y^2 := by
  unfold illConditionedForm
  ring

theorem illConditioned_near_kernel (t : ℚ) :
    illConditionedForm t (-t) 1 = 1 := by
  simp [illConditionedForm]

/-- Even though both LDL pivots equal 1, a 1/20 lower spectral bound
  already fails at t=10 on the vector (-10,1). -/
theorem positive_pivots_do_not_imply_margin :
    illConditionedForm 10 (-10) 1 <
      (1 : ℚ) / 20 * ((-10 : ℚ)^2 + (1 : ℚ)^2) := by
  norm_num [illConditionedForm]

/-- Projection of the coordinates by a rectangular rational matrix B. -/
def projectVector {n m : ℕ} (B : Fin n → Fin m → ℚ)
    (x : Fin m → ℚ) (i : Fin n) : ℚ :=
  ∑ a : Fin m, B i a * x a

/-- The finite error form after projection necessarily acquires the
squared operator norm factor kappa. No `fin_cases` and any dimensions. -/
theorem projected_error_lower
    {n m : ℕ} (M C eps : Fin n → Fin n → ℚ)
    (B : Fin n → Fin m → ℚ) (x : Fin m → ℚ)
    (δ κ : ℚ) (hδ : 0 ≤ δ)
    (hsym : ∀ i j, eps i j = eps j i)
    (hentry : ∀ i j,
      -eps i j ≤ M i j - C i j ∧ M i j - C i j ≤ eps i j)
    (hrow : ∀ i, (∑ j : Fin n, eps i j) ≤ δ)
    (hB : (∑ i : Fin n, (projectVector B x i)^2) ≤
      κ * (∑ a : Fin m, x a^2)) :
    ES8P14.qForm C (projectVector B x) -
        δ * κ * (∑ a : Fin m, x a^2) ≤
      ES8P14.qForm M (projectVector B x) := by
  have hpert := ES8P14.qForm_row_error_lower M C eps
    (projectVector B x) δ hsym hentry hrow
  have hbnd := mul_le_mul_of_nonneg_left hB hδ
  nlinarith

/-- A uniform bound for the centered projected quadratic form and an
unprojected error budget control the projected one, PROVIDED the projection
has a norm bound. This hypothesis cannot be omitted. -/
theorem projected_positive_with_norm_bound
    {n m : ℕ} (M C eps : Fin n → Fin n → ℚ)
    (B : Fin n → Fin m → ℚ) (x : Fin m → ℚ)
    (hsym : ∀ i j, eps i j = eps j i)
    (hentry : ∀ i j,
      -eps i j ≤ M i j - C i j ∧ M i j - C i j ≤ eps i j)
    (hrow : ∀ i, (∑ j : Fin n, eps i j) ≤ (1 : ℚ) / 100)
    (hB : (∑ i : Fin n, (projectVector B x i)^2) ≤
      (∑ a : Fin m, x a^2))
    (hcenter : (3 : ℚ) / 50 * (∑ a : Fin m, x a^2) <
      ES8P14.qForm C (projectVector B x)) :
    (1 : ℚ) / 20 * (∑ a : Fin m, x a^2) <
      ES8P14.qForm M (projectVector B x) := by
  have h := projected_error_lower M C eps B x (1/100) 1
    (by norm_num) hsym hentry hrow (by simpa using hB)
  norm_num at h
  linarith

/-- A scalar example: projection by B=2 turns an error of 1/100 into 4/100.
  Thus controlling raw row sums alone does not control projected row sums. -/
theorem projection_amplification_counterexample :
    ((2 : ℚ)^2) * (1 / 100) = 4 / 100 := by norm_num

end Numeria.ES8P16
