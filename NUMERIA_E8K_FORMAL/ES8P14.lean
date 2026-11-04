import Mathlib.Tactic

/-! NUMERIA ES8-P14: algebraic matrix perturbation lemma.
All theorems below concern finite rational matrices. There is no new axiom.
-/

namespace Numeria.ES8P14

def qForm {n : ℕ} (M : Fin n → Fin n → ℚ) (x : Fin n → ℚ) : ℚ :=
  ∑ i : Fin n, ∑ j : Fin n, x i * M i j * x j

/-- If an entry error lies in [-e,e], its cross-term can be bounded
    using only squares, regardless of the signs of the coordinates. -/
theorem error_pair_lower (e t u v : ℚ) (hlo : -e ≤ t) (hhi : t ≤ e) :
    -(e * (u ^ 2 + v ^ 2)) / 2 ≤ u * t * v := by
  have hp : 0 ≤ (e + t) * (u + v) ^ 2 :=
    mul_nonneg (by linarith) (sq_nonneg _)
  have hm : 0 ≤ (e - t) * (u - v) ^ 2 :=
    mul_nonneg (by linarith) (sq_nonneg _)
  nlinarith [hp, hm]

private theorem symmetric_envelope {n : ℕ}
    (eps : Fin n → Fin n → ℚ) (x : Fin n → ℚ)
    (hsym : ∀ i j, eps i j = eps j i) :
    (∑ i : Fin n, ∑ j : Fin n,
      eps i j * (x i ^ 2 + x j ^ 2) / 2) =
    ∑ i : Fin n, (∑ j : Fin n, eps i j) * x i ^ 2 := by
  have hs :
      (∑ i : Fin n, ∑ j : Fin n, eps i j * x j ^ 2) =
      ∑ i : Fin n, ∑ j : Fin n, eps i j * x i ^ 2 := by
    calc
      _ = ∑ j : Fin n, ∑ i : Fin n, eps i j * x j ^ 2 := by
        rw [Finset.sum_comm]
      _ = ∑ j : Fin n, ∑ i : Fin n, eps j i * x j ^ 2 := by
        apply Finset.sum_congr rfl
        intro j _
        apply Finset.sum_congr rfl
        intro i _
        rw [hsym]
      _ = _ := rfl
  calc
    (∑ i : Fin n, ∑ j : Fin n,
        eps i j * (x i ^ 2 + x j ^ 2) / 2) =
      ((∑ i : Fin n, ∑ j : Fin n, eps i j * x i ^ 2) +
       (∑ i : Fin n, ∑ j : Fin n, eps i j * x j ^ 2)) / 2 := by
        simp only [mul_add, add_div, Finset.sum_add_distrib, Finset.sum_div]
    _ = ∑ i : Fin n, ∑ j : Fin n, eps i j * x i ^ 2 := by
      rw [hs]
      ring
    _ = _ := by
      simp only [Finset.sum_mul]

/-- A rational symmetric entrywise error with row sums at most δ
    changes any quadratic form by no more than δ * Σ x_i². -/
theorem qForm_row_error_lower {n : ℕ}
    (M C eps : Fin n → Fin n → ℚ) (x : Fin n → ℚ) (δ : ℚ)
    (hsym : ∀ i j, eps i j = eps j i)
    (hentry : ∀ i j, -eps i j ≤ M i j - C i j ∧ M i j - C i j ≤ eps i j)
    (hrow : ∀ i, (∑ j : Fin n, eps i j) ≤ δ) :
    qForm C x - δ * (∑ i : Fin n, x i ^ 2) ≤ qForm M x := by
  have hb :
      -(∑ i : Fin n, ∑ j : Fin n,
        eps i j * (x i ^ 2 + x j ^ 2) / 2) ≤
      ∑ i : Fin n, ∑ j : Fin n,
        x i * (M i j - C i j) * x j := by
    simp only [← Finset.sum_neg_distrib]
    apply Finset.sum_le_sum
    intro i _
    apply Finset.sum_le_sum
    intro j _
    exact error_pair_lower _ _ _ _ (hentry i j).1 (hentry i j).2
  have hr :
      (∑ i : Fin n, (∑ j : Fin n, eps i j) * x i ^ 2) ≤
      δ * (∑ i : Fin n, x i ^ 2) := by
    calc
      _ ≤ ∑ i : Fin n, δ * x i ^ 2 := by
        apply Finset.sum_le_sum
        intro i _
        exact mul_le_mul_of_nonneg_right (hrow i) (sq_nonneg _)
      _ = _ := by rw [Finset.mul_sum]
  have hd :
      (∑ i : Fin n, ∑ j : Fin n,
        x i * (M i j - C i j) * x j) = qForm M x - qForm C x := by
    simp only [qForm, mul_sub, sub_mul, Finset.sum_sub_distrib]
  rw [symmetric_envelope eps x hsym] at hb
  rw [hd] at hb
  linarith

/-- Adding a scalar multiple of the identity adds δ times the squared norm. -/
theorem diagonal_shift {n : ℕ} (C : Fin n → Fin n → ℚ)
    (x : Fin n → ℚ) (s : ℚ) :
    qForm (fun i j => C i j + if i = j then s else 0) x =
      qForm C x + s * (∑ i : Fin n, x i ^ 2) := by
  unfold qForm
  simp only [mul_add, add_mul, Finset.sum_add_distrib]
  congr 1
  calc
    (∑ i : Fin n, ∑ j : Fin n, x i * (if i = j then s else 0) * x j) =
      ∑ i : Fin n, s * x i ^ 2 := by
        apply Finset.sum_congr rfl
        intro i _
        simp [mul_assoc, pow_two, mul_left_comm, mul_comm]
    _ = _ := by rw [Finset.mul_sum]

/-- Abstract theorem used to transfer the P13 LDL witness to box-enclosed matrices.
    The matrix C-3/50*I must be positive. Error row budget is 1/100. -/
theorem positive_from_LDL_and_box {n : ℕ}
    (M C eps : Fin n → Fin n → ℚ) (x : Fin n → ℚ)
    (hcoarse : 0 < qForm (fun i j => C i j - if i = j then (3 : ℚ) / 50 else 0) x)
    (hsym : ∀ i j, eps i j = eps j i)
    (hentry : ∀ i j, -eps i j ≤ M i j - C i j ∧ M i j - C i j ≤ eps i j)
    (hrow : ∀ i, (∑ j : Fin n, eps i j) ≤ (1 : ℚ) / 100) :
    (1 : ℚ) / 20 * (∑ i : Fin n, x i ^ 2) < qForm M x := by
  let D : Fin n → Fin n → ℚ := fun i j =>
    C i j - if i = j then (3 : ℚ) / 50 else 0
  have hsh : qForm C x = qForm D x + (3 : ℚ) / 50 * (∑ i : Fin n, x i ^ 2) := by
    have hid : ∀ i j, C i j = D i j + if i = j then (3 : ℚ) / 50 else 0 := by
      intro i j
      simp only [D]
      split_ifs <;> ring
    conv_lhs => arg 1; ext i j; rw [hid i j]
    exact diagonal_shift D x (3 / 50)
  have hpert := qForm_row_error_lower M C eps x (1/100) hsym hentry hrow
  change 0 < qForm D x at hcoarse
  linarith

end Numeria.ES8P14
