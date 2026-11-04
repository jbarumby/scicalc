import Mathlib.Tactic

/-!
NUMERIA ES8-P12: formal exact Newton-Cotes moments and abstract LDL
positivity. No extra axioms and no admitted proof.
Analytic digamma, Cauchy, input node enclosures and concrete historical
matrix LDL factorizations are not claimed as formalized.
-/

namespace Numeria.ES8P12

def weightNum (j : Fin 16) : ℤ :=
  match j.val with
  | 0 => 35310023
  | 1 => 265553865
  | 2 => -232936065
  | 3 => 1047777585
  | 4 => -1562840685
  | 5 => 2461884669
  | 6 => -2000332805
  | 7 => 1018807605
  | 8 => 1018807605
  | 9 => -2000332805
  | 10 => 2461884669
  | 11 => -1562840685
  | 12 => 1047777585
  | 13 => -232936065
  | 14 => 265553865
  | _ => 35310023

def weight (j : Fin 16) : ℚ := (weightNum j : ℚ) / 2066448384
def node (j : Fin 16) : ℚ := (j.val : ℚ) / 15

set_option maxHeartbeats 8000000 in
theorem newton_cotes_exact_moments :
    ∀ k : Fin 16,
      (∑ j : Fin 16, weight j * (node j) ^ k.val) =
        (1 : ℚ) / (k.val + 1) := by
  intro k
  fin_cases k <;> norm_num [weight, weightNum, node, Fin.sum_univ_succ]

theorem has_negative_weight : weight ⟨2, by decide⟩ < 0 := by
  norm_num [weight, weightNum]

theorem quadrature_polynomial_exact (c : Fin 16 → ℚ) :
    (∑ j : Fin 16, weight j * (∑ k : Fin 16, c k * node j ^ k.val)) =
      ∑ k : Fin 16, c k / (k.val + 1) := by
  calc
    (∑ j : Fin 16, weight j * (∑ k : Fin 16, c k * node j ^ k.val)) =
        ∑ j : Fin 16, ∑ k : Fin 16, c k * (weight j * node j ^ k.val) := by
          congr 1
          funext j
          rw [Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro k _
          ring
    _ = ∑ k : Fin 16, ∑ j : Fin 16, c k * (weight j * node j ^ k.val) := by
      rw [Finset.sum_comm]
    _ = ∑ k : Fin 16, c k * (∑ j : Fin 16, weight j * node j ^ k.val) := by
      congr 1
      funext k
      rw [Finset.mul_sum]
    _ = ∑ k : Fin 16, c k / (k.val + 1) := by
      apply Finset.sum_congr rfl
      intro k _
      rw [newton_cotes_exact_moments k]
      ring

def factoredQ {n : ℕ} (d : Fin n → ℚ) (U : Fin n → Fin n → ℚ)
    (x : Fin n → ℚ) : ℚ :=
  ∑ i : Fin n, d i * (∑ j : Fin n, U i j * x j) ^ 2

theorem factoredQ_nonneg {n : ℕ} (d : Fin n → ℚ)
    (U : Fin n → Fin n → ℚ) (hd : ∀ i, 0 ≤ d i)
    (x : Fin n → ℚ) : 0 ≤ factoredQ d U x := by
  unfold factoredQ
  apply Finset.sum_nonneg
  intro i _
  exact mul_nonneg (hd i) (sq_nonneg _)

theorem factoredQ_pos {n : ℕ} (d : Fin n → ℚ)
    (U : Fin n → Fin n → ℚ) (hd : ∀ i, 0 < d i)
    (hinj : ∀ x : Fin n → ℚ,
       (∀ i, (∑ j : Fin n, U i j * x j) = 0) → ∀ i, x i = 0)
    (x : Fin n → ℚ) (hx : ∃ i, x i ≠ 0) :
    0 < factoredQ d U x := by
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

theorem positive_quadratic_two {a b c : ℚ}
    (ha : 0 < a) (hdet : 0 < a * c - b ^ 2)
    (x y : ℚ) (hnonzero : x ≠ 0 ∨ y ≠ 0) :
    0 < a * x ^ 2 + 2 * b * x * y + c * y ^ 2 := by
  let q : ℚ := a * x ^ 2 + 2 * b * x * y + c * y ^ 2
  have hid : a * q = (a * x + b * y) ^ 2 + (a * c - b ^ 2) * y ^ 2 := by
    dsimp [q]
    ring
  by_cases hy : y = 0
  · have hx : x ≠ 0 := by
      rcases hnonzero with hx | hy'; exact hx; exact False.elim (hy' hy)
    have hxx : 0 < x ^ 2 := sq_pos_of_ne_zero hx
    have hpos : 0 < a * x ^ 2 := mul_pos ha hxx
    simp [hy, q] at *
    nlinarith
  · have hyy : 0 < y ^ 2 := sq_pos_of_ne_zero hy
    have hp : 0 < (a * c - b ^ 2) * y ^ 2 := mul_pos hdet hyy
    have hsq : 0 ≤ (a * x + b * y) ^ 2 := sq_nonneg _
    have haQ : 0 < a * q := by linarith
    by_contra hnq
    have hq : q ≤ 0 := le_of_not_gt hnq
    have hc : a * q ≤ 0 := mul_nonpos_of_nonneg_of_nonpos (le_of_lt ha) hq
    linarith

end Numeria.ES8P12
