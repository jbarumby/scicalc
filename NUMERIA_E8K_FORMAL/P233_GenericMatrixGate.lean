import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-! NUMERIA ES8-P233: arbitrary finite-dimension real positivity interface.
Analytic enclosure, Gram factor, error row bounds are hypotheses; RH is NOT proved. -/
namespace Numeria.ES8P233

noncomputable def qForm {n : ℕ} (M : Fin n → Fin n → ℝ) (x : Fin n → ℝ) : ℝ :=
  ∑ i : Fin n, ∑ j : Fin n, x i * M i j * x j

noncomputable def factoredQ {n : ℕ} (d : Fin n → ℝ) (U : Fin n → Fin n → ℝ)
    (x : Fin n → ℝ) : ℝ :=
  ∑ k : Fin n, d k * (∑ j : Fin n, U k j * x j) ^ 2

theorem qForm_eq_factored {n : ℕ} (C : Fin n → Fin n → ℝ)
    (d : Fin n → ℝ) (U : Fin n → Fin n → ℝ)
    (h : ∀ i j : Fin n, C i j = ∑ k : Fin n, d k * U k i * U k j)
    (x : Fin n → ℝ) : qForm C x = factoredQ d U x := by
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

theorem factoredQ_nonneg {n : ℕ} (d : Fin n → ℝ) (U : Fin n → Fin n → ℝ)
    (hdiag : ∀ k, 0 ≤ d k) (x : Fin n → ℝ) : 0 ≤ factoredQ d U x := by
  unfold factoredQ
  apply Finset.sum_nonneg
  intro k hk
  exact mul_nonneg (hdiag k) (sq_nonneg _)

theorem qForm_eq_dotProduct {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) (x : Fin n → ℝ) :
    qForm M x = dotProduct x (M.mulVec x) := by
  classical
  change (∑ i : Fin n, ∑ j : Fin n, x i * M i j * x j) =
    ∑ i : Fin n, x i * (∑ j : Fin n, M i j * x j)
  calc
    _ = ∑ i : Fin n, ∑ j : Fin n, x i * (M i j * x j) := by
      apply Finset.sum_congr rfl
      intro i _
      apply Finset.sum_congr rfl
      intro j _
      ring
    _ = _ := by
      apply Finset.sum_congr rfl
      intro i _
      rw [Finset.mul_sum]

theorem diagonal_shift {n : ℕ} (C : Fin n → Fin n → ℝ) (x : Fin n → ℝ) (s : ℝ) :
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
        simp [pow_two, mul_left_comm, mul_comm]
    _ = _ := by rw [Finset.mul_sum]

/-- A bound for one entry error, converted into a quadratic-form bound. -/
theorem error_pair_lower (e t u v : ℝ) (hlo : -e ≤ t) (hhi : t ≤ e) :
    -(e * (u ^ 2 + v ^ 2)) / 2 ≤ u * t * v := by
  have hp : 0 ≤ (e + t) * (u + v) ^ 2 :=
    mul_nonneg (by linarith) (sq_nonneg _)
  have hm : 0 ≤ (e - t) * (u - v) ^ 2 :=
    mul_nonneg (by linarith) (sq_nonneg _)
  nlinarith [hp, hm]

private theorem sum_div_two {n : ℕ} (f : Fin n → ℝ) :
    (∑ i : Fin n, f i) / 2 = ∑ i : Fin n, f i / 2 := by
  rw [div_eq_mul_inv, Finset.sum_mul]
  simp only [div_eq_mul_inv]

private theorem double_sum_div_two {n : ℕ} (f : Fin n → Fin n → ℝ) :
    (∑ i : Fin n, ∑ j : Fin n, f i j) / 2 =
      ∑ i : Fin n, ∑ j : Fin n, f i j / 2 := by
  rw [sum_div_two]
  apply Finset.sum_congr rfl
  intro i _
  exact sum_div_two (f i)

private theorem symmetric_envelope {n : ℕ} (eps : Fin n → Fin n → ℝ) (x : Fin n → ℝ)
    (hsym : ∀ i j, eps i j = eps j i) :
    (∑ i : Fin n, ∑ j : Fin n, eps i j * (x i ^ 2 + x j ^ 2) / 2) =
    ∑ i : Fin n, (∑ j : Fin n, eps i j) * x i ^ 2 := by
  have hs :
      (∑ i : Fin n, ∑ j : Fin n, eps i j * x j ^ 2) =
      ∑ i : Fin n, ∑ j : Fin n, eps i j * x i ^ 2 := by
    calc
      _ = ∑ j : Fin n, ∑ i : Fin n, eps i j * x j ^ 2 := by rw [Finset.sum_comm]
      _ = ∑ j : Fin n, ∑ i : Fin n, eps j i * x j ^ 2 := by
        apply Finset.sum_congr rfl
        intro j _
        apply Finset.sum_congr rfl
        intro i _
        rw [hsym]
      _ = _ := rfl
  calc
    (∑ i : Fin n, ∑ j : Fin n, eps i j * (x i ^ 2 + x j ^ 2) / 2) =
      ∑ i : Fin n, ∑ j : Fin n,
        (eps i j * x i ^ 2 / 2 + eps i j * x j ^ 2 / 2) := by
          apply Finset.sum_congr rfl
          intro i _
          apply Finset.sum_congr rfl
          intro j _
          ring
    _ = (∑ i : Fin n, ∑ j : Fin n, eps i j * x i ^ 2 / 2) +
        (∑ i : Fin n, ∑ j : Fin n, eps i j * x j ^ 2 / 2) := by
          simp only [Finset.sum_add_distrib]
    _ = ((∑ i : Fin n, ∑ j : Fin n, eps i j * x i ^ 2) +
        (∑ i : Fin n, ∑ j : Fin n, eps i j * x j ^ 2)) / 2 := by
          rw [← double_sum_div_two, ← double_sum_div_two]
          ring
    _ = ∑ i : Fin n, ∑ j : Fin n, eps i j * x i ^ 2 := by rw [hs]; ring
    _ = _ := by simp only [Finset.sum_mul]

theorem qForm_row_error_lower {n : ℕ} (M C eps : Fin n → Fin n → ℝ)
    (x : Fin n → ℝ) (delta : ℝ)
    (hsym : ∀ i j, eps i j = eps j i)
    (hentry : ∀ i j, -eps i j ≤ M i j - C i j ∧ M i j - C i j ≤ eps i j)
    (hrow : ∀ i, (∑ j : Fin n, eps i j) ≤ delta) :
    qForm C x - delta * (∑ i : Fin n, x i ^ 2) ≤ qForm M x := by
  have hb :
      -(∑ i : Fin n, ∑ j : Fin n, eps i j * (x i ^ 2 + x j ^ 2) / 2) ≤
      ∑ i : Fin n, ∑ j : Fin n, x i * (M i j - C i j) * x j := by
    simp only [← Finset.sum_neg_distrib]
    apply Finset.sum_le_sum
    intro i _
    apply Finset.sum_le_sum
    intro j _
    simpa only [neg_div] using
      (error_pair_lower (eps i j) (M i j - C i j) (x i) (x j)
        (hentry i j).1 (hentry i j).2)
  have hr :
      (∑ i : Fin n, (∑ j : Fin n, eps i j) * x i ^ 2) ≤
      delta * (∑ i : Fin n, x i ^ 2) := by
    calc
      _ ≤ ∑ i : Fin n, delta * x i ^ 2 := by
        apply Finset.sum_le_sum
        intro i _
        exact mul_le_mul_of_nonneg_right (hrow i) (sq_nonneg _)
      _ = _ := by rw [Finset.mul_sum]
  have hd :
      (∑ i : Fin n, ∑ j : Fin n, x i * (M i j - C i j) * x j) =
      qForm M x - qForm C x := by
    simp only [qForm, mul_sub, sub_mul, Finset.sum_sub_distrib]
  rw [symmetric_envelope eps x hsym] at hb
  rw [hd] at hb
  linarith

/-- A factorization of a shifted centre gives the spectral-gap inequality. -/
theorem center_lower_of_gram {n : ℕ} (C : Matrix (Fin n) (Fin n) ℝ)
    (d : Fin n → ℝ) (U : Fin n → Fin n → ℝ) (shift : ℝ)
    (hfac : ∀ i j : Fin n,
      C i j = (∑ k : Fin n, d k * U k i * U k j) + (if i = j then shift else 0))
    (hpos : ∀ k, 0 ≤ d k) (x : Fin n → ℝ) :
    shift * (∑ i : Fin n, x i ^ 2) ≤ dotProduct x (C.mulVec x) := by
  let G : Matrix (Fin n) (Fin n) ℝ :=
    fun i j => ∑ k : Fin n, d k * U k i * U k j
  have hgf : ∀ i j : Fin n, G i j = ∑ k : Fin n, d k * U k i * U k j := by
    intro i j
    rfl
  have hfact : qForm G x = factoredQ d U x :=
    qForm_eq_factored G d U hgf x
  have hrewrite : qForm C x = qForm G x + shift * (∑ i : Fin n, x i ^ 2) := by
    calc
      qForm C x = qForm (fun i j => G i j + if i = j then shift else 0) x := by
        congr 1
        funext i
        funext j
        exact hfac i j
      _ = _ := diagonal_shift G x shift
  have hge : 0 ≤ factoredQ d U x := factoredQ_nonneg d U hpos x
  rw [← qForm_eq_dotProduct C x, hrewrite, hfact]
  linarith

/-- Universal finite-dimensional positivity theorem under a centre gap and error budget. -/
theorem box_posDef {n : ℕ} (W C : Matrix (Fin n) (Fin n) ℝ)
    (eps : Fin n → Fin n → ℝ) (rho shift : ℝ)
    (hHerm : W.IsHermitian)
    (hbox : ∀ i j, |W i j - C i j| ≤ eps i j)
    (hsym : ∀ i j, eps i j = eps j i)
    (hrow : ∀ i, (∑ j : Fin n, eps i j) ≤ rho)
    (hcenter : ∀ x : Fin n → ℝ,
      shift * (∑ i : Fin n, x i ^ 2) ≤ dotProduct x (C.mulVec x))
    (hmargin : rho < shift) : W.PosDef := by
  classical
  apply Matrix.PosDef.of_dotProduct_mulVec_pos hHerm
  intro x hx
  have hq : 0 < ∑ i : Fin n, x i ^ 2 := by
    by_contra h
    have hz : ∀ i, x i = 0 := by
      intro i
      by_contra hi
      have hp : 0 < x i ^ 2 := sq_pos_of_ne_zero hi
      have hle : x i ^ 2 ≤ ∑ j : Fin n, x j ^ 2 :=
        Finset.single_le_sum (fun j _ => sq_nonneg (x j)) (Finset.mem_univ i)
      exact h (lt_of_lt_of_le hp hle)
    exact hx (funext hz)
  have hentries : ∀ i j,
      -eps i j ≤ W i j - C i j ∧ W i j - C i j ≤ eps i j := by
    intro i j
    exact abs_le.mp (hbox i j)
  have hpert := qForm_row_error_lower W C eps x rho hsym hentries hrow
  rw [qForm_eq_dotProduct C x, qForm_eq_dotProduct W x] at hpert
  have hc := hcenter x
  have hm : 0 < (shift - rho) * (∑ i : Fin n, x i ^ 2) :=
    mul_pos (sub_pos.mpr hmargin) hq
  have hs : (shift - rho) * (∑ i : Fin n, x i ^ 2) ≤
      dotProduct x (W.mulVec x) := by
    nlinarith [hc, hpert]
  exact lt_of_lt_of_le hm hs

/-- Finite arithmetic certificate: Gram diagonal, boxes and gap must be proven. -/
theorem gram_box_posDef {n : ℕ} (W C : Matrix (Fin n) (Fin n) ℝ)
    (d : Fin n → ℝ) (U : Fin n → Fin n → ℝ)
    (eps : Fin n → Fin n → ℝ) (rho shift : ℝ)
    (hHerm : W.IsHermitian)
    (hfac : ∀ i j : Fin n,
      C i j = (∑ k : Fin n, d k * U k i * U k j) + (if i = j then shift else 0))
    (hpos : ∀ k, 0 ≤ d k)
    (hbox : ∀ i j, |W i j - C i j| ≤ eps i j)
    (hsym : ∀ i j, eps i j = eps j i)
    (hrow : ∀ i, (∑ j : Fin n, eps i j) ≤ rho)
    (hmargin : rho < shift) : W.PosDef := by
  apply box_posDef W C eps rho shift hHerm hbox hsym hrow _ hmargin
  intro x
  exact center_lower_of_gram C d U shift hfac hpos x

/-- Explicit interface to analytic intervals. Inclusion of the TRUE entries in
the intervals lower/upper is a stated HYPOTHESIS and is not asserted by this file. -/
theorem analytic_enclosure_to_posDef {n : ℕ}
    (W C lower upper : Matrix (Fin n) (Fin n) ℝ)
    (d : Fin n → ℝ) (U : Fin n → Fin n → ℝ)
    (eps : Fin n → Fin n → ℝ) (rho shift : ℝ)
    (hHerm : W.IsHermitian)
    (hfac : ∀ i j : Fin n,
      C i j = (∑ k : Fin n, d k * U k i * U k j) + (if i = j then shift else 0))
    (hpos : ∀ k, 0 ≤ d k)
    (hAnalyticInclusion : ∀ i j, lower i j ≤ W i j ∧ W i j ≤ upper i j)
    (hLower : ∀ i j, C i j - eps i j ≤ lower i j)
    (hUpper : ∀ i j, upper i j ≤ C i j + eps i j)
    (hsym : ∀ i j, eps i j = eps j i)
    (hrow : ∀ i, (∑ j : Fin n, eps i j) ≤ rho)
    (hmargin : rho < shift) : W.PosDef := by
  have hbox : ∀ i j, |W i j - C i j| ≤ eps i j := by
    intro i j
    apply abs_le.mpr
    constructor
    · have h0 := (hAnalyticInclusion i j).1
      have h1 := hLower i j
      linarith
    · have h0 := (hAnalyticInclusion i j).2
      have h1 := hUpper i j
      linarith
  exact gram_box_posDef W C d U eps rho shift hHerm hfac hpos hbox hsym hrow hmargin

/-- n=200 merely specializes a conditional theorem; no 200x200 input is certified here. -/
theorem P233_200_modes_conditional
    (W C lower upper : Matrix (Fin 200) (Fin 200) ℝ)
    (d : Fin 200 → ℝ) (U : Fin 200 → Fin 200 → ℝ)
    (eps : Fin 200 → Fin 200 → ℝ) (rho shift : ℝ)
    (hHerm : W.IsHermitian)
    (hfac : ∀ i j : Fin 200,
      C i j = (∑ k : Fin 200, d k * U k i * U k j) + (if i = j then shift else 0))
    (hpos : ∀ k, 0 ≤ d k)
    (hAnalyticInclusion : ∀ i j, lower i j ≤ W i j ∧ W i j ≤ upper i j)
    (hLower : ∀ i j, C i j - eps i j ≤ lower i j)
    (hUpper : ∀ i j, upper i j ≤ C i j + eps i j)
    (hsym : ∀ i j, eps i j = eps j i)
    (hrow : ∀ i, (∑ j : Fin 200, eps i j) ≤ rho)
    (hmargin : rho < shift) : W.PosDef := by
  exact analytic_enclosure_to_posDef W C lower upper d U eps rho shift
    hHerm hfac hpos hAnalyticInclusion hLower hUpper hsym hrow hmargin

#print axioms P233_200_modes_conditional
end Numeria.ES8P233
