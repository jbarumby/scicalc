import Mathlib.Tactic

/-!
# NUMERIA ES8-P18 — rational projection tail and uniformity obligations

No `sorry`, `admit` or new axioms. We establish exact rational constants,
a generic Cauchy–Schwarz-to-projection inequality, and a telescoping square
series budget. The analytic E7 integral inequality and the identification of
the actual Fourier coefficients with an infinite Hilbert operator are NOT
proved here. The latter is justified separately in the mathematical report.
-/

namespace Numeria.ES8P18

/-- An elementary reciprocal-square estimate which telescopes. -/
theorem reciprocal_sq_step (q : ℚ) (hq : 0 < q) :
    (1 : ℚ) / (q + 1) ^ 2 ≤ 1 / q - 1 / (q + 1) := by
  have hq1 : 0 < q + 1 := by linarith
  have hr : (1 : ℚ) / q - 1 / (q + 1) = 1 / (q * (q + 1)) := by
    field_simp [ne_of_gt hq, ne_of_gt hq1]
    <;> ring
  rw [hr]
  apply (div_le_div_iff₀ (pow_pos hq1 2) (mul_pos hq hq1)).2
  nlinarith

/-- For every finite number of terms, the tail beginning at m=4 has
    rational sum at most 1/3; no limiting axiom is used. -/
theorem inverse_sq_tail_four (n : ℕ) :
    (∑ j ∈ Finset.range n, (1 : ℚ) / ((j : ℚ) + 4) ^ 2) ≤
      1 / 3 - 1 / ((n : ℚ) + 3) := by
  induction n with
  | zero => norm_num
  | succ n ih =>
    rw [Finset.sum_range_succ]
    have hs := reciprocal_sq_step ((n : ℚ) + 3) (by positivity)
    have hcast : (((n + 1 : ℕ) : ℚ) + 3) = ((n : ℚ) + 4) := by
      push_cast
      ring
    rw [hcast]
    linarith

/-- Separating m=2 and m=3, and bounding all m>=4 by telescoping,
    the EVEN reciprocal-square coefficient budget is 25/144. -/
theorem even_inverse_square_budget (n : ℕ) :
    (1 : ℚ) / 16 + 1 / 36 +
      (1 / 4) * (∑ j ∈ Finset.range n,
        (1 : ℚ) / ((j : ℚ) + 4) ^ 2) ≤ 25 / 144 := by
  have h := inverse_sq_tail_four n
  have hpos : 0 ≤ (1 : ℚ) / ((n : ℚ) + 3) := by positivity
  nlinarith

/-- Coefficient and spectral constant calculations are exact rationals. -/
theorem improved_parity_constants :
    ((409 : ℚ) / 400) ^ 2 / 4 = 167281 / 640000 ∧
    ((1609 : ℚ) / 800) ^ 2 * (25 / 144) = 2588881 / 3686400 := by
  constructor <;> norm_num

/-- Even mode controls the larger of the two parity block bounds. -/
theorem odd_less_even :
    (167281 : ℚ) / 640000 < 2588881 / 3686400 := by norm_num

/-- Optimal within these elementary coefficient-square estimates. -/
theorem improved_projection_constant :
    (1 : ℚ) + 2588881 / 3686400 = 6275281 / 3686400 ∧
    (6275281 : ℚ) / 3686400 < 171 / 100 ∧
    (171 : ℚ) / 100 < 5148881 / 2560000 := by
  constructor
  · norm_num
  constructor <;> norm_num

/-- If the odd and even pivot evaluations obey their weighted Cauchy–Schwarz
    estimates, the full projected quadratic norm obeys a dimension-free bound.
    Here `o` and `e` are the squared norms of the disjoint parity blocks. -/
theorem projection_bound_from_parity_CauchySchwarz
    (o e pivotOdd pivotEven : ℚ)
    (ho : 0 ≤ o) (he : 0 ≤ e)
    (hodd : pivotOdd ^ 2 ≤ ((167281 : ℚ) / 640000) * o)
    (heven : pivotEven ^ 2 ≤ ((2588881 : ℚ) / 3686400) * e) :
    o + e + pivotOdd ^ 2 + pivotEven ^ 2 ≤
      ((6275281 : ℚ) / 3686400) * (o + e) := by
  have habs : ((167281 : ℚ) / 640000) * o ≤
    ((2588881 : ℚ) / 3686400) * o := by
      apply mul_le_mul_of_nonneg_right
      · exact le_of_lt odd_less_even
      · exact ho
  nlinarith

/-- Uniform bound with a convenient small rational constant. -/
theorem projection_bound_convenient
    (o e pivotOdd pivotEven : ℚ)
    (ho : 0 ≤ o) (he : 0 ≤ e)
    (hodd : pivotOdd ^ 2 ≤ ((167281 : ℚ) / 640000) * o)
    (heven : pivotEven ^ 2 ≤ ((2588881 : ℚ) / 3686400) * e) :
    o + e + pivotOdd ^ 2 + pivotEven ^ 2 ≤
      (171 : ℚ) / 100 * (o + e) := by
  have h := projection_bound_from_parity_CauchySchwarz o e pivotOdd pivotEven ho he hodd heven
  have hsum : 0 ≤ o + e := by linarith
  have hk : (6275281 : ℚ) / 3686400 * (o + e) ≤
      (171 : ℚ) / 100 * (o + e) := by
    apply mul_le_mul_of_nonneg_right
    · norm_num
    · exact hsum
  exact le_trans h hk

/-- Better uniform projected tail budget than ES8-P17, conditional on E7
    and the true projection coefficient estimates. -/
theorem projected_arch_budget_improved :
    (171 : ℚ) / 100 * (406 / 91125) = 3857 / 506250 ∧
    (3857 : ℚ) / 506250 < 1 / 100 ∧
    (1 : ℚ) / 100 - 3857 / 506250 = 2411 / 1012500 := by
  constructor
  · norm_num
  constructor <;> norm_num

/-- The analytic premises are *explicit arguments*, not hidden Lean axioms. -/
theorem arch_tail_after_projection
    (raw projected : ℚ)
    (hraw : raw ≤ (406 : ℚ) / 91125)
    (hproject : projected ≤ (171 : ℚ) / 100 * raw) :
    projected < (1 : ℚ) / 100 := by
  have h : projected ≤ (171 : ℚ) / 100 * (406 / 91125) :=
    le_trans hproject (mul_le_mul_of_nonneg_left hraw (by norm_num))
  nlinarith [projected_arch_budget_improved.2.1]

end Numeria.ES8P18
