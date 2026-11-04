import Mathlib.Tactic

/-!
NUMERIA–RIEMANN ES8-P17: dimension-free arithmetic budget for projected
ARCHIMEDEAN tail. No axiom, `sorry`, or assertion about RH.

This module establishes the rational implication from the hypotheses:
  raw_tail ≤ 16/(3645*n) + 2/(30375*n³),
  projected_tail ≤ (21/10)*raw_tail, n ≥ 1.
The analytic E7 tail bound and the actual operator-norm bound for the
Fourier pole-cancellation matrix are NOT derived in this Lean file.
Those remain separate analytic obligations, proved algebraically on paper
and tested using exact rational code in ES8-P17.
-/
namespace Numeria.ES8P17

/-- Dimension-dependent raw ARCHIMEDEAN tail majorant from ES8-P15. -/
noncomputable def rawBudget (n : ℚ) : ℚ :=
  16 / (3645 * n) + 2 / (30375 * n ^ 3)

/-- Exact, machine-checkable global constant for the row-sum bound. -/
theorem rational_budget_identity :
    (16 : ℚ) / 3645 + 2 / 30375 = 406 / 91125 := by norm_num

/-- The rational coefficient-squared majorants for the two parity rows. -/
theorem parity_coefficient_budget :
    ((409 : ℚ) / 400) ^ 2 / 2 = 167281 / 320000 ∧
    ((1609 : ℚ) / 800) ^ 2 / 4 = 2588881 / 2560000 := by
  constructor <;> norm_num

/-- The convenient projection amplification bound dominates the sharp
    rational envelope obtained from parity sums. -/
theorem projection_kappa_budget :
    (5148881 : ℚ) / 2560000 < 21 / 10 := by norm_num

/-- The universal projected ARCHIMEDEAN tail budget has positive slack. -/
theorem projected_budget_strict :
    (21 : ℚ) / 10 * (406 / 91125) < 1 / 100 := by norm_num

theorem projected_budget_slack :
    (1 : ℚ) / 100 - (21 / 10) * (406 / 91125) = 391 / 607500 := by norm_num

/-- For n≥1, the rational majorant decreases below its n=1 value.
    This is a pure rational inequality: no statement about Weil integrals. -/
theorem rawBudget_le_at_one (n : ℚ) (hn : 1 ≤ n) :
    rawBudget n ≤ 406 / 91125 := by
  have hpos : 0 < n := by linarith
  have hpowpos : 0 < n ^ 3 := pow_pos hpos 3
  have hsq : 1 ≤ n ^ 2 := by
    nlinarith [mul_nonneg (show 0 ≤ n - 1 by linarith)
      (show 0 ≤ n + 1 by linarith)]
  have hcube : 1 ≤ n ^ 3 := by
    nlinarith [mul_nonneg (show 0 ≤ n - 1 by linarith) (sq_nonneg n)]
  have h1 : (1 : ℚ) / n ≤ 1 :=
    (div_le_iff₀ hpos).2 (by linarith)
  have h3 : (1 : ℚ) / n ^ 3 ≤ 1 :=
    (div_le_iff₀ hpowpos).2 (by linarith)
  calc
    rawBudget n = (16 / 3645) * (1 / n) + (2 / 30375) * (1 / n ^ 3) := by
      unfold rawBudget
      field_simp [ne_of_gt hpos]
      <;> ring
    _ ≤ (16 / 3645 : ℚ) * 1 + (2 / 30375 : ℚ) * 1 := by
      apply add_le_add
      · exact mul_le_mul_of_nonneg_left h1 (by norm_num)
      · exact mul_le_mul_of_nonneg_left h3 (by norm_num)
    _ = 406 / 91125 := by norm_num

/-- Dimension-independent tail budget after any projection satisfying the
    norm-squared envelope 21/10. The analytic and projection estimates are
    hypotheses, NOT consequences of this theorem. -/
theorem projectedTail_within_budget
    (n rawError projectedError : ℚ)
    (hn : 1 ≤ n)
    (hraw : rawError ≤ rawBudget n)
    (hprojection : projectedError ≤ (21 : ℚ) / 10 * rawError) :
    projectedError < 1 / 100 := by
  have hr := le_trans hraw (rawBudget_le_at_one n hn)
  have hp : projectedError ≤ (21 : ℚ) / 10 * (406 / 91125) :=
    le_trans hprojection
      (mul_le_mul_of_nonneg_left hr (by norm_num))
  exact lt_of_le_of_lt hp projected_budget_strict

end Numeria.ES8P17
