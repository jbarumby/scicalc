import Mathlib.Tactic

/-!
NUMERIA–E8K V16 — prerequisite algebraic lemma for the Proposition 6 tail bound.

This module proves only the rational domain condition used by the tail-bound
formula. It does NOT prove the analytic integral estimate or Proposition 6.
-/

namespace E8K.V16

/-- If 0 ≤ w and 2w ≤ R with R > 0, then 1 - (w/R)^2 is strictly positive.
This isolates the denominator/domain obligation before any integral estimate. -/
theorem tail_denominator_pos {w R : ℚ}
    (hw : 0 ≤ w) (hR : 0 < R) (h2w : 2 * w ≤ R) :
    0 < 1 - (w / R) ^ 2 := by
  have hratio : 0 ≤ w / R := div_nonneg hw (le_of_lt hR)
  have hhalf : w / R ≤ (1 : ℚ) / 2 := by
    apply (div_le_iff₀ hR).2
    nlinarith
  have hproduct : 0 ≤ (w / R) * ((1 : ℚ) / 2 - w / R) :=
    mul_nonneg hratio (by linarith)
  nlinarith

/-- Reusable rational corollary: the denominator is at least 3/4. -/
theorem tail_denominator_lower_bound {w R : ℚ}
    (hw : 0 ≤ w) (hR : 0 < R) (h2w : 2 * w ≤ R) :
    (3 : ℚ) / 4 ≤ 1 - (w / R) ^ 2 := by
  have hratio : 0 ≤ w / R := div_nonneg hw (le_of_lt hR)
  have hhalf : w / R ≤ (1 : ℚ) / 2 := by
    apply (div_le_iff₀ hR).2
    nlinarith
  have hproduct : 0 ≤ (w / R) * ((1 : ℚ) / 2 - w / R) :=
    mul_nonneg hratio (by linarith)
  nlinarith

end E8K.V16
