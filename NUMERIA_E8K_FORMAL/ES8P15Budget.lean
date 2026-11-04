import Mathlib.Tactic

/-!
NUMERIA ES8-P15: strictly algebraic budget for the UNPROJECTED
archimedean tail under the E7 Proposition 5 pointwise bound.
For L=9/20, modes i,j<=n, pi<4 and R(n)=200 n^2.
Does NOT formalize E7's analytic integral or projected matrices.
-/

namespace Numeria.ES8P15Budget

theorem raw_tail_row_budget (n : ℕ) (hn : 1 ≤ n) :
    ((16 : ℚ) / 3645) / (n : ℚ) +
    ((2 : ℚ) / 30375) / ((n : ℚ) ^ 3) < 1 / 100 := by
  have hq : (1 : ℚ) ≤ (n : ℚ) := by exact_mod_cast hn
  have hn0 : (0 : ℚ) < (n : ℚ) := by linarith
  have hpow : (1 : ℚ) ≤ (n : ℚ) ^ 3 := one_le_pow₀ hq
  have hfirst : ((16 : ℚ) / 3645) / (n : ℚ) ≤ (16 : ℚ) / 3645 := by
    apply (div_le_iff₀ hn0).2
    nlinarith [hq]
  have hsecond : ((2 : ℚ) / 30375) / ((n : ℚ) ^ 3) ≤
      (2 : ℚ) / 30375 := by
    apply (div_le_iff₀ (pow_pos hn0 3)).2
    nlinarith [hpow]
  linarith

theorem raw_tail_scaling_identity (n : ℕ) (hn : 1 ≤ n) :
    ((128000 : ℚ) / 729) * (n : ℚ) ^ 3 *
      ((1 : ℚ) / (200 * (n : ℚ) ^ 2) ^ 2 +
        (3 : ℚ) / (200 * (n : ℚ) ^ 2) ^ 3) =
    ((16 : ℚ) / 3645) / (n : ℚ) +
      ((2 : ℚ) / 30375) / ((n : ℚ) ^ 3) := by
  have hq : (1 : ℚ) ≤ (n : ℚ) := by exact_mod_cast hn
  have hn0 : (n : ℚ) ≠ 0 := ne_of_gt (by linarith : (0 : ℚ) < (n : ℚ))
  field_simp [hn0]
  <;> ring

theorem unprojected_e7_tail_budget (n : ℕ) (hn : 1 ≤ n) :
    ((128000 : ℚ) / 729) * (n : ℚ) ^ 3 *
      ((1 : ℚ) / (200 * (n : ℚ) ^ 2) ^ 2 +
        (3 : ℚ) / (200 * (n : ℚ) ^ 2) ^ 3) < 1 / 100 := by
  rw [raw_tail_scaling_identity n hn]
  exact raw_tail_row_budget n hn

end Numeria.ES8P15Budget
