import Mathlib.Tactic
import Mathlib.Basic.Real.Basic

/-!
# NUMERIA P92-T — Bellman ramp inequalities on a fixed support cell

This file builds on P92-S and removes the abstract Bellman-step premise
for piecewise-affine ramps *when the same integer n counts the positive
terms at both mesh sizes*. This is a genuine, universally quantified
theorem about arbitrary real L,h,b under cell constraints.

The global R1 theorem still requires gluing across integer support
transitions, then transferring to positive mixtures. NO proof of full
Preissmann Lemmas 5 and 6 or RH is asserted here.
-/

namespace NUMERIA.P92T

/-- Finite arithmetic ramp sum with exactly n active summands:
    h * Σ_{j=0}^{n-1} (L - j*h). -/
def affineRampSum (L h : ℝ) : ℕ → ℝ
  | 0 => 0
  | n + 1 => affineRampSum L h n + h * (L - (n : ℝ) * h)

/-- Its rational-polynomial closed form, valid algebraically for every n. -/
def rampPoly0 (L h : ℝ) (n : ℕ) : ℝ :=
  h * (n : ℝ) * L - h ^ 2 * ((n : ℝ) * ((n : ℝ) - 1)) / 2

/-- The first term h*L is removed to obtain the shifted tail U₁. -/
def rampPoly1 (L h : ℝ) (n : ℕ) : ℝ :=
  rampPoly0 L h n - h * L

theorem affineRampSum_closed (L h : ℝ) (n : ℕ) :
    affineRampSum L h n = rampPoly0 L h n := by
  induction n with
  | zero =>
      simp [affineRampSum, rampPoly0]
  | succ n ih =>
      change affineRampSum L h n + h * (L - (n : ℝ) * h) =
        rampPoly0 L h (n + 1)
      rw [ih]
      simp only [rampPoly0, Nat.cast_add, Nat.cast_one]
      ring

/-- Monotonicity of U₀ on any interval where exactly n≥1 ramp terms
    can be active at both ends of the interval. -/
theorem rampPoly0_mono_same_cell (L h₁ h₂ : ℝ) (n : ℕ)
    (hn : 1 ≤ n) (hle : h₁ ≤ h₂)
    (hlower : ((n : ℝ) - 1) * h₂ ≤ L) :
    rampPoly0 L h₁ n ≤ rampPoly0 L h₂ n := by
  have hn0 : (0 : ℝ) ≤ n := Nat.cast_nonneg _
  have hn1 : (0 : ℝ) ≤ (n : ℝ) - 1 := by
    have : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
    linarith
  have htwice : h₁ + h₂ ≤ 2 * h₂ := by linarith
  have hmul :
      ((n : ℝ) - 1) * (h₁ + h₂) ≤
        ((n : ℝ) - 1) * (2 * h₂) :=
    mul_le_mul_of_nonneg_left htwice hn1
  have hfactor : 0 ≤ L - ((n : ℝ) - 1) * (h₁ + h₂) / 2 := by
    linarith
  have hdiff :
      rampPoly0 L h₂ n - rampPoly0 L h₁ n =
        (h₂ - h₁) * (n : ℝ) *
          (L - ((n : ℝ) - 1) * (h₁ + h₂) / 2) := by
    unfold rampPoly0
    ring
  have hpos :
      0 ≤ (h₂ - h₁) * (n : ℝ) *
        (L - ((n : ℝ) - 1) * (h₁ + h₂) / 2) :=
    mul_nonneg (mul_nonneg (sub_nonneg.mpr hle) hn0) hfactor
  linarith [hdiff]

/-- Antitonicity of U₁ on a common positive-support cell. -/
theorem rampPoly1_antitone_same_cell (L h₁ h₂ : ℝ) (n : ℕ)
    (hn : 1 ≤ n) (hle : h₁ ≤ h₂)
    (hupper : L ≤ (n : ℝ) * h₁) :
    rampPoly1 L h₂ n ≤ rampPoly1 L h₁ n := by
  have hn0 : (0 : ℝ) ≤ n := Nat.cast_nonneg _
  have hn1 : (0 : ℝ) ≤ (n : ℝ) - 1 := by
    have : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
    linarith
  have htwice : 2 * h₁ ≤ h₁ + h₂ := by linarith
  have hmul :
      (n : ℝ) * (2 * h₁) ≤
        (n : ℝ) * (h₁ + h₂) :=
    mul_le_mul_of_nonneg_left htwice hn0
  have hfactor : L - (n : ℝ) * (h₁ + h₂) / 2 ≤ 0 := by
    linarith
  have hdiff :
      rampPoly1 L h₂ n - rampPoly1 L h₁ n =
        (h₂ - h₁) * ((n : ℝ) - 1) *
          (L - (n : ℝ) * (h₁ + h₂) / 2) := by
    unfold rampPoly1 rampPoly0
    ring
  have hneg :
      (h₂ - h₁) * ((n : ℝ) - 1) *
        (L - (n : ℝ) * (h₁ + h₂) / 2) ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos
      (mul_nonneg (sub_nonneg.mpr hle) hn1) hfactor
  linarith [hdiff]

/-- Bellman inequality for b≤h on a single ramp-support cell:
    the next-step coefficient is b, and monotonicity of U₀ applies. -/
theorem bellman_ramp_small_step_same_cell
    (L h b : ℝ) (n : ℕ)
    (hn : 1 ≤ n) (hb : b ≤ h)
    (hlower : ((n : ℝ) - 1) * h ≤ L) :
    min h b * L + rampPoly1 L b n ≤ rampPoly0 L h n := by
  have hm := rampPoly0_mono_same_cell L b h n hn hb hlower
  rw [min_eq_right hb]
  unfold rampPoly1
  linarith

/-- Bellman inequality for h≤b on a single ramp-support cell:
    the current coefficient is h, and antitonicity of U₁ applies. -/
theorem bellman_ramp_large_step_same_cell
    (L h b : ℝ) (n : ℕ)
    (hn : 1 ≤ n) (hb : h ≤ b)
    (hupper : L ≤ (n : ℝ) * h) :
    min h b * L + rampPoly1 L b n ≤ rampPoly0 L h n := by
  have hm := rampPoly1_antitone_same_cell L h b n hn hb hupper
  rw [min_eq_left hb]
  unfold rampPoly1 at *
  linarith

/-- All spacing regimes are now covered *provided the active-term count
    n remains the same*. Crossing the cell interfaces is a separate gate. -/
theorem bellman_ramp_same_cell (L h b : ℝ) (n : ℕ)
    (hn : 1 ≤ n)
    (hlower : ((n : ℝ) - 1) * max h b ≤ L)
    (hupper : L ≤ (n : ℝ) * min h b) :
    min h b * L + rampPoly1 L b n ≤ rampPoly0 L h n := by
  rcases le_total b h with hb | hb
  · exact bellman_ramp_small_step_same_cell L h b n hn hb
      (le_trans (mul_le_mul_of_nonneg_left
        (le_max_left h b) (by
          have : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
          linarith)) hlower)
  · exact bellman_ramp_large_step_same_cell L h b n hn hb
      (le_trans hupper (mul_le_mul_of_nonneg_left
        (min_le_left h b) (Nat.cast_nonneg n)))

end NUMERIA.P92T
