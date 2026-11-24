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
noncomputable def rampPoly0 (L h : ℝ) (n : ℕ) : ℝ :=
  h * (n : ℝ) * L - h ^ 2 * ((n : ℝ) * ((n : ℝ) - 1)) / 2

/-- The first term h*L is removed to obtain the shifted tail U₁. -/
noncomputable def rampPoly1 (L h : ℝ) (n : ℕ) : ℝ :=
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


/-- The actual finite positive-part ramp sum. -/
noncomputable def finiteRamp0 (L h : ℝ) (n : ℕ) : ℝ :=
  h * ∑ j ∈ Finset.range n, max (L - (j : ℝ) * h) 0

/-- Exact identification of a positive-part ramp sum with an affine sum
    when all n included summands are nonnegative. -/
theorem finiteRamp0_eq_affine (L h : ℝ) :
    ∀ n : ℕ, (∀ j : ℕ, j < n → (j : ℝ) * h ≤ L) →
      finiteRamp0 L h n = affineRampSum L h n := by
  intro n
  induction n with
  | zero =>
      intro _
      simp [finiteRamp0, affineRampSum]
  | succ n ih =>
      intro hsupp
      have hsuppn : ∀ j : ℕ, j < n → (j : ℝ) * h ≤ L := by
        intro j hj
        exact hsupp j (Nat.lt_succ_of_lt hj)
      have hlast : (n : ℝ) * h ≤ L :=
        hsupp n (Nat.lt_succ_self _)
      have hmax :
          max (L - (n : ℝ) * h) 0 = L - (n : ℝ) * h :=
        max_eq_left (sub_nonneg.mpr hlast)
      calc
        finiteRamp0 L h (n + 1) =
            finiteRamp0 L h n + h * (L - (n : ℝ) * h) := by
              simp only [finiteRamp0, Finset.sum_range_succ, hmax]
              ring
        _ = affineRampSum L h (n + 1) := by
              rw [ih hsuppn]
              rfl

/-- The positive-part ramp sum really equals P₀ on its support cell. -/
theorem finiteRamp0_closed (L h : ℝ) (n : ℕ)
    (hsupp : ∀ j : ℕ, j < n → (j : ℝ) * h ≤ L) :
    finiteRamp0 L h n = rampPoly0 L h n :=
  (finiteRamp0_eq_affine L h n hsupp).trans (affineRampSum_closed L h n)

/-- Support condition follows from the last positive ramp being nonnegative. -/
theorem ramp_support_of_cell (L h : ℝ) (n : ℕ)
    (hh : 0 ≤ h)
    (hlower : ((n : ℝ) - 1) * h ≤ L) :
    ∀ j : ℕ, j < n → (j : ℝ) * h ≤ L := by
  intro j hj
  have hjN : j + 1 ≤ n := Nat.succ_le_iff.mpr hj
  have hjR : (j : ℝ) + 1 ≤ (n : ℝ) := by exact_mod_cast hjN
  have hbound : (j : ℝ) ≤ (n : ℝ) - 1 := by linarith
  exact (mul_le_mul_of_nonneg_right hbound hh).trans hlower

/-- Real Bellman inequality for *actual finite sums of positive-part ramps*,
    not merely auxiliary polynomials, on every fixed-support cell. -/
theorem finiteRamp_bellman_same_cell
    (L h b : ℝ) (n : ℕ)
    (hn : 1 ≤ n) (hh : 0 ≤ h) (hb : 0 ≤ b)
    (hlower : ((n : ℝ) - 1) * max h b ≤ L)
    (hupper : L ≤ (n : ℝ) * min h b) :
    min h b * L + (finiteRamp0 L b n - b * L) ≤
      finiteRamp0 L h n := by
  have hn1 : (0 : ℝ) ≤ (n : ℝ) - 1 := by
    have hnR : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
    linarith
  have hlh : ((n : ℝ) - 1) * h ≤ L :=
    (mul_le_mul_of_nonneg_left (le_max_left h b) hn1).trans hlower
  have hlb : ((n : ℝ) - 1) * b ≤ L :=
    (mul_le_mul_of_nonneg_left (le_max_right h b) hn1).trans hlower
  have hsuph := ramp_support_of_cell L h n hh hlh
  have hsupb := ramp_support_of_cell L b n hb hlb
  rw [finiteRamp0_closed L b n hsupb,
      finiteRamp0_closed L h n hsuph]
  change min h b * L + rampPoly1 L b n ≤ rampPoly0 L h n
  exact bellman_ramp_same_cell L h b n hn hlower hupper

/-- All terms beyond the support count vanish, with no hidden convergence
    assumption. -/
theorem ramp_terms_zero_after_support (L h : ℝ) (n j : ℕ)
    (hh : 0 ≤ h) (hupper : L ≤ (n : ℝ) * h) (hnj : n ≤ j) :
    max (L - (j : ℝ) * h) 0 = 0 := by
  have hnR : (n : ℝ) ≤ (j : ℝ) := by exact_mod_cast hnj
  have hmul : (n : ℝ) * h ≤ (j : ℝ) * h :=
    mul_le_mul_of_nonneg_right hnR hh
  have hnonpos : L - (j : ℝ) * h ≤ 0 := by linarith
  exact max_eq_right hnonpos


/-- The support switch at L=n*h is continuous: the newly included
    affine summand has value h*(L-n*h)=0. -/
theorem rampPoly0_boundary_agree (L h : ℝ) (n : ℕ)
    (hL : L = (n : ℝ) * h) :
    rampPoly0 L h n = rampPoly0 L h (n + 1) := by
  calc
    rampPoly0 L h n = affineRampSum L h n :=
      (affineRampSum_closed L h n).symm
    _ = affineRampSum L h (n + 1) := by
      change affineRampSum L h n =
        affineRampSum L h n + h * (L - (n : ℝ) * h)
      rw [hL]
      ring
    _ = rampPoly0 L h (n + 1) := affineRampSum_closed L h (n + 1)

/-- The shifted U₁ polynomial has the same boundary continuity. -/
theorem rampPoly1_boundary_agree (L h : ℝ) (n : ℕ)
    (hL : L = (n : ℝ) * h) :
    rampPoly1 L h n = rampPoly1 L h (n + 1) := by
  unfold rampPoly1
  rw [rampPoly0_boundary_agree L h n hL]

/-- The actual positive-part ramp sums also agree at support switches;
    this equality is independent of the analytic monotonicity issue. -/
theorem finiteRamp0_boundary_agree (L h : ℝ) (n : ℕ)
    (hL : L = (n : ℝ) * h) :
    finiteRamp0 L h n = finiteRamp0 L h (n + 1) := by
  unfold finiteRamp0
  rw [Finset.sum_range_succ]
  rw [hL]
  simp

end NUMERIA.P92T
