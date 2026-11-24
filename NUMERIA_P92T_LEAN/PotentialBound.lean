import MeasureKernel
import Mathlib.Tactic

/-!
# NUMERIA P92-W — A genuine analytic upper bound for ramp potentials

For L≥0 and h>0, the potential with Nat.ceil support satisfies
  0 ≤ U₀(L,h) ≤ L(L+h).
For T≤X this gives the uniform bound
  U_T(w,h) ≤ (X-w)_+((X-w)_+ + h).
These are exact Lean theorems, with no integrability assumptions.

One cannot integrate the uniform bound against the *full* μ_X:
μ_X has infinite mass near zero. Only the restriction T>w is usable.
Finite local mass/integrability is the NEXT formal gate.
-/

namespace NUMERIA.P92W

open NUMERIA.P92T

theorem fullRamp0_nonneg (L h : ℝ) (hh : 0 ≤ h) :
    0 ≤ fullRamp0 L h := by
  unfold fullRamp0 finiteRamp0
  apply mul_nonneg hh
  apply Finset.sum_nonneg
  intro j hj
  exact le_max_right _ _

/-- An exact polynomial bound valid for all nonnegative ramp heights. -/
theorem fullRamp0_le_quadratic (L h : ℝ)
    (hL : 0 ≤ L) (hh : 0 < h) :
    fullRamp0 L h ≤ L * (L + h) := by
  rcases eq_or_lt_of_le hL with hzero | hpos
  · subst L
    simp [fullRamp0, finiteRamp0]
  · let n : ℕ := Nat.ceil (L / h)
    have hcell := ceil_ramp_support L h hpos hh
    have hsupport :
        ∀ j : ℕ, j < n → (j : ℝ) * h ≤ L :=
      ramp_support_of_cell L h n (le_of_lt hh) hcell.2.1
    have heq : fullRamp0 L h = rampPoly0 L h n := by
      change finiteRamp0 L h n = rampPoly0 L h n
      exact finiteRamp0_closed L h n hsupport
    have hn1 : (1 : ℝ) ≤ (n : ℝ) := by
      exact_mod_cast hcell.1
    have hn0 : 0 ≤ (n : ℝ) - 1 := by linarith
    have hquad : 0 ≤ h ^ 2 * ((n : ℝ) * ((n : ℝ) - 1)) / 2 := by
      apply div_nonneg
      · exact mul_nonneg (sq_nonneg _) (mul_nonneg (by positivity) hn0)
      · norm_num
    have hpoly : rampPoly0 L h n ≤ (h * (n : ℝ)) * L := by
      unfold rampPoly0
      nlinarith
    have hstep : h * (n : ℝ) ≤ L + h := by
      dsimp [n] at hcell
      nlinarith [hcell.2.1]
    have hmul : (h * (n : ℝ)) * L ≤ (L + h) * L :=
      mul_le_mul_of_nonneg_right hstep hL
    rw [heq]
    nlinarith [hpoly, hmul]

theorem rampPotential_bound (T w h : ℝ) (hh : 0 < h) :
    0 ≤ rampPotential T w h ∧
      rampPotential T w h ≤
        max (T - w) 0 * (max (T - w) 0 + h) := by
  unfold rampPotential
  exact ⟨fullRamp0_nonneg _ _ (le_of_lt hh),
    fullRamp0_le_quadratic _ _ (le_max_right _ _) hh⟩

/-- A uniform integrand bound on the actual support of μ_X. -/
theorem rampPotential_bound_on_support
    (T X w h : ℝ) (hTX : T ≤ X) (hh : 0 < h) :
    0 ≤ rampPotential T w h ∧
      rampPotential T w h ≤
        max (X - w) 0 * (max (X - w) 0 + h) := by
  obtain ⟨hpos, hbase⟩ := rampPotential_bound T w h hh
  let L : ℝ := max (T - w) 0
  let M : ℝ := max (X - w) 0
  have hLM : L ≤ M := by
    dsimp [L, M]
    exact max_le_max (sub_le_sub_right hTX w) le_rfl
  have hL : 0 ≤ L := le_max_right _ _
  have hM : 0 ≤ M := le_max_right _ _
  have hprod : L * (L + h) ≤ M * (M + h) := by
    have hp : 0 ≤ (M - L) * (M + L + h) :=
      mul_nonneg (sub_nonneg.mpr hLM) (by linarith)
    nlinarith
  exact ⟨hpos, hbase.trans hprod⟩

/-- Every potential vanishes below w, including the singular region T≈0. -/
theorem rampPotential_zero_below (T w h : ℝ) (hT : T ≤ w) :
    rampPotential T w h = 0 := by
  have hmax : max (T - w) 0 = 0 :=
    max_eq_right (sub_nonpos.mpr hT)
  simp [rampPotential, fullRamp0, finiteRamp0, hmax]

end NUMERIA.P92W
