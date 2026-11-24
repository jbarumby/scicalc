import TruncatedKernel
import Mathlib.Tactic

/-!
# NUMERIA P92-U — Finite interchange and uniform ramp support

Once μ is a positive measure on ℝ and the individual ramps are
integrable, Bochner linearity gives interchange of a FINITE sum and the
integral. This is *not* an infinite Tonelli theorem.

If all support parameters T satisfy T≤X and h>0, every term with
j≥ceil((X-w)/h) vanishes. This provides a uniform finite support
bound when μ is restricted to T≤X, as occurs for f_X.

No f_X density integral or global Preissmann 5/6 proof is assumed.
-/

namespace NUMERIA.P92U

open MeasureTheory

/-- Exact interchange of a finite sum of positive-part ramp functions
with a positive-measure Bochner integral, under explicit integrability
of each summand. -/
theorem integral_finite_ramp_sum
    (μ : Measure ℝ) (S : Finset ℕ) (w h : ℝ)
    (hint : ∀ j ∈ S,
      Integrable (fun T : ℝ =>
        max (T - (w + (j : ℝ) * h)) 0) μ) :
    (∫ T : ℝ, h *
        (∑ j ∈ S, max (T - (w + (j : ℝ) * h)) 0) ∂μ) =
      h * (∑ j ∈ S,
        (∫ T : ℝ, max (T - (w + (j : ℝ) * h)) 0 ∂μ)) := by
  calc
    (∫ T : ℝ, h *
        (∑ j ∈ S, max (T - (w + (j : ℝ) * h)) 0) ∂μ) =
      h * (∫ T : ℝ,
        (∑ j ∈ S, max (T - (w + (j : ℝ) * h)) 0) ∂μ) := by
      rw [integral_const_mul]
    _ = h * (∑ j ∈ S,
        (∫ T : ℝ, max (T - (w + (j : ℝ) * h)) 0 ∂μ)) := by
      rw [integral_finsetSum S hint]

/-- Pointwise tail vanishing whenever T has a uniform upper support
bound X and the mesh index is past a compatible cutoff N. -/
theorem ramp_tail_zero_after
    (T X w h : ℝ) (N j : ℕ)
    (hh : 0 < h) (hTX : T ≤ X)
    (hcut : X - w ≤ (N : ℝ) * h) (hNj : N ≤ j) :
    max (T - (w + (j : ℝ) * h)) 0 = 0 := by
  have hcast : (N : ℝ) ≤ (j : ℝ) := by exact_mod_cast hNj
  have hmul : (N : ℝ) * h ≤ (j : ℝ) * h :=
    mul_le_mul_of_nonneg_right hcast (le_of_lt hh)
  have hnonpos : T - (w + (j : ℝ) * h) ≤ 0 := by
    linarith
  exact max_eq_right hnonpos

/-- Nat.ceil supplies a uniform finite cutoff for all T≤X.
No extra assumption X≥w is needed; if X<w the full ramp vanishes. -/
theorem ramp_tail_zero_after_ceil
    (T X w h : ℝ) (j : ℕ)
    (hh : 0 < h) (hTX : T ≤ X)
    (hNj : Nat.ceil ((X - w) / h) ≤ j) :
    max (T - (w + (j : ℝ) * h)) 0 = 0 := by
  have hc : (X - w) / h ≤
      ((Nat.ceil ((X - w) / h) : ℕ) : ℝ) :=
    Nat.le_ceil _
  have hm := mul_le_mul_of_nonneg_right hc (le_of_lt hh)
  have hcancel : ((X - w) / h) * h = X - w := by
    field_simp [ne_of_gt hh]
  have hcut : X - w ≤
      ((Nat.ceil ((X - w) / h) : ℕ) : ℝ) * h := by
    linarith
  exact ramp_tail_zero_after T X w h
    (Nat.ceil ((X - w) / h)) j hh hTX hcut hNj

end NUMERIA.P92U
