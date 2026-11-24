import ShiftRamp
import Mathlib.MeasureTheory.Integral.Bochner.Basic
import Mathlib.Tactic

/-!
# NUMERIA P92-U — Integral lift of global convex-ramp Bellman

R1 for individual ramps and finite positive combinations was checked in
P92-T. Here we prove the measure-theoretic transfer for an arbitrary
positive measure μ on ℝ, using Bochner integral monotonicity,
linearity and explicit integrability hypotheses.

No R1 is assumed as an axiom: the pointwise Bellman inequality is
imported as the P92-T theorem bellman_for_all_ramps.

Scope: this is the R1 Bellman inequality for any integrable mixture
of positive-part ramps. Establishing integrability and exact integral
representations for f_X, u⁻² and φ_h is a separate task. Preissmann
Lemmas 5 and 6 and the Riemann hypothesis are not proved here.
-/

namespace NUMERIA.P92U

open MeasureTheory
open NUMERIA.P92T

/-- Positive-measure integration preserves sharp Bellman inequalities
provided the three integrands are integrable. The measure is arbitrary. -/
theorem bellman_integral_ramps
    (μ : Measure ℝ) (w h b : ℝ)
    (hh : 0 < h) (hb : 0 < b)
    (hw : Integrable (fun T : ℝ => max (T - w) 0) μ)
    (hshift : Integrable (fun T : ℝ =>
      rampPotential T (w + b) b) μ)
    (hbase : Integrable (fun T : ℝ =>
      rampPotential T w h) μ) :
    min h b * (∫ T : ℝ, max (T - w) 0 ∂μ) +
      (∫ T : ℝ, rampPotential T (w + b) b ∂μ) ≤
    ∫ T : ℝ, rampPotential T w h ∂μ := by
  have hscaled :
      Integrable (fun T : ℝ => min h b * max (T - w) 0) μ :=
    hw.const_mul (min h b)
  have hsum :
      Integrable (fun T : ℝ =>
        min h b * max (T - w) 0 + rampPotential T (w + b) b) μ :=
    hscaled.add hshift
  have hpointwise : ∀ T : ℝ,
      min h b * max (T - w) 0 + rampPotential T (w + b) b ≤
        rampPotential T w h := by
    intro T
    exact bellman_for_all_ramps T w h b hh hb
  have hbound :
      (∫ T : ℝ, min h b * max (T - w) 0 +
           rampPotential T (w + b) b ∂μ) ≤
      (∫ T : ℝ, rampPotential T w h ∂μ) :=
    integral_mono hsum hbase hpointwise
  rw [integral_add hscaled hshift, integral_const_mul] at hbound
  exact hbound

/-- A reusable Bellman statement with mixed function and potential. -/
noncomputable def mixedRamp (μ : Measure ℝ) (u : ℝ) : ℝ :=
  ∫ T : ℝ, max (T - u) 0 ∂μ

noncomputable def mixedPotential (μ : Measure ℝ) (u h : ℝ) : ℝ :=
  ∫ T : ℝ, rampPotential T u h ∂μ

theorem mixedRamp_bellman
    (μ : Measure ℝ) (w h b : ℝ)
    (hh : 0 < h) (hb : 0 < b)
    (hw : Integrable (fun T : ℝ => max (T - w) 0) μ)
    (hshift : Integrable (fun T : ℝ =>
      rampPotential T (w + b) b) μ)
    (hbase : Integrable (fun T : ℝ =>
      rampPotential T w h) μ) :
    min h b * mixedRamp μ w +
      mixedPotential μ (w + b) b ≤ mixedPotential μ w h := by
  exact bellman_integral_ramps μ w h b hh hb hw hshift hbase

end NUMERIA.P92U
