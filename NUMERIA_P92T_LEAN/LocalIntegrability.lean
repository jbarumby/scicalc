import PotentialBound
import Mathlib.MeasureTheory.Integral.IntegrableOn
import Mathlib.Tactic

/-!
# NUMERIA P92-W — Local integrability criterion for Bellman potentials

This theorem correctly avoids the false global finite-mass condition.
For w>0, the potential vanishes at T≤w, and is uniformly bounded
when the measure is supported in T≤X. Thus it suffices that:
  (1) μ((w,+∞)) < ∞,
  (2) the potential is strongly measurable,
  (3) μ-almost all T>w also satisfy T≤X.

This is a genuine analytic reduction, not a proof of those three
conditions for the specific μ_X. That specialization requires
work with withDensity and the finite-support measurability lemma.
-/

namespace NUMERIA.P92W

open MeasureTheory NUMERIA.P92T

/-- A precise local-mass-to-integrability implication. -/
theorem rampPotential_integrable_of_local_tail
    (μ : Measure ℝ) (X w h : ℝ) (hh : 0 < h)
    (hfinite : μ (Set.Ioi w) < ⊤)
    (hmeas : AEStronglyMeasurable (fun T : ℝ =>
      rampPotential T w h) μ)
    (hcut : ∀ᵐ T : ℝ ∂(μ.restrict (Set.Ioi w)), T ≤ X) :
    Integrable (fun T : ℝ => rampPotential T w h) μ := by
  have hlocal : IntegrableOn
      (fun T : ℝ => rampPotential T w h) (Set.Ioi w) μ := by
    apply Measure.integrableOn_of_bounded hfinite.ne hmeas
    filter_upwards [hcut] with T hTX
    obtain ⟨hpos, hbound⟩ :=
      rampPotential_bound_on_support T X w h hTX hh
    rw [Real.norm_eq_abs, abs_of_nonneg hpos]
    exact hbound
  have hsupport : Function.support (fun T : ℝ => rampPotential T w h) ⊆
      Set.Ioi w := by
    intro T hT
    by_contra hn
    have hTle : T ≤ w := le_of_not_gt hn
    exact hT (rampPotential_zero_below T w h hTle)
  exact (integrableOn_iff_integrable_of_support_subset hsupport).mp hlocal

end NUMERIA.P92W
