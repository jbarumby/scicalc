import LocalIntegrability
import Mathlib.MeasureTheory.Measure.WithDensity
import Mathlib.MeasureTheory.Measure.Dirac.Basic
import Mathlib.Tactic

/-!
# P92-W — μ_X is supported on T≤X

Unlike a claim that μ_X has finite total mass, this support property
is true for the infinite-mass density 2/T³ near zero. It uses actual
measure-theoretic almost-everywhere lemmas (withDensity + dirac).
-/

namespace NUMERIA.P92W

open MeasureTheory

theorem measurable_truncatedDensity (X : ℝ) :
    Measurable (truncatedDensity X) := by
  unfold truncatedDensity
  fun_prop

theorem reciprocalMeasure_ae_le_cutoff (X : ℝ) :
    ∀ᵐ T : ℝ ∂reciprocalMeasure X, T ≤ X := by
  unfold reciprocalMeasure
  rw [ae_add_measure_iff]
  constructor
  · apply (ae_withDensity_iff' (measurable_truncatedDensity X).aemeasurable).2
    apply Filter.Eventually.of_forall
    intro T hT
    by_contra hbad
    have hXT : X ≤ T := le_of_not_ge hbad
    exact hT (density_zero_of_cutoff X T hXT)
  · have hDirac : ∀ᵐ T : ℝ ∂Measure.dirac X, T ≤ X := by
      simpa only [ae_dirac_eq, Filter.eventually_pure] using
        (show X ≤ X from le_rfl)
    exact Measure.ae_smul_measure hDirac _

end NUMERIA.P92W
