import TruncatedFTC
import Mathlib.MeasureTheory.Measure.WithDensity
import Mathlib.Tactic

/-!
# NUMERIA P92-W — The actual positive truncated reciprocal measure

This is the explicit mathematical object, not a hypothetical integrable
measure: μ_X = [2 1_{0<T<X} T^{-3}] • volume + X^{-2} δ_X.
The density is an ℝ≥0∞-valued weight and the atom is a genuine Dirac
measure. It is not falsely asserted that μ_X has finite total mass:
the continuous part diverges at the origin.

Remaining obligations: show (i) the weighted Bochner integral agrees
with the previously checked interval integral, (ii) the relevant
ramps/potentials are integrable for w>0, (iii) identify mixedRamp μ_X
with f_X. Those claims are NOT axioms or sorry-proofs in this file.
-/

namespace NUMERIA.P92W

open MeasureTheory

/-- The positive density 2/T³ supported on 0<T<X. -/
noncomputable def truncatedDensity (X T : ℝ) : ℝ≥0∞ :=
  ENNReal.ofReal (if 0 < T ∧ T < X then 2 / T ^ 3 else 0)

/-- The complete P92-R measure: continuous density plus boundary atom. -/
noncomputable def reciprocalMeasure (X : ℝ) : Measure ℝ :=
  (volume : Measure ℝ).withDensity (truncatedDensity X) +
    (ENNReal.ofReal (1 / X ^ 2)) • Measure.dirac X

theorem density_of_between (X T : ℝ) (h0 : 0 < T) (hTX : T < X) :
    truncatedDensity X T = ENNReal.ofReal (2 / T ^ 3) := by
  simp [truncatedDensity, h0, hTX]

theorem density_zero_outside (X T : ℝ)
    (hout : ¬(0 < T ∧ T < X)) :
    truncatedDensity X T = 0 := by
  simp [truncatedDensity, hout]

theorem density_zero_of_nonpos (X T : ℝ) (hT : T ≤ 0) :
    truncatedDensity X T = 0 := by
  apply density_zero_outside
  intro h
  exact (not_lt.mpr hT) h.1

theorem density_zero_of_cutoff (X T : ℝ) (hT : X ≤ T) :
    truncatedDensity X T = 0 := by
  apply density_zero_outside
  intro h
  exact (not_lt.mpr hT) h.2

theorem density_zero_at_endpoint (X : ℝ) :
    truncatedDensity X X = 0 := by
  exact density_zero_of_cutoff X X le_rfl

end NUMERIA.P92W
