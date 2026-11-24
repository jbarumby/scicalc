import TruncatedKernel
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Tactic

/-!
# NUMERIA P92-V — Fundamental theorem of calculus for the f_X density

Target: for 0<u≤X, the interval integral of 2(T-u)/T³ equals
the difference of the explicit primitive F_u(T)=-2/T+u/T².

Combining this with the independently certified P92-U boundary atom
should establish the full elementary identity
  ∫_u^X 2(T-u)/T³ dT + (X-u)/X² = 1/u - 1/X.

This is a theorem about an INTERVAL integral, not yet the identification
of the density measure μ_X or integrability of Bellman potentials.
-/

namespace NUMERIA.P92V

open MeasureTheory intervalIntegral

private noncomputable def reciprocalPrimitive (u t : ℝ) : ℝ :=
  -2 * t⁻¹ + u * (t⁻¹) ^ 2

private theorem reciprocalPrimitive_hasDerivAt
    (u t : ℝ) (ht : 0 < t) :
    HasDerivAt (reciprocalPrimitive u) (2 * (t - u) / t ^ 3) t := by
  have hne : t ≠ 0 := ne_of_gt ht
  have hinv : HasDerivAt (fun x : ℝ => x⁻¹) (-(t ^ 2)⁻¹) t :=
    hasDerivAt_inv hne
  have hsq : HasDerivAt (fun x : ℝ => (x⁻¹) ^ 2)
      (2 * t⁻¹ * (-(t ^ 2)⁻¹)) t := by
    convert hinv.pow 2 using 1 <;> ring
  have hsum := (hinv.const_mul (-2)).add (hsq.const_mul u)
  unfold reciprocalPrimitive
  convert hsum using 1
  · ring
  · field_simp [hne]
    ring

private theorem reciprocalDensity_continuousOn
    (u X : ℝ) (hu : 0 < u) (hUX : u ≤ X) :
    ContinuousOn (fun t : ℝ => 2 * (t - u) / t ^ 3) (Set.uIcc u X) := by
  intro t ht
  have htIcc : t ∈ Set.Icc u X := by
    simpa only [Set.uIcc_of_le hUX] using ht
  have htpos : 0 < t := lt_of_lt_of_le hu htIcc.1
  have hne : t ^ 3 ≠ 0 := pow_ne_zero _ (ne_of_gt htpos)
  have hn : ContinuousAt (fun t : ℝ => 2 * (t - u)) t := by fun_prop
  have hd : ContinuousAt (fun t : ℝ => t ^ 3) t := by fun_prop
  exact (hn.div hd hne).continuousWithinAt

theorem truncated_density_FTC
    (u X : ℝ) (hu : 0 < u) (hUX : u ≤ X) :
    (∫ T : ℝ in u..X, 2 * (T - u) / T ^ 3) =
      (-2 / X + u / X ^ 2) - (-2 / u + u / u ^ 2) := by
  have hder : ∀ T ∈ Set.uIcc u X,
      HasDerivAt (reciprocalPrimitive u) (2 * (T - u) / T ^ 3) T := by
    intro T hT
    have hTcc : T ∈ Set.Icc u X := by
      simpa only [Set.uIcc_of_le hUX] using hT
    exact reciprocalPrimitive_hasDerivAt u T (lt_of_lt_of_le hu hTcc.1)
  have hint : IntervalIntegrable
      (fun T : ℝ => 2 * (T - u) / T ^ 3) volume u X :=
    (reciprocalDensity_continuousOn u X hu hUX).intervalIntegrable
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt hder hint
  unfold reciprocalPrimitive at hFTC
  convert hFTC using 1
  · field_simp [ne_of_gt hu, ne_of_gt (lt_of_lt_of_le hu hUX)] <;> ring

/-- The complete truncated reciprocal formula: continuous density
plus the nonnegative endpoint atomic mass X⁻² at T=X. -/
theorem truncated_reciprocal_integral_with_atom
    (u X : ℝ) (hu : 0 < u) (hUX : u ≤ X) :
    (∫ T : ℝ in u..X, 2 * (T - u) / T ^ 3) +
      (1 / X ^ 2) *
        (∫ T : ℝ, max (T - u) 0 ∂Measure.dirac X) =
        1 / u - 1 / X := by
  have hX : 0 < X := lt_of_lt_of_le hu hUX
  rw [truncated_density_FTC u X hu hUX,
      NUMERIA.P92U.dirac_ramp_boundary_atom u X hUX hX]
  exact NUMERIA.P92U.reciprocal_trunc_primitive_plus_atom u X hu hX

end NUMERIA.P92V
