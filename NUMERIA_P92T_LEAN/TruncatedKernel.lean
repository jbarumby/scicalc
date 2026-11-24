import IntegralRamp
import Mathlib.Tactic

/-!
P92-U: boundary atom in the truncated reciprocal kernel

For 0 < u ≤ X, the representation
(1/u - 1/X) = ∫[u,X] 2(T-u)/T³ dT + (X-u)/X²
has a boundary atom at T=X with weight 1/X².

This module proves the *algebraic* antiderivative accounting and the
Dirac-mass contribution. It does NOT yet formalize the definite
integral of the continuous density; that is a separate calculus gate.
-/

namespace NUMERIA.P92U

open MeasureTheory

/-- The integrand 2(T-u)/T³ has the elementary antiderivative
    F(T)=-2/T+u/T². This identity is the exact endpoint
    evaluation plus atom, with no integral theorem asserted. -/
theorem reciprocal_trunc_primitive_plus_atom
    (u X : ℝ) (hu : 0 < u) (hX : 0 < X) :
    ((-2 / X + u / X ^ 2) - (-2 / u + u / u ^ 2)) +
      (X - u) / X ^ 2 = 1 / u - 1 / X := by
  field_simp
  ring

/-- Omitting the atom changes the primitive difference by exactly
    (X-u)/X²; hence it cannot be silently dropped. -/
theorem reciprocal_trunc_missing_atom
    (u X : ℝ) (hu : 0 < u) (hX : 0 < X) :
    (1 / u - 1 / X) -
      ((-2 / X + u / X ^ 2) - (-2 / u + u / u ^ 2))
      = (X - u) / X ^ 2 := by
  field_simp
  ring

/-- The atom's contribution is nonnegative in the intended domain. -/
theorem reciprocal_trunc_atom_nonneg
    (u X : ℝ) (hUX : u ≤ X) (hX : 0 < X) :
    0 ≤ (X - u) / X ^ 2 := by
  exact div_nonneg (sub_nonneg.mpr hUX) (sq_nonneg X)

/-- An actual Dirac measure at X integrates the ramp to its
    pointwise value at X. This verifies the atomic interpretation,
    rather than relying on an informal delta symbol. -/
theorem dirac_ramp_at_endpoint (u X : ℝ) :
    (∫ T : ℝ, max (T - u) 0 ∂Measure.dirac X) =
      max (X - u) 0 := by
  exact integral_dirac (fun T : ℝ => max (T - u) 0) X

/-- With the positive mass X⁻², the Dirac term has precisely the
    correction required by the truncated reciprocal representation. -/
theorem dirac_ramp_boundary_atom
    (u X : ℝ) (hUX : u ≤ X) (hX : 0 < X) :
    (1 / X ^ 2) *
       (∫ T : ℝ, max (T - u) 0 ∂Measure.dirac X) =
      (X - u) / X ^ 2 := by
  rw [dirac_ramp_at_endpoint]
  rw [max_eq_left (sub_nonneg.mpr hUX)]
  ring

end NUMERIA.P92U
