import P233_GenericMatrixGate
import Mathlib.LinearAlgebra.Matrix.PosDef

/-!
NUMERIA ES8-P234: a finite-dimensional Lean bridge for an exact diagonal-
congruence certificate. The 200x200 rational inequalities are checked by the
accompanying independent Python exact-integer auditor, not by this Lean file.
The true analytic Weil entry inclusion remains an explicit hypothesis.
-/
namespace Numeria.ES8P234

/-- A symmetric matrix with small symmetric entrywise errors around the
identity, with every row error sum bounded by 1/4, is positive definite. -/
theorem near_identity_posDef {n : ℕ} (G : Matrix (Fin n) (Fin n) ℝ)
    (eps : Fin n → Fin n → ℝ)
    (hG : G.IsHermitian)
    (hbox : ∀ i j : Fin n, |G i j - (1 : Matrix (Fin n) (Fin n) ℝ) i j| ≤ eps i j)
    (hsym : ∀ i j, eps i j = eps j i)
    (hrow : ∀ i, (∑ j : Fin n, eps i j) ≤ (1 / 4 : ℝ)) : G.PosDef := by
  apply Numeria.ES8P233.box_posDef G 1 eps (1 / 4 : ℝ) 1 hG hbox hsym hrow
  · intro x
    simp
  · norm_num

/-- Transfer positivity across an invertible real congruence. The unit
hypothesis encodes invertibility of the triangular witness. -/
theorem invertible_congruence_posDef {n : ℕ}
    (W P : Matrix (Fin n) (Fin n) ℝ)
    (eps : Fin n → Fin n → ℝ)
    (hP : IsUnit P)
    (hG : (star P * W * P).IsHermitian)
    (hbox : ∀ i j : Fin n,
      |(star P * W * P) i j - (1 : Matrix (Fin n) (Fin n) ℝ) i j| ≤ eps i j)
    (hsym : ∀ i j, eps i j = eps j i)
    (hrow : ∀ i, (∑ j : Fin n, eps i j) ≤ (1 / 4 : ℝ)) : W.PosDef := by
  have hcert : (star P * W * P).PosDef :=
    near_identity_posDef (star P * W * P) eps hG hbox hsym hrow
  exact (Matrix.IsUnit.posDef_star_left_conjugate_iff hP).mp hcert

/-- The 200-mode case is an INSTANCE of a finite theorem: it requires
true analytic enclosures, interval error and integer P196 verification. -/
theorem P234_200_modes_conditional
    (W P : Matrix (Fin 200) (Fin 200) ℝ)
    (eps : Fin 200 → Fin 200 → ℝ)
    (hP : IsUnit P)
    (hG : (star P * W * P).IsHermitian)
    (hbox : ∀ i j : Fin 200,
      |(star P * W * P) i j - (1 : Matrix (Fin 200) (Fin 200) ℝ) i j| ≤ eps i j)
    (hsym : ∀ i j, eps i j = eps j i)
    (hrow : ∀ i, (∑ j : Fin 200, eps i j) ≤ (1 / 4 : ℝ)) : W.PosDef := by
  exact invertible_congruence_posDef W P eps hP hG hbox hsym hrow

#print axioms P234_200_modes_conditional
end Numeria.ES8P234
