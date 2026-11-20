import Development.AutomorphicArithmetic.CharacterBoxes.EveryCharacter
import Formalism.GammaFactorNonvanishing.RhOfRight

open Complex DirichletCharacter

theorem Formalism.VW.gammaFactor_ne_zero {N : ℕ} (χ : DirichletCharacter ℂ N) {s : ℂ} (hs : 0 < s.re) : gammaFactor χ s ≠ 0 := by
  rcases χ.even_or_odd with h | h
  · rw [h.gammaFactor_def]; exact Gammaℝ_ne_zero_of_re_pos hs
  · rw [h.gammaFactor_def]; exact Gammaℝ_ne_zero_of_re_pos (by rw [add_re, one_re]; linarith)

theorem Formalism.VW.critLine (h : UpperTheory.EveryCharacter) (N : ℕ) [NeZero N] (χ : DirichletCharacter ℂ N) (hχ : χ.IsPrimitive) (s : ℂ)
    (h0 : 0 < s.re) (h1 : s.re < 1) (hz : LFunction χ s = 0) : s.re = 1 / 2 := by
  have hs0 : s ≠ 0 := fun e => by simp [e] at h0
  have hs1 : s ≠ 1 := fun e => by simp [e] at h1
  by_contra hne
  rcases lt_or_gt_of_ne hne with hlt | hgt
  · have hinv : χ⁻¹.IsPrimitive := by
      rw [isPrimitive_def, conductor_inv]; exact hχ
    have hc : completedLFunction χ s = 0 := by
      rw [LFunction_eq_completed_div_gammaFactor χ s (Or.inl hs0)] at hz
      exact (div_eq_zero_iff.1 hz).resolve_right (Formalism.VW.gammaFactor_ne_zero χ h0)
    have hfe := hinv.completedLFunction_one_sub s
    rw [inv_inv, hc, mul_zero] at hfe
    have h2 : (1 - s) ≠ 0 := sub_ne_zero.2 hs1.symm
    refine h N χ⁻¹ (1 - s) (by rw [sub_re, one_re]; linarith) (Or.inr fun e => hs0 (by linear_combination -e)) ?_
    rw [LFunction_eq_completed_div_gammaFactor χ⁻¹ (1 - s) (Or.inl h2), hfe, zero_div]
  · exact h N χ s hgt (Or.inr hs1) hz

theorem Formalism.VW.rh_of_every (h : UpperTheory.EveryCharacter) : RiemannHypothesis :=
  Formalism.VW.rh_of_right fun s hs hs1 => by
    have := h 1 (1 : DirichletCharacter ℂ 1) s hs (Or.inr hs1)
    rwa [LFunction_modOne_eq] at this
