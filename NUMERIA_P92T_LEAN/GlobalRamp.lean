import NumeriaP92T
import Mathlib.Tactic
import Mathlib.Algebra.Order.Floor.Semiring

/-!
P92-T3: global ramp-potential monotonicity via exact Nat.ceil support count.

This module uses the arbitrary-cell induction proved in NumeriaP92T.lean.
It introduces no assumptions about the number of cells between h and b:
the number is chosen as Nat.ceil (L/h), and the ordering of those counts
follows from Nat.ceil_mono. The only analytic prerequisites are positivity
of the real step and of the (truncated) ramp height.

Full Preissmann Lemmas 5/6 and RH are NOT formalized by this module.
-/

namespace NUMERIA.P92T

/-- A positive ramp has finitely many active terms, and Nat.ceil(L/h)
    is the exact support count (with boundary equality allowed). -/
theorem ceil_ramp_support (L h : ℝ) (hL : 0 < L) (hh : 0 < h) :
    1 ≤ Nat.ceil (L / h) ∧
      (((Nat.ceil (L / h) : ℕ) : ℝ) - 1) * h ≤ L ∧
      L ≤ ((Nat.ceil (L / h) : ℕ) : ℝ) * h := by
  have hpos : 0 < L / h := div_pos hL hh
  have hc : (L / h) * h = L := by
    field_simp [ne_of_gt hh]
  have hceil : L / h ≤ (Nat.ceil (L / h) : ℝ) := Nat.le_ceil _
  have hceilupper :
      (Nat.ceil (L / h) : ℝ) < L / h + 1 :=
    Nat.ceil_lt_add_one (le_of_lt hpos)
  have hleft :
      ((Nat.ceil (L / h) : ℝ) - 1) * h ≤ L := by
    have ha : (Nat.ceil (L / h) : ℝ) - 1 < L / h := by
      linarith
    have hm := mul_lt_mul_of_pos_right ha hh
    linarith [hc]
  have hright :
      L ≤ (Nat.ceil (L / h) : ℝ) * h := by
    have hm := mul_le_mul_of_nonneg_right hceil (le_of_lt hh)
    linarith [hc]
  exact ⟨(Nat.one_le_ceil_iff).2 hpos, hleft, hright⟩

/-- The full (finite) potential of the positive-part ramp, with the
    support count inferred automatically from h. -/
noncomputable def fullRamp0 (L h : ℝ) : ℝ :=
  finiteRamp0 L h (Nat.ceil (L / h))

/-- The full shifted potential U1=U0-hL. -/
noncomputable def fullRamp1 (L h : ℝ) : ℝ :=
  fullRamp0 L h - h * L

/-- Global, unconditionally quantified monotonicity of both ramp
    potentials. The proof now crosses an arbitrary number of support
    switches, and each chosen support count comes from Nat.ceil. -/
theorem fullRamp_both_mono (L a b : ℝ)
    (hL : 0 ≤ L) (ha : 0 < a) (hb : 0 < b) (hab : a ≤ b) :
    fullRamp0 L a ≤ fullRamp0 L b ∧
      fullRamp1 L b ≤ fullRamp1 L a := by
  rcases eq_or_lt_of_le hL with hzero | hpos
  · subst L
    simp [fullRamp0, fullRamp1, finiteRamp0]
  · have hA := ceil_ramp_support L a hpos ha
    have hB := ceil_ramp_support L b hpos hb
    have hfrac : L / b ≤ L / a := by
      apply (div_le_div_iff₀ hb ha).2
      exact mul_le_mul_of_nonneg_left hab (le_of_lt hL)
    have hn : Nat.ceil (L / b) ≤ Nat.ceil (L / a) :=
      Nat.ceil_mono hfrac
    have hcompare :=
      finiteRamp_mono_across_cells L a b
        (Nat.ceil (L / b)) (Nat.ceil (L / a))
        hB.1 hn ha hb hab hA.2.1 hA.2.2 hB.2.1 hB.2.2
    simpa [fullRamp0, fullRamp1] using hcompare

/-- Bellman's inequality for every positive L≥0 and every two positive
    real step lengths, with NO same-support hypothesis or abstract
    Bellman-step assumption. -/
theorem fullRamp_bellman (L h b : ℝ)
    (hL : 0 ≤ L) (hh : 0 < h) (hb : 0 < b) :
    min h b * L + fullRamp1 L b ≤ fullRamp0 L h := by
  rcases le_total b h with hbh | hhb
  · have hmono := (fullRamp_both_mono L b h hL hb hh hbh).1
    rw [min_eq_right hbh]
    unfold fullRamp1
    linarith
  · have hmono := (fullRamp_both_mono L h b hL hh hb hhb).2
    rw [min_eq_left hhb]
    unfold fullRamp1 at *
    linarith

end NUMERIA.P92T
