import GlobalRamp
import Mathlib.Tactic

/-!
P92-T4: identify the shifted ramp potential exactly. This turns
the fullRamp_bellman height formulation into the genuine Bellman R1
for f_T(w)=max(T-w,0) and U0(w,h) defined by the finite ceiling sum.

This is only the ramp case. Positive mixtures over an arbitrary measure
and the infinite two-index compensator still require formalization.
-/

namespace NUMERIA.P92T

/-- Shifting the positive-part finite sum deletes its first term
    and shifts the remaining mesh indices by one. -/
theorem finiteRamp0_shift (L b : ℝ) (n : ℕ) (hL : 0 ≤ L) :
    finiteRamp0 L b (n + 1) - b * L =
      finiteRamp0 (L - b) b n := by
  unfold finiteRamp0
  rw [Finset.sum_range_succ']
  simp only [Nat.cast_zero, zero_mul, sub_zero, max_eq_left hL]
  have hsum :
      (∑ j ∈ Finset.range n,
          max (L - ((j + 1 : ℕ) : ℝ) * b) 0) =
        ∑ j ∈ Finset.range n,
          max ((L - b) - (j : ℝ) * b) 0 := by
    apply Finset.sum_congr rfl
    intro j hj
    have harg : L - ((j + 1 : ℕ) : ℝ) * b =
        (L - b) - (j : ℝ) * b := by
      push_cast
      ring
    rw [harg]
  rw [hsum]
  ring

/-- The shifted full potential U₁(L,b) is exactly U₀((L-b)₊,b).
    Every support-size transition is accounted for by Nat.ceil_add_one. -/
theorem fullRamp_shift (L b : ℝ) (hL : 0 ≤ L) (hb : 0 < b) :
    fullRamp1 L b = fullRamp0 (max (L - b) 0) b := by
  rcases le_total L b with hsmall | hlarge
  · by_cases hz : L = 0
    · subst L
      have hmax : max (-b) (0 : ℝ) = 0 :=
        max_eq_right (neg_nonpos.mpr (le_of_lt hb))
      simp [fullRamp0, fullRamp1, finiteRamp0, hmax]
    · have hpos : 0 < L := lt_of_le_of_ne hL (Ne.symm hz)
      have hfrac : L / b ≤ 1 := by
        apply (div_le_iff₀ hb).2
        simpa using hsmall
      have hceil_le : Nat.ceil (L / b) ≤ 1 := by
        apply (Nat.ceil_le).2
        simpa only [Nat.cast_one] using hfrac
      have hceil_pos : 1 ≤ Nat.ceil (L / b) :=
        (Nat.one_le_ceil_iff).2 (div_pos hpos hb)
      have hceil : Nat.ceil (L / b) = 1 := by omega
      have hmax : max (L - b) 0 = 0 :=
        max_eq_right (sub_nonpos.mpr hsmall)
      calc
        fullRamp1 L b = b * L - b * L := by
          simp [fullRamp1, fullRamp0, finiteRamp0, hceil,
            max_eq_left hL]
        _ = 0 := by ring
        _ = fullRamp0 (max (L - b) 0) b := by
          simp [hmax, fullRamp0, finiteRamp0]
  · have hnew : 0 ≤ (L - b) / b :=
      div_nonneg (sub_nonneg.mpr hlarge) (le_of_lt hb)
    have hdiv : L / b = (L - b) / b + 1 := by
      field_simp [ne_of_gt hb]
      ring
    have hcount :
        Nat.ceil (L / b) = Nat.ceil ((L - b) / b) + 1 := by
      rw [hdiv, Nat.ceil_add_one hnew]
    have hmax : max (L - b) 0 = L - b :=
      max_eq_left (sub_nonneg.mpr hlarge)
    calc
      fullRamp1 L b =
          finiteRamp0 L b (Nat.ceil (L / b)) - b * L := rfl
      _ = finiteRamp0 L b (Nat.ceil ((L - b) / b) + 1) - b * L := by
        rw [hcount]
      _ = finiteRamp0 (L - b) b (Nat.ceil ((L - b) / b)) :=
        finiteRamp0_shift L b _ hL
      _ = fullRamp0 (max (L - b) 0) b := by
        rw [hmax]
        rfl

/-- The elementary positive-part identity for translating the origin
    of a ramp by b≥0. -/
theorem positive_part_sub_shift (A b : ℝ) (hb : 0 ≤ b) :
    max (A - b) 0 = max (max A 0 - b) 0 := by
  rcases le_total A 0 with hneg | hpos
  · have hab : A - b ≤ 0 := by linarith
    rw [max_eq_right hab, max_eq_right hneg]
    have hnb : -b ≤ 0 := neg_nonpos.mpr hb
    simpa using (max_eq_right hnb).symm
  · rw [max_eq_left hpos]

/-- True Bellman for the ramp f_T(w)=(T-w)_+, with the potential
    as a positive-part sum and the shifted state w+b on the LHS.
    No support-count hypothesis and no assumed Bellman step. -/
noncomputable def rampPotential (T w h : ℝ) : ℝ :=
  fullRamp0 (max (T - w) 0) h

theorem bellman_for_all_ramps (T w h b : ℝ)
    (hh : 0 < h) (hb : 0 < b) :
    min h b * max (T - w) 0 +
        rampPotential T (w + b) b ≤ rampPotential T w h := by
  let L : ℝ := max (T - w) 0
  have hL : 0 ≤ L := le_max_right _ _
  have hshiftHeight : max (T - (w + b)) 0 =
      max (L - b) 0 := by
    have harg : T - (w + b) = (T - w) - b := by ring
    rw [harg]
    exact positive_part_sub_shift (T - w) b (le_of_lt hb)
  change min h b * L + fullRamp0 (max (T - (w + b)) 0) b ≤
      fullRamp0 L h
  rw [hshiftHeight, ← fullRamp_shift L b hL hb]
  exact fullRamp_bellman L h b hL hh hb


/-- Unconditional Bellman for every FINITE nonnegative mixture of
    shifted ramps, with arbitrary real support parameters Tᵢ.
    Infinite positive mixtures still require Tonelli in Lean. -/
theorem bellman_finite_ramp_mixture
    {ι : Type*} (S : Finset ι) (T c : ι → ℝ)
    (hc : ∀ i ∈ S, 0 ≤ c i)
    (w h b : ℝ) (hh : 0 < h) (hb : 0 < b) :
    min h b * (∑ i ∈ S, c i * max (T i - w) 0) +
        (∑ i ∈ S, c i * rampPotential (T i) (w + b) b) ≤
      ∑ i ∈ S, c i * rampPotential (T i) w h := by
  calc
    min h b * (∑ i ∈ S, c i * max (T i - w) 0) +
        (∑ i ∈ S, c i * rampPotential (T i) (w + b) b) =
      ∑ i ∈ S, c i *
        (min h b * max (T i - w) 0 +
          rampPotential (T i) (w + b) b) := by
            rw [Finset.mul_sum, ← Finset.sum_add_distrib]
            apply Finset.sum_congr rfl
            intro i hi
            ring
    _ ≤ ∑ i ∈ S, c i * rampPotential (T i) w h := by
      apply Finset.sum_le_sum
      intro i hi
      exact mul_le_mul_of_nonneg_left
        (bellman_for_all_ramps (T i) w h b hh hb) (hc i hi)

end NUMERIA.P92T
