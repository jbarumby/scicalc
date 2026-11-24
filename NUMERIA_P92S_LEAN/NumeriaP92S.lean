import Mathlib.Tactic
import Mathlib.Basic.Real.Basic

/-!
# NUMERIA P92-S — First machine-checked kernel of the P92-R reconstruction

This file formalizes elementary components of P92-R in Lean 4 / Mathlib:
1. the two-pole partial-fraction identity (the algebraic heart of R7);
2. its termwise finite-sum lifting;
3. backward Bellman telescoping, *conditional on the Bellman inequality*;
4. the exact truncated-harmonic boundary accounting, *conditional on R2*;
5. the implication R6 + the two energy estimates => the sharp R7;
6. a strictly positive sharp right-hand side for gap sizes <= h and K > 3.

Not yet formalized: the universal analytic R1 inequality for ramps, its
integration against positive measures, R2 for every convex kernel, the
R5 limit to R6, and the infinite-series transfer of exterior tails.
No theorem here is a Lean proof of the full Preissmann Lemma 5/6 or RH.
-/

namespace NUMERIA.P92S

/-- The exact signed two-pole decomposition in the variable
    u = y_r - y_s and h = y_s - y_t. -/
theorem two_pole_partial_fractions (u h : ℝ)
    (hu : u ≠ 0) (hh : h ≠ 0) (huh : u + h ≠ 0) :
    1 / (u ^ 2 * (u + h) ^ 2) =
      (1 / h ^ 2) * (1 / u ^ 2 + 1 / (u + h) ^ 2) +
        (2 / h ^ 3) * (1 / (u + h) - 1 / u) := by
  field_simp
  ring

/-- Finite-sum version: each summand in the energy identity is rewritten
    exactly. The sum is finite; no convergence theorem is assumed. -/
theorem finite_two_pole_identity
    {ι : Type*} (S : Finset ι) (δ u : ι → ℝ) (h : ℝ)
    (hh : h ≠ 0)
    (hu : ∀ i ∈ S, u i ≠ 0)
    (huh : ∀ i ∈ S, u i + h ≠ 0) :
    (∑ i ∈ S, δ i / (u i ^ 2 * (u i + h) ^ 2)) =
      ∑ i ∈ S, δ i *
        ((1 / h ^ 2) * (1 / u i ^ 2 + 1 / (u i + h) ^ 2) +
          (2 / h ^ 3) * (1 / (u i + h) - 1 / u i)) := by
  apply Finset.sum_congr rfl
  intro i hi
  calc
    δ i / (u i ^ 2 * (u i + h) ^ 2)
        = δ i * (1 / (u i ^ 2 * (u i + h) ^ 2)) := by ring
    _ = _ := by rw [two_pole_partial_fractions (u i) h (hu i hi) hh (huh i hi)]

/-- A backward Bellman cost: the final point has weight h, while each
    nonfinal point is weighted by the smaller adjacent gap min(h,b). -/
def pathCost (f : ℝ → ℝ) : ℝ → ℝ → List ℝ → ℝ
  | w, h, [] => h * f w
  | w, h, b :: bs => min h b * f w + pathCost f (w + b) b bs

/-- The exact *discrete* induction R1 => R2.
    The assumptions are the two Bellman inequalities; this theorem does
    NOT pretend to derive them from convexity yet. -/
theorem pathCost_le_of_bellman (f : ℝ → ℝ) (U : ℝ → ℝ → ℝ)
    (base : ∀ w h, h * f w ≤ U w h)
    (step : ∀ w h b, min h b * f w + U (w + b) b ≤ U w h)
    (w h : ℝ) (remainingGaps : List ℝ) :
    pathCost f w h remainingGaps ≤ U w h := by
  induction remainingGaps generalizing w h with
  | nil =>
      simpa [pathCost] using base w h
  | cons b bs ih =>
      calc
        pathCost f w h (b :: bs) =
            min h b * f w + pathCost f (w + b) b bs := rfl
        _ ≤ min h b * f w + U (w + b) b :=
          add_le_add_right (ih (w + b) b) _
        _ ≤ U w h := step w h b

/-- Exact algebra that converts a truncated R2 bound into R3.
    Here B stands for H_K, so this theorem deliberately does not
    assume the harmonic estimate without a proof. -/
theorem boundary_accounting
    {ι : Type*} (S : Finset ι) (δ z : ι → ℝ)
    (X a B : ℝ) (K : ℕ)
    (hX : 0 < X)
    (hdelta : (∑ i ∈ S, δ i) ≤ X)
    (htruncated :
      (∑ i ∈ S, δ i * (1 / z i - 1 / X)) ≤ B - (K : ℝ) * a / X) :
    (∑ i ∈ S, δ i / z i) ≤ B + (X - (K : ℝ) * a) / X := by
  have hid :
      (∑ i ∈ S, δ i / z i) =
        (∑ i ∈ S, δ i * (1 / z i - 1 / X)) +
          (∑ i ∈ S, δ i) / X := by
    rw [Finset.sum_div, ← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro i hi
    ring
  have hquot : (∑ i ∈ S, δ i) / X ≤ 1 :=
    (div_le_iff₀ hX).2 (by simpa using hdelta)
  rw [hid]
  have heq : B - (K : ℝ) * a / X + 1 =
      B + (X - (K : ℝ) * a) / X := by
    field_simp
    ring
  calc
    _ ≤ B - (K : ℝ) * a / X + 1 := add_le_add htruncated hquot
    _ = _ := heq

/-- R6 (compensator) + energy bounds + exact algebra imply R7.
    The assumed equality M = ... is supplied by the finite two-pole
    identity (or an absolutely convergent infinite-sum version).
    No hidden limiting or analytic assumptions are inserted. -/
theorem sharp_bound_of_energy_and_compensator
    (K h ds dt Es Et C M : ℝ)
    (hh : 0 < h)
    (henergy_s : Es ≤ K / ds) (henergy_t : Et ≤ K / dt)
    (hcomp : C ≤ -(ds + dt) / h)
    (hexact :
      M = (Es + Et) / h ^ 2 - (ds + dt) / h ^ 4 + 2 * C / h ^ 3) :
    M ≤ K / h ^ 2 * (1 / ds + 1 / dt) -
      3 * (ds + dt) / h ^ 4 := by
  have hh2 : 0 ≤ h ^ 2 := le_of_lt (pow_pos hh 2)
  have hh3 : 0 ≤ h ^ 3 := le_of_lt (pow_pos hh 3)
  have he :
      (Es + Et) / h ^ 2 ≤ (K / ds + K / dt) / h ^ 2 :=
    div_le_div_of_nonneg_right (add_le_add henergy_s henergy_t) hh2
  have hc :
      2 * C / h ^ 3 ≤ 2 * (-(ds + dt) / h) / h ^ 3 :=
    div_le_div_of_nonneg_right
      (mul_le_mul_of_nonneg_left hcomp (by norm_num)) hh3
  rw [hexact]
  calc
    (Es + Et) / h ^ 2 - (ds + dt) / h ^ 4 + 2 * C / h ^ 3
       ≤ (K / ds + K / dt) / h ^ 2 - (ds + dt) / h ^ 4 +
           2 * (-(ds + dt) / h) / h ^ 3 := by linarith
    _ = K / h ^ 2 * (1 / ds + 1 / dt) -
          3 * (ds + dt) / h ^ 4 := by
          field_simp
          ring

/-- Algebraic positivity of the sharp two-point bound, dimensionally
    correct: the natural scale in the proof is h^3 times RHS. -/
theorem sharp_rhs_positive
    (K h ds dt : ℝ) (hK : 3 < K)
    (hh : 0 < h) (hds : 0 < ds) (hdt : 0 < dt)
    (hds_le : ds ≤ h) (hdt_le : dt ≤ h) :
    0 < K / h ^ 2 * (1 / ds + 1 / dt) -
      3 * (ds + dt) / h ^ 4 := by
  have hprod : ds * dt ≤ h ^ 2 := by
    calc
      ds * dt ≤ h * dt := mul_le_mul_of_nonneg_right hds_le (le_of_lt hdt)
      _ ≤ h * h := mul_le_mul_of_nonneg_left hdt_le (le_of_lt hh)
      _ = h ^ 2 := by ring
  have hsq : 0 < h ^ 2 := pow_pos hh 2
  have hcoeff : 0 < K * h ^ 2 - 3 * (ds * dt) := by
    nlinarith [mul_lt_mul_of_pos_right hK hsq]
  have hnumer : 0 < (K * h ^ 2 - 3 * (ds * dt)) * (ds + dt) :=
    mul_pos hcoeff (add_pos hds hdt)
  have hden : 0 < h ^ 4 * ds * dt := by positivity
  have heq :
      K / h ^ 2 * (1 / ds + 1 / dt) -
          3 * (ds + dt) / h ^ 4 =
        ((K * h ^ 2 - 3 * (ds * dt)) * (ds + dt)) /
          (h ^ 4 * ds * dt) := by
    field_simp
    ring
  rw [heq]
  exact div_pos hnumer hden


/-- The nearest-neighbor weights of a finite chain of consecutive gaps.
For an internal point the weight is the minimum of its adjacent gaps;
the right endpoint carries the final gap, exactly as in P92-R/R3. -/
def nearestMass : List ℝ → ℝ
  | [] => 0
  | [a] => a
  | a :: b :: rest => min a b + nearestMass (b :: rest)

/-- The total nearest-neighbor mass is at most the sum of consecutive gaps.
This is the missing **unconditional finite** input used for boundary accounting
in R3: the total gap sum is the interval length X. -/
theorem nearestMass_le_sum (gaps : List ℝ) :
    nearestMass gaps ≤ gaps.sum := by
  induction gaps with
  | nil => simp [nearestMass]
  | cons a rest ih =>
      cases rest with
      | nil => simp [nearestMass]
      | cons b tail =>
          change min a b + nearestMass (b :: tail) ≤ a + (b :: tail).sum
          exact add_le_add (min_le_left a b) ih

end NUMERIA.P92S
