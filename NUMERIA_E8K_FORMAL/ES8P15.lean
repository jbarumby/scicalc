import Mathlib.Tactic
import ES8P14

/-!
NUMERIA ES8-P15 — conditional tail, tsum and dimension-uniform lemmas.

The analytic E7 estimate for the actual digamma/Fourier integral is
a hypothesis, not established here. The spectral gap and row-error
bounds are also hypotheses for arbitrary n, not verified for NUMERIA.
-/

open Filter
open scoped Topology BigOperators

namespace Numeria.ES8P15

noncomputable def rationalTail (A B R : ℝ) : ℝ := A / R ^ 2 + B / R ^ 3

theorem log_tail_bounded (R : ℝ) (hR : 1 ≤ R) :
    (8 + Real.log (3 * R) + (1 : ℝ) / 3) / (3 * R ^ 3) ≤
      rationalTail 1 3 R := by
  have hpos : 0 < R := by linarith
  have hlog : Real.log (3 * R) ≤ 3 * R :=
    Real.log_le_self (by positivity)
  have heq : (3 * R + 9) / (3 * R ^ 3) = rationalTail 1 3 R := by
    unfold rationalTail
    field_simp [ne_of_gt hpos]
    <;> ring
  calc
    _ ≤ (3 * R + 9) / (3 * R ^ 3) := by
      have hnum : 8 + Real.log (3 * R) + (1 : ℝ) / 3 ≤ 3 * R + 9 := by
        linarith [hlog]
      exact (div_le_div_iff₀ (by positivity) (by positivity)).2
        (mul_le_mul_of_nonneg_right hnum (by positivity))
    _ = _ := heq

theorem e7_bound_implies_rational_enclosure
    (K R T : ℝ) (hK : 0 ≤ K) (hR : 1 ≤ R)
    (hE7 : |T| ≤ K * ((8 + Real.log (3 * R) + (1 : ℝ) / 3) /
      (3 * R ^ 3))) :
    |T| ≤ K * rationalTail 1 3 R := by
  exact hE7.trans (mul_le_mul_of_nonneg_left (log_tail_bounded R hR) hK)

theorem rationalTail_tendsto (A B : ℝ) :
    Tendsto (fun n : ℕ => rationalTail A B ((n : ℝ) + 1)) atTop (𝓝 0) := by
  have hi : Tendsto (fun n : ℕ => (1 : ℝ) / ((n : ℝ) + 1))
      atTop (𝓝 0) := tendsto_one_div_add_atTop_nhds_zero_nat
  have h2 := (tendsto_const_nhds (x := A)).mul (hi.pow 2)
  have h3 := (tendsto_const_nhds (x := B)).mul (hi.pow 3)
  have hsum := h2.add h3
  convert hsum using 1 <;> simp [rationalTail, div_eq_mul_inv, one_div]

theorem bounded_archTail_tendsto (T : ℕ → ℝ) (A B : ℝ)
    (hbound : ∀ n : ℕ, |T n| ≤ rationalTail A B ((n : ℝ) + 1)) :
    Tendsto T atTop (𝓝 0) := by
  have hlim := rationalTail_tendsto A B
  have hneg : Tendsto (fun n : ℕ => -rationalTail A B ((n : ℝ) + 1))
      atTop (𝓝 0) := by simpa using hlim.neg
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le hneg hlim
    (fun n => (abs_le.mp (hbound n)).1)
    (fun n => (abs_le.mp (hbound n)).2)

theorem tsum_margin_of_partial_sums (f g : ℕ → ℝ)
    (hf : Summable f) (hg : Summable g)
    (hfinite : ∀ n : ℕ,
      (1 : ℝ) / 20 * (∑ k ∈ Finset.range n, g k) ≤
        (∑ k ∈ Finset.range n, f k)) :
    (1 : ℝ) / 20 * (∑' k, g k) ≤ ∑' k, f k := by
  have hlim : Tendsto
      (fun n : ℕ => (∑ k ∈ Finset.range n, f k) -
        (1 : ℝ) / 20 * (∑ k ∈ Finset.range n, g k))
      atTop (𝓝 ((∑' k, f k) - (1 : ℝ) / 20 * (∑' k, g k))) :=
    hf.tendsto_sum_tsum_nat.sub
      ((tendsto_const_nhds (x := (1 : ℝ) / 20)).mul hg.tendsto_sum_tsum_nat)
  have hn : 0 ≤ (∑' k, f k) - (1 : ℝ) / 20 * (∑' k, g k) :=
    ge_of_tendsto hlim
      (Eventually.of_forall (fun n => sub_nonneg.mpr (hfinite n)))
  linarith

theorem quadratic_margin_at_limit (q s : ℕ → ℝ) (Q S : ℝ)
    (hq : Tendsto q atTop (𝓝 Q))
    (hs : Tendsto s atTop (𝓝 S))
    (hfinite : ∀ n : ℕ, (1 : ℝ) / 20 * s n ≤ q n) :
    (1 : ℝ) / 20 * S ≤ Q := by
  have hlim : Tendsto (fun n : ℕ => q n - (1 : ℝ) / 20 * s n)
      atTop (𝓝 (Q - (1 : ℝ) / 20 * S)) :=
    hq.sub ((tendsto_const_nhds (x := (1 : ℝ) / 20)).mul hs)
  have hn : 0 ≤ Q - (1 : ℝ) / 20 * S :=
    ge_of_tendsto hlim
      (Eventually.of_forall (fun n => sub_nonneg.mpr (hfinite n)))
  linarith

theorem uniform_margin_all_dimensions
    (M C eps : (n : ℕ) → Fin n → Fin n → ℚ)
    (hcenter : ∀ (n : ℕ) (x : Fin n → ℚ),
      (∃ i, x i ≠ 0) →
        0 < ES8P14.qForm
          (fun i j => C n i j - if i = j then (3 : ℚ) / 50 else 0) x)
    (hsym : ∀ (n : ℕ) (i j : Fin n), eps n i j = eps n j i)
    (hentry : ∀ (n : ℕ) (i j : Fin n),
      -eps n i j ≤ M n i j - C n i j ∧
        M n i j - C n i j ≤ eps n i j)
    (hrow : ∀ (n : ℕ) (i : Fin n),
      (∑ j : Fin n, eps n i j) ≤ (1 : ℚ) / 100)
    (n : ℕ) (x : Fin n → ℚ) (hx : ∃ i, x i ≠ 0) :
    (1 : ℚ) / 20 * (∑ i : Fin n, x i ^ 2) < ES8P14.qForm (M n) x := by
  exact ES8P14.positive_from_LDL_and_box (M n) (C n) (eps n) x
    (hcenter n x hx) (hsym n) (hentry n) (hrow n)

end Numeria.ES8P15
