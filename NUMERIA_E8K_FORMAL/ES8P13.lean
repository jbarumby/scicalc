import Mathlib.Tactic

/-!
NUMERIA–RIEMANN ES8-P13: concrete exact rational LDL certificates.

This proves only the POSITIVITY OF THE EXPLICIT RATIONAL COARSE MATRICES Cn - 3/50 * I.
The exact 90-digit original arithmetic endpoints are checked by a separate Python
reference replay: max absolute row sum of A_n - Cn <= 1/100.
The bridge from the special-function integral to those endpoints, and the
general perturbation lemma, are not formalized by this file.
Thus there is NO Lean proof of the complete numerical-to-analytic chain or RH.
-/

namespace Numeria.ES8P13


/-- For the actual *finite rational* matrix represented by its LDL witness. -/
def qForm {n : ℕ} (C : Fin n → Fin n → ℚ) (x : Fin n → ℚ) : ℚ :=
  ∑ i : Fin n, ∑ j : Fin n, x i * C i j * x j

def factoredQ {n : ℕ} (d : Fin n → ℚ) (U : Fin n → Fin n → ℚ)
    (x : Fin n → ℚ) : ℚ :=
  ∑ k : Fin n, d k * (∑ i : Fin n, U k i * x i) ^ 2

/-- Pure algebra: a coefficientwise LDL identity gives equality of quadratic forms. -/
theorem qForm_eq_factored {n : ℕ} (C : Fin n → Fin n → ℚ)
    (d : Fin n → ℚ) (U : Fin n → Fin n → ℚ)
    (h : ∀ i j : Fin n, C i j = ∑ k : Fin n, d k * U k i * U k j)
    (x : Fin n → ℚ) : qForm C x = factoredQ d U x := by
  unfold qForm factoredQ
  simp_rw [h]
  calc
    (∑ i : Fin n, ∑ j : Fin n,
        x i * (∑ k : Fin n, d k * U k i * U k j) * x j) =
      ∑ i : Fin n, ∑ j : Fin n, ∑ k : Fin n,
        d k * U k i * x i * U k j * x j := by
          apply Finset.sum_congr rfl
          intro i _
          apply Finset.sum_congr rfl
          intro j _
          rw [Finset.mul_sum, Finset.sum_mul]
          exact Finset.sum_congr rfl (by intro k _; ring)
    _ = ∑ k : Fin n, ∑ i : Fin n, ∑ j : Fin n,
          d k * U k i * x i * U k j * x j := by
            calc
              _ = ∑ i : Fin n, ∑ k : Fin n, ∑ j : Fin n,
                  d k * U k i * x i * U k j * x j := by
                    apply Finset.sum_congr rfl
                    intro i _
                    rw [Finset.sum_comm]
              _ = _ := by rw [Finset.sum_comm]
    _ = ∑ k : Fin n, d k * (∑ i : Fin n, U k i * x i) ^ 2 := by
        apply Finset.sum_congr rfl
        intro k _
        simp only [pow_two, Finset.mul_sum, Finset.sum_mul]
        apply Finset.sum_congr rfl
        intro i _
        apply Finset.sum_congr rfl
        intro j _
        ring

/-- The standard finite SOS strict positivity criterion. -/
theorem factoredQ_pos {n : ℕ} (d : Fin n → ℚ)
    (U : Fin n → Fin n → ℚ) (hd : ∀ i, 0 < d i)
    (hinj : ∀ x : Fin n → ℚ,
       (∀ i, (∑ j : Fin n, U i j * x j) = 0) → ∀ i, x i = 0)
    (x : Fin n → ℚ) (hx : ∃ i, x i ≠ 0) : 0 < factoredQ d U x := by
  have hy : ∃ i : Fin n, (∑ j : Fin n, U i j * x j) ≠ 0 := by
    by_contra h
    push Not at h
    obtain ⟨j, hj⟩ := hx
    exact hj (hinj x h j)
  obtain ⟨i, hi⟩ := hy
  unfold factoredQ
  apply Finset.sum_pos'
  · intro k _
    exact mul_nonneg (le_of_lt (hd k)) (sq_nonneg _)
  · refine ⟨i, Finset.mem_univ _, ?_⟩
    exact mul_pos (hd i) (sq_pos_of_ne_zero hi)


/-! ## Concrete projected matrix N=6, projected dimension 4 -/

def matrix4 (i j : Fin 4) : ℚ :=
  match i.val, j.val with
  | 0, 0 => (429 : ℚ) / 2500
  | 0, 2 => (-647 : ℚ) / 2500
  | 1, 1 => (2081 : ℚ) / 2500
  | 1, 3 => (413 : ℚ) / 10000
  | 2, 0 => (-647 : ℚ) / 2500
  | 2, 2 => (1873 : ℚ) / 2500
  | 3, 1 => (413 : ℚ) / 10000
  | 3, 3 => (10753 : ℚ) / 10000
  | _, _ => 0
def pivots4 (i : Fin 4) : ℚ :=
  match i.val with
  | 0 => (429 : ℚ) / 2500
  | 1 => (2081 : ℚ) / 2500
  | 2 => (96227 : ℚ) / 268125
  | 3 => (89337403 : ℚ) / 83240000
  | _ => 0
def upper4 (i j : Fin 4) : ℚ :=
  match i.val, j.val with
  | 0, 0 => 1
  | 0, 2 => (-647 : ℚ) / 429
  | 1, 1 => 1
  | 1, 3 => (413 : ℚ) / 8324
  | 2, 2 => 1
  | 3, 3 => 1
  | _, _ => 0

set_option maxHeartbeats 6000000 in
theorem factor_entries4 :
    ∀ i j : Fin 4, matrix4 i j =
      ∑ k : Fin 4, pivots4 k * upper4 k i * upper4 k j := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    norm_num [matrix4, pivots4, upper4, Fin.sum_univ_succ]

theorem pivots_pos4 : ∀ i : Fin 4, 0 < pivots4 i := by
  intro i
  fin_cases i <;> norm_num [pivots4]


theorem upper_inj4 (x : Fin 4 → ℚ)
    (h : ∀ i : Fin 4, (∑ j : Fin 4, upper4 i j * x j) = 0) :
    ∀ i : Fin 4, x i = 0 := by

  have h3 : x (3 : Fin 4) = 0 := by

    have hi := h (3 : Fin 4)

    simp [upper4, Fin.sum_univ_succ] at hi

    linarith

  have h2 : x (2 : Fin 4) = 0 := by

    have hi := h (2 : Fin 4)

    simp [upper4, Fin.sum_univ_succ, h3] at hi

    linarith

  have h1 : x (1 : Fin 4) = 0 := by

    have hi := h (1 : Fin 4)

    simp [upper4, Fin.sum_univ_succ, h3, h2] at hi

    linarith

  have h0 : x (0 : Fin 4) = 0 := by

    have hi := h (0 : Fin 4)

    simp [upper4, Fin.sum_univ_succ, h3, h2, h1] at hi

    linarith

  intro i
  fin_cases i

  all_goals assumption


theorem concrete_positive4 (x : Fin 4 → ℚ)
    (hx : ∃ i, x i ≠ 0) : 0 < qForm matrix4 x := by
  rw [qForm_eq_factored matrix4 pivots4 upper4 factor_entries4]
  exact factoredQ_pos pivots4 upper4 pivots_pos4 upper_inj4 x hx


/-! ## Concrete projected matrix N=8, projected dimension 6 -/

def matrix6 (i j : Fin 6) : ℚ :=
  match i.val, j.val with
  | 0, 0 => (429 : ℚ) / 2500
  | 0, 2 => (-647 : ℚ) / 2500
  | 0, 4 => (-1741 : ℚ) / 10000
  | 1, 1 => (2081 : ℚ) / 2500
  | 1, 3 => (413 : ℚ) / 10000
  | 1, 5 => (-343 : ℚ) / 10000
  | 2, 0 => (-647 : ℚ) / 2500
  | 2, 2 => (1873 : ℚ) / 2500
  | 2, 4 => (-787 : ℚ) / 10000
  | 3, 1 => (413 : ℚ) / 10000
  | 3, 3 => (10753 : ℚ) / 10000
  | 3, 5 => (-1783 : ℚ) / 10000
  | 4, 0 => (-1741 : ℚ) / 10000
  | 4, 2 => (-787 : ℚ) / 10000
  | 4, 4 => (6887 : ℚ) / 5000
  | 5, 1 => (-343 : ℚ) / 10000
  | 5, 3 => (-1783 : ℚ) / 10000
  | 5, 5 => (5777 : ℚ) / 5000
  | _, _ => 0
def pivots6 (i : Fin 6) : ℚ :=
  match i.val with
  | 0 => (429 : ℚ) / 2500
  | 1 => (2081 : ℚ) / 2500
  | 2 => (96227 : ℚ) / 268125
  | 3 => (89337403 : ℚ) / 83240000
  | 4 => (843185691 : ℚ) / 962270000
  | 5 => (1004981693723 : ℚ) / 893374030000
  | _ => 0
def upper6 (i j : Fin 6) : ℚ :=
  match i.val, j.val with
  | 0, 0 => 1
  | 0, 2 => (-647 : ℚ) / 429
  | 0, 4 => (-1741 : ℚ) / 1716
  | 1, 1 => 1
  | 1, 3 => (413 : ℚ) / 8324
  | 1, 5 => (-343 : ℚ) / 8324
  | 2, 2 => 1
  | 2, 4 => (-732025 : ℚ) / 769816
  | 3, 3 => 1
  | 3, 5 => (-14700033 : ℚ) / 89337403
  | 4, 4 => 1
  | 5, 5 => 1
  | _, _ => 0

set_option maxHeartbeats 6000000 in
theorem factor_entries6 :
    ∀ i j : Fin 6, matrix6 i j =
      ∑ k : Fin 6, pivots6 k * upper6 k i * upper6 k j := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    norm_num [matrix6, pivots6, upper6, Fin.sum_univ_succ]

theorem pivots_pos6 : ∀ i : Fin 6, 0 < pivots6 i := by
  intro i
  fin_cases i <;> norm_num [pivots6]


theorem upper_inj6 (x : Fin 6 → ℚ)
    (h : ∀ i : Fin 6, (∑ j : Fin 6, upper6 i j * x j) = 0) :
    ∀ i : Fin 6, x i = 0 := by

  have h5 : x (5 : Fin 6) = 0 := by

    have hi := h (5 : Fin 6)

    simp [upper6, Fin.sum_univ_succ] at hi

    linarith

  have h4 : x (4 : Fin 6) = 0 := by

    have hi := h (4 : Fin 6)

    simp [upper6, Fin.sum_univ_succ, h5] at hi

    linarith

  have h3 : x (3 : Fin 6) = 0 := by

    have hi := h (3 : Fin 6)

    simp [upper6, Fin.sum_univ_succ, h5, h4] at hi

    linarith

  have h2 : x (2 : Fin 6) = 0 := by

    have hi := h (2 : Fin 6)

    simp [upper6, Fin.sum_univ_succ, h5, h4, h3] at hi

    linarith

  have h1 : x (1 : Fin 6) = 0 := by

    have hi := h (1 : Fin 6)

    simp [upper6, Fin.sum_univ_succ, h5, h4, h3, h2] at hi

    linarith

  have h0 : x (0 : Fin 6) = 0 := by

    have hi := h (0 : Fin 6)

    simp [upper6, Fin.sum_univ_succ, h5, h4, h3, h2, h1] at hi

    linarith

  intro i
  fin_cases i

  all_goals assumption


theorem concrete_positive6 (x : Fin 6 → ℚ)
    (hx : ∃ i, x i ≠ 0) : 0 < qForm matrix6 x := by
  rw [qForm_eq_factored matrix6 pivots6 upper6 factor_entries6]
  exact factoredQ_pos pivots6 upper6 pivots_pos6 upper_inj6 x hx


/-! ## Concrete projected matrix N=10, projected dimension 8 -/

def matrix8 (i j : Fin 8) : ℚ :=
  match i.val, j.val with
  | 0, 0 => (429 : ℚ) / 2500
  | 0, 2 => (-647 : ℚ) / 2500
  | 0, 4 => (-1741 : ℚ) / 10000
  | 0, 6 => (-197 : ℚ) / 2000
  | 1, 1 => (2081 : ℚ) / 2500
  | 1, 3 => (413 : ℚ) / 10000
  | 1, 5 => (-343 : ℚ) / 10000
  | 1, 7 => (-167 : ℚ) / 2500
  | 2, 0 => (-647 : ℚ) / 2500
  | 2, 2 => (1873 : ℚ) / 2500
  | 2, 4 => (-787 : ℚ) / 10000
  | 2, 6 => (47 : ℚ) / 2500
  | 3, 1 => (413 : ℚ) / 10000
  | 3, 3 => (10753 : ℚ) / 10000
  | 3, 5 => (-1783 : ℚ) / 10000
  | 3, 7 => (-987 : ℚ) / 5000
  | 4, 0 => (-1741 : ℚ) / 10000
  | 4, 2 => (-787 : ℚ) / 10000
  | 4, 4 => (6887 : ℚ) / 5000
  | 4, 6 => (391 : ℚ) / 2500
  | 5, 1 => (-343 : ℚ) / 10000
  | 5, 3 => (-1783 : ℚ) / 10000
  | 5, 5 => (5777 : ℚ) / 5000
  | 5, 7 => (-2439 : ℚ) / 10000
  | 6, 0 => (-197 : ℚ) / 2000
  | 6, 2 => (47 : ℚ) / 2500
  | 6, 4 => (391 : ℚ) / 2500
  | 6, 6 => (17243 : ℚ) / 10000
  | 7, 1 => (-167 : ℚ) / 2500
  | 7, 3 => (-987 : ℚ) / 5000
  | 7, 5 => (-2439 : ℚ) / 10000
  | 7, 7 => (15071 : ℚ) / 10000
  | _, _ => 0
def pivots8 (i : Fin 8) : ℚ :=
  match i.val with
  | 0 => (429 : ℚ) / 2500
  | 1 => (2081 : ℚ) / 2500
  | 2 => (96227 : ℚ) / 268125
  | 3 => (89337403 : ℚ) / 83240000
  | 4 => (843185691 : ℚ) / 962270000
  | 5 => (1004981693723 : ℚ) / 893374030000
  | 6 => (96879442105713 : ℚ) / 59959871360000
  | 7 => (7023056433067541 : ℚ) / 5024908468615000
  | _ => 0
def upper8 (i j : Fin 8) : ℚ :=
  match i.val, j.val with
  | 0, 0 => 1
  | 0, 2 => (-647 : ℚ) / 429
  | 0, 4 => (-1741 : ℚ) / 1716
  | 0, 6 => (-985 : ℚ) / 1716
  | 1, 1 => 1
  | 1, 3 => (413 : ℚ) / 8324
  | 1, 5 => (-343 : ℚ) / 8324
  | 1, 7 => (-167 : ℚ) / 2081
  | 2, 2 => 1
  | 2, 4 => (-732025 : ℚ) / 769816
  | 2, 6 => (-556643 : ℚ) / 1539632
  | 3, 3 => 1
  | 3, 5 => (-14700033 : ℚ) / 89337403
  | 3, 7 => (-16155692 : ℚ) / 89337403
  | 4, 4 => 1
  | 4, 6 => (-171717487 : ℚ) / 2248495176
  | 5, 5 => 1
  | 5, 7 => (-248883660259 : ℚ) / 1004981693723
  | 6, 6 => 1
  | 7, 7 => 1
  | _, _ => 0

set_option maxHeartbeats 6000000 in
theorem factor_entries8 :
    ∀ i j : Fin 8, matrix8 i j =
      ∑ k : Fin 8, pivots8 k * upper8 k i * upper8 k j := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    norm_num [matrix8, pivots8, upper8, Fin.sum_univ_succ]

theorem pivots_pos8 : ∀ i : Fin 8, 0 < pivots8 i := by
  intro i
  fin_cases i <;> norm_num [pivots8]


theorem upper_inj8 (x : Fin 8 → ℚ)
    (h : ∀ i : Fin 8, (∑ j : Fin 8, upper8 i j * x j) = 0) :
    ∀ i : Fin 8, x i = 0 := by

  have h7 : x (7 : Fin 8) = 0 := by

    have hi := h (7 : Fin 8)

    simp [upper8, Fin.sum_univ_succ] at hi

    linarith

  have h6 : x (6 : Fin 8) = 0 := by

    have hi := h (6 : Fin 8)

    simp [upper8, Fin.sum_univ_succ, h7] at hi

    linarith

  have h5 : x (5 : Fin 8) = 0 := by

    have hi := h (5 : Fin 8)

    simp [upper8, Fin.sum_univ_succ, h7, h6] at hi

    linarith

  have h4 : x (4 : Fin 8) = 0 := by

    have hi := h (4 : Fin 8)

    simp [upper8, Fin.sum_univ_succ, h7, h6, h5] at hi

    linarith

  have h3 : x (3 : Fin 8) = 0 := by

    have hi := h (3 : Fin 8)

    simp [upper8, Fin.sum_univ_succ, h7, h6, h5, h4] at hi

    linarith

  have h2 : x (2 : Fin 8) = 0 := by

    have hi := h (2 : Fin 8)

    simp [upper8, Fin.sum_univ_succ, h7, h6, h5, h4, h3] at hi

    linarith

  have h1 : x (1 : Fin 8) = 0 := by

    have hi := h (1 : Fin 8)

    simp [upper8, Fin.sum_univ_succ, h7, h6, h5, h4, h3, h2] at hi

    linarith

  have h0 : x (0 : Fin 8) = 0 := by

    have hi := h (0 : Fin 8)

    simp [upper8, Fin.sum_univ_succ, h7, h6, h5, h4, h3, h2, h1] at hi

    linarith

  intro i
  fin_cases i

  all_goals assumption


theorem concrete_positive8 (x : Fin 8 → ℚ)
    (hx : ∃ i, x i ≠ 0) : 0 < qForm matrix8 x := by
  rw [qForm_eq_factored matrix8 pivots8 upper8 factor_entries8]
  exact factoredQ_pos pivots8 upper8 pivots_pos8 upper_inj8 x hx


end Numeria.ES8P13
