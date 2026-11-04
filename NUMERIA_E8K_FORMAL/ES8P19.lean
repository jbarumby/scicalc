import Mathlib.Tactic

/-!
# NUMERIA–RIEMANN ES8-P19 — finite-support quadratic forms and pole blocks

There are no added axioms, `sorry` or `admit`.
The statements are algebraic and conditional on their explicit assumptions.
They do not identify an analytic Weil form, prove that it is bounded/closed,
or establish positivity beyond the supplied finite-section hypotheses.

For a two-dimensional pole block, `projectedBlock - Ahh` is the sum of
FOUR outer products (two for each pole coordinate). This is the genuine
finite-rank mechanism. The unrelated diagonal shift (3/50) I is NOT
finite-rank in infinite-dimensional ℓ².
-/

namespace Numeria.ES8P19

variable {n : ℕ}

/-- Four blocks of an arbitrary rational matrix; symmetry is not assumed. -/
def projectedBlock
    (A00 : Fin 2 → Fin 2 → ℚ)
    (A0h : Fin 2 → Fin n → ℚ)
    (Ah0 : Fin n → Fin 2 → ℚ)
    (Ahh : Fin n → Fin n → ℚ)
    (C : Fin 2 → Fin n → ℚ)
    (i j : Fin n) : ℚ :=
  Ahh i j
  - (∑ a : Fin 2, C a i * A0h a j)
  - (∑ a : Fin 2, Ah0 i a * C a j)
  + (∑ a : Fin 2, ∑ b : Fin 2, C a i * A00 a b * C b j)

/-- The difference from the H-H block is a sum of four outer products:
    columns C₀,C₁,-Ah0(·,0),-Ah0(·,1), with corresponding row factors.
    No dependence on the size `n` occurs in the number of terms. -/
def fourOuterProducts
    (A00 : Fin 2 → Fin 2 → ℚ)
    (A0h : Fin 2 → Fin n → ℚ)
    (Ah0 : Fin n → Fin 2 → ℚ)
    (C : Fin 2 → Fin n → ℚ)
    (i j : Fin n) : ℚ :=
  (∑ a : Fin 2,
     C a i * ((∑ b : Fin 2, A00 a b * C b j) - A0h a j))
  + (∑ a : Fin 2, (-Ah0 i a) * C a j)

/-- Exact block-compression identity, with arbitrary finite dimension n.
    This does NOT use fin_cases on n. -/
theorem projected_minus_reference_is_four_outer_products
    (A00 : Fin 2 → Fin 2 → ℚ)
    (A0h : Fin 2 → Fin n → ℚ)
    (Ah0 : Fin n → Fin 2 → ℚ)
    (Ahh : Fin n → Fin n → ℚ)
    (C : Fin 2 → Fin n → ℚ)
    (i j : Fin n) :
    projectedBlock A00 A0h Ah0 Ahh C i j - Ahh i j =
      fourOuterProducts A00 A0h Ah0 C i j := by
  simp only [projectedBlock, fourOuterProducts, Fin.sum_univ_two]
  ring

/-- A quadratic expression for finite-support vectors, indexed by all N. -/
def finiteSupportForm (A : ℕ → ℕ → ℚ) (x : ℕ →₀ ℚ) : ℚ :=
  ∑ i ∈ x.support, ∑ j ∈ x.support, x i * A i j * x j

def finiteSupportNormSq (x : ℕ →₀ ℚ) : ℚ :=
  ∑ i ∈ x.support, (x i) ^ 2

/-- Every finite section is nonnegative -> the induced algebraic form is
    nonnegative on c00. The hypothesis MUST cover every finite section.
    Checking three finite sizes does not satisfy it. -/
theorem all_finite_sections_imply_c00_nonnegative
    (A : ℕ → ℕ → ℚ)
    (h : ∀ (s : Finset ℕ) (x : ℕ → ℚ),
      0 ≤ ∑ i ∈ s, ∑ j ∈ s, x i * A i j * x j)
    (x : ℕ →₀ ℚ) :
    0 ≤ finiteSupportForm A x := by
  exact h x.support (fun i => x i)

/-- A uniform quantitative finite-section margin passes to c00.
    Again, the uniform margin is an EXPLICIT assumption, not an inference
    from the P13/P14 finite certificates. -/
theorem uniform_finite_margin_implies_c00_margin
    (A : ℕ → ℕ → ℚ) (γ : ℚ)
    (h : ∀ (s : Finset ℕ) (x : ℕ → ℚ),
       γ * (∑ i ∈ s, (x i)^2) ≤
          ∑ i ∈ s, ∑ j ∈ s, x i * A i j * x j)
    (x : ℕ →₀ ℚ) :
    γ * finiteSupportNormSq x ≤ finiteSupportForm A x := by
  exact h x.support (fun i => x i)

/-- Elementary diagnostic: a finite-rank shift in each finite section
    is not a dimension-free rank bound, since the rank equals n for I_n.
    This lemma records only the strictly positive diagonal coefficient. -/
theorem diagonal_shift_is_nonzero : (0 : ℚ) < 3 / 50 := by norm_num

end Numeria.ES8P19
