import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Tactic

/-!
NUMERIA P120 — Ponte condicional entre uma cota quadrática e positividade.
Não importa os certificados P111/P115 nem prova a inclusão Arb.
O enunciado exige explicitamente a simetria hermitiana e as duas cotas.
-/

namespace Numeria.ES8P120

/-- A soma dos quadrados de um vetor racional não nulo é positiva. -/
theorem normSq_pos_of_ne_zero {ι : Type*} [Fintype ι]
    (x : ι → ℚ) (hx : x ≠ 0) :
    0 < ∑ i, x i ^ 2 := by
  classical
  by_contra h
  have hx0 : ∀ i, x i = 0 := by
    intro i
    by_contra hi
    have hp : 0 < x i ^ 2 := sq_pos_of_ne_zero hi
    have hle : x i ^ 2 ≤ ∑ j, x j ^ 2 :=
      Finset.single_le_sum (fun j _ => sq_nonneg (x j)) (Finset.mem_univ i)
    exact h (lt_of_lt_of_le hp hle)
  apply hx
  funext i
  exact hx0 i

/-- A ponte escalar é uma desigualdade elementar, não uma nova hipótese analítica. -/
theorem scalar_bridge (q center err shift rho : ℚ)
    (hq : 0 < q)
    (hc : shift * q ≤ center)
    (he : -(rho * q) ≤ err)
    (hm : rho < shift) :
    0 < center + err := by
  have hmargin : 0 < (shift - rho) * q :=
    mul_pos (sub_pos.mpr hm) hq
  have hsum : (shift - rho) * q ≤ center + err := by
    nlinarith [hc, he]
  exact lt_of_lt_of_le hmargin hsum

/--
Prova FINITA CONDICIONAL: a matriz completa é positiva se a forma central
superar shift*||x||² e o erro não cair abaixo de -rho*||x||².

Este teorema não demonstra h_center, h_error nem hW_hermitian a partir das
bolas Arb; essas obrigações devem ser demonstradas separadamente.
-/
theorem matrix_bridge {ι : Type*} [Fintype ι] [DecidableEq ι]
    (W C E : Matrix ι ι ℚ)
    (hW_hermitian : W.IsHermitian)
    (h_decomp : W = C + E)
    (shift rho : ℚ)
    (h_margin : rho < shift)
    (h_center : ∀ x : ι → ℚ, x ≠ 0 →
      shift * (∑ i, x i ^ 2) ≤ Matrix.dotProduct x (C.mulVec x))
    (h_error : ∀ x : ι → ℚ, x ≠ 0 →
      -(rho * (∑ i, x i ^ 2)) ≤ Matrix.dotProduct x (E.mulVec x)) :
    W.PosDef := by
  classical
  apply Matrix.PosDef.of_dotProduct_mulVec_pos hW_hermitian
  intro x hx
  have hq : 0 < ∑ i, x i ^ 2 := normSq_pos_of_ne_zero x hx
  have hsplit : Matrix.dotProduct x (W.mulVec x) =
      Matrix.dotProduct x (C.mulVec x) + Matrix.dotProduct x (E.mulVec x) := by
    rw [h_decomp, Matrix.add_mulVec, Matrix.dotProduct_add]
  have hp : 0 < Matrix.dotProduct x (W.mulVec x) := by
    rw [hsplit]
    exact scalar_bridge _ _ _ shift rho hq (h_center x hx) (h_error x hx) h_margin
  simpa using hp

end Numeria.ES8P120
