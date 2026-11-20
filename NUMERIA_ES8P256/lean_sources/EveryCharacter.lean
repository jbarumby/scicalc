import Development.AbelianVarietyTheory.Proofs.GlobalInterpolation.ShapeOf
import Development.AutomorphicArithmetic.CharacterBoxes.ShapeWin
import Development.AutomorphicArithmetic.CharacterBoxes.SupplyChi
import Development.AutomorphicArithmetic.CharacterBoxes.Imprimitive

open Formalism.BasicNotions Formalism.MainNotions GlobalTheory.Fc

theorem UpperTheory.Top.plugWin (hG : GlobalTheory.GLq.Green) (hV : UpperTheory.Top.VWin)
    (h : Nonempty (Hyps (fun 𝒞 => Closed 𝒞 fun Fm _ => LinkN Fm) DRS)) : Formalism.HC.ZeroFreeAll :=
  h.elim fun h => UpperTheory.Top.zeroFreeAll_of_topG (la := 3 / 2) (by norm_num)
    (GlobalTheory.Fc.shape_of _ _ _ _ (UpperTheory.Top.shapeWin hG) trivial trivial h) fun ψ => UpperTheory.Top.supplyChi ψ (hV ψ)

def UpperTheory.EveryCharacter : Prop :=
  ∀ (N : ℕ) [NeZero N] (χ : DirichletCharacter ℂ N) (s : ℂ), 1 / 2 < s.re → (χ ≠ 1 ∨ s ≠ 1) → DirichletCharacter.LFunction χ s ≠ 0

theorem UpperTheory.Top.plugAll (hG : GlobalTheory.GLq.Green) (hV : UpperTheory.Top.VWin)
    (h : Nonempty (Hyps (fun 𝒞 => Closed 𝒞 fun Fm _ => LinkN Fm) DRS)) : UpperTheory.EveryCharacter :=
  UpperTheory.Top.imprimitive (UpperTheory.Top.plugWin hG hV h)
