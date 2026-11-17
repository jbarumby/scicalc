import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
NUMERIA ES8-P232: real-field transfer of the exact rational prefix4 certificate; cloned proof over ℝ. Original P231 is unchanged. for the first four even modes.
Same coefficients as P231_PREFIX4_DATA.json. Generated from P231_PREFIX4_DATA.json (SHA-256 51f6144c096ccc28658f66a940c3bf675e3c053141aa1994f7e3030257fc94f1).
P193 center gzip SHA-256: c105996b75a5782de338dc5b9711b97839f7093f190367bfb67e43cc418b0c65.
P230 full 200-mode finite lower bound: 5363043335327603757336221739842322080321267587773192483503406329186517435022496372140870254782681904878003585055869617350776827783501151810276163439946821247017919939007371493499440837685554440356724880113314567741875041813634432099334796699272215493569434940450633119251278132274032569862149622536695586531477279443858852674426180097211690248449622452035837992988283119360192401720306558705719770292351633542204600960171527657166292586342191671707380179415087053061226519254272313342124345766680318883905046627824418173967039725810713125530490766193404707503963746889283098503702811769802196301209849998935862685764810668665622982049713353640722197109997920154986382171729393136022527950404495423631939722239757488433094454752516827475536402635882443498139607995797316758362792426525697833353470471388665734628458314957944949268434634899648742314935464803093732031911826702076056316678952771288404615876200694096257/536304333532760375733622174016066079868894616642419189056627215906078178518261775569883173378711915438573043572066591941897434240901554700807134628886314835247888949139830047251815745570747434603165223466058353864073929646375958548831095190134486384510983662387157381103078767483068158114190581760770181375234579340099229143573623149235889565371059802098076176404345464328580842350266334178456721074778172040807334621947373299470882944050806547999846588794005962150122539208291622171958374709356866309946578742726572079680115900696831771790678826125489820603917668247611388962826065072830808919702395878024758452594303798195423607903253050300149068414655008216170326640935613540448556545177693443804749852170599324725366791453518964569669460172914392588970527920268547750013790473682944000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.
The rational boxes are validated by the accompanying Python extractor against
P227/P229 and the P230 source hash. This file proves a finite conditional
statement; it does not prove the analytic inclusion of Weil entries, any tail
or global positivity, or RH.
-/
namespace Numeria.ES8P232Real

noncomputable def center4 (i j : Fin 4) : ℝ :=
  match i.val, j.val with
  | 0, 0 => (180997957162753355382313656816260521747 : ℝ) / 2500000000000000000000000000000000000000
  | 0, 1 => (24795717273920126521880074950497363967 : ℝ) / 156250000000000000000000000000000000000
  | 0, 2 => (892188608018018443865950218511544379369 : ℝ) / 5000000000000000000000000000000000000000
  | 0, 3 => (487518916522187238945885447966477997509 : ℝ) / 5000000000000000000000000000000000000000
  | 1, 0 => (24795717273920126521880074950497363967 : ℝ) / 156250000000000000000000000000000000000
  | 1, 1 => (3499621805400713562994264372774225507247 : ℝ) / 10000000000000000000000000000000000000000
  | 1, 2 => (2013756796477639638044753278542819254829 : ℝ) / 5000000000000000000000000000000000000000
  | 1, 3 => (489798784477703620803781150023313708063 : ℝ) / 2000000000000000000000000000000000000000
  | 2, 0 => (892188608018018443865950218511544379369 : ℝ) / 5000000000000000000000000000000000000000
  | 2, 1 => (2013756796477639638044753278542819254829 : ℝ) / 5000000000000000000000000000000000000000
  | 2, 2 => (507482835474212042244003043378820269007 : ℝ) / 1000000000000000000000000000000000000000
  | 2, 3 => (1078123359871745904704511643332422181261 : ℝ) / 2500000000000000000000000000000000000000
  | 3, 0 => (487518916522187238945885447966477997509 : ℝ) / 5000000000000000000000000000000000000000
  | 3, 1 => (489798784477703620803781150023313708063 : ℝ) / 2000000000000000000000000000000000000000
  | 3, 2 => (1078123359871745904704511643332422181261 : ℝ) / 2500000000000000000000000000000000000000
  | 3, 3 => (6929731715611809008325159334768709448893 : ℝ) / 10000000000000000000000000000000000000000
  | _, _ => 0

noncomputable def pivots4 (i : Fin 4) : ℝ :=
  match i.val with
  | 0 => (180997957162753355357313656816260521747 : ℝ) / 2500000000000000000000000000000000000000
  | 1 => (426771134277624333978771584804134827639745616940800256141381179380654784597 : ℝ) / 201108841291948172619237396462511690830000000000000000000000000000000000000000
  | 2 => (16825260112150767959595309229260926805880858922289799342069964342928943718023358433389399518167683341629477748379 : ℝ) / 4267711342776243339787715848041348276397456169408002561413811793806547845970000000000000000000000000000000000000000
  | 3 => (72821767271843988687368008817099593939223417729067092176363639094051171326385689359601462220090776331928040977440326960924005068410704609610757036261 : ℝ) / 42063150280376919898988273073152317014702147305724498355174910857322359295058396083473498795419208354073694370947500000000000000000000000000000000000000
  | _ => 0

noncomputable def upper4 (i j : Fin 4) : ℝ :=
  match i.val, j.val with
  | 0, 0 => (1 : ℝ)
  | 0, 1 => (132243825460907341450027066402652607824 : ℝ) / 60332652387584451785771218938753507249
  | 0, 2 => (99132067557557604873994468723504931041 : ℝ) / 40221768258389634523847479292502338166
  | 0, 3 => (162506305507395746315295149322159332503 : ℝ) / 120665304775168903571542437877507014498
  | 1, 0 => (0 : ℝ)
  | 1, 1 => (1 : ℝ)
  | 1, 2 => (2339236158717583150167105554034241021756279018027322792093811906916664114910 : ℝ) / 426771134277624333978771584804134827639745616940800256141381179380654784597
  | 1, 3 => (6270522002623955756882689703121889194512218824344929977963585809523416474201 : ℝ) / 426771134277624333978771584804134827639745616940800256141381179380654784597
  | 2, 0 => (0 : ℝ)
  | 2, 1 => (0 : ℝ)
  | 2, 2 => (1 : ℝ)
  | 2, 3 => (85499972143170398719049413484272430358723035250972715765060937424790410999889142751448011049242347032544910146127 : ℝ) / 16825260112150767959595309229260926805880858922289799342069964342928943718023358433389399518167683341629477748379
  | 3, 0 => (0 : ℝ)
  | 3, 1 => (0 : ℝ)
  | 3, 2 => (0 : ℝ)
  | 3, 3 => (1 : ℝ)
  | _, _ => 0

noncomputable def shift : ℝ := (1 : ℝ) / 100000000000000000000
noncomputable def epsilon : ℝ := (1 : ℝ) / 10000000000000000000000000000000000000000
noncomputable def rho : ℝ := (1 : ℝ) / 2500000000000000000000000000000000000000

noncomputable def shifted4 (i j : Fin 4) : ℝ :=
  center4 i j - if i = j then shift else 0

noncomputable def err4 (_i _j : Fin 4) : ℝ := epsilon

noncomputable def qForm {n : ℕ} (M : Fin n → Fin n → ℝ) (x : Fin n → ℝ) : ℝ :=
  ∑ i : Fin n, ∑ j : Fin n, x i * M i j * x j

noncomputable def factoredQ {n : ℕ} (d : Fin n → ℝ) (U : Fin n → Fin n → ℝ)
    (x : Fin n → ℝ) : ℝ :=
  ∑ k : Fin n, d k * (∑ j : Fin n, U k j * x j) ^ 2

set_option maxHeartbeats 6000000 in
theorem factor_entries4 : ∀ i j : Fin 4,
    shifted4 i j = ∑ k : Fin 4, pivots4 k * upper4 k i * upper4 k j := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    norm_num [shifted4, center4, shift, pivots4, upper4, Fin.sum_univ_succ]

theorem pivots_pos4 : ∀ i : Fin 4, 0 < pivots4 i := by
  intro i
  fin_cases i <;> norm_num [pivots4]

theorem upper_inj4 (x : Fin 4 → ℝ)
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
  fin_cases i <;> assumption

theorem qForm_eq_factored {n : ℕ} (C : Fin n → Fin n → ℝ)
    (d : Fin n → ℝ) (U : Fin n → Fin n → ℝ)
    (h : ∀ i j : Fin n, C i j = ∑ k : Fin n, d k * U k i * U k j)
    (x : Fin n → ℝ) : qForm C x = factoredQ d U x := by
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

theorem factoredQ_nonneg4 (x : Fin 4 → ℝ) : 0 ≤ factoredQ pivots4 upper4 x := by
  unfold factoredQ
  apply Finset.sum_nonneg
  intro k hk
  exact mul_nonneg (le_of_lt (pivots_pos4 k)) (sq_nonneg _)

theorem qForm_eq_dotProduct (M : Matrix (Fin 4) (Fin 4) ℝ) (x : Fin 4 → ℝ) :
    qForm M x = dotProduct x (M.mulVec x) := by
  classical
  change (∑ i : Fin 4, ∑ j : Fin 4, x i * M i j * x j) =
    ∑ i : Fin 4, x i * (∑ j : Fin 4, M i j * x j)
  calc
    _ = ∑ i : Fin 4, ∑ j : Fin 4, x i * (M i j * x j) := by
      apply Finset.sum_congr rfl
      intro i _
      apply Finset.sum_congr rfl
      intro j _
      ring
    _ = _ := by
      apply Finset.sum_congr rfl
      intro i _
      rw [Finset.mul_sum]

theorem diagonal_shift4 (C : Fin 4 → Fin 4 → ℝ) (x : Fin 4 → ℝ) (s : ℝ) :
    qForm (fun i j => C i j + if i = j then s else 0) x =
      qForm C x + s * (∑ i : Fin 4, x i ^ 2) := by
  unfold qForm
  simp only [mul_add, add_mul, Finset.sum_add_distrib]
  congr 1
  calc
    (∑ i : Fin 4, ∑ j : Fin 4, x i * (if i = j then s else 0) * x j) =
      ∑ i : Fin 4, s * x i ^ 2 := by
        apply Finset.sum_congr rfl
        intro i _
        simp [pow_two, mul_left_comm, mul_comm]
    _ = _ := by rw [Finset.mul_sum]

theorem center_shift_id4 (x : Fin 4 → ℝ) :
    qForm center4 x = qForm shifted4 x + shift * (∑ i : Fin 4, x i ^ 2) := by
  have hid : ∀ i j : Fin 4,
      center4 i j = shifted4 i j + if i = j then shift else 0 := by
    intro i j
    simp [shifted4]
  calc
    qForm center4 x = qForm (fun i j => shifted4 i j + if i = j then shift else 0) x := by
      congr 1
      funext i
      funext j
      exact hid i j
    _ = _ := diagonal_shift4 shifted4 x shift

theorem center_lower4 (x : Fin 4 → ℝ) :
    shift * (∑ i : Fin 4, x i ^ 2) ≤ qForm center4 x := by
  have hf := qForm_eq_factored shifted4 pivots4 upper4 factor_entries4 x
  rw [center_shift_id4, hf]
  have hn := factoredQ_nonneg4 x
  linarith

/-- A bound for one entry error, converted into a quadratic-form bound. -/
theorem error_pair_lower (e t u v : ℝ) (hlo : -e ≤ t) (hhi : t ≤ e) :
    -(e * (u ^ 2 + v ^ 2)) / 2 ≤ u * t * v := by
  have hp : 0 ≤ (e + t) * (u + v) ^ 2 :=
    mul_nonneg (by linarith) (sq_nonneg _)
  have hm : 0 ≤ (e - t) * (u - v) ^ 2 :=
    mul_nonneg (by linarith) (sq_nonneg _)
  nlinarith [hp, hm]

private theorem sum_div_two (f : Fin 4 → ℝ) :
    (∑ i : Fin 4, f i) / 2 = ∑ i : Fin 4, f i / 2 := by
  rw [div_eq_mul_inv, Finset.sum_mul]
  simp only [div_eq_mul_inv]

private theorem double_sum_div_two (f : Fin 4 → Fin 4 → ℝ) :
    (∑ i : Fin 4, ∑ j : Fin 4, f i j) / 2 =
      ∑ i : Fin 4, ∑ j : Fin 4, f i j / 2 := by
  rw [sum_div_two]
  apply Finset.sum_congr rfl
  intro i _
  exact sum_div_two (f i)

private theorem symmetric_envelope (eps : Fin 4 → Fin 4 → ℝ) (x : Fin 4 → ℝ)
    (hsym : ∀ i j, eps i j = eps j i) :
    (∑ i : Fin 4, ∑ j : Fin 4, eps i j * (x i ^ 2 + x j ^ 2) / 2) =
    ∑ i : Fin 4, (∑ j : Fin 4, eps i j) * x i ^ 2 := by
  have hs :
      (∑ i : Fin 4, ∑ j : Fin 4, eps i j * x j ^ 2) =
      ∑ i : Fin 4, ∑ j : Fin 4, eps i j * x i ^ 2 := by
    calc
      _ = ∑ j : Fin 4, ∑ i : Fin 4, eps i j * x j ^ 2 := by rw [Finset.sum_comm]
      _ = ∑ j : Fin 4, ∑ i : Fin 4, eps j i * x j ^ 2 := by
        apply Finset.sum_congr rfl
        intro j _
        apply Finset.sum_congr rfl
        intro i _
        rw [hsym]
      _ = _ := rfl
  calc
    (∑ i : Fin 4, ∑ j : Fin 4, eps i j * (x i ^ 2 + x j ^ 2) / 2) =
      ∑ i : Fin 4, ∑ j : Fin 4,
        (eps i j * x i ^ 2 / 2 + eps i j * x j ^ 2 / 2) := by
          apply Finset.sum_congr rfl
          intro i _
          apply Finset.sum_congr rfl
          intro j _
          ring
    _ = (∑ i : Fin 4, ∑ j : Fin 4, eps i j * x i ^ 2 / 2) +
        (∑ i : Fin 4, ∑ j : Fin 4, eps i j * x j ^ 2 / 2) := by
          simp only [Finset.sum_add_distrib]
    _ = ((∑ i : Fin 4, ∑ j : Fin 4, eps i j * x i ^ 2) +
        (∑ i : Fin 4, ∑ j : Fin 4, eps i j * x j ^ 2)) / 2 := by
          rw [← double_sum_div_two, ← double_sum_div_two]
          ring
    _ = ∑ i : Fin 4, ∑ j : Fin 4, eps i j * x i ^ 2 := by rw [hs]; ring
    _ = _ := by simp only [Finset.sum_mul]

theorem qForm_row_error_lower (M C eps : Fin 4 → Fin 4 → ℝ)
    (x : Fin 4 → ℝ) (delta : ℝ)
    (hsym : ∀ i j, eps i j = eps j i)
    (hentry : ∀ i j, -eps i j ≤ M i j - C i j ∧ M i j - C i j ≤ eps i j)
    (hrow : ∀ i, (∑ j : Fin 4, eps i j) ≤ delta) :
    qForm C x - delta * (∑ i : Fin 4, x i ^ 2) ≤ qForm M x := by
  have hb :
      -(∑ i : Fin 4, ∑ j : Fin 4, eps i j * (x i ^ 2 + x j ^ 2) / 2) ≤
      ∑ i : Fin 4, ∑ j : Fin 4, x i * (M i j - C i j) * x j := by
    simp only [← Finset.sum_neg_distrib]
    apply Finset.sum_le_sum
    intro i _
    apply Finset.sum_le_sum
    intro j _
    simpa only [neg_div] using
      (error_pair_lower (eps i j) (M i j - C i j) (x i) (x j)
        (hentry i j).1 (hentry i j).2)
  have hr :
      (∑ i : Fin 4, (∑ j : Fin 4, eps i j) * x i ^ 2) ≤
      delta * (∑ i : Fin 4, x i ^ 2) := by
    calc
      _ ≤ ∑ i : Fin 4, delta * x i ^ 2 := by
        apply Finset.sum_le_sum
        intro i _
        exact mul_le_mul_of_nonneg_right (hrow i) (sq_nonneg _)
      _ = _ := by rw [Finset.mul_sum]
  have hd :
      (∑ i : Fin 4, ∑ j : Fin 4, x i * (M i j - C i j) * x j) =
      qForm M x - qForm C x := by
    simp only [qForm, mul_sub, sub_mul, Finset.sum_sub_distrib]
  rw [symmetric_envelope eps x hsym] at hb
  rw [hd] at hb
  linarith

/-- This is P120's finite conditional bridge: hypotheses must be proved elsewhere. -/
theorem matrix_bridge (W C E : Matrix (Fin 4) (Fin 4) ℝ)
    (hW_hermitian : W.IsHermitian)
    (h_decomp : W = C + E)
    (h_margin : rho < shift)
    (h_center : ∀ x : Fin 4 → ℝ, x ≠ 0 →
      shift * (∑ i, x i ^ 2) ≤ dotProduct x (C.mulVec x))
    (h_error : ∀ x : Fin 4 → ℝ, x ≠ 0 →
      -(rho * (∑ i, x i ^ 2)) ≤ dotProduct x (E.mulVec x)) : W.PosDef := by
  classical
  apply Matrix.PosDef.of_dotProduct_mulVec_pos hW_hermitian
  intro x hx
  have hq : 0 < ∑ i, x i ^ 2 := by
    by_contra h
    have hz : ∀ i, x i = 0 := by
      intro i
      by_contra hi
      have hp : 0 < x i ^ 2 := sq_pos_of_ne_zero hi
      have hle : x i ^ 2 ≤ ∑ j, x j ^ 2 :=
        Finset.single_le_sum (fun j _ => sq_nonneg (x j)) (Finset.mem_univ i)
      exact h (lt_of_lt_of_le hp hle)
    exact hx (funext hz)
  have hsplit : dotProduct x (W.mulVec x) =
      dotProduct x (C.mulVec x) + dotProduct x (E.mulVec x) := by
    rw [h_decomp, Matrix.add_mulVec, dotProduct_add]
  have hpos : 0 < dotProduct x (W.mulVec x) := by
    rw [hsplit]
    have hm := h_center x hx
    have he := h_error x hx
    have hmargin : 0 < (shift - rho) * (∑ i, x i ^ 2) :=
      mul_pos (sub_pos.mpr h_margin) hq
    have hsum : (shift - rho) * (∑ i, x i ^ 2) ≤
        dotProduct x (C.mulVec x) + dotProduct x (E.mulVec x) := by
      nlinarith [hm, he]
    exact lt_of_lt_of_le hmargin hsum
  simpa using hpos

/-- Every rational matrix in the P227 prefix box is positive definite,
provided it is Hermitian. P227/P229-to-Weil analytic inclusion is external. -/
theorem P232_box_positive_real (W : Matrix (Fin 4) (Fin 4) ℝ)
    (hW_hermitian : W.IsHermitian)
    (hbox : ∀ i j : Fin 4, |W i j - center4 i j| ≤ epsilon) : W.PosDef := by
  let C : Matrix (Fin 4) (Fin 4) ℝ := center4
  let E : Matrix (Fin 4) (Fin 4) ℝ := W - C
  have hdecomp : W = C + E := by
    ext i j
    change W i j = center4 i j + (W i j - center4 i j)
    ring
  have hentry : ∀ i j : Fin 4,
      -err4 i j ≤ W i j - center4 i j ∧ W i j - center4 i j ≤ err4 i j := by
    intro i j
    have h := abs_le.mp (hbox i j)
    simpa [err4, epsilon] using h
  have hrows : ∀ i : Fin 4, (∑ j : Fin 4, err4 i j) ≤ rho := by
    intro i
    norm_num [err4, epsilon, rho, Fin.sum_univ_succ]
  have hcenter : ∀ x : Fin 4 → ℝ, x ≠ 0 →
      shift * (∑ i, x i ^ 2) ≤ dotProduct x (C.mulVec x) := by
    intro x hx
    rw [← qForm_eq_dotProduct C x]
    simpa [C] using center_lower4 x
  have herror : ∀ x : Fin 4 → ℝ, x ≠ 0 →
      -(rho * (∑ i, x i ^ 2)) ≤ dotProduct x (E.mulVec x) := by
    intro x hx
    have hpert := qForm_row_error_lower W center4 err4 x rho
      (by intro i j; rfl) hentry hrows
    have hdiff : qForm E x = qForm W x - qForm center4 x := by
      change qForm (fun i j => W i j - center4 i j) x = _
      simp only [qForm, mul_sub, sub_mul, Finset.sum_sub_distrib]
    rw [← qForm_eq_dotProduct E x, hdiff]
    nlinarith [hpert]
  have hm : rho < shift := by norm_num [rho, shift]
  exact matrix_bridge W C E hW_hermitian hdecomp hm hcenter herror

end Numeria.ES8P232Real
