import Mathlib.Tactic

/-!
ES8-P11: rational interval, signed-weight and algebraic correction lemmas.
All proofs are intended to be kernel checked.
No analytic digamma/Cauchy identity, finite LDL replay or RH proof asserted.
-/
namespace Numeria.ES8P11

structure RatInterval where
  lo : ℚ
  hi : ℚ

def Contains (I : RatInterval) (x : ℚ) : Prop := I.lo ≤ x ∧ x ≤ I.hi
def plus (I J : RatInterval) : RatInterval := ⟨I.lo + J.lo, I.hi + J.hi⟩
def minus (I J : RatInterval) : RatInterval := ⟨I.lo - J.hi, I.hi - J.lo⟩
def times (w : ℚ) (I : RatInterval) : RatInterval :=
  if 0 ≤ w then ⟨w * I.lo, w * I.hi⟩ else ⟨w * I.hi, w * I.lo⟩

theorem contains_plus {I J : RatInterval} {x y : ℚ}
    (hx : Contains I x) (hy : Contains J y) : Contains (plus I J) (x + y) := by
  unfold Contains plus at *
  constructor <;> linarith [hx.1, hx.2, hy.1, hy.2]

theorem contains_minus {I J : RatInterval} {x y : ℚ}
    (hx : Contains I x) (hy : Contains J y) : Contains (minus I J) (x - y) := by
  unfold Contains minus at *
  constructor <;> linarith [hx.1, hx.2, hy.1, hy.2]

theorem contains_times_nonneg {I : RatInterval} {x w : ℚ}
    (hx : Contains I x) (hw : 0 ≤ w) :
    w * I.lo ≤ w * x ∧ w * x ≤ w * I.hi := by
  exact ⟨mul_le_mul_of_nonneg_left hx.1 hw,
         mul_le_mul_of_nonneg_left hx.2 hw⟩

theorem contains_times_nonpos {I : RatInterval} {x w : ℚ}
    (hx : Contains I x) (hw : w ≤ 0) :
    w * I.hi ≤ w * x ∧ w * x ≤ w * I.lo := by
  constructor
  · exact mul_le_mul_of_nonpos_left hx.2 hw
  · exact mul_le_mul_of_nonpos_left hx.1 hw

theorem contains_times {I : RatInterval} {x w : ℚ}
    (hx : Contains I x) : Contains (times w I) (w * x) := by
  by_cases hw : 0 ≤ w
  · simpa [times, hw, Contains] using contains_times_nonneg hx hw
  · have hn : w ≤ 0 := le_of_lt (lt_of_not_ge hw)
    simpa [times, hw, Contains] using contains_times_nonpos hx hn

theorem weighted_sum_included {ι : Type*} [Fintype ι]
    (w x : ι → ℚ) (I : ι → RatInterval)
    (hx : ∀ i, Contains (I i) (x i)) :
    (∑ i, (times (w i) (I i)).lo) ≤ (∑ i, w i * x i) ∧
    (∑ i, w i * x i) ≤ (∑ i, (times (w i) (I i)).hi) := by
  constructor
  · apply Finset.sum_le_sum
    intro i hi
    exact (contains_times (hx i)).1
  · apply Finset.sum_le_sum
    intro i hi
    exact (contains_times (hx i)).2

theorem enlarge {I : RatInterval} {x eps : ℚ}
    (hx : Contains I x) (heps : 0 ≤ eps) :
    I.lo - eps ≤ x ∧ x ≤ I.hi + eps := by
  constructor <;> linarith [hx.1, hx.2]

theorem midpoint_radius {x mid rad : ℚ}
    (h : |x - mid| ≤ rad) : mid - rad ≤ x ∧ x ≤ mid + rad := by
  rcases (abs_le.mp h) with ⟨hl, hu⟩
  constructor <;> linarith

def oldK : ℚ := 16416784384 / 1601806640625
def integratedK : ℚ := 250582204083791 / 150744971763610839843750
def oldE : ℚ := 4154891126177792 / 209922127267236328125
def newE : ℚ := 31709674433579248304 / 9877816817861795858001708984375

theorem constants_positive :
    0 < oldK ∧ 0 < integratedK ∧ 0 < oldE ∧ 0 < newE := by
  norm_num [oldK, integratedK, oldE, newE]

theorem exact_error_scaling : newE * oldK = oldE * integratedK := by
  norm_num [newE, oldK, oldE, integratedK]

theorem improvement_over_six_million : 6000000 * newE < oldE := by
  norm_num [newE, oldE]

theorem strictly_sharper : newE < oldE := by
  norm_num [newE, oldE]

theorem exact_ratio :
    oldE / newE = 1544972804863488000000 / 250582204083791 := by
  norm_num [oldE, newE]

def bTwo (u b : ℚ) : ℚ :=
  b * (3 * u ^ 2 + b) / (12 * u ^ 2 * (u ^ 2 + b) ^ 2)

def bFour (u b : ℚ) : ℚ :=
  b * (10 * u ^ 6 + 5 * b * u ^ 4 + 4 * b ^ 2 * u ^ 2 + b ^ 3) /
    (120 * u ^ 4 * (u ^ 2 + b) ^ 4)

theorem bTwo_nonneg {u b : ℚ} (hu : 0 < u) (hb : 0 ≤ b) :
    0 ≤ bTwo u b := by
  unfold bTwo
  have h_u : 0 ≤ u ^ 2 := sq_nonneg u
  have h_den : 0 < 12 * u ^ 2 * (u ^ 2 + b) ^ 2 := by positivity
  apply div_nonneg
  · positivity
  · exact le_of_lt h_den

theorem bFour_nonneg {u b : ℚ} (hu : 0 < u) (hb : 0 ≤ b) :
    0 ≤ bFour u b := by
  unfold bFour
  have h_den : 0 < 120 * u ^ 4 * (u ^ 2 + b) ^ 4 := by positivity
  apply div_nonneg
  · positivity
  · exact le_of_lt h_den

theorem disk_geometry :
    (0 : ℚ) < 3/10 ∧ (0 : ℚ) < 7/20 ∧
    (0 : ℚ) < 3/40 ∧ (0 : ℚ) < 137/200 ∧
    7/20 = 3/10 + 1/20 := by norm_num


/-! ## Generic directed fixed-point rounding (including scale 10^90). -/

def floorAt (scale x : ℚ) : ℚ :=
  (Int.floor (x * scale) : ℚ) / scale
def ceilAt (scale x : ℚ) : ℚ :=
  (Int.ceil (x * scale) : ℚ) / scale
def rounded (scale x : ℚ) : RatInterval :=
  ⟨floorAt scale x, ceilAt scale x⟩

theorem floorAt_le {scale : ℚ} (hs : 0 < scale) (x : ℚ) :
    floorAt scale x ≤ x := by
  have hf : (Int.floor (x * scale) : ℚ) ≤ x * scale := Int.floor_le _
  show (Int.floor (x * scale) : ℚ) / scale ≤ x
  rw [div_le_iff₀ hs]
  exact hf

theorem ceilAt_ge {scale : ℚ} (hs : 0 < scale) (x : ℚ) :
    x ≤ ceilAt scale x := by
  have hc : x * scale ≤ (Int.ceil (x * scale) : ℚ) := Int.le_ceil _
  show x ≤ (Int.ceil (x * scale) : ℚ) / scale
  rw [le_div_iff₀ hs]
  exact hc

theorem rounded_contains {scale : ℚ} (hs : 0 < scale) (x : ℚ) :
    Contains (rounded scale x) x := by
  exact ⟨floorAt_le hs x, ceilAt_ge hs x⟩

theorem floorAt_mono {scale : ℚ} (hs : 0 < scale) {x y : ℚ}
    (hxy : x ≤ y) : floorAt scale x ≤ floorAt scale y := by
  have hm : x * scale ≤ y * scale :=
    mul_le_mul_of_nonneg_right hxy (le_of_lt hs)
  have hf : (Int.floor (x * scale) : ℚ) ≤ (Int.floor (y * scale) : ℚ) := by
    exact_mod_cast Int.floor_mono hm
  show (Int.floor (x * scale) : ℚ) / scale ≤
       (Int.floor (y * scale) : ℚ) / scale
  exact div_le_div_of_nonneg_right hf (le_of_lt hs)

theorem ceilAt_mono {scale : ℚ} (hs : 0 < scale) {x y : ℚ}
    (hxy : x ≤ y) : ceilAt scale x ≤ ceilAt scale y := by
  have hm : x * scale ≤ y * scale :=
    mul_le_mul_of_nonneg_right hxy (le_of_lt hs)
  have hc : (Int.ceil (x * scale) : ℚ) ≤ (Int.ceil (y * scale) : ℚ) := by
    exact_mod_cast Int.ceil_mono hm
  show (Int.ceil (x * scale) : ℚ) / scale ≤
       (Int.ceil (y * scale) : ℚ) / scale
  exact div_le_div_of_nonneg_right hc (le_of_lt hs)

theorem scale90_pos : (0 : ℚ) < (10 : ℚ) ^ 90 := by positivity

theorem rounded90_contains (x : ℚ) :
    Contains (rounded ((10 : ℚ) ^ 90) x) x :=
  rounded_contains scale90_pos x

end Numeria.ES8P11
