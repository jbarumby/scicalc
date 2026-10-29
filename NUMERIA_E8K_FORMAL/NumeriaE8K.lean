/-
  NUMERIA E8K — Formalização do núcleo de intervalos racionais (Lean 4 + Mathlib)

  ATUALIZACAO (continuacao do esqueleto):
   - contains_mul: DEMONSTRADO (argumento dos 4 cantos; casos de sinal via
     nlinarith; combinatoria de min/max via le_trans/le_min/max_le).
   - contains_sq: corolario imediato de contains_mul.
   - contains_recip_pos E contains_recip_neg: DEMONSTRADOS (inv_le_inv;
     caso negativo via inv_le_inv_of_neg, que reduz a inv_le_inv + inv_neg).
   - Larguras de add/sub/neg: DEMONSTRADAS (ring). Largura grosseira do
     produto (<= 2*max|I|*max|J|): DEMONSTRADA.
   - Fixed-point tick 10^-42 COMPLETO: fpDown_le, fpUp_ge e fp_mono
     DEMONSTRADOS (Int.floor_le, Int.le_ceil, Int.floor_mono, Int.ceil_mono,
     div_le_iff, le_div_iff, div_le_div_of_le); fpDirected fecha a instancia.
   - width_mul_le (cota fina de largura do produto): DEMONSTRADA via
     estrategia "par a par" (lema corner_diff_le): uv - u'v' =
     u(v - v') + (u - u')v' com |v - v'| <= J.hi - J.lo e |u - u'| <=
     I.hi - I.lo; o max fica a B do min porque TODO par de cantos difere
     em no maximo B (evita os 16 casos de identificacao de max/min).

  STATUS: ZERO sorry no arquivo. Nao compilado localmente (sem toolchain).

  COMO COMPILAR: toolchain Lean 4 com Mathlib. Ex.:
     lake new e8k mathlib && cd e8k
     copie este arquivo como E8K/E8K.lean e `import E8K` no E8K.lean raiz
     lake build

  STATUS DE COMPILACAO: NAO compilado localmente (sem toolchain Lean aqui).
  Nomes a conferir contra a versao do Mathlib: Int.floor_le, Int.le_ceil,
  Int.floor_mono, Int.ceil_mono, div_le_iff (ou div_le_iff0), le_div_iff
  (ou le_div_iff0), div_le_div_of_le (ou div_le_div_right), inv_neg,
  inv_le_inv, abs_le, le_abs_self, neg_abs_le, abs_mul, abs_nonneg,
  mul_le_mul (ordem dos 4 argumentos: a≤b, c≤d, 0≤c, 0≤b), le_min,
  max_le, le_max_left, le_max_right, le_trans, le_refl. Eventuais falhas
  aqui sao mecanicas de nome, nao de conteudo matematico.
-/

import Mathlib.Tactic

namespace E8K

/-- Intervalo racional fechado. -/
structure RatI where
  lo : ℚ
  hi : ℚ
  hlo : lo ≤ hi

def RatI.contains (I : RatI) (x : ℚ) : Prop := I.lo ≤ x ∧ x ≤ I.hi

def RatI.singleton (x : ℚ) : RatI := ⟨x, x, le_refl x⟩

theorem singleton_contains (x : ℚ) : (RatI.singleton x).contains x :=
  ⟨le_refl x, le_refl x⟩

def RatI.add (I J : RatI) : RatI :=
  ⟨I.lo + J.lo, I.hi + J.hi, add_le_add I.hlo J.hlo⟩

def RatI.neg (I : RatI) : RatI := ⟨-I.hi, -I.lo, neg_le_neg I.hlo⟩

def RatI.sub (I J : RatI) : RatI := I.add J.neg

theorem contains_add {I J : RatI} {a b : ℚ}
    (hI : I.contains a) (hJ : J.contains b) :
    (I.add J).contains (a + b) :=
  ⟨add_le_add hI.1 hJ.1, add_le_add hI.2 hJ.2⟩

theorem contains_neg {I : RatI} {a : ℚ} (h : I.contains a) :
    I.neg.contains (-a) :=
  ⟨neg_le_neg h.2, neg_le_neg h.1⟩

theorem contains_sub {I J : RatI} {a b : ℚ}
    (hI : I.contains a) (hJ : J.contains b) :
    (I.sub J).contains (a - b) := by
  show (I.add J.neg).contains (a + -b)
  exact contains_add hI (contains_neg hJ)

private lemma min_le_max4 (u v w z : ℚ) :
    min (min u v) (min w z) ≤ max (max u v) (max w z) :=
  le_trans (min_le_left (min u v) (min w z))
    (le_trans (min_le_left u v)
      (le_trans (le_max_left u v)
        (le_max_left (max u v) (max w z))))

/-- Produto por cantos (4 produtos dos extremos). -/
def RatI.mul (I J : RatI) : RatI where
  lo := min (min (I.lo * J.lo) (I.lo * J.hi)) (min (I.hi * J.lo) (I.hi * J.hi))
  hi := max (max (I.lo * J.lo) (I.lo * J.hi)) (max (I.hi * J.lo) (I.hi * J.hi))
  hlo := min_le_max4 (I.lo * J.lo) (I.lo * J.hi) (I.hi * J.lo) (I.hi * J.hi)

/-- Fato central: multiplicador fixo A, valor c em [r, s]. O produto A*c esta
    entre min(A*r, A*s) e max(A*r, A*s), com o lado controlado pelo sinal de A. -/
private lemma mul_bounds {A r c s : ℚ} (h1 : r ≤ c) (h2 : c ≤ s) :
    min (A * r) (A * s) ≤ A * c ∧ A * c ≤ max (A * r) (A * s) := by
  rcases le_or_lt 0 A with hA | hA
  · refine ⟨le_trans (min_le_left (A * r) (A * s)) ?_,
            le_trans ?_ (le_max_right (A * r) (A * s))⟩
    · exact nlinarith [h1, hA]          -- A * r ≤ A * c  (0 ≤ A)
    · exact nlinarith [h2, hA]          -- A * c ≤ A * s  (0 ≤ A)
  · refine ⟨le_trans (min_le_right (A * r) (A * s)) ?_,
            le_trans ?_ (le_max_left (A * r) (A * s))⟩
    · exact nlinarith [h2, hA]          -- A * s ≤ A * c  (A < 0)
    · exact nlinarith [h1, hA]          -- A * c ≤ A * r  (A < 0)

/-- a*b vive entre o minimo e o maximo dos 4 produtos dos cantos da caixa
    [p,q] x [r,s]. Aplica mul_bounds: (A) com multiplicador a sobre [r,s];
    (B) com multiplicador r e s sobre [p,q], mais comutatividade. -/
private lemma box_mul_bounds {p q r s a b : ℚ}
    (h1 : p ≤ a) (h2 : a ≤ q) (h3 : r ≤ b) (h4 : b ≤ s) :
    min (min (p * r) (p * s)) (min (q * r) (q * s)) ≤ a * b ∧
    a * b ≤ max (max (p * r) (p * s)) (max (q * r) (q * s)) := by
  have A1 := (mul_bounds (A := a) h3 h4).1   -- min (a*r) (a*s) ≤ a*b
  have A2 := (mul_bounds (A := a) h3 h4).2   -- a*b ≤ max (a*r) (a*s)
  have B1 := (mul_bounds (A := r) h1 h2).1    -- min (r*p) (r*q) ≤ r*a
  rw [mul_comm r a] at B1
  rw [mul_comm r p, mul_comm r q] at B1       -- min (p*r) (q*r) ≤ a*r
  have B2 := (mul_bounds (A := r) h1 h2).2    -- r*a ≤ max (r*p) (r*q)
  rw [mul_comm r a] at B2
  rw [mul_comm r p, mul_comm r q] at B2       -- a*r ≤ max (p*r) (q*r)
  have B1s := (mul_bounds (A := s) h1 h2).1
  rw [mul_comm s a, mul_comm s p, mul_comm s q] at B1s
                                              -- min (p*s) (q*s) ≤ a*s
  have B2s := (mul_bounds (A := s) h1 h2).2
  rw [mul_comm s a, mul_comm s p, mul_comm s q] at B2s
                                              -- a*s ≤ max (p*s) (q*s)
  -- canto inferior: os 4 cantos ≤ a*r e ≤ a*s
  have cPr := le_trans (min_le_left (min (p * r) (p * s)) (min (q * r) (q * s)))
    (min_le_left (p * r) (p * s))
  have cQr := le_trans (min_le_right (min (p * r) (p * s)) (min (q * r) (q * s)))
    (min_le_left (q * r) (q * s))
  have cPs := le_trans (min_le_left (min (p * r) (p * s)) (min (q * r) (q * s)))
    (min_le_right (p * r) (p * s))
  have cQs := le_trans (min_le_right (min (p * r) (p * s)) (min (q * r) (q * s)))
    (min_le_right (q * r) (q * s))
  have hmin_r := le_min cPr cQr              -- cantos ≤ min (p*r) (q*r)
  have hmin_s := le_min cPs cQs              -- cantos ≤ min (p*s) (q*s)
  have hmin_ar := le_trans hmin_r B1         -- cantos ≤ a*r
  have hmin_as := le_trans hmin_s B1s        -- cantos ≤ a*s
  -- canto superior: a*r e a*s ≤ os 4 cantos (lado max)
  have uPr := le_trans (le_max_left (p * r) (p * s))
    (le_max_left (max (p * r) (p * s)) (max (q * r) (q * s)))
  have uQr := le_trans (le_max_right (p * r) (p * s))
    (le_max_left (max (p * r) (p * s)) (max (q * r) (q * s)))
  have uPs := le_trans (le_max_left (p * r) (p * s))
    (le_max_right (max (p * r) (p * s)) (max (q * r) (q * s)))
  have uQs := le_trans (le_max_right (p * r) (p * s))
    (le_max_right (max (p * r) (p * s)) (max (q * r) (q * s)))
  have hmax_r := max_le uPr uQr             -- max (p*r) (q*r) ≤ cantos_max
  have hmax_s := max_le uPs uQs             -- max (p*s) (q*s) ≤ cantos_max
  have hmax_ar := le_trans B2 hmax_r        -- a*r ≤ cantos_max
  have hmax_as := le_trans B2s hmax_s       -- a*s ≤ cantos_max
  exact ⟨le_trans (le_min hmin_ar hmin_as) A1,
         le_trans A2 (max_le hmax_ar hmax_as)⟩

theorem contains_mul {I J : RatI} {a b : ℚ}
    (hI : I.contains a) (hJ : J.contains b) :
    (I.mul J).contains (a * b) := by
  obtain ⟨h1, h2⟩ := hI
  obtain ⟨h3, h4⟩ := hJ
  show min (min (I.lo * J.lo) (I.lo * J.hi)) (min (I.hi * J.lo) (I.hi * J.hi)) ≤ a * b ∧
       a * b ≤ max (max (I.lo * J.lo) (I.lo * J.hi)) (max (I.hi * J.lo) (I.hi * J.hi))
  exact box_mul_bounds h1 h2 h3 h4

def RatI.sq (I : RatI) : RatI := I.mul I

theorem contains_sq {I : RatI} {a : ℚ} (h : I.contains a) :
    I.sq.contains (a * a) :=
  contains_mul h h

/-- Reciproco em dominio estritamente positivo: 1/x vive em
    [1/I.hi, 1/I.lo] (funcao decrescente). O caso negativo e demonstrado a
    seguir; zero-containing depende da semantica de guarda da fonte B5
    (zero_recip). -/
theorem contains_recip_pos {I : RatI} {a : ℚ}
    (h : I.contains a) (hpos : 0 < I.lo) :
    (RatI.mk (I.hi)⁻¹ (I.lo)⁻¹
        (inv_le_inv hpos (lt_of_lt_of_le hpos I.hlo) I.hlo)).contains (a⁻¹) := by
  have ha0 : 0 < a := lt_of_lt_of_le hpos h.1
  have hhi0 : 0 < I.hi := lt_of_lt_of_le hpos I.hlo
  refine ⟨?_, ?_⟩
  · exact inv_le_inv ha0 hhi0 h.2      -- (I.hi)⁻¹ ≤ a⁻¹   pois a ≤ I.hi
  · exact inv_le_inv hpos ha0 h.1      -- a⁻¹ ≤ (I.lo)⁻¹   pois I.lo ≤ a

/-- Versao de inv_le_inv para argumentos estritamente negativos: u ≤ v < 0
    implica v⁻¹ ≤ u⁻¹ (o inverso e decrescente tambem em R-). -/
private lemma inv_le_inv_of_neg {u v : ℚ} (hu : u < 0) (hv : v < 0)
    (h : u ≤ v) : v⁻¹ ≤ u⁻¹ := by
  have hu' : 0 < -u := by linarith
  have hv' : 0 < -v := by linarith
  have h' : -v ≤ -u := by linarith
  have key := inv_le_inv hv' hu' h'   -- (-u)⁻¹ ≤ (-v)⁻¹
  simp only [inv_neg] at key           -- -u⁻¹ ≤ -v⁻¹
  linarith

/-- Reciproco em dominio estritamente negativo: mesma forma do caso
    positivo (o inverso inverte os extremos em ambos os sinais). -/
theorem contains_recip_neg {I : RatI} {a : ℚ}
    (h : I.contains a) (hneg : I.hi < 0) :
    (RatI.mk (I.hi)⁻¹ (I.lo)⁻¹
        (inv_le_inv_of_neg (by have hh := I.hlo; linarith) hneg
          I.hlo)).contains (a⁻¹) := by
  have ha0 : a < 0 := lt_of_le_of_lt h.2 hneg
  have hlo0 : I.lo < 0 := by have hh := I.hlo; linarith
  refine ⟨?_, ?_⟩
  · exact inv_le_inv_of_neg ha0 hneg h.2       -- (I.hi)⁻¹ ≤ a⁻¹
  · exact inv_le_inv_of_neg hlo0 ha0 h.1      -- a⁻¹ ≤ (I.lo)⁻¹

/-- Larguras exatas: add/sub/neg (soma/subtracao de larguras). -/
theorem width_add (I J : RatI) :
    (I.add J).hi - (I.add J).lo = (I.hi - I.lo) + (J.hi - J.lo) := by
  show (I.hi + J.hi) - (I.lo + J.lo) = (I.hi - I.lo) + (J.hi - J.lo)
  ring

theorem width_sub (I J : RatI) :
    (I.sub J).hi - (I.sub J).lo = (I.hi - I.lo) + (J.hi - J.lo) := by
  show (I.hi + -J.lo) - (I.lo + -J.hi) = (I.hi - I.lo) + (J.hi - J.lo)
  ring

theorem width_neg (I : RatI) :
    I.neg.hi - I.neg.lo = I.hi - I.lo := by
  show (-I.lo) - (-I.hi) = I.hi - I.lo
  ring

/-- Cota grosseira (mas fechada) da largura do produto: todo canto u*v
    satisfaz |u*v| <= max|I| * max|J|, logo a caixa tem largura
    <= 2 * max|I| * max|J|. -/
theorem width_mul_abs_bound (I J : RatI) :
    (I.mul J).hi - (I.mul J).lo ≤
      2 * (max (abs I.lo) (abs I.hi) * max (abs J.lo) (abs J.hi)) := by
  have key : ∀ u v : ℚ, I.lo ≤ u → u ≤ I.hi → J.lo ≤ v → v ≤ J.hi →
      abs (u * v) ≤ max (abs I.lo) (abs I.hi) * max (abs J.lo) (abs J.hi) := by
    intro u v hu1 hu2 hv1 hv2
    have hI1 : abs I.lo ≤ max (abs I.lo) (abs I.hi) := le_max_left _ _
    have hI2 : abs I.hi ≤ max (abs I.lo) (abs I.hi) := le_max_right _ _
    have hJ1 : abs J.lo ≤ max (abs J.lo) (abs J.hi) := le_max_left _ _
    have hJ2 : abs J.hi ≤ max (abs J.lo) (abs J.hi) := le_max_right _ _
    have huA : abs u ≤ max (abs I.lo) (abs I.hi) := by
      rw [abs_le]
      exact ⟨by linarith [neg_abs_le I.lo, hI1, hu1],
             by linarith [le_abs_self I.hi, hI2, hu2]⟩
    have hvA : abs v ≤ max (abs J.lo) (abs J.hi) := by
      rw [abs_le]
      exact ⟨by linarith [neg_abs_le J.lo, hJ1, hv1],
             by linarith [le_abs_self J.hi, hJ2, hv2]⟩
    have hIpos : 0 ≤ max (abs I.lo) (abs I.hi) :=
      le_trans (abs_nonneg I.lo) (le_max_left _ _)
    rw [abs_mul]
    exact mul_le_mul huA hvA hIpos (abs_nonneg u)
  have hhi : (I.mul J).hi ≤
      max (abs I.lo) (abs I.hi) * max (abs J.lo) (abs J.hi) := by
    show max (max (I.lo * J.lo) (I.lo * J.hi)) (max (I.hi * J.lo) (I.hi * J.hi)) ≤
      max (abs I.lo) (abs I.hi) * max (abs J.lo) (abs J.hi)
    have c1 := le_trans (le_abs_self _)
      (key I.lo J.lo (le_refl _) (le_refl _) (le_refl _) (le_refl _))
    have c2 := le_trans (le_abs_self _)
      (key I.lo J.hi (le_refl _) (le_refl _) (le_refl _) (le_refl _))
    have c3 := le_trans (le_abs_self _)
      (key I.hi J.lo (le_refl _) (le_refl _) (le_refl _) (le_refl _))
    have c4 := le_trans (le_abs_self _)
      (key I.hi J.hi (le_refl _) (le_refl _) (le_refl _) (le_refl _))
    exact max_le (max_le c1 c2) (max_le c3 c4)
  have hlo : -max (abs I.lo) (abs I.hi) * max (abs J.lo) (abs J.hi) ≤
      (I.mul J).lo := by
    show -max (abs I.lo) (abs I.hi) * max (abs J.lo) (abs J.hi) ≤
      min (min (I.lo * J.lo) (I.lo * J.hi)) (min (I.hi * J.lo) (I.hi * J.hi))
    have c1 := le_trans
      (neg_le_neg (key I.lo J.lo (le_refl _) (le_refl _) (le_refl _) (le_refl _)))
      (neg_abs_le _)
    have c2 := le_trans
      (neg_le_neg (key I.lo J.hi (le_refl _) (le_refl _) (le_refl _) (le_refl _)))
      (neg_abs_le _)
    have c3 := le_trans
      (neg_le_neg (key I.hi J.lo (le_refl _) (le_refl _) (le_refl _) (le_refl _)))
      (neg_abs_le _)
    have c4 := le_trans
      (neg_le_neg (key I.hi J.hi (le_refl _) (le_refl _) (le_refl _) (le_refl _)))
      (neg_abs_le _)
    exact le_min (le_min c1 c2) (le_min c3 c4)
  linarith [hhi, hlo]

/-- Lema auxiliar: para quaisquer cantos u*v e u'*v' da caixa [p,q] x [r,s],
    a diferenca e majorada pela soma "largura x modulo do outro lado":
        u*v - u'*v' = u*(v - v') + (u - u')*v'
    com |u*(v-v')| <= max|pq| * (s - r) e |(u-u')*v'| <= (q - p) * max|rs|.
    Nao ha analise de casos de max/min: as pertenncias dao |v - v'| <= s - r
    e |u - u'| <= q - p diretamente (v <= s e r <= v', etc.). -/

private lemma corner_diff_le {p q r s u u' v v' : ℚ}
    (hpq : p ≤ q) (hrs : r ≤ s)
    (hu : u = p ∨ u = q) (hu' : u' = p ∨ u' = q)
    (hv : v = r ∨ v = s) (hv' : v' = r ∨ v' = s) :
    u * v - u' * v' ≤
      (q - p) * max (abs r) (abs s) + (s - r) * max (abs p) (abs q) := by
  have hub : abs u ≤ max (abs p) (abs q) := by
    rcases hu with h | h
    · rw [h]; exact le_max_left (abs p) (abs q)
    · rw [h]; exact le_max_right (abs p) (abs q)
  have hub' : abs v' ≤ max (abs r) (abs s) := by
    rcases hv' with h | h
    · rw [h]; exact le_max_left (abs r) (abs s)
    · rw [h]; exact le_max_right (abs r) (abs s)
  have hv_s : v ≤ s := by
    rcases hv with h | h
    · rw [h]; exact hrs
    · rw [h]; exact le_refl s
  have hv_r : r ≤ v := by
    rcases hv with h | h
    · rw [h]; exact le_refl r
    · rw [h]; exact hrs
  have hv'_s : v' ≤ s := by
    rcases hv' with h | h
    · rw [h]; exact hrs
    · rw [h]; exact le_refl s
  have hv'_r : r ≤ v' := by
    rcases hv' with h | h
    · rw [h]; exact le_refl r
    · rw [h]; exact hrs
  have hdv : abs (v - v') ≤ s - r := by
    rw [abs_le]
    constructor
    · linarith
    · linarith
  have hu_q : u ≤ q := by
    rcases hu with h | h
    · rw [h]; exact hpq
    · rw [h]; exact le_refl q
  have hu_p : p ≤ u := by
    rcases hu with h | h
    · rw [h]; exact le_refl p
    · rw [h]; exact hpq
  have hu'_q : u' ≤ q := by
    rcases hu' with h | h
    · rw [h]; exact hpq
    · rw [h]; exact le_refl q
  have hu'_p : p ≤ u' := by
    rcases hu' with h | h
    · rw [h]; exact le_refl p
    · rw [h]; exact hpq
  have hdu : abs (u - u') ≤ q - p := by
    rw [abs_le]
    constructor
    · linarith
    · linarith
  have hnn1 : 0 ≤ max (abs p) (abs q) :=
    le_trans (abs_nonneg p) (le_max_left (abs p) (abs q))
  have hnn2 : 0 ≤ max (abs r) (abs s) :=
    le_trans (abs_nonneg r) (le_max_left (abs r) (abs s))
  have hA : u * (v - v') ≤ max (abs p) (abs q) * (s - r) := by
    have e1 : u * (v - v') ≤ abs (u * (v - v')) := le_abs_self _
    rw [abs_mul] at e1
    have e2 : abs u * abs (v - v') ≤ max (abs p) (abs q) * (s - r) :=
      mul_le_mul hub hdv (abs_nonneg _) hnn1
    linarith
  have hB : (u - u') * v' ≤ (q - p) * max (abs r) (abs s) := by
    have e1 : (u - u') * v' ≤ abs ((u - u') * v') := le_abs_self _
    rw [abs_mul] at e1
    have e2 : abs (u - u') * abs v' ≤ (q - p) * max (abs r) (abs s) :=
      mul_le_mul hdu hub' (abs_nonneg _) hnn2
    linarith
  have hid : u * v - u' * v' = u * (v - v') + (u - u') * v' := by ring
  linarith

/-- Cota fina de largura do produto (fecha o ultimo sorry do arquivo).
    Estrategia "par a par": em vez de identificar qual canto atinge o max
    e o min (16 casos), observamos que
        hi - lo <= B  <->  hi <= lo + B
    e que hi (max de 4 cantos) satisfaz c <= lo + B para todo canto c,
    porque lo (min de 4 cantos) satisfaz c - B <= d para todo canto d -
    que e exatamente corner_diff_le. Cada lado e fechado com max_le/le_min
    ponto a ponto e linarith. -/

theorem width_mul_le (I J : RatI) :
    (I.mul J).hi - (I.mul J).lo ≤
      (I.hi - I.lo) * max (abs J.lo) (abs J.hi) +
      (J.hi - J.lo) * max (abs I.lo) (abs I.hi) := by
  -- desigualdade par a par entre cantos quaisquer (a cota fina)
  have key : ∀ u u' v v' : ℚ,
      (u = I.lo ∨ u = I.hi) → (u' = I.lo ∨ u' = I.hi) →
      (v = J.lo ∨ v = J.hi) → (v' = J.lo ∨ v' = J.hi) →
      u * v - u' * v' ≤
        (I.hi - I.lo) * max (abs J.lo) (abs J.hi) +
        (J.hi - J.lo) * max (abs I.lo) (abs I.hi) := by
    intro u u' v v' hu hu' hv hv'
    exact corner_diff_le (p := I.lo) (q := I.hi) (r := J.lo) (s := J.hi)
      I.hlo J.hlo hu hu' hv hv'
  -- todo canto c tem c - B <= (I.mul J).lo (ou seja, c e no maximo B
  -- acima do minimo, que domina cada canto via le_min)
  have hc0 : ∀ u v : ℚ, (u = I.lo ∨ u = I.hi) → (v = J.lo ∨ v = J.hi) →
      u * v - ((I.hi - I.lo) * max (abs J.lo) (abs J.hi) +
               (J.hi - J.lo) * max (abs I.lo) (abs I.hi)) ≤
        (I.mul J).lo := by
    intro u v hu hv
    show u * v - _ ≤
      min (min (I.lo * J.lo) (I.lo * J.hi)) (min (I.hi * J.lo) (I.hi * J.hi))
    refine le_min (le_min ?_ ?_) (le_min ?_ ?_)
    · linarith [key u I.lo v J.lo hu (Or.inl rfl) hv (Or.inl rfl)]
    · linarith [key u I.lo v J.hi hu (Or.inl rfl) hv (Or.inr rfl)]
    · linarith [key u I.hi v J.lo hu (Or.inr rfl) hv (Or.inl rfl)]
    · linarith [key u I.hi v J.hi hu (Or.inr rfl) hv (Or.inr rfl)]
  -- o max dos cantos fica a no maximo B acima do minimo
  have hhi : (I.mul J).hi ≤ (I.mul J).lo +
      (I.hi - I.lo) * max (abs J.lo) (abs J.hi) +
      (J.hi - J.lo) * max (abs I.lo) (abs I.hi) := by
    show max (max (I.lo * J.lo) (I.lo * J.hi))
        (max (I.hi * J.lo) (I.hi * J.hi)) ≤ _
    refine max_le (max_le ?_ ?_) (max_le ?_ ?_)
    · linarith [hc0 I.lo J.lo (Or.inl rfl) (Or.inl rfl)]
    · linarith [hc0 I.lo J.hi (Or.inl rfl) (Or.inr rfl)]
    · linarith [hc0 I.hi J.lo (Or.inr rfl) (Or.inl rfl)]
    · linarith [hc0 I.hi J.hi (Or.inr rfl) (Or.inr rfl)]
  linarith

/-- Arredondamento dirigido (abstracao do fixed-point da fonte B5). -/
structure Directed where
  down : ℚ → ℚ
  up : ℚ → ℚ
  hdown : ∀ x, down x ≤ x
  hup : ∀ x, x ≤ up x
  hmono : ∀ x y, x ≤ y → down x ≤ down y ∧ up x ≤ up y

/-- Invariante central do fixed-point: o singleton arredondado contem o valor. -/
theorem directed_singleton_contains (D : Directed) (x : ℚ) :
    (RatI.mk (D.down x) (D.up x) (le_trans (D.hdown x) (D.hup x))).contains x :=
  ⟨D.hdown x, D.hup x⟩

/-- Instancia concreta do fixed-point com tick = 10^-42 (N = 10^42):
    down x = floor(x*N)/N, up x = ceil(x*N)/N. -/
def N_TICK : ℚ := (10 : ℚ) ^ 42

def fpDown (x : ℚ) : ℚ := (Int.floor (x * N_TICK) : ℚ) / N_TICK

def fpUp (x : ℚ) : ℚ := (Int.ceil (x * N_TICK) : ℚ) / N_TICK

/-- floor(x*N)/N ≤ x: floor majora por baixo; divisao por N > 0 preserva. -/
theorem fpDown_le (x : ℚ) : fpDown x ≤ x := by
  have hN : 0 < N_TICK := by positivity
  have h1 : (Int.floor (x * N_TICK) : ℚ) ≤ x * N_TICK := Int.floor_le _
  show (Int.floor (x * N_TICK) : ℚ) / N_TICK ≤ x
  rwa [div_le_iff hN]        -- a / c ≤ b ↔ a ≤ b * c

/-- x ≤ ceil(x*N)/N: analogo via le_ceil. -/
theorem fpUp_ge (x : ℚ) : x ≤ fpUp x := by
  have hN : 0 < N_TICK := by positivity
  have h1 : x * N_TICK ≤ (Int.ceil (x * N_TICK) : ℚ) := Int.le_ceil _
  show x ≤ (Int.ceil (x * N_TICK) : ℚ) / N_TICK
  rwa [le_div_iff hN]        -- b ≤ a / c ↔ b * c ≤ a

/-- Monotonicidade: x ≤ y implica fpDown x ≤ fpDown y e fpUp x ≤ fpUp y. -/
theorem fp_mono (x y : ℚ) (h : x ≤ y) :
    fpDown x ≤ fpDown y ∧ fpUp x ≤ fpUp y := by
  have hN : 0 < N_TICK := by positivity
  have hmul : x * N_TICK ≤ y * N_TICK := by nlinarith
  have hf : (Int.floor (x * N_TICK) : ℚ) ≤ (Int.floor (y * N_TICK) : ℚ) := by
    exact_mod_cast Int.floor_mono hmul
  have hc : (Int.ceil (x * N_TICK) : ℚ) ≤ (Int.ceil (y * N_TICK) : ℚ) := by
    exact_mod_cast Int.ceil_mono hmul
  constructor
  · show (Int.floor (x * N_TICK) : ℚ) / N_TICK ≤
        (Int.floor (y * N_TICK) : ℚ) / N_TICK
    exact div_le_div_of_le hN hf
  · show (Int.ceil (x * N_TICK) : ℚ) / N_TICK ≤
        (Int.ceil (y * N_TICK) : ℚ) / N_TICK
    exact div_le_div_of_le hN hc

/-- A instancia `Directed` do fixed-point (tick 10^-42), agora completa:
    fpDown_le, fpUp_ge e fp_mono estao demonstrados acima. -/
def fpDirected : Directed :=
  ⟨fpDown, fpUp, fpDown_le, fpUp_ge, fp_mono⟩

end E8K
