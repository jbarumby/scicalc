# NUMERIA P92-T — Formalização das rampas de Bellman em Lean 4

**Base:** [PR #108, P92-S](https://github.com/jbarumby/scicalc/pull/108), sobre [PR #107, P92-R](https://github.com/jbarumby/scicalc/pull/107).

**Alvo final:** a desigualdade geral de Bellman para a função convexa de rampa
\[
f_T(u)=(T-u)_+,\qquad
\min(h,b)f_T(w)+U_0(w+b,b)\le U_0(w,h),
\]
com `U_0(w,h)=h\sum_{j\ge0}f_T(w+jh)`. Isso permitirá eliminar a hipótese analítica `step` do teorema de telescopagem P92-S.

## Monotonicidade e Bellman em cada faixa de suporte fixo

Defina `L=T-w`, e, para um inteiro `n>=1`, as formas fechadas
\[
P_0(L,h,n)=hnL-\frac{h^2n(n-1)}{2},
\qquad
P_1(L,h,n)=P_0(L,h,n)-hL.
\]
São exatamente as somas das rampas se `(n-1)h\le L\le nh`. A P92-T formaliza em Lean:

- `affineRampSum_closed`: fórmula fechada da soma afim de `n` termos.
- `rampPoly0_mono_same_cell`: monotonicidade de `P_0` dentro de uma célula de suporte.
- `rampPoly1_antitone_same_cell`: antitonicidade de `P_1` dentro da célula.
- `bellman_ramp_small_step_same_cell`: desigualdade de Bellman no ramo `b<=h` sob suporte comum.
- `bellman_ramp_large_step_same_cell`: desigualdade de Bellman no ramo `h<=b` sob suporte comum.
- `bellman_ramp_same_cell`: combinação dos ramos quando `(n-1)max(h,b)<=L<=n min(h,b)`.

**Resultados ampliados já compilados no Lean:** o arquivo também formaliza o núcleo real `finiteRamp0` como soma de `max(L-jh,0)` e demonstra `finiteRamp0_eq_affine`, `finiteRamp0_closed`, `ramp_support_of_cell` e `finiteRamp_bellman_same_cell`. Assim, **a identificação entre soma de rampas e polinômio está provada em toda célula de suporte fixo**. Os lemas `ramp_terms_zero_after_support`, `rampPoly0_boundary_agree`, `rampPoly1_boundary_agree` e `finiteRamp0_boundary_agree` comprovam que as parcelas posteriores ao corte são zero e que os dois lados de cada nó `L=nh` concordam exatamente.

**Misturas positivas finitas:** `bellman_finite_positive_mixture` prova, por soma de desigualdades, que Bellman é estável sob combinações **finitas** de núcleos já válidos, com coeficientes não negativos. Esse resultado é condicional à desigualdade para cada núcleo componente. Sua compilação está registrada no GitHub Actions.

**Atenção ao escopo:** apesar dessas provas sem `sorry` ou `admit`, a **desigualdade global** para passos quaisquer ainda **não está provada no Lean**, pois falta formalizar a colagem de diferentes números de parcelas. A colagem tem demonstração matemática no arquivo `GLUING_RAMPAS.md`, que não substitui um teorema Lean. A integração de misturas positivas arbitrárias por Tonelli e as etapas R5–R6 de séries infinitas também continuam pendentes.

## Porta matemática para o próximo avanço

1. Definir `U_0(L,h)` via soma finita `Finset.range (Nat.ceil (L/h))`, com o suporte exato.
2. Usar o suporte exato `n=ceil(L/h)` para passar dos lemas por célula ao potencial global; a igualdade nos nós já foi formalizada.
3. Colar as desigualdades de células sucessivas; como há apenas um número finito de nós entre dois passos positivos `h<=b`, concluir monotonicidade global.
4. Integrar misturas positivas de rampas: prova para `f_X`, para `u^{-2}` e para o núcleo `phi_h`; eliminar os parâmetros condicionais da P92-S.
5. Formalizar convergência pareada R5→R6.

## Reprodução

Da raiz:
```bash
cd NUMERIA_P92T_LEAN
lake update
lake exe cache get
lake build
```

A CI `.github/workflows/numeria-p92t-lean-bellman.yml` compila em Lean 4.35.0 / Mathlib v4.35.0, sem `sorry` ou `admit`.

**Escopo:** formalização parcial das rampas; a P92-R ainda não é uma prova Lean completa dos Lemas 5/6 de Preissmann. A Hipótese de Riemann permanece em aberto. Manter PR draft sem merge automático.
