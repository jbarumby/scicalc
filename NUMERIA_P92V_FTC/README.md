# NUMERIA P92-V — Teorema Fundamental do Cálculo para o núcleo truncado

**PR:** [#111](https://github.com/jbarumby/scicalc/pull/111), baseada na [P92-U #110](https://github.com/jbarumby/scicalc/pull/110). Escopo: integral real com a densidade do núcleo truncado de Preissmann.

## Problema específico

A P92-U formalizou Bellman para misturas de rampas com respeito a medidas positivas, assumindo integrabilidade, e provou a contribuição do átomo X⁻² δ_X. Faltava verificar que a integral da densidade contínua equivale à diferença dos extremos da primitiva.

Para 0<u≤X e
\[
F_u(T)=-\frac2T+\frac u{T^2},
\]
temos
\[
F_u'(T)=\frac{2(T-u)}{T^3} \quad(T>0).
\]

O módulo Lean [TruncatedFTC.lean](../NUMERIA_P92T_LEAN/TruncatedFTC.lean) demonstra e verifica no kernel Lean
\[
\boxed{\int_u^X\frac{2(T-u)}{T^3}\,dT
=F_u(X)-F_u(u)}
\]
pelo teorema fundamental do cálculo, com integrabilidade decorrente da continuidade em intervalo positivo. Combinado à P92-U:
\[
\boxed{
\int_u^X\frac{2(T-u)}{T^3}\,dT+
\frac1{X^2}\int(T-u)_+\,d\delta_X(T)
=\frac1u-\frac1X.}
\]

O arquivo também demonstra a equivalência com a parte positiva
\[
\int_u^X\frac{2(T-u)}{T^3}\,dT + X^{-2}\int(T-u)_+\,d\delta_X(T)
=\left(\frac1u-\frac1X\right)_+
\quad (0<u\le X),
\]
bem como a nulidade de \((1/u-1/X)_+\) para \(u\ge X>0\).

## O que isso não demonstrará sozinho

A avaliação de uma integral de intervalo não identifica ainda no Lean a medida concreta
\[
\mu_X(dT)=\mathbf1_{(0,X)}(T)\frac2{T^3}dT + X^{-2}\delta_X(dT)
\]
e a função (1/u-1/X)₊ para todo u>0. É necessário formalizar a densidade da medida, provar a integrabilidade dos três potenciais de Bellman e associar a representação integral aos telescopamentos quantitativos R2/R3.

A compensação pareada R5→R6 e a transferência das caudas R4 permanecem abertas. Não há prova dos Lemas 5/6 completos em Lean, nem da Hipótese de Riemann.

## Verificação

Lean 4.35.0, Mathlib v4.35.0. O novo módulo está incluído como alvo obrigatório Lake na pasta NUMERIA_P92T_LEAN; workflow: .github/workflows/numeria-p92v-truncated-ftc.yml, sem sorry ou admit.

**Estado certificado:** [GitHub Actions #45312336677](https://github.com/jbarumby/scicalc/actions/runs/45312336677) — `completed/success` no commit `1a304bec763eca15a9736868e1b842f8eb8fc2ae`. Compilaram os sete alvos Lean obrigatórios; os controles sem `sorry` e `admit` também passaram. Os teoremas principais incluem `truncated_density_FTC`, `truncated_reciprocal_integral_with_atom`, `truncated_reciprocal_as_positive_part` e `truncated_reciprocal_zero_above_cutoff`.

**Próximo gate P92-W:** identificar a medida concreta \(\mu_X\) pelo `withDensity` e pela massa de Dirac, provar a integrabilidade dos três potenciais para \(w,h,b>0\) e então conectar a identidade provada ao teorema `mixedRamp_bellman` da P92-U e ao R2/R3 da P92-S. A PR continua draft, sem merge automático.
