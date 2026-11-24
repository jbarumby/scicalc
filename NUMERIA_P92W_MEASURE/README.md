# P92-W — medida concreta, controle do potencial e gate de integrabilidade

**Base:** [PR #111, P92-V](https://github.com/jbarumby/scicalc/pull/111). A PR #112 contém arquivos Lean separados para medida, potencial e redução local de integrabilidade.

## 1. Definição correta da medida

Para `X>0`, tomar

\[
\boxed{\mu_X=\bigl[\mathbf1_{\{0<T<X\}}\,2T^{-3}\bigr]\,dT+X^{-2}\delta_X.}
\]

A densidade de `MeasureKernel.lean` usa `ENNReal.ofReal` para produzir um peso não negativo, `volume.withDensity` e a medida `Measure.dirac X` com coeficiente não negativo. Portanto, define uma **medida positiva genuína**.

**Não é uma medida finita:** `∫_0^ε 2T^{-3} dT=∞`. Essa divergência não invalida as integrais de Bellman, porque seus integrandos são zero em `T≤w` para `w>0`.

## 2. Massa na cauda relevante

Para `w>0`, a massa estritamente acima de `w` é

\[
\mu_X((w,\infty))=
\begin{cases}
\displaystyle\int_w^X\frac2{T^3}\,dT+\frac1{X^2}
=\frac1{w^2},&w<X,\\
0,&w\ge X.
\end{cases}
\]

Em particular,

\[
\boxed{\mu_X((w,\infty))\le w^{-2}<\infty.}
\]

**Atenção ao ponto de fronteira:** para `w=X`, o átomo em `T=X` **não pertence** ao intervalo aberto `(w,\infty)`, logo a massa é `0`. Já na cauda fechada `[X,\infty)`, a massa seria `X^{-2}`. Essa distinção deve ser preservada no Lean.

A prova da igualdade acima é uma dedução **analítica em papel**, combinando a primitiva `-1/T²` com o átomo. **Ainda não é um teorema Lean** até que o cálculo da medida `withDensity` seja formalizado.

## 3. Majorante real do potencial

Para `L=(T-w)_+`, `h>0` e `n=\lceil L/h\rceil`,

\[
0\le U_T(w,h)
=h\sum_{j=0}^{n-1}(L-jh)_+
\le hnL
\le (L+h)L.
\]

A última desigualdade segue de `n< L/h+1` quando `L>0` (e é trivial quando `L=0`). Se `T≤X`,

\[
\boxed{0\le U_T(w,h)
\le (X-w)_+\bigl((X-w)_++h\bigr).}
\]

O potencial **se anula para `T≤w`**. Os arquivos `PotentialBound.lean` e `LocalIntegrability.lean` tratam dessas afirmações, e somente devem ser classificados como prova kernel após CI aprovada.

## 4. Integrabilidade de Bellman — argumento matemático correto

Para `w>0, h>0`, a função `U_T(w,h)` é contínua em `T`: em cada faixa de suporte tem fórmula polinomial e, nos pontos de mudança, a nova parcela se anula (P92-T). Portanto, é mensurável.

Como ela é zero em `T≤w` e uniformemente limitada no suporte `T≤X`,

\[
0\le\int U_T(w,h)\,d\mu_X(T)
\le\frac{(X-w)_+\bigl((X-w)_++h\bigr)}{w^2}<\infty.
\]

O mesmo raciocínio vale para `T\mapsto(T-w)_+`, usando a cota `(X-w)_+`, e para o potencial deslocado `U_T(w+b,b)` quando `b>0` (substituir `w` por `w+b>0`).

**Assim, a hipótese de integrabilidade dos três termos de** `bellman_integral_ramps` **é satisfeita matematicamente para a medida concreta** `μ_X` e `w,h,b>0`.

Para que isso seja um **resultado formal Lean completo**, restam as seguintes obrigações separadas:
1. demonstrar `μ_X((w,\infty))<∞` e a concentração `T≤X` por `withDensity` e `Measure.dirac`;
2. demonstrar a mensurabilidade forte da função `T↦rampPotential T w h` (por suporte local finito ou por colagem contínua);
3. aplicar `LocalIntegrability.rampPotential_integrable_of_local_tail` à medida `μ_X` e seus três termos;
4. estabelecer a igualdade de `∫ (T-w)_+ dμ_X` com a integral de intervalo mais o átomo já provados na P92-V, e então com `max(1/w-1/X,0)`.

Não declarar (1)–(4) formalizados sem compilação do Lean. O intervalo `w=0` não faz parte do enunciado, pois a massa singular exige `w>0`.

## 5. Próximas conexões

Ao fechar (1)–(4), instanciar `mixedRamp_bellman` da P92-U com `μ_X`. Depois formalizar a desigualdade R2 por telescopagem em listas finitas e a contabilidade truncada do R3 na P92-S; posteriormente R4 e o compensador pareado R5–R6. Lemas 5/6 e RH não estão integralmente formalizados.

**Governança:** PR #112 rascunho, sem merge automático; status técnico da CI deve ser registrado por commit exato.
