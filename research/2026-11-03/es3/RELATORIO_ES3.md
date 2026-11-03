# NUMERIA–RIEMANN — ES3: auditoria Spin–Zeta

**Data:** 3 de novembro de 2026  
**Fonte primária da hipótese:** `Qwen_latex_20261103_m86z47xnk.txt` (documento especulativo submetido pelo usuário).  
**Status:** auditoria matemática negativa para a proposta quadrática simples; controle espectral e testes numéricos, **sem prova de RH**.

## 1. Verificação do operador de Sanctuary

Sejam as matrizes usuais de Pauli e, para orientações `n_x,n_y,n_z∈{−1,1}`,

\[
s=\tfrac12(I+n_z\sigma_z+n_x\sigma_x+i n_y\sigma_y),\qquad A=2s-I.
\]

O cálculo **exato**, para as oito orientações, encontra:

\[
A^2=I,\quad s^2=s,\quad \operatorname{tr}s=1,\quad
\operatorname{spec}(s)=\{0,1\},
\]
\[
\operatorname{spec}(s+s^\dagger)=\{1-\sqrt2,1+\sqrt2\},\quad
\operatorname{spec}(s^\dagger s)=\{0,2\}.
\]

Assim, o texto de origem troca indevidamente o espectro de `A` pelo de `s` e confunde a matriz de ressonância com a sua parte sem traço. A matriz `s^†s` é positiva semidefinida e **singular**. Essa positividade não identifica os zeros da função zeta.

## 2. Auditoria dos sítios finitos

No modelo de spins independentes em `m` fatores de dimensão dois, usamos

\[
\mathcal A_m=\sum_{v=1}^{m}I^{\otimes(v-1)}\otimes A\otimes I^{\otimes(m-v)}.
\]

Os termos comutam; como cada `A` tem espectro `{-1,1}` e é diagonalizável, o espectro da soma é

\[
\{-m,-m+2,\ldots,m-2,m\},\quad
\operatorname{mult}(-m+2k)=\binom{m}{k}.
\]

**Verificado simbolicamente para `m=1,2,3,4`**. O espectro não é `±√(2^m)`; o número `2^{3m}` de escolhas de orientação não é a dimensão do espaço tensorial e não determina multiplicidade espectral. O argumento de cardinalidade vale para o **setor finito de spin**; se se tensoriza com espaço adélico infinito, será preciso especificar a ação efetiva nesse fator. Nenhum operador adélico é definido ou testado aqui.

## 3. Auditoria do Hamiltoniano quadrático

Do texto de origem:

\[
H=\alpha(a+a^\dagger)+\beta a^\dagger a+\gamma(a^2+(a^\dagger)^2),\quad
 a=(u+ip)/\sqrt2,\quad p=-i\,d/du.
\]

Para parâmetros reais, no regime **confinante** `β>2|γ|`, a reescrita exata é

\[
H=\tfrac{\beta+2\gamma}{2}u^2+\tfrac{\beta-2\gamma}{2}p^2
+\sqrt2\alpha u-\tfrac\beta2.
\]

Completar o quadrado e escalar canonicamente `u,p` fornece

\[
\boxed{E_n=\underbrace{\frac{\Omega-\beta}{2}
-\frac{\alpha^2}{\beta+2\gamma}}_{E_0}+n\Omega,\qquad
\Omega=\sqrt{\beta^2-4\gamma^2},\ n=0,1,2,\ldots}
\]

Os níveis possuem espaçamento **constante**, e o número de níveis até `T` é `O(T)`. Já a contagem real de zeros não triviais da zeta satisfaz, incondicionalmente,

\[
N(T)=\frac{T}{2\pi}\log\frac{T}{2\pi}-\frac{T}{2\pi}+O(\log T).
\]

**Resultado negativo rigoroso para a família quadrática no regime acima:** ela não pode ter todo o conjunto das ordenadas positivas dos zeros como espectro, mesmo que alguns níveis sejam ajustados. Acoplar uma quantidade finita de estados de spin por perturbação limitada não altera o termo dominante da contagem. Regimes não confinantes não fornecem o espectro discreto adequado desse modelo simples. Esse argumento NÃO exclui operadores espectrais mais gerais de Hilbert–Pólya ou Connes.

## 4. Experimentos numéricos independentes de uma demonstração de RH

Parâmetros ilustrativos: `α=0.4`, `β=3.0`, `γ=0.5`, `Ω≈2.82842712474619`; diagonalizamos truncamentos da base de Hermite com `M=12,24,48` e comparamos os seis menores autovalores com os níveis exatos. O erro máximo cai de aproximadamente `1.02×10^−2` em `M=12` para `1.24×10^−14` em `M=48` (precisão de ponto flutuante, NÃO certificação intervalar).

Os primeiros **20** valores `γ_j=Im(zetazero(j))` vieram de `mpmath` com 65 dígitos. Para `j=0,…`, ajustamos o modelo **mais flexível** `γ_j≈c+Ω j` aos oito primeiros e medimos erro nos seis seguintes. Permitir `c` arbitrário favorece o modelo: a família quadrática original tem menos liberdade espectral. Resultado: `RMSE` de treinamento `≈1.232`, de validação `≈5.328`. Tal resultado numérico é apenas ilustrativo; a incompatibilidade assintótica acima é a razão matemática decisiva.

## 5. Coeficientes de Li: armadilha de positividade trivial

Para um zero conhecido na reta crítica `ρ=1/2+iγ`, definimos `θ(γ)=π−2 arctan(2γ)` de modo que `1−1/ρ=e^{iθ}`. Pareando `ρ` e `\bar ρ`, a contribuição ao coeficiente de Li é

\[
\bigl[1-(1-1/\rho)^n\bigr]
+\bigl[1-(1-1/\bar\rho)^n\bigr]
=2(1-\cos(n\theta))\ge0.
\]

Por isso a soma **parcial** dos primeiros 20 pares sobre a reta crítica já é não negativa **por construção**. Não prova que os demais zeros estejam nessa reta. Por exemplo, para `n=1`:

\[
\lambda_1=1+\gamma_{\rm Euler}/2-\tfrac12\log(4\pi)
\approx0.02309570896612,
\]

mas a soma dos 20 pares conhecidos é aproximadamente `0.01595936775678`. Não é a mesma quantidade, e o resíduo não é uma cota certificada da cauda completa. A definição geral de `λ_n` exige a soma de **todos** os zeros, contados com multiplicidade e limite simétrico; RH equivale a `λ_n≥0` para todo `n≥1`.

A fórmula candidata `Tr(1−cos(nH))` só teria relação com `λ_n` sob uma construção que identificasse **todos os ângulos** `θ_ρ`, com multiplicidades e condições de classe de traço. **Se** `H` tiver espectro das ordenadas `γ`, será necessário primeiro aplicar a transformação espectral `θ(H)=π−2 arctan(2H)` (para o setor `γ>0`) e justificar a normalização/fator 2. A simples positividade de `1−cos(nH)` para `H` autoadjunto é tautológica, não o critério de Li.

## 6. Relação com NUMERIA E7, V24, NB1, NL2

- **Útil**: teste adversarial para rejeitar falsos operadores de Hilbert–Pólya, verificação de autoadjunção, normalização e contagem espectral; tecnologia de testes complementa NB1/NL2.
- **Não há** transporte da forma de Weil, identificação de matrizes projetadas E7/V24 nem certificado de todos os zeros.
- Para continuar além do contraexemplo é necessário definir um operador explícito não quadrático em espaço funcional bem especificado e provar a correspondência de traços/fórmula explícita sem pressupor RH.

## 7. Referências

1. Sanctuary, B. C. *Structure of a spin 1/2*, arXiv:0908.3219 — artigo de base (não estabelece conexão com zeta).
2. Li, X.-J. (1997), *The positivity of a sequence of numbers and the Riemann hypothesis*, *J. Number Theory* 65, 325–333.
3. Bombieri, E.; Lagarias, J. C. (1999), *Complements to Li's Criterion for the Riemann Hypothesis*, *J. Number Theory* 77, 274–287, DOI 10.1006/jnth.1999.2392.
4. Fórmula de Riemann–von Mangoldt (teorema de contagem dos zeros).
5. `Qwen_latex_20261103_m86z47xnk.txt`: hipótese recebida, **não** é uma fonte de resultados provados.

**Conclusão:** o modelo de spin simples ajuda como teste adversarial; a conjectura Spin–Zeta apresentada não gera um candidato espectral adequado sem uma nova construção analítica.
