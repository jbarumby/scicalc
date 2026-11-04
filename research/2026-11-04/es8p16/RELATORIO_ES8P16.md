# ES8-P16 — Revisão matemática de `deepseek_latex_20261104_525471.tex`

## Conclusão

A rota de extensão a uma forma sobre `ℓ²` é matematicamente natural **se** forem provadas a definição coerente da família infinita, a convergência das formas e uma estimativa uniforme da parte central e do erro. O arquivo enviado identifica corretamente essas dificuldades, mas contém saltos lógicos que impedem sua utilização como demonstração.

## 1. O que os certificados de fato dizem

O código ES8-P14 contém centros `center4`, `center6`, `center8`; extremos intervalares `lower/upper` e envelopes `error4/6/8`. Ao comparar as matrizes por `fractions.Fraction`, obtivemos **208 igualdades exatas** entre cantos nas oito transições: 4→6 e 6→8 para cada uma das quatro famílias. São certificados de nestedness **finita**. Para dizer que existe *uma única* matriz infinita, precisaríamos fornecer as entradas para todos os pares `(i,j)` e provar compatibilidade das truncagens para todo `N`. Os mesmos oito primeiros índices permitem extensões distintas; isso foi formalizado no módulo Lean.

## 2. Pivô LDL e autovalor não são a mesma coisa

Para todo `t∈ℚ`, considere

\[
U_t=\begin{pmatrix}1&t\\0&1\end{pmatrix},\quad
A_t=U_t^\top U_t=\begin{pmatrix}1&t\\t&t^2+1\end{pmatrix}.
\]

Os pivôs LDL são `1,1`, mas no vetor `v=(-t,1)`:

\[
v^\top A_t v=1,\qquad \|v\|^2=t^2+1,
\qquad \frac{v^\top A_t v}{\|v\|^2}=\frac1{t^2+1}.
\]

Assim, nenhum limite espectral uniforme estritamente positivo decorre só de pivôs `≥δ`. Seria necessário controlar a menor singularidade de `U_t`, ou trabalhar diretamente com a forma. Em particular, **adicionar `3/50` aos pivôs de `C−3/50 I` não recupera os pivôs de `C`**. O segundo pivô do caso 4×4 já fornece um contraexemplo racional calculado no JSON.

## 3. O que significa orçamento de erro depois de projetar

Seja `y=B x`. Um limite bruto da forma `|yᵀ E y|≤δ||y||²` implica somente

\[
|x^\top B^\top E Bx|\le\delta\,\|Bx\|^2
\le\delta\,\|B\|_{op}^2\|x\|^2.
\]

A ES8-P15 controla apenas a cauda arquimediana **não projetada** sob hipótese E7. Portanto ela não permite concluir, sem controle de `B_N`, o orçamento `1/100` da forma projetada. O novo arquivo Lean contém o teorema geral de perturbação por projeção, com a hipótese `||Bx||²≤κ||x||²`. O exemplo 1×1 `B=[2]` aumenta o erro `1/100` para `4/100`.

## 4. Rota Hilbert–Schmidt: evidência concreta, não prova de infinito

Para a matriz de envoltórias racional da ES8-P14, a soma \(\sum_{ij}e_{ij}^2\) é:

| Dimensão | Soma exata | Menor que `1/10000`? |
|---|---|---|
| 4×4 | `3913513/1000000000000` | Sim |
| 6×6 | `18220307/1000000000000` | Sim |
| 8×8 | `12382243/200000000000` | Sim |

Esses valores tornam plausível investigar um limite de **norma de Hilbert–Schmidt**. Mas três termos finitos não estabelecem \(\sum_{i,j=0}^{\infty}e_{ij}^2<\infty\), tampouco `||E∞||≤1/100`. As envoltórias da ES8-P14 são limites de erro, não necessariamente a matriz de erro real. Além disso, o operador completo `C∞` pode ser não limitado, o que exigiria uma rota por formas fechadas e definição de domínio.

## 5. Alterações necessárias no LaTeX anexado

- `ℓ²(ℕ;ℚ)` não é o espaço de Hilbert completo padrão; use `ℓ²(ℕ;ℝ)` ou `ℓ²(ℕ;ℂ)`, mantendo ℚ apenas como subespaço denso para certificados.
- A rota por resolventes fortes não resulta **automaticamente** de nestedness sem hipóteses adequadas de domínio e operador; no caso limitado é melhor provar a extensão contínua da forma por densidade diretamente.
- Na rota KLMN, não basta declarar a forma fechada: é necessário verificar a relativa limitação da perturbação e seu domínio; a terminologia KLMN não substitui uma prova de convergência.
- O enunciado proposto `infinite_positive` mistura matrizes centrais, `+3/50 I` e erro e toma como hipótese a própria margem desejada `UniformFinitePositivity Minf (1/20)`. Isso não produz a estimativa espectral nova desejada. A organização correta é provar `C_N−3/50 I≥0`, `||E_N||≤1/100`, daí `C_N+E_N≥1/20 I`.
- Uma margem de positividade para **L fixo** não demonstra o critério universal de Weil: suportes arbitrários e identificação analítica completa continuam necessários.
- O texto chama `ErrorFormNorm` de “falso com os dados atuais”, mas só temos 4, 6, 8 e nenhuma cota para os demais índices. O correto é **não demonstrado**.
- O “esqueleto Lean” com `sorry` é uma proposta; não foi testado nem deve ser tratado como prova.

## 6. Rota de continuidade recomendada

1. Construir e autenticar uma definição global da família `C_N`, com uma prova de compatibilidade para todo `N`. Testar se a projeção utilizada nas matrizes históricas preserva esse aninhamento em qualquer dimensão.
2. Derivar as normas de `B_N` (ou seus coeficientes) e aplicar o novo lema de projeção para obter o orçamento após o projetor.
3. Investigar cotas estruturais para `C_N−3/50 I` válidas para **todo** `N`. Controle de pivôs isolados não basta: verificar o condicionamento da fatoração.
4. Fechar a desigualdade analítica E7 para a integral de fato, a avaliação rigorosa dos termos arquimedianos e de primos e os limites de interpolação.
5. Se a forma inteira for definida e continuar limitada ou fechável, aplicar a passagem à dimensão infinita já formalizada condicionalmente em ES8-P15.

**Status científico:** nenhuma etapa demonstra RH ou positividade universal de Weil. Esta auditoria identifica condições necessárias, formaliza relações algébricas e testa os dados de baixa dimensão.
