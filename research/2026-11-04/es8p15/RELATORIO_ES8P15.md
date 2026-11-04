# NUMERIA–RIEMANN — ES8-P15: caudas, séries infinitas e margens uniformes

**Data:** 04/11/2026. **Status:** prova parcial e condicional; consulte o resultado da compilação Lean no CI. **Base:** ES8-P14, `e687605389d7baf97bae7d59c6bd09c567a696c4`.

## 1. Cauda arquimediana — o que foi e o que não foi provado

A estimativa E7, Proposição 5, para modos de mesma paridade, pressupõe `R ≥ max(1,2w_max)` e dá:

\[
 |T_{ij}(R)| \le K_{ij}\frac{8+\log(3R)+1/3}{3R^3},\qquad
 K_{ij}=\frac{16w_iw_j}{\pi L}\ge0,
 \quad w_k=\frac{k\pi}{2L}.
\]

**Novo passo formalizável:** como `log(3R) ≤ 3R` para `R ≥ 1`, segue

\[
 |T_{ij}(R)|\le K_{ij}\left(\frac1{R^2}+\frac3{R^3}\right).
\]

Em Lean `log_tail_bounded` prova a desigualdade logarítmica sem recorrer a aproximações; `e7_bound_implies_rational_enclosure` prova a implicação para qualquer resto que já satisfaça a cota E7. `rationalTail_tendsto` e `bounded_archTail_tendsto` formalizam o limite da envoltória e a passagem pelo princípio do confronto. **A identificação do verdadeiro integral `T_ij` com a hipótese E7 continua externa à formalização.**

### 1.1. Proposta quantitativa uniforme para a cauda **não projetada**

Para a família de modos `1≤i,j≤n`, com `L=9/20` e `π<4`,
`K_ij = 4πij/L³ ≤ (128000/729)n²`. Escolhendo `R(n)=200n²`,
a condição `R(n)≥2w_n` vale para todo `n≥1`. A soma por linha dos
majorantes (antes da projeção) satisfaz, com aritmética racional exata:

\[

\sum_{j=1}^{n} E^{\rm cauda}_{ij}(200n^2)
\le\frac{16}{3645n}+\frac{2}{30375n^3}
\le\frac{406}{91125}<\frac1{100}.
\]

Esta identidade de majorantes está verificada por testes Python (`Fraction`).
A sua formalização Lean como teorema separado está **pendente**. Mesmo
uma prova mecanizada desta desigualdade controlaria apenas **a cauda bruta**,
não a cauda após os coeficientes da projeção e nem os demais erros de entrada.

A passagem `R → ∞` aplica-se inicialmente a **cada par fixo** `(i,j)`. A constante `K_ij` depende dos índices. Não há, portanto, uniformidade em `N` sem estimativa adicional. O truncamento histórico `R = 160` só pode ser usado quando `R ≥ 2w_max`; como `w_max ∝ N/L`, essa condição eventualmente falha para `N` crescente.

## 2. `Finset.sum → tsum` e formas quadráticas

O teorema `tsum_margin_of_partial_sums` afirma, para funções `f,g : ℕ → ℝ`, ambas somáveis:

\[
  \forall n,\quad \frac1{20}\sum_{k<n}g(k)\le\sum_{k<n}f(k)
  \quad\Longrightarrow\quad
  \frac1{20}\sum_{k\ge0}g(k)\le\sum_{k\ge0}f(k).
\]

Usa a convergência `Summable.tendsto_sum_tsum_nat`, sem `fin_cases`. Já `quadratic_margin_at_limit` serve para as somas *duplas* de entradas matriciais: se as formas e os quadrados das normas dos truncamentos convergem e a desigualdade vale para cada truncamento, vale também no limite.

A passagem ao limite preserva a desigualdade **não estrita**, mesmo se
cada caso finito for estrito. Para obter estrita positividade no limite,
é necessária hipótese adicional, como margem positiva em norma não nula.

Esses teoremas **não** provam a somabilidade da expressão de Weil, nem a continuidade da forma de Weil na topologia escolhida: isso exigirá uma estimativa apropriada e/ou um operador limitado (ou domínio fechado controlado).

## 3. Margens uniformes: teorema condicional e obstruções

`uniform_margin_all_dimensions` reaproveita `ES8P14.positive_from_LDL_and_box`, de maneira polimórfica em `n`. Se, **para todo `n`**, forem provados:

- `C_n-(3/50)I` positiva definida;
- erro de entrada simétrico, limitado por `ε_ij^(n)`;
- `max_i ∑_j ε_ij^(n) ≤ 1/100`;

então `M_n` satisfaz a margem `1/20` em todos os vetores racionais não nulos, em todas as dimensões.

**Nenhuma das três condições universais foi derivada dos certificados ES8-P13 ou ES8-P14.** O Lean verifica aqui a implicação, e não suas hipóteses no NUMERIA. Uma extensão posterior para vetores reais ainda dependerá de apropriação da estrutura contínua da forma.

### Contraexemplos que impedem inferências indevidas

1. `C_n=(1/(n+1))I_n` é positiva definida para cada `n≥1`, mas não satisfaz a margem `1/20` para `n≥20`; isso prova que positividade dimensão a dimensão **não** garante uniformidade.
2. Um erro por entrada `ε_ij = 1/10000` pode parecer pequeno, mas sua soma por linha é `n/10000`, acima de `1/100` quando `n>100`. Portanto, um limite de erro pontual que independe de `N` não basta.
3. No NUMERIA, `K_ij ∝ i j/L^3` na cota conservadora e o número de entradas por linha cresce. Sem compensações adicionais, é arriscado extrapolar o orçamento fixo `0,01` para `N → ∞`.

## 4. Critérios de conclusão da fase

- Lean 4.35.0 + Mathlib: compilar `ES8P14.lean` e `ES8P15.lean`, sem `sorry`, `admit` ou axiomas adicionais.
- Testes locais: executar `python -m unittest -v test_es8p15.py`, conferindo a fórmula racional e contraexemplos.
- `RESULTADOS_ES8P15.json` registra parâmetros e verificações numéricas de alta precisão apenas como testes — não como substituto da prova Lean.
- Não alegar positividade universal ou prova da Hipótese de Riemann.

## 5. Próxima obrigação científica

Provar uma **estimativa de núcleo e de projeção em índices arbitrários**, capaz de controlar simultaneamente o menor autovalor de `C_n` e as somas de erro por linha, com `R=R(n)` adequado. Estabelecer a convergência das formas truncadas e a cobertura/densidade das funções-teste para suportes arbitrários. Uma cota uniforme `1/20` pode revelar-se falsa para a família escolhida; nesse caso deve-se buscar uma margem `c(n)>0` com passagem ao limite apropriada, sem pressupor a meta numérica atual.