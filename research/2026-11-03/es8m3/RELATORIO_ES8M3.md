# NUMERIA–RIEMANN — ES8-M3: certificação de Hankel e Stieltjes 3×3

**Data:** 03/11/2026  
**Objeto:** calcular os momentos completos normalizados `m₀,...,m₅` pela fórmula explícita de Guinand–Weil, e testar a positividade definida das primeiras matrizes 3×3 de Hankel e Stieltjes.  
**Resultado:** **certificação racional condicional aos lemas analíticos declarados**, ambas positivas definidas, sem conclusão sobre a Hipótese de Riemann nem sobre E7/V24.

## 1. Definições e identidade independente de RH

Para `a=1/50`, `zρ=(ρ-1/2)/i` (complexo para um zero eventualmente fora da reta crítica), a fórmula explícita clássica aplicada à gaussiana admite

\[
F(a)=\sum_{\rho}e^{-a z_\rho^2}=2e^{a/4}+A(a)-P(a),
\quad A(a)=\frac1\pi\int_0^\infty e^{-at^2}\Omega(t)dt,
\]
\[
\Omega(t)=\Re\psi(\tfrac14+it/2)-\log\pi,
\quad P(a)=\frac1{\sqrt{\pi a}}\sum_{n\ge2}
\frac{\Lambda(n)}{\sqrt n}\exp\left(-\frac{(\log n)^2}{4a}\right).
\]

**Sem supor RH**, define-se `μ_k(a)=(-1)^k F^(k)(a)` e `m_k=μ_k/400^k`. A diferenciação em `a>0` é permitida pela convergência apropriada dos lados da fórmula. Os três termos de `m_k` são

\[
m_k=\frac{2(-1)^ke^{a/4}}{(1600)^k}
+\frac1{\pi}\int_0^\infty\left(\frac{t^2}{400}\right)^ke^{-at^2}\Omega(t)dt
+\frac{(-1)^{k+1}}{400^k}\,P^{(k)}(a).
\]

Os **seis** momentos aqui resultantes são do **lado analítico completo** da fórmula, não da lista dos 20 ou 25 primeiros zeros críticos usados nos experimentos anteriores.

## 2. Melhoria decisiva: cotas de quarta derivada

A quadratura de ponto médio de segunda ordem da ES8-M2 fornecia intervalos suficientemente estreitos para 2×2, mas não para 3×3. Na ES8-M3 implementamos Simpson composto sobre blocos de largura `2h` (com `h=T/M`). Seu erro por bloco é limitado por

\[
|E_{\mathrm{Simpson},[l,r]}|\le\frac{h^5}{90}\sup_{t\in[l,r]}|g_k^{(4)}(t)|,
\qquad g_k(t)=\left(\frac{t^2}{400}\right)^k e^{-at^2}\Omega(t).
\]

Da expansão em frações parciais da digama (DLMF §5.7.6), escreve-se

\[
\Omega(t)=\Omega(0)+\sum_{n\ge0} b_n(t),\qquad
b_n(t)=\frac{t^2}{(n+1/4)[4(n+1/4)^2+t^2]}.
\]

Para `j=1,2,3,4`, uma decomposição em frações simples dá

\[
|b_n^{(j)}(t)|\le\frac{j!}{2^j(n+1/4)^{j+1}}.
\]

A soma é majorada pela integral de funções positivas decrescentes. Escolhemos os majorantes globais inteiros

\[
|\Omega'(t)|\le10,\quad |\Omega''(t)|\le36,
\quad |\Omega'''(t)|\le208,\quad |\Omega^{(4)}(t)|\le1632,
\]

além de `|Ω(t)|≤7+18t²` do módulo anterior. O programa aplica Leibniz explicitamente à quarta derivada de `g_k`, limita cada fator positivo pelo extremo direito do bloco e o fator exponencial pelo extremo esquerdo. Todos os coeficientes são racionais e os erros por bloco são somados **para cima**.

A cauda `t≥48` usa `|Ω(t)|≤7+18t²` e a recorrência por partes

\[
\int_T^\infty t^{2j}e^{-at^2}dt\le e^{-aT^2}C_{2j},\quad
C_0=\frac1{2aT},\quad C_{2j}=\frac{T^{2j-1}}{2a}+\frac{2j-1}{2a}C_{2j-2}.
\]

## 3. Cauda da série de digama: Euler–Maclaurin e constante de Euler

Na expansão de `Ω`, a cauda `n≥q` é uma soma de valores da função

\[
f(x)=\frac{t^2}{x(4x^2+t^2)},\qquad x=n+\tfrac14.
\]

As derivadas `f''` e `f^(4)` são não negativas e `f'''` é não positiva. Com `x=q+1/4`, a integral da cauda é

\[
I=\int_x^\infty f(u)du=\tfrac12\log\left(1+\frac{t^2}{4x^2}\right).
\]

Usamos Euler–Maclaurin com correções até `B₄`, e o limitante de resto da Bernoulli periódica de ordem 4, para obter a inclusão **racional**

\[
I+\frac{f(x)}2-\frac{f'(x)}{12}+\frac{f'''(x)}{360}
\le \sum_{n\ge q}b_n(t)
\le I+\frac{f(x)}2-\frac{f'(x)}{12}.
\]

A constante de Euler é cercada via o teorema de Euler–Maclaurin para números harmônicos:

\[
\frac1{2N}-\frac1{12N^2}
< H_N-\log N-\gamma
<\frac1{2N}-\frac1{12N^2}+\frac1{120N^4}
\quad(N=4096).
\]

Essas correções reduziram substancialmente a largura dos momentos. Séries racionais com restos explícitos continuam a fornecer os limites para `π`, logaritmos e `exp(-x)`.

## 4. Derivadas dos primos até ordem 5

Para cada potência de primo `n=p^j`, denote `u=log n`, `x=u²/(4a)` e `P_n(a)=(πa)^(-1/2) Λ(n)n^(-1/2)e^(-x)`. Então

\[
P_n^{(k)}(a)=P_n(a)a^{-k}p_k(x),\quad p_0=1,
\quad p_{k+1}(x)=(x-k-\tfrac12)p_k(x)-x p_k'(x).
\]

A recorrência é implementada com coeficientes `Fraction` exatos, evitando diferenciação numérica. Para `0≤k≤5`, a soma dos módulos dos coeficientes, normalizada por `x^k`, é inferior a 2 sempre que `x≥200` (teste racional explícito no código). Como `log 64>4`, há um majorante decrescente a partir do corte `n≥64`, e a cauda da soma de cada derivada satisfaz

\[
|\text{cauda}_k|\le\frac{40\,625^5}{91}e^{-182}
<\frac{40\,625^5}{91}(10/27)^{182}<10^{-64}.
\]

As constantes usam `Λ(n)≤log n`, `1/sqrt(πa)<5`, `2k log u≤u²` em `u≥4`, e o termo `e^{-u/2}` proveniente de `n^{-1/2}`. É uma cota conservadora, mas inclusiva. A parcela finita `n≤64` é calculada com intervalos racionais e arredondamento externo.

## 5. Resultados da execução refinada

Parâmetros finais: **`a=1/50`, `M=16384`, `terms=96`, `T=48`, `N=64`**. Para rastreabilidade também fornecemos a execução com `M=8192`, cujos intervalos contêm os da execução refinada.

| `k` | limite inferior de `m_k` | limite superior de `m_k` |
|---:|---:|---:|
| 0 | 0.037082581983342 | 0.037082584941936 |
| 1 | 0.018705340578467 | 0.018705340842767 |
| 2 | 0.009549222496999 | 0.009549222571496 |
| 3 | 0.005003322178959 | 0.005003322217620 |
| 4 | 0.002766294357188 | 0.002766294385774 |
| 5 | 0.001691292143400 | 0.001691292170359 |

Extremos decimais arredondados para fora; os extremos como **frações racionais exatas** estão em `RESULTADOS_ES8M3_REFINED.json`.

Formamos as duas matrizes reais simétricas:

\[
H_2=(m_{i+j})_{0\le i,j\le2},\qquad
H_2^{[1]}=(m_{i+j+1})_{0\le i,j\le2}.
\]

Para matrizes reais simétricas, a positividade definida é equivalente à positividade dos menores principais líderes de ordem 1, 2 e 3 (critério de Sylvester). A propagação intervalar desses menores obteve:

| Matriz | determinante de ordem 2 | determinante 3×3 | Status |
|---|---|---|---|
| Hankel `H₂` | `[0.000004220050078101, 0.000004220090980475]` | `[7.220715, 7.328818]·10⁻¹²` | CERTIFIED_POSITIVE_DEFINITE |
| Stieltjes `H₂^[1]` | `[0.000002401193661275, 0.000002401197129545]` | `[6.447066, 6.464591]·10⁻¹²` | CERTIFIED_POSITIVE_DEFINITE |

Os menores de ordem 1 estão, respectivamente, dentro dos intervalos de `m₀` e `m₁` (ambos positivos). O refinamento de `M=8192` para `16384` preservou e estreitou todos os seis intervalos. Os resultados foram confrontados também com cálculos numéricos `mpmath`, usados **somente como diagnósticos**, não como certificados.

## 6. Cotejo com E7/V24 e alcance do resultado

A ES7 certifica o **integral arquimediano escalar gaussiano**. A V24 trabalha com produtos `F_i(r) conj(F_j(r))` de transformadas de Fourier de modos suavizados, além de polos e potências de primos, seguidos da projeção `Bᵀ(G+P)B`. **Não foram provadas as igualdades entre os dois sistemas de funções-teste, entre as matrizes ou entre seus pivôs.** A nova matriz `H₂` é uma matriz de momentos para a variável `z²/400`, não uma matriz projetada da V24.

A fórmula explícita de Weil e a equivalência do critério de positividade para todas as funções-teste admissíveis são conhecidas. Ter encontrado positividade para **dois casos finitos** da família gaussiana (mesmo com intervalos rigorosos) não autoriza estender a conclusão a todas as dimensões, a todos os parâmetros ou à Hipótese de Riemann.

## 7. Fontes e reprodução

- NIST DLMF, §5.7.6: https://dlmf.nist.gov/5.7.E6 (frações parciais da digama).
- NIST DLMF, §2.10(i), §24.17(i): https://dlmf.nist.gov/2.10 e https://dlmf.nist.gov/24.17 (Euler–Maclaurin e restos).
- Guinand–Weil e critérios: https://openmaths.org/courses/NT-EXPF/NT-EXPF-03.html; Bombieri–Lagarias, https://websites.umich.edu/~lagarias/doc/bombieri.pdf.
- Fase anterior ES8-M2 (branch do projeto): `research/numeria-es8m2-moments-2026-11-03`.

**Status científico final:** computação racional verificável para duas matrizes de momentos finitas, baseada em desigualdades analíticas apresentadas e testadas. Não há prova de positividade universal, não há determinação de todos os zeros, não há equivalência estabelecida com E7/V24, e a Hipótese de Riemann permanece sem demonstração.
