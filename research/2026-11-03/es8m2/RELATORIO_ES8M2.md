# NUMERIA–RIEMANN — ES8-M2: primeiros determinantes certificados de Hankel e Stieltjes

**Data:** 03/11/2026  
**Objetivo:** estender a ES8-M de momentos completos μ₀ e μ₁ até μ₂ e μ₃, certificar por aritmética racional as duas matrizes 2×2 de momentos e incorporar os arquivos `References_LeapSpace (7).csv` e `(8).csv`.  
**Status científico:** **certificação de positividade definida para duas matrizes FINITAS de momentos gaussianos** (a=1/50), sob fórmulas analíticas conhecidas; **nenhuma prova de positividade universal de Weil ou da Hipótese de Riemann**. Não se reivindica novidade da identidade global explícita. A equivalência com as matrizes E7/V24 permanece **não estabelecida**.

## 1. Identidade analítica e convenções

Na fórmula explícita de Guinand–Weil, com zρ=(ρ−1/2)/i (que pode ser **complexo**, sem pressupor RH), adotamos a função par e inteira hₐ(t)=exp(−at²), a>0. A fórmula clássica determina a soma de todos os zeros (com multiplicidades, inclusive eventuais zeros fora da reta crítica):

\[
F(a)=\sum_\rho e^{-a z_\rho^2}
=2e^{a/4}+A(a)-P(a),
\]
\[
A(a)=\frac1\pi\int_0^\infty e^{-at^2}\Omega(t)\,dt,\quad
\Omega(t)=\operatorname{Re}\psi(\tfrac14+\tfrac{it}{2})-\log\pi,
\]
\[
P(a)=\frac1{\sqrt{\pi a}}\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}
\exp\!\left(-\frac{(\log n)^2}{4a}\right).
\]

A extensão por diferenciação é justificada pela convergência gaussiana uniforme em parâmetros a>0 em conjuntos compactos:

\[
\mu_k(a)=\sum_\rho z_\rho^{2k}e^{-az_\rho^2}
=(-1)^k\frac{d^k}{da^k}F(a).
\]

Na ES8-M2 certificamos **somente** a=1/50 e k=0,1,2,3; para melhorar condicionamento definimos **mₖ=μₖ/400ᵏ**, que equivale à variável real x=z²/400 quando RH é verdadeira.

Os termos de primos diferenciados usam xₙ=(log n)²/(4a) e os polinômios pₖ tais que P⁽ᵏ⁾ₙ=Pₙ a^(−k) pₖ(xₙ):

\[
\begin{aligned}
p_0(x)&=1,\\
p_1(x)&=x-\tfrac12,\\
p_2(x)&=x^2-3x+\tfrac34,\\
p_3(x)&=x^3-\tfrac{15}2x^2+\tfrac{45}4x-\tfrac{15}8.
\end{aligned}
\]

Os sinais da recomposição são importantes:

\[
m_k=\underbrace{\frac{2(-1)^k e^{a/4}}{(4\cdot400)^k}}_{\text{polo}}
+\underbrace{\frac1{\pi\cdot400^k}\int_0^\infty t^{2k}e^{-at^2}\Omega(t)\,dt}_{\text{arquimediano}}
+\underbrace{\frac{(-1)^{k+1}}{400^k}P^{(k)}(a)}_{\text{primos}}.
\]

**Não** substituir a soma completa dos zeros pelas primeiras ordenadas conhecidas para estabelecer a validade desta fórmula. Valores de zeros críticos são usados apenas como diagnóstico adicional.

## 2. Certificado da quadratura arquimediana

O código `es8m2_certify.py` reaproveita `es7_archimedean.py` e suas desigualdades de digama, pi, gamma, logaritmos e exp(-x), com todas as operações arredondadas externamente em inteiros (escala S=10²⁸). Usa ponto médio composto com M=16.384, T=48, 192 termos de digama, a=1/50.

Para f(t)=e^(−at²) Ω(t), ES7 fornece as cotas analíticas explícitas

\[
|f|\le e^{-at^2}(7+18t^2),\quad
|f'|\le e^{-at^2}\{36t+2at(7+18t^2)\},
\]
\[
|f''|\le e^{-at^2}\{36+180at^2+72a^2t^4+28a^2t^2+14a\}.
\]

Com gₖ(t)=t^(2k)f(t)/400ᵏ, para k≥1,

\[
|g_k''(t)|\le\frac{e^{-at^2}}{400^k}
\bigl(t^{2k}C(t)+4kt^{2k-1}D(t)+2k(2k-1)t^{2k-2}B(t)\bigr),
\]

onde B=7+18t², D=36t+2atB, C=36+180at²+72a²t⁴+28a²t²+14a. Para k=0, vale apenas o termo C. Em cada subintervalo [l,r], substituímos o fator exponencial por exp(−al²), avaliado por excesso; os polinômios não negativos são majorados em r. A soma das cotas com h³/24 dá erro **rigoroso** da quadratura do ponto médio.

Para a cauda em t≥T, usa-se o limite para |Ω(t)| e integrações por partes: I₀≤exp(−aT²)/(2aT) e I₂j≤exp(−aT²)T^(2j−1)/(2a)+(2j−1)I₂j−2/(2a). Dessa forma a cauda é englobada com aritmética dirigida, sem assumir localização de zeros.

## 3. Certificados dos termos dos primos após três derivadas

O código utiliza os intervalos racionalmente validados para ln n, π, √n e exp(−x) da implementação ES8-M corrigida (escala 10⁴⁰). A soma em n≤64 é uma soma de potências de primos com intervalos externos.

Para u=ln n≥ln 64>4 e k≤3, x=(25/2)u²≥200, vale |pₖ(x)|≤xᵏ. A função majorante de cada termo omitido é decrescente em n≥64. Após a mudança y=ln t, as desigualdades u/2≤u²/8 e 2k ln u≤u² implicam

\[
\sum_{n>64}|P^{(k)}_n(1/50)|
\le 5\,625^3\int_4^\infty u e^{-(91/8)u^2}\,du
=\frac{20\,625^3}{91}e^{-182}.
\]

Como e>27/10, o lado direito é menor que a fração racional `1/10^70` (comparação exata entre inteiros), aplicada **simetricamente** para cada ordem k=0,1,2,3. Os termos incluídos possuem arredondamentos próprios que já estão incorporados ao certificado.

Os intervalos resultantes para P⁽ᵏ⁾(1/50) são fornecidos em `RESULTADOS_ES8M2.json` sem substituir as frações pelos decimais exibidos.

## 4. Resultado efetivo — momentos completos, normalizados

| k | Limite inferior de mₖ | Limite superior de mₖ |
|:--:|--:|--:|
| 0 | 0,03695273542288 | 0,03721240109963 |
| 1 | 0,01866165694507 | 0,01874902372677 |
| 2 | 0,00952953721872 | 0,00956890792479 |
| 3 | 0,00499206970653 | 0,00501457478077 |

Valores decimais arredondados externamente; `RESULTADOS_ES8M2.json` conserva as frações exatas.

### Duas matrizes 2×2 completas

\[
H_1=\begin{pmatrix}m_0&m_1\\m_1&m_2\end{pmatrix},\qquad
H_1^{[1]}=\begin{pmatrix}m_1&m_2\\m_2&m_3\end{pmatrix}.
\]

Propagando intervalos sem supor independência estatística (aritmética intervalar inclusiva), os determinantes ficam:

| Matriz | Limite inferior det | Limite superior det | Conclusão |
|---|--:|--:|---|
| Hankel H₁ | **0,00000061657683** | 0,00000782459985 | Positiva definida |
| Stieltjes H₁^[1] | **0,00000159629343** | 0,00000320630195 | Positiva definida |

As diagonais também têm limites inferiores positivos; para uma matriz real simétrica 2×2, isso e o determinante positivo implicam **positividade definida**. Esses são resultados **finitos**, relativos à funcional linear gaussiana e um a fixo; **não** se estendem a todo n e todas as funções de teste, tampouco atestam os pivôs projetados E7/V24.

### Controle de refinamento obrigatório

| M | Determinante Hankel | Determinante Stieltjes |
|---:|---|---|
| 4.096 | Inconclusivo | Inconclusivo |
| 8.192 | Inconclusivo | Inconclusivo |
| 16.384 | Certificado positivo | Certificado positivo |

Os quatro momentos refinados permanecem dentro dos intervalos mais largos de M=8.192. A melhora é de **largura do intervalo**, não de dedução a partir dos zeros conhecidos.

### Conferência em ponto flutuante (NÃO certificado)

Os valores aproximados `mpmath` para m₀…m₃ são, respectivamente, `0.0370825834626389645`, `0.0187053407106170003`, `0.0095492225342477344`, `0.0050033221982893095`. Eles caem dentro dos intervalos rigorosos correspondentes. Esse diagnóstico não é usado na demonstração intervalar.

## 5. O que os dois CSVs realmente acrescentam

Os arquivos trazem 22 linhas bibliográficas, incluindo duplicações por DOI. Após eliminar duplicações, **19 referências distintas** foram registradas em `REFERENCIAS_TRIADAS_ES8M2.csv`, e os CSVs originais foram preservados em `sources/`.

**Pertinência direta ao contexto:**

1. **Jeffrey C. Lagarias**, *Li coefficients for automorphic L-functions*, *Annales de l'Institut Fourier* 57 (2007), DOI `10.5802/aif.2311`. Relaciona coeficientes de Li a uma forma quadrática de Weil em contexto de funções L automórficas; útil como comparação conceitual, não fornece nossa certificação numérica. Confirmado pela página editorial.
2. **Alain Connes**, *Trace formula in noncommutative geometry and the zeros of the Riemann zeta function*, *Selecta Mathematica* 5 (1999), DOI `10.1007/s000290050042`. Relaciona fórmulas explícitas e traços em espaços adélicos; não prova a identificação do espectro da ES4.
3. **M. V. Berry e J. P. Keating**, *The Riemann Zeros and Eigenvalue Asymptotics*, *SIAM Review* 41 (1999), DOI `10.1137/S0036144598347497`. Contexto para comparação de contagem espectral.
4. **Artigo sobre espaçamento de zeros**, DOI `10.1090/S0025-5718-1987-0866115-0`. Útil à avaliação espectral, mas não fundamenta os limites de integral desta etapa.

Os artigos sobre cordas, amplitudes adélicas e outros modelos quânticos foram mantidos como **contexto secundário ou não diretamente pertinente**; registros que proclamam provas de RH não foram aceitos sem auditoria. A aritmética certificada da presente etapa depende das fórmulas clássicas e das desigualdades implementadas, não desses resumos CSV.

## 6. Limites científicos e porta seguinte

- **Weil como critério de equivalência** RH ↔ positividade da forma em *todas* as funções admissíveis é estabelecido na literatura. A **positividade universal exigida pelo critério não foi demonstrada no NUMERIA**.
- Positividade dessas duas pequenas Hankel/Stieltjes não significa que uma medida de suporte real tenha sido identificada para os zeros complexos. A demonstração de recíproca via Carleman proposta no arquivo Qwen permanece não estabelecida.
- A gaussianização em parâmetros e a variável x=z²/400 não identificam as matrizes projetadas E7/V24: esse cotejo requer demonstrar a equivalência entre bases, projeções, polo, primes e arquimediano **entrada por entrada**.
- O passo ES8-M3 é construir os momentos μ₄ e μ₅, testar a positividade das duas matrizes 3×3 e obter estimativas de cauda mais apertadas, talvez com uma implementação independente em Arb/FLINT.
- Outra frente paralela é retornar ao alvo E7/V24 e certificar as **entradas dos produtos** F_i(r)overline{F_j(r)}, que são distintas dessas gaussianas.

## 7. Reprodutibilidade e fontes

- `es8m2_certify.py`, `test_es8m2_certify.py`, `RESULTADOS_ES8M2.json`, `PARAMETROS_ES8M2.json`, `MANIFESTO_SHA256.json`.
- `es7_archimedean.py`, `es6_weil_gaussian_baseline.py` são cópias versionadas das implementações corrigidas, não modificações da V24 nem de suas provas.
- Ver `README.md` para comandos, testes e como recomputar a execução com os parâmetros certificados.
- Fonte da normalização de Guinand: Paul Garrett, *Guinand's explicit formula*, notas (2021), `https://www-users.cse.umn.edu/~garrett/m/complex/notes_2020-21/12g_guinand_explicit_fml.pdf`; fonte complementar acessível: `https://openmaths.org/courses/NT-EXPF/NT-EXPF-03.html`.
- Condição de Stieltjes clássica: Schmüdgen, *The Stieltjes condition and multidimensional K-moment problems*, `https://link.springer.com/article/10.1007/s00013-021-01599-9`.
- Referências principais de CSV confirmadas via editores e site do autor: DOI acima.

**Conclusão:** ES8-M2 amplia a certificação para quatro momentos completos e duas matrizes 2×2, sem resolver a condição universal de Weil, sem identificar espectro da zeta e sem provar RH.
