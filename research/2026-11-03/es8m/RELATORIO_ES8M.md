# NUMERIA–ES8-M — Auditoria de momentos gaussianos e positividade de Stieltjes

**03/11/2026 — relatório científico reproduzível.**

## 1. Resposta sobre Weil e escopo

O **critério de positividade de Weil** (equivalência entre RH e uma positividade universal para uma classe adequada de autocorrelações) é um teorema conhecido. Não está provada a positividade universal exigida nesse critério e, portanto, **a Hipótese de Riemann não foi demonstrada**. O NUMERIA possui certificações finitas com hipóteses e normalizações específicas; nenhuma delas implica a positividade universal. A ES8-M não deve ser apresentada como comprovação de Weil.

Referência: Jeffrey Lagarias, *The Riemann Hypothesis: Arithmetic and Geometry*, Teorema 3.2, https://dept.math.lsa.umich.edu/~lagarias/doc/mt-holyoke-rev.pdf ; Clay Mathematics Institute, https://www.claymath.org/millennium-problems/ .

## 2. Extração crítica dos dois textos recebidos

- DeepSeek, `sources/deepseek_latex_20261103_74b66a.tex`: propõe momentos gaussianos de Guinand–Weil, Hankel e certificação racional. O corolário **RH ⇒ Hankel PSD** está correto. A declaração de que a recíproca é falsa não foi provada no manuscrito; dificuldades de uma verificação finita não provam a falsidade de uma equivalência matemática infinita. A ES7 anterior não tinha implementado matrizes de Hankel, ao contrário do que sugere o enquadramento do texto.
- Qwen, `sources/Qwen_latex_20261103_pag4a82de.txt`: propõe reforçar Hankel com a condição deslocada de Stieltjes, momentos altos e funções-teste polinomiais. O critério de Stieltjes é clássico e correto para sequências reais de momentos. O manuscrito não demonstra adequadamente a conclusão de que estas condições particulares da zeta são equivalentes à RH; **Carleman sozinho não faz essa identificação**. Suas expressões com `Ω^(j)` não resultam da mera diferenciação da gaussiana, e a afirmação de que `F=G/ζ` tem zeros nos zeros de ζ está invertida (na ausência de cancelamentos, esses pontos são **polos**).
- Ambos os textos são preservados integralmente por rastreabilidade e são **fontes de hipóteses**, não teoremas do projeto.

Referência matemática da caracterização de Stieltjes: https://www.degruyterbrill.com/document/doi/10.1515/math-2022-0001/html , Teorema 4.2; condicional de Carleman: Teorema 4.4.

## 3. Fórmula usada, diferenciação e duas famílias de matrizes

Para `a>0` e `z_ρ=(ρ−1/2)/i`, definimos

\[
F(a)=\sum_{\rho}e^{-az_\rho^2}=2e^{a/4}+A(a)-P(a),
\quad A(a)=\frac1\pi\int_0^\infty e^{-at^2}\Omega(t)dt,
\]
\[
P(a)=\frac1{\sqrt{\pi a}}\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}
  e^{-(\log n)^2/(4a)},\quad
\Omega(t)=\operatorname{Re}\psi(\tfrac14+it/2)-\log\pi.
\]

Essa é a convenção ES6/ES7 da fórmula de Guinand–Weil. Os zeros são **todos** os não triviais, sem pressupor RH. A soma é absolutamente convergente para `a>0` (as partes imaginárias de `z_ρ` são limitadas na faixa crítica e as ordenadas crescem). Para todo `k` fixo, a diferenciação local em `a>0` é legitimada pelo decaimento gaussiano e pelas cotas clássicas de contagem de zeros.

\[
\mu_k(a)=(-1)^k F^{(k)}(a)=\sum_\rho z_\rho^{2k}e^{-az_\rho^2},
\quad\mu_0=2e^{a/4}+A-P,
\quad\mu_1=P'(a)-A'(a)-\tfrac12e^{a/4}.
\]

Se RH for verdadeira, cada `x_ρ=z_ρ²` pertence a `[0,∞)` e temos a medida de Stieltjes positiva `ν_a=Σ_ρ e^(−a x_ρ)δ_(x_ρ)`. Nesse caso, para **qualquer** `d`:

\[
H_d(a)=(\mu_{i+j})_{0\le i,j\le d}\succeq0,
\quad H_d^{[1]}(a)=(\mu_{i+j+1})_{0\le i,j\le d}\succeq0.
\]

O deslocamento satisfaz `dH_d/da = −H_d^{[1]}`; a derivada **não** é a matriz deslocada com mesmo sinal. Não inferimos recíproca, unicidade da medida associada aos zeros nem equivalência com RH nesta auditoria.

## 4. Certificado racional novo para `P'(1/50)`

Termo a termo no intervalo de convergência:

\[
P'(a)=\sum_{n\ge2}P_n(a)\left(-\frac1{2a}+\frac{(\log n)^2}{4a^2}\right).
\]

No ponto `a=1/50`, o fator é `c_n=−25+625(log n)^2`, estritamente positivo para todo `n≥2`, pois `log 2 > 1/2` (cota racional disponível no módulo). Os logaritmos, a raiz quadrada de `π`, os exponenciais e todos os produtos de `n≤64` são cercados com inteiros escalonados, frações exatas e arredondamento externo.

A cauda é controlada **sem RH**: para `x=log t≥log 64>4`, usamos `Λ(t)≤log t` (como majorante pontual para potências primas), e as desigualdades

\[
x^2e^{-12.5x^2}\le e^{-6.25x^2},\qquad e^{x/2}\le e^{x^2/8}\quad(x\ge4).
\]

A função integranda majorante é decrescente a partir de 64. Após a substituição `x=log t`, obtemos

\[
0\le P'(1/50)-P'_{64}(1/50)
\le\frac{625}{\sqrt{\pi/50}}\int_4^\infty xe^{-(49/8)x^2}dx
\le\frac{12500}{49}e^{-98}<10^{-37}.
\]

A última desigualdade é verificada por cotas racionais `exp(-98)<10^-40`; `1/sqrt(π/50)<5`. Os extremos íntegros e a fórmula de cauda estão no JSON. **Esse certificado vale para `P'`, não para uma integral arquimediana nova.**

### Correção de proveniência importante

A versão original ES6 copiada como insumo continha uma escolha invertida do lado do resto alternado em `arctan_alt` para `m` par. A ES8-M corrigiu `return (s,s+term) if m%2 == 0 else (s-term,s)` em sua **cópia isolada**. Verificações de alta precisão confirmam o enquadramento de `arctan(1/5)`, `arctan(1/239)` e `π`; nenhum código antigo do repositório foi sobrescrito. A largura do arredondamento escalonado em ES6 era mais larga que o erro de Taylor, mas isso não autoriza ignorar o erro lógico. Reavaliação e comparação numérica das novas cotas foram executadas.

## 5. Resultados realmente certificados

Para `a=1/50`, combinando `A` e `A'` com os **extremos inteiros** da ES7, `P` e `P'` com cotas racionais da ES6 corrigida/ES8-M e `exp(a/4)` com a expansão racional, os intervalos são:

\[
0.03656551942685 < \mu_0(1/50) < 0.03759961709565,
\]
\[
7.41235693525210 < \mu_1(1/50) < 7.55191533348501.
\]

As desigualdades decimais exibidas são **afrouxamentos** dos limites racionais do arquivo `RESULTADOS_ES8M.json`. A positividade é rigorosa **condicionada à validade dos lemas analíticos e dos certificados ES6/ES7 reusados**, não uma certificação de um novo quadrado universal de Weil.

Assim, apenas os dois testes de dimensão `1×1` para as matrizes **completas** `H_0=[μ0]` e `H_0^[1]=[μ1]` são positivos certificados.

## 6. Matrizes-modelo exatas e diagnósticos reais não certificados

Na **bancada artificial** racional (nós `1,2,3,4,5` e pesos `1,2,3,5,8`), os blocos de Hankel e deslocados com dimensões 1×1, 2×2, 3×3 e 4×4 foram checados por LDL em frações exatas. Todos os pivôs são estritamente positivos. Isso testa o programa e a condição clássica de Stieltjes, **não a função zeta**.

Como controle adversarial, os nós complexos `1+i` e `1−i`, com pesos unitários, dão momentos `2,2,0,−4`; a Hankel `[[2,2],[2,0]]` tem pivôs LDL `2,−2` e é indefinida, demonstração exata **somente neste exemplo artificial**.

Na bancada numérica, foram usados 32 zeros conhecidos da **reta crítica** calculados por `mpmath` com 45 dígitos e a variável reescalada `y=γ²/400`. Para `a=0.02`, os menores autovalores de `H_d` e `H_d^[1]` (d=0,1,2) são positivos. Isso é praticamente tautológico: esse subconjunto real gera uma medida positiva por construção. Os valores não são certificados por Arb e **não verificam ausência de zeros fora da linha crítica**.

Momentos até `μ5` foram calculados numericamente também pelo lado analítico da fórmula de Guinand–Weil usando a identidade de Laguerre

\[
(-\partial_a)^k\{a^{-1/2}e^{-b/a}\}
=k!a^{-k-1/2}e^{-b/a}L_k^{-1/2}(b/a).
\]

A concordância numérica com a soma finita dos zeros conhecidos é um teste de implementação, não uma certificação de `μ2..μ5` nem de matrices completas de dimensão 2 ou 3.

## 7. Gates para ES8-M2 / E7–V24

1. **Aberto:** certificar `μ2` e `μ3` com quadratura intervalar das integrais `t^4 Ω(t)` e `t^6 Ω(t)` e dos termos primos derivados, para testar `H_1` e `H_1^[1]` reais, por LDL intervalar.
2. **Aberto:** obter `μ4..μ7` e trabalhar com monômios ortogonalizados e escalados, para combater o condicionamento de Hankel 3×3 e 4×4.
3. **Aberto:** provar (não apenas supor) a passagem de qualquer resultado de momentos gaussianos para produtos de transformadas das funções compactamente suportadas da V24, incluindo projeção `B* (G+P) B`.
4. **Aberto:** um argumento global que estabeleça positividade universal para a classe admissível de Weil, sem circularidade em RH. Este último é a dificuldade fundamental.

**Status final:** `CERTIFIED_MU0_MU1_ONLY`; `EXACT_RATIONAL_TOY_MATRICES`; `UNVERIFIED_HIGHER_ZETA_HANKEL`; `NO_RH_PROOF`.
