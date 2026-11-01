# NUMERIA–RIEMANN — E8-K V17

## Especialização da fórmula explícita de Guinand–Weil às funções seno do NUMERIA

**Situação:** conferência das normalizações e dos resíduos completada; **fórmula explícita usada como teorema publicado**, não provada novamente; inclusões numéricas arquimedianas da Fase4 continuam condicionais. Não é demonstração da Hipótese de Riemann.

### 1. Fonte do teorema e hipóteses

Empregamos a formulação publicada em Chirre–Gonçalves et al., *Fourier Interpolation with Zeros of Zeta and L-Functions*, *Constructive Approximation* (2022), equação (1.1):
https://link.springer.com/article/10.1007/s00365-022-09599-w

O teorema afirma, sob hipóteses de holomorfia em uma faixa horizontal e decaimento polinomial uniforme nela, uma igualdade entre: soma dos zeros não triviais; avaliações nos pontos imaginários ±i/2; integral do termo gama/digama; e soma de von Mangoldt sobre potências primas. Ele vale **sem pressupor RH**. Usamos a convenção

\[F_k(r)=\int_{-L}^{L}f_k(t)e^{irt}\,dt,\qquad C_{ij}(u)=\frac1{2\pi}\int_{\mathbb R}H_{ij}(r)e^{-iru}\,dr.\]

A normalização de Fourier é indispensável; converter uma convenção com `2π` para outra sem ajustar o argumento de Fourier alteraria o fator do termo primo.

Referência adicional: E. Bombieri, *Remarks on Weil’s quadratic functional in the theory of prime numbers, I* (2000), seções 2 e 12: https://www.bdim.eu/item?id=RLIN_2000_9_11_3_183_0 .

### 2. As funções do NUMERIA satisfazem a classe do teorema

Para `L=9/20` e inteiro `k>=1`, defina

\[f_k(t)=\mathbf{1}_{[-L,L]}(t)\frac{1}{\sqrt L}\sin\left(\frac{k\pi(t+L)}{2L}\right),\quad w_k=\frac{k\pi}{2L}.\]

As funções são reais, compactamente suportadas e nulas em ambos os extremos; são `C⁰` e suaves no interior, mas **não são C∞ na borda**. A fórmula explicitada por Chirre–Gonçalves, porém, é enunciada em termos das **propriedades da transformada**, e não exige `f_k∈C_c^∞`.

A integral que define `F_k(z)` é inteira em `z`. Fazendo duas integrações por partes, usando `f_k(±L)=0`, os termos de fronteira da segunda integração são finitos, e obtemos `F_k(x+iy)=O_y(|x|^-2)` para `y` em qualquer intervalo compacto. Logo

\[H_{ij}(z)=\tfrac12[F_i(z)F_j(-z)+F_j(z)F_i(-z)]\]

é par, inteira e `O(|Re z|^-4)` uniformemente em toda faixa horizontal limitada. Portanto, as condições suficientes de holomorfia e decaimento exigidas pela fórmula explícita estão satisfeitas. Isso trata a questão de admissibilidade que ainda estava em aberto na V16, **tomando como premissa o teorema publicado**.

### 3. Transformada exata e amplitude arquimediana

Integração direta de `sin(w_k y)e^{iry}` dá, inicialmente para `r≠±w_k`,

\[F_k(r)=\frac{w_k}{\sqrt L}\frac{e^{-irL}-(-1)^ke^{irL}}{w_k^2-r^2}.\]

As singularidades em `±w_k` são removíveis, por sua representação integral inteira. Para `r` real:

\[|F_k(r)|^2=\frac{k^2\pi^2}{L^3(r^2-w_k^2)^2}\begin{cases}\cos^2(Lr),&k\text{ ímpar},\\\sin^2(Lr),&k\text{ par}.\end{cases}\]

Perto de `r=w_k`, colocando `x=L(r-w_k)` e `sinc(x)=sin(x)/x`, resulta a forma regularizada

\[|F_k(r)|^2=\frac{4Lw_k^2\operatorname{sinc}^2(L(r-w_k))}{(r+w_k)^2},\qquad |F_k(w_k)|^2=L.\]

As fórmulas conferem com os dois ramos de `amp2` do B5 no domínio `r≥0`. Como `H_kk=|F_k|²` é par e `Ω(r)=Re ψ(1/4+ir/2)−log π` é par,

\[\frac1{2\pi}\int_{-\infty}^\infty H_{kk}(r)\Omega(r)dr = \frac1\pi\int_0^\infty |F_k(r)|^2\Omega(r)dr.\]

O corpo AST de `eval_point` foi inspecionado: `amp2(r,k)*omega(r)/PI`. A verificação da forma **não substitui** a prova independente de todas as cotas de truncamento de `amp2`/`omega`.

### 4. Termos primos e paridade

Para `u≥0`, a inversa de Fourier de `H_ij` é a autocorrelação simétrica

\[C_{ij}(u)=\tfrac12\left[\int f_i(t+u)f_j(t)dt + \int f_j(t+u)f_i(t)dt\right].\]

Como `supp f_k⊆[-L,L]`, temos `C_ij(u)=0` se `|u|≥2L`. O programa produz limites racionais exatos, pela série de `atanh`, demonstrando

\[\log 2 < 2L=9/10 < \log3.\]

Portanto, somente `n=2` contribui à soma de potências primas, com a convenção do teorema:

\[P_{ij}=-\frac{2\log2}{\sqrt2}\,C_{ij}(\log2).\]

Pela paridade centrada, `f_k(-t)=(-1)^{k+1}f_k(t)`; se `i-j` é ímpar, `H_ij=0`, e assim `C_ij=0` para todo `u`. Além disso, `C_ij(0)=δ_ij` pela ortonormalidade das funções seno.

### 5. Momentos dos polos: duas condições anuladas

Considere `M_{k,s}=∫_{-L}^L f_k(t)e^{s t/2}dt`, `s∈{−1,+1}`. Cálculo por integração elementar fornece

\[M_{k,s}=\frac{2\pi\sqrt L\,k}{k^2\pi^2+L^2}\left(e^{-sL/2}-(-1)^k e^{sL/2}\right).\]

O termo entre parênteses depende apenas da **paridade** de `k`, para cada sinal `s`. Consequentemente, para `p=1` (k ímpar) ou `p=2` (k par), o coeficiente

\[\beta_{p,k}=-\frac{k}{p}\frac{p^2\pi^2+L^2}{k^2\pi^2+L^2}\]

faz `M_{k,+}+β_{p,k}M_{p,+}=0` **e** `M_{k,-}+β_{p,k}M_{p,-}=0`, simultaneamente. O verificador checa os dois coeficientes da identidade polinomial em `X=π²`, de forma racional exata, para `k=3,…,10` (oito colunas). Como os dois momentos desaparecem em qualquer função projetada, as avaliações `H(i/2)` e `H(-i/2)` também desaparecem.

Este passo demonstra a **relação entre a projeção B e a anulação dos polos**, o elo estrutural pendente da V16. Não implica positividade para dimensão arbitrária.

### 6. Execução V17 e evidências

- Origens E7, Fase4, B5, Fase3 e relatório/código V16 verificados contra SHA-256 fixados.
- Todos os 30 registros da Fase4 reconferidos por hashes individuais; não se reutiliza a palavra `CONCLUIDO` como certificado analítico.
- 9/9 controles estruturais e 6/6 mutações negativas isoladas aprovados.
- 8/8 cancelamentos de polos demonstrados por coeficientes exatos.
- 25 verificações numéricas independentes de transformada, momentos e ortogonalidade aprovadas, com erro máximo observado de aproximadamente `2,86e−13`; **essas 25 verificações são apenas testes auxiliares, não majorantes rigorosos**.
- As matrizes N=6,8,10 positivas condicionalmente na V16 permanecem **herdadas**, não são recalculadas aqui.

### 7. Pendências científicas

1. Demonstrar **sem dependências circulares** todas as cotas de erro nodal da Fase4, o erro global de interpolação e a cauda arquimediana além de R=160. Sem isso, as conclusões sobre os sinais da matriz são condicionais.
2. Formalizar os passos de admissibilidade e integração em Lean 4, fixando versão Mathlib e registrando `lake build` bem-sucedido. A formalização herdada ainda não foi aceita pelo compilador nesta execução.
3. Estender os resultados dos três N testados e da família seno a todo o conjunto admissível necessário ao critério de Weil, se isso for matematicamente possível.
4. Uma auditoria alternativa das normalizações pode ser feita a partir da fórmula de Bombieri, mas **este relatório utiliza a fórmula de Guinand–Weil publicada como teorema externo**, não a demonstra por contorno/resíduos do zero.

**Classificação:** `FOURIER_AND_POLE_SPECIALIZATION_AUDITED_FULL_EXPLICIT_FORMULA_EXTERNAL`.

**RH permanece não demonstrada.**