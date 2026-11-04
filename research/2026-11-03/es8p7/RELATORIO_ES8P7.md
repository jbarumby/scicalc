# NUMERIA–RIEMANN — ES8-P7: auditoria nodal e admissibilidade da fórmula de Weil

Data: 03/11/2026. **Status: auditoria amostral; não prova RH.**

## 1. Origem e objetivo

A ES8-P6 rederivou majorantes de Cauchy, pesos de interpolação e cauda da Fase 4, mas avaliou apenas nove controles nodais. A ES8-P7 combina uma nova implementação em `fractions.Fraction` com os avaliadores históricos do ZIP V16/E7 recuperado, preservados sem alterações. O objetivo é detectar divergências de sinais, fases, amplitudes e digama nos nós **efetivamente usados na quadratura**.

A Fase 4 contém **30 pares de mesma paridade**, 1.600 células por par e 16 nós por célula: **768.000 avaliações originais**. Não foi recalculado todo esse conjunto. A ES8-P7 executa uma **amostra determinística de 15 abscissas × 30 pares = 450 verificações intervalares**, nas células 0, 40, 100, 400 e 1599, com índices dos nós 0, 8 e 15 de cada célula. Todos os valores dos nós são frações exatas `r=c/10+j/150`.

## 2. Método independente para os nós

Escreva `L=9/20`, `w_k=kπ/(2L)` e, para os modos senoidais, `A_k(r)=sqrt(L)[sinc(L(r−w_k))−(−1)^k sinc(L(r+w_k))]`. A parcela real da entrada entre modos de mesma paridade tem a forma

\[
q_{ij}(r)=\frac{s_{ij}}{\pi}A_i(r)A_j(r)\Omega(r),\quad
s_{ij}=(-1)^{|i-j|/2},\quad
\Omega(r)=\Re\psi(1/4+ir/2)-\log\pi.
\]

Pares de paridades opostas têm parcela real exatamente nula, por simetria.

A implementação independente utiliza `Fraction` para todo cálculo certificado, com:

1. `π` por Machin e séries alternadas de arctangente, com orientação correta do primeiro termo omitido; `log` por série de `atanh` e resto geométrico; `γ` por Euler–Maclaurin de `H_256−log 256` com correção dos termos em `n^-2` e `n^-4`; `sqrt L` por `isqrt` com arredondamento externo;
2. para `b=r/2` e `u=N+1/4`, soma exata de `g(n+1/4)=b²/[u(u²+b²)]` até `N−1`, com arredondamento racional externo das parcelas à escala `10^34`; a cauda satisfaz
   \[
   \frac12\log(1+b^2/u^2)+\frac{g(u)}2 \le \sum_{n=N}^{\infty}g(n+1/4) \le
   \frac12\log(1+b^2/u^2)+\frac{g(u)}2-\frac{g'(u)}{12}.
   \]
   A desigualdade decorre da representação positiva `g(u)=∫₀∞e^{-ut}(1−cos bt)dt` e de `1/t+1/2≤1/(1−e^{-t})≤1/t+1/2+t/12`;
3. `sinc` por Taylor de `sin` com resto limitado pela primeira derivada omitida e redução modular a `2π`; todos os arredondamentos do centro racional para `10^-18` são compensados explicitamente usando `|sin'|≤1` e `|sinc'|≤1/2`. Para `x` contendo zero, usa a série de `sinc` obtida da integral de `cos(tx)`.

Esses intervalos foram confrontados com a implementação histórica de inteiros fixos (40 dígitos, 64 termos), autenticamente copiada de Fase 4. A implementação histórica **não usa** o mesmo código desta etapa; ambas compartilham identidades clássicas de funções especiais, portanto concordância não é prova formal independente de todos os lemas.

A avaliação `mpmath` de `Re digamma` a 65 dígitos constitui diagnóstico numérico, não *enclosure* rigoroso; seus resultados são usados apenas como controle de plausibilidade.

## 3. Critério de admissibilidade das funções-teste

As funções da base histórica são

\[
\phi_k(x)=\frac{1}{\sqrt L}\sin\!\left(\frac{k\pi(x+L)}{2L}\right)\mathbf1_{[-L,L]}(x),
\quad F_k(z)=\int_{-L}^L\phi_k(x)e^{izx}dx.
\]

As funções `φ_k` são contínuas, compactamente suportadas e se anulam em `±L`. Cada `F_k` é **inteira**, por diferenciação sob uma integral compacta. A extensão holomorfa correta da entrada real simetrizada é

\[
H_{ij}(z)=\tfrac12\bigl(F_i(z)F_j(-z)+F_j(z)F_i(-z)\bigr).
\]

Ela é **inteira e par**. Não se deve utilizar `Re(F_i(z)conj(F_j(z)))` para `z` complexo: essa expressão não é holomorfa.

Fixada qualquer `η>1/2`, usando que `φ_k(±L)=0`, duas integrações por partes no interior de `[-L,L]` dão, para `z≠0`,

\[
|F_k(z)|\le\frac{e^{Lη}}{|z|^2}\left(
\frac{2\omega_k}{\sqrt L}+\frac{2L\omega_k^2}{\sqrt L}\right),\qquad |\Im z|\le η.
\]

O limite trivial é `|F_k(z)|≤2sqrt(L) exp(Lη)`. Combinando as duas desigualdades, uniformemente na faixa `|Im z|≤η`, segue `F_k(x+iy)=O((1+|x|)^{-2})`, portanto **`H_ij=O((1+|x|)^-4)`**. Isso satisfaz a classe clássica de Guinand–Weil: holomorfia para `|Im z|≤1/2+ε` e decaimento pelo menos `(1+|Re z|)^-(1+δ)`. A formulação publicada por Carneiro e Milinovich é utilizada como teorema **externo**; não é provada novamente neste pacote.

Ainda por Fubini, a transformada inversa de `H_ij` é uma combinação de convoluções compactamente suportadas em `[-2L,2L]=[-0,9;0,9]`. Entre as potências de primos `log n` com `n≥2`, somente `log 2` pode estar nesse suporte. A classe suavizada da V24 (`χ_δ φ_k`, `χ_δ∈C∞_c`) possui decaimento ainda mais rápido; para a fórmula explícita, não é necessário exigir que a base não suavizada já seja `C∞` porque seus produtos espectrais satisfazem a hipótese de holomorfia e decaimento.

**Atenção:** a função arquimediana complexa `Ω_C(z)` apresenta polos em `z=±i/2`; ela é usada no integrando da quadratura em discos de raio menor que `1/2`, mas não é a função-teste `H_ij` da fórmula explícita. A admissibilidade de `H_ij` não torna `H_ij Ω_C` inteira.

## 4. Artigos incorporados

- **K. Boyadzhiev (2007)**, *Evaluation of series with Hurwitz and Lerch zeta function coefficients by using Hankel contour integrals*, DOI `10.1016/j.amc.2006.08.061`: fornece representações integrais alternativas de digama/log-gama e será útil para uma terceira implementação. **Não é uma prova pronta de inclusões intervalares**, pois é necessário derivar majorantes de contorno e arredondamento.
- **C. Ferreira e J. López (2004)**, *Asymptotic expansions of the Hurwitz–Lerch zeta function*, DOI `10.1016/j.jmaa.2004.05.040`: relevante para regimes assintóticos com parâmetros complexos; as restrições de argumento e do corte em `z=1` exigem cuidado antes de especializar a digama.
- **G. Nemes (2017)**, *Error bounds for the asymptotic expansion of the Hurwitz zeta function*, DOI `10.1098/rspa.2017.0363`: **fonte externa adicional**, fora do CSV recebido, especialmente pertinente a restos explícitos para poligama em uma próxima auditoria.
- As outras cinco referências do inventário (incluindo Nemes) são mantidas no JSON para rastreabilidade, sem inferir aplicação direta da ES8-P7 onde não foi demonstrada.

Referência de admissibilidade: Carneiro e Milinovich, lema da fórmula explícita de Guinand–Weil, https://par.nsf.gov/servlets/purl/10637681, §2.2. Referências de funções especiais: https://dlmf.nist.gov/5.7 e https://dlmf.nist.gov/5.5.

## 5. Limites e próximo gate

Esta fase **não executa as 768 mil avaliações** e não substitui todos os limites de erro da quadratura da Fase 4; tampouco formaliza os lemas em Lean. Não afirma a positividade para todas as funções-teste, todos os modos ou todas as dimensões. A identidade com as matrizes históricas foi examinada na ES8-P5 no nível de compatibilidade intervalar, não por igualdade funcional universal.

**ES8-P8:** executar o replay nodal exaustivo em lotes verificáveis; conferir cada nó e rastrear o erro certificado total, especialmente os modos da fronteira e os cancelamentos, com proveniência independente. A positividade universal de Weil e RH permanecem questões abertas.
