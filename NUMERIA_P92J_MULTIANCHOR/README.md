# NUMERIA–RIEMANN P92-J — Compensação espectral para múltiplas âncoras e redução exata de Schur

**Base:** P92-I, [PR #98](https://github.com/jbarumby/scicalc/pull/98). Esta nota desenvolve **a trilha espectral da matriz Δ** da P88, sem afirmar uma prova geral dos Lemas 5–6 de Preissmann nem da Hipótese de Riemann.

## 1. Teorema J.1 — núcleo de dimensão arbitrária com folhas de uma escala

Sejam \(x_1<\dots<x_N\) pontos reais distintos, \(\delta_i=\min_{j\neq i}|x_i-x_j|\), e defina a matriz real simétrica de entradas não negativas

\[
\Delta_{ij}=\frac{\sqrt{\delta_i\delta_j}(\sqrt{\delta_i}-\sqrt{\delta_j})^2}{2(x_i-x_j)^2}\quad(i\neq j),\qquad \Delta_{ii}=0.
\]

Escolha **qualquer** partição dos índices em um conjunto núcleo \(U\), de cardinal \(r\), e um conjunto de folhas \(V\), de cardinal \(M\), sujeito apenas a \(\delta_j=h>0\) **igual para todo \(j\in V\)**. Não se exige que as folhas sejam índices consecutivos nem uniformemente espaçadas. Ordenando os índices de \(U\) antes dos de \(V\), temos exatamente

\[
\Delta=\begin{pmatrix}C&W\\W^{\mathsf T}&0_{M\times M}\end{pmatrix},\quad C=\Delta_{UU},\ W=\Delta_{UV}.
\tag{J1}
\]

Para qualquer número real \(t>0\), **equivalência exata de Schur**:

\[
\boxed{tI_N-\Delta\succ0\quad\Longleftrightarrow\quad
 tI_r-C-\frac1tWW^{\mathsf T}\succ0.}\tag{J2}
\]

A demonstração é a identidade de complemento de Schur em relação ao bloco \(tI_M\succ0\). Como Δ é simétrica e tem entradas não negativas, \(\|\Delta\|_2=\lambda_{\max}(\Delta)\): para qualquer vetor \(v\), \(|v^{\mathsf T}\Delta v|\le|v|^{\mathsf T}\Delta|v|\le\lambda_{\max}(\Delta)\|v\|_2^2\). Portanto **uma fatoração LDL exata da matriz reduzida de ordem \(r\) certifica \(\|\Delta\|_2<t\)**. O tamanho \(M\) das folhas não determina a ordem da fatoração; para dados racionais com raízes de distâncias racionais, o certificado pode ser verificado com `fractions.Fraction` sem floats.

**Limitação:** J2 é um *critério exato a posteriori* válido para qualquer \(r,M\), NÃO uma afirmação de que todas as configurações o satisfaçam para um \(t\) fixo. Em particular, configurações sem um conjunto útil de folhas com \(\delta\) comum não são resolvidas por J2.

## 2. Teorema J.2 — duas populações uniformes, âncoras ilimitadas

Sejam \(L\ge h>0\) reais e \(r\ge1\), \(M\ge2\) inteiros arbitrários. Tome **r âncoras**, espaçadas entre si por \(L\), seguidas por **M pontos densos** espaçados por \(h\):

\[
\begin{aligned}
x_i&=iL, &&0\le i\le r-1,\\
x_{r+j}&=rL+jh, &&0\le j\le M-1.
\end{aligned}\tag{J3}
\]

Logo \(\delta_i=L\) nas âncoras e \(\delta_{r+j}=h\) no grupo denso (\(M\ge2\) garante a última igualdade, inclusive no ponto inicial). Assim \(\Delta_{UU}=0\), \(\Delta_{VV}=0\), e Δ é bipartida:

\[
\Delta=\begin{pmatrix}0&W\\W^{\mathsf T}&0\end{pmatrix}.
\]

Defina \(\lambda=L/h\ge1\), \(s=\sqrt\lambda\) e indexe as âncoras por sua distância ao grupo (\(m=1,\ldots,r\)). Para \(j=0,\ldots,M-1\),

\[
\boxed{W_{m,j}=\frac{\sqrt\lambda(\sqrt\lambda-1)^2}{2(m\lambda+j)^2}.}\tag{J4}
\]

**Teorema J.2 (prova analítica para TODOS os parâmetros reais, e TODAS as dimensões):**

\[
\boxed{\|\Delta\|_2^2<\frac5{48},\qquad
\|\Delta\|_2<\sqrt{\frac5{48}}<\frac13.}\tag{J5}
\]

O número de âncoras \(r\) pode crescer sem limite; o tamanho do grupo \(M\) e o contraste \(L/h\) também. A cota de J5 não é declarada ótima.

### Demonstração analítica de J5

Pela forma bipartida \(\|\Delta\|_2=\|W\|_2\le\|W\|_{HS}\). Assim,

\[
\|\Delta\|_2^2\le\frac{\lambda(\sqrt\lambda-1)^4}{4}
\sum_{m=1}^{r}\sum_{j=0}^{M-1}(m\lambda+j)^{-4}.
\]

Usando que \(f(y)=(a+y)^{-4}\) decresce, para todo \(a>0\),

\[
\sum_{j=0}^{M-1}(a+j)^{-4}\le a^{-4}+\int_0^\infty(a+y)^{-4}\,dy
=a^{-4}+\frac{1}{3a^3}.
\]

Somando em todas as âncoras, e majorando a soma finita pela infinita,

\[
\sum_{m=1}^{r}\sum_{j=0}^{M-1}(m\lambda+j)^{-4}
\le\frac{\zeta(4)}{\lambda^4}+\frac{\zeta(3)}{3\lambda^3}.
\]

Com \(u=1/\sqrt\lambda\in(0,1]\), resulta

\[
\|\Delta\|_2^2\le\frac{(1-u)^4}{4}
\left(\zeta(4)u^2+\frac{\zeta(3)}3\right).
\]

A comparação elementar \(\zeta(4)<\zeta(3)<5/4\) decorre de \(\sum_{n\ge3}n^{-3}<\int_2^\infty y^{-3}\,dy=1/8\), portanto

\[
\|\Delta\|_2^2<\frac5{48}(1-u)^4(1+3u^2)\le\frac5{48}.
\]

A última desigualdade é exata: \(f(u)=(1-u)^4(1+3u^2)\) é decrescente em \([0,1]\), pois

\[
f'(u)=-(1-u)^3(4-6u+18u^2)\le0,
\quad 4-6u+18u^2=18(u-1/6)^2+7/2>0.
\]

Por fim \(5/48<1/9\), provando J5. Nenhuma etapa exige que as âncoras ou os pontos densos sejam poucos, nem invoca resultados gerais de Preissmann.

## 3. J.3 — contraexemplo exato à suficiência da soma de linha

Tome \(h=1,\ L=256,\ r=4,\ M=280\), ou seja, **284 pontos**. A linha correspondente à última âncora tem entradas positivas

\[
w_j=\frac{1800}{(256+j)^2},\quad j=0,\ldots,279.
\]

Como \(y\mapsto(256+y)^{-2}\) decresce,

\[
\sum_{j=0}^{279}w_j>
1800\int_0^{280}\frac{dy}{(256+y)^2}
=\frac{7875}{2144}
>\frac{2725}{768}.
\tag{J6}
\]

A última fração \(2725/768\) é o **limiar inferior racional rigoroso** usado na P92-I para a cota residual requerida pela P88. O critério ingênuo da soma máxima de linha, portanto, **não certifica** esse exemplo. Contudo, J5 certifica \(\|\Delta\|_2<1/3\), com margem enorme. A matriz tem 284 índices, mas J2 reduz um certificado LDL de dimensão 284 a dimensão **quatro**, sem tolerância ou autovalores float.

A falha de um teste suficiente não constitui falha da desigualdade espectral: precisamente o oposto é demonstrado aqui.

## 4. Verificador e certificado racional independente

- [certificate.py](certificate.py) reproduz a matriz original \(\Delta\) diretamente a partir das coordenadas e das distâncias aos vizinhos e compara suas entradas com J4;
- reconstrói 43 configurações com \(s=\sqrt{L/h}\) racional, incluindo o exemplo 4+280;
- certifica \(\|W\|_{HS}^2\) por frações, frente à envoltória racional \(\frac5{48}(1-1/s)^4(1+3/s^2)\);
- usa **43 fatorações LDL racionais** de \((1/9)I_r-WW^{\mathsf T}\), ou equivalente à redução J2 em \(t=1/3\);
- para 28 casos menores, faz **duas fatorações LDL independentes da matriz integral**, de \((1/3)I\pm\Delta\) — 56 verificações completas;
- preserva os valores das somas de linha como frações exatas no JSON e demonstra J6 por comparação integral **racional**, não por aproximação decimal.

Reproduzir na raiz do repositório:

    python -m unittest -v NUMERIA_P92J_MULTIANCHOR.test_p92j
    python -O -m unittest -v NUMERIA_P92J_MULTIANCHOR.test_p92j
    python -m NUMERIA_P92J_MULTIANCHOR.certificate --out p92j-multianchor.json

Os 43 casos e os testes LDL constituem regressões. Os Teoremas J.1 e J.2 são demonstrados **analiticamente** nas seções 1–2, independentemente dos testes.

## 5. Limite e próximo alvo P92-K

**Demonstrado:** redução exata para configurações com muitas folhas de mesma escala e uma cota uniforme \(<1/3\) para a família J3, com dois grupos e ambas as cardinalidades arbitrárias; exemplo exato em que a soma de linha é inadequada.

**Não demonstrado:** uma partição canônica de configurações arbitrárias em grupos desse tipo mantendo a cota residual; controle universal de \(\|\Delta\|_2\) para todos os pontos; Lema 5 ou Lema 6 gerais de Preissmann; extensão da P93 geométrica de cinco lacunas; Hipótese de Riemann.

**Prioridade P92-K:** permitir **múltiplos grupos densos de escalas diferentes** intercalados entre âncoras, usando complemento de Schur com blocos de folha de norma pequena. O desafio permanece obter uma estimativa que não cresça com o número de grupos sem hipóteses injustificadas de separação.

Manter a PR em rascunho, sem merge automático.
