# NUMERIA–RIEMANN P92-L — Estimativa espectral adaptativa e separação sublinear entre agrupamentos

**Continuação da P92-K, [PR #100](https://github.com/jbarumby/scicalc/pull/100).**

Esta fase trabalha exclusivamente com a matriz simétrica não negativa `Delta` definida na P88. Demonstra novas **condições suficientes**, com provas para todos os parâmetros reais das classes indicadas. **Não** demonstra a cota espectral global para todos os conjuntos de pontos, o Lema 6 geral de Preissmann ou RH.

## L.1 — Lema de cardinalidade de dois blocos

Sejam K>=1 agrupamentos ordenados de pontos. O agrupamento r tem m_r>=2 pontos igualmente espaçados por h_r>0, e a distância entre as extremidades de dois agrupamentos consecutivos é pelo menos o máximo dos seus passos internos. Consequentemente, todos os pontos do agrupamento r têm distância ao vizinho mais próximo `delta_i=h_r`. Denote por D_rs a distância entre a extremidade direita do grupo r e a esquerda do grupo s (r<s), e ponha `H_rs=max(h_r,h_s)`, `rho_rs=D_rs/H_rs`.

Entre dois grupos de passos a,b, o bloco da matriz P88 é

\[
 W_{jk}=\frac{\beta}{(D+ja+kb)^2},\quad
 \beta=\frac{\sqrt{ab}(\sqrt a-\sqrt b)^2}{2},
 \qquad 0\le j<m_r,\ 0\le k<m_s,
\]

com índices contados a partir das extremidades enfrentadas. Para `H=max(a,b)` e `u=sqrt(min(a,b)/H)` em [0,1],

\[
 \frac{\beta}{H^2}=\frac{u(1-u)^2}{2}\le\frac{2}{27}
\]

porque `u(1-u)^2` atinge máximo `4/27` em u=1/3. Como `D+ja+kb >= D`, cada entrada satisfaz `0<=W_jk <=2 H²/(27D²)` e

\[
\boxed{\|W_{rs}\|_2 \le \|W_{rs}\|_{HS}
\le\frac{2\sqrt{m_rm_s}}{27\rho_{rs}^2}
\le\frac{m_r+m_s}{27\rho_{rs}^2}.}\tag{L1}
\]

Se h_r=h_s, o bloco W_rs é exatamente zero. **L1 independe das magnitudes absolutas das escalas.**

A P92-K provou, independentemente da cardinalidade, a segunda cota

\[
\|W_{rs}\|_2\le\frac{1}{2\sqrt2\,\rho_{rs}}
\qquad(\rho_{rs}\ge1).\tag{K-pair}
\]

Como `5/14>1/(2sqrt2)` (verificar quadrados `25/196 > 1/8`), obtém-se um majorante **racional** seguro, adequado para verificação exata:

\[
 A_{rs}=\begin{cases}
 0,&r=s\text{ ou }h_r=h_s,\\
 \displaystyle\min\left\{\frac5{14\rho_{rs}},\frac{m_r+m_s}{27\rho_{rs}^2}\right\},&r\ne s,
 \end{cases}\tag{L2}
\]

onde D_rs é interpretada com a ordem adequada dos grupos. A matriz A é simétrica e não negativa, e cada entrada domina a norma do bloco correspondente de Delta.

## L.2 — Teorema para separação sublinear e número ilimitado de agrupamentos

**Teorema P92-L.1.** Suponha que `2<=m_r<=M` e, para algum c>=1 e p>1/2,

\[
\boxed{D_{rs}\ge c(s-r)^p\max(h_r,h_s),\qquad r<s.}\tag{L-sep}
\]

Então, para **todo** K e todos os passos h_r>0 reais sob essas condições,

\[
\boxed{\|\Delta\|_2\le\frac{4M}{27c^2}\,\zeta(2p).}\tag{L-global}
\]

*Prova.* Os blocos diagonais de Delta são nulos. Pela L1, para r≠s,

\[
 \|W_{rs}\|_2\le\frac{2M}{27c^2|r-s|^{2p}}.
\]

Dado um vetor em blocos z=(z_1,...,z_K), com eta_r=||z_r||_2, a desigualdade triangular e Cauchy–Schwarz fornecem

\[
 \|\Delta z\|_2\le\|A_{\mathrm{count}}\eta\|_2,
\]

onde a matriz não negativa `A_count` tem entradas `2M/[27c²|r-s|^{2p}]` fora da diagonal. Assim,

\[
 \|\Delta\|_2 \le\|A_{\mathrm{count}}\|_2
 \le\max_r\sum_{s\ne r}(A_{\mathrm{count}})_{rs}
 \le\frac{4M}{27c^2}\sum_{k=1}^{\infty}k^{-2p}.
\]

A série converge pois p>1/2. Não existe dependência em K, no total de pontos nem nas razões entre escalas, exceto pelo máximo M especificado. QED.

### Corolário L.1 — Separação **linear** em vez de quadrática

Com `p=1`, `c=1`, `zeta(2)=pi²/6`:

\[
\boxed{D_{rs}\ge|r-s|\max(h_r,h_s),\quad m_r\le M
\implies \|\Delta\|_2\le\frac{2M\pi^2}{81}.}\tag{L-linear}
\]

Quando `M<=14`, a desigualdade `2M/27<=28/27<sqrt(6/5)` (verificar exatamente `5*28²=3920<4374=6*27²`) assegura

\[
 \|\Delta\|_2 <\frac{\pi^2}{3}\sqrt{\frac65},
\]

que é precisamente o **resíduo espectral procurado na P88**, **nessa classe**, não no caso geral.

### Corolário L.2 — Separação de expoente **3/4**

Com `p=3/4` e c=1, a comparação monótona da série com integral dá

\[
\zeta(3/2)=1+\sum_{n=2}^\infty n^{-3/2}
<1+\int_1^\infty x^{-3/2}\,dx=3.
\]

Logo:

\[
\boxed{\|\Delta\|_2<\frac{4M}{9}.}\tag{L-3/4}
\]

Para `M<=4`, a constante é <16/9, abaixo da cota residual da P88. A condição `D_rs>=|r-s|^(3/4)H_rs` é **sublinear** na diferença entre índices. Não requer progressão geométrica nem igualdade dos passos de diferentes grupos.

A hipótese `c>=1` na L-sep assegura, em particular, que os agrupamentos adjacentes estão isolados por pelo menos seus maiores passos internos, de modo que as δ_i utilizadas são realmente os h_r.

## L.3 — Critério espectral adaptativo, sem limite prévio de cardinalidade

**Teorema P92-L.2.** Para agrupamentos uniformes de tamanhos *arbitrários*, separados pelo menos por seus passos entre grupos adjacentes, construa A pela fórmula racional L2. Se para algum racional T>0 existir um certificado exato de positividade definida

\[
\boxed{T I_K-A\succ0,}\tag{L-LDL}
\]

então

\[
\boxed{\|\Delta\|_2< T.}
\]

*Prova.* Para cada vetor z=(z_1,...,z_K), a desigualdade de operadores por blocos dá `||Delta z||_2 <= ||A eta||_2`, eta_r=||z_r||_2. Como A é real simétrica e tem entradas não negativas, seu raio espectral, e portanto sua norma espectral, é seu maior autovalor. A condição L-LDL implica `lambda_max(A)<T`, concluindo.

**Importante:** L-LDL é um critério de certificação finito que pode ser aplicado a qualquer K; não significa que toda configuração o satisfaça. Se os dados forem reais não racionais, uma implementação deve majorar cada A_rs por uma fração racional exata antes de executar LDL. O verificador P92-L usa configurações racionais com `h_r=(raiz_r)^2`.

### Exemplo independente do critério de somas de linha

Construa uma **âncora longa** com 64 pontos e espaçamento 1, seguida por 12 agrupamentos de 2 pontos cada e espaçamento `10^-6`. A primeira separação entre agrupamentos é 1; as demais, `10^-6`, de forma que cada ponto tem δ correto e a matriz de comparação tem estrutura de estrela:

\[
A=\begin{pmatrix}0&w^\top\\w&0\end{pmatrix},\quad K=13,\ N=88.
\]

Usando os coeficientes racionais (L2), o verificador demonstra **exatamente**:

\[
\boxed{\sum_{s} A_{0s}>4,\qquad
 \sum_s A_{0s}^2<(3/2)^2,\qquad (3/2)I-A\succ0.}\tag{L-star}
\]

Como uma estrela simétrica tem autovalores não nulos ±sqrt(sum_s A_0s²), isso também prova por uma **segunda via algébrica** `||Delta||<3/2`, embora a maior soma de linha de A ultrapasse 4. Esse exemplo **não satisfaz a L-sep linear** e não é apresentado como se satisfizesse. Serve especificamente para demonstrar a melhoria fornecida pelo gate espectral adaptativo, mesmo fora da família de separação linear.

O alvo da P88 é aproximadamente 3,604 (<4); logo o teste de soma de linha com A seria inconclusivo, mas L-LDL aprova com margem substancial. A conclusão é apenas sobre **esta configuração concreta** e sobre a validade condicional do método em geral.

## L.4 — Certificado e verificações

O pacote `NUMERIA_P92L_ADAPTIVE_BLOCKS` inclui:

- `certificate.py`: configurações racionais, distâncias aos vizinhos recalculadas, comparação das duas estimativas de cada bloco contra o Hilbert–Schmidt **exato**, matrizes A exatas e LDL racional.
- `test_p92l.py`: **9 testes**, executados em Python normal e com `python -O`, incluindo reconstrução independente da matriz Delta e caso estrela.
- `certificate.json`: **8 configurações com separação linear**, **249 interações de grupos** e **4 LDL completos sobre Delta**, além do gate exato da estrela com **78 interações**.
- Workflow GitHub Actions: gera JSON e manifesto SHA-256, audita as margens em Fraction e proíbe flags de teorema global sem prova.

Comandos, desde a raiz do repositório:

    python -m unittest -v NUMERIA_P92L_ADAPTIVE_BLOCKS.test_p92l
    python -O -m unittest -v NUMERIA_P92L_ADAPTIVE_BLOCKS.test_p92l
    python -m NUMERIA_P92L_ADAPTIVE_BLOCKS.certificate --out p92l-certificate.json

Os testes finitos **não demonstram sozinhos** os teoremas de parâmetros reais. As provas são as desigualdades analíticas L1–L-global e a dominação por matrizes de blocos.

## L.5 — Fronteiras e próximo alvo

**Demonstrado nesta fase:** controle espectral uniforme para qualquer número de agrupamentos com cardinalidade máxima M e separação `c|r-s|^p` com p>1/2; casos lineares com M<=14 e sublineares com M<=4 abaixo da meta residual; gate racional LDL de blocos com cardinalidades arbitrárias e um exemplo em que a soma de linha falha.

**Não demonstrado:** decomposição universal de qualquer conjunto de pontos nesses agrupamentos, controle universal de todas as configurações de Δ, Lema 5/6 gerais de Preissmann, ou Hipótese de Riemann. Nenhum desses resultados é implicado pela passagem dos exemplos finitos.

**Próximo passo P92-M:** reduzir a dependência do tamanho máximo M em L-global usando a forma real de cada acoplamento (núcleos integrais, amortecimento de blocos grandes, sparsificação controlada e composição de certificados LDL). Evitar tratar a convergência de soma de linhas como condição necessária.

**Política:** PR em modo *draft*; sem merge automático.