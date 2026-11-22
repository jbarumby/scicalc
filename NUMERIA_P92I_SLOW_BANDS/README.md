# NUMERIA–RIEMANN P92-I — Faixas lentas, escalas contínuas e composição de blocos

**Base:** [P92-H, PR #97](https://github.com/jbarumby/scicalc/pull/97). Este é um avanço **na trilha espectral** da matriz Δ da P88; não resolve a P93 geométrica de cinco lacunas e não demonstra a Hipótese de Riemann.

## 1. Teorema I.1 — controle por uma faixa de vizinhos próximos

Sejam `N>=2` e `x_0<...<x_(N-1)` números reais distintos. Ponha

\[
\delta_i=\min_{j\ne i}|x_i-x_j|>0,
\qquad
\Delta_{ij}=\frac{\sqrt{\delta_i\delta_j}(\sqrt{\delta_i}-\sqrt{\delta_j})^2}{2(x_i-x_j)^2}\quad(i\ne j),
\qquad\Delta_{ii}=0.
\]

Seja \(h=\min_i\delta_i\) e \(H=\max_i\delta_i\), com \(R=H/h\). **Para quaisquer pontos reais e todo \(1\le R\le9\)**, temos a desigualdade universal nesta classe:

\[
\boxed{\|\Delta\|_2\le \zeta(2)\sqrt R\,(\sqrt R-1)^2,\quad \zeta(2)=\pi^2/6.}\tag{I1}
\]

Não há hipótese de monotonicidade nem de progressão aritmética/geométrica; os \(\delta_i\) podem assumir uma quantidade arbitrária de valores distintos, incluindo valores reais não racionais.

**Corolário I.1a (faixa de fator quatro).** Se \(H/h\le4\), segue

\[
\boxed{\|\Delta\|_2\le2\zeta(2)=\pi^2/3<\frac{\pi^2}3\sqrt{6/5}.}\tag{I2}
\]

**Corolário I.1b (faixa um pouco mais ampla).** Ponha \(s_*=2037/1000\). Como a função \(s(s-1)^2\) cresce para \(s\ge1\), e

\[
[s_*(s_*-1)^2]^2-\frac{24}{5}
= -\frac{1592982496617591}{1000000000000000000}<0,
\]

se \(H/h\le s_*^2=4149369/1000000\approx4.149369\), então

\[
\boxed{\|\Delta\|_2<\frac{\pi^2}3\sqrt{6/5}.}\tag{I3}
\]

Essa é uma classe nova de validade **analítica** da estimativa residual da P88, independentemente da dimensão. O valor \(s_*\) é uma escolha racional certificada, não uma otimização exata do limiar.

### Demonstração do Teorema I.1

Para \(a,b\in[h,H]\), tome \(u=\max(\sqrt a,\sqrt b)/\sqrt h\) e \(v=\min(\sqrt a,\sqrt b)/\sqrt h\). Então \(1\le v\le u\le\sqrt R\le3\) e

\[
\frac{\sqrt{ab}(\sqrt a-\sqrt b)^2}{h^2}=uv(u-v)^2.
\]

Fixando \(u\), a derivada de \(v\mapsto uv(u-v)^2\) é \(u(u-v)(u-3v)\le0\) quando \(v\ge1\) e \(u\le3\). Portanto o máximo ocorre em \(v=1\). A função \(u\mapsto u(u-1)^2\) é crescente para \(u\ge1\). Daí

\[
\sqrt{ab}(\sqrt a-\sqrt b)^2\le h^2\sqrt R(\sqrt R-1)^2.
\]

Como todos os intervalos consecutivos \(x_{k+1}-x_k\ge h\), para \(i\ne j\) vale \(|x_i-x_j|\ge|i-j|h\). Logo

\[
0\le\Delta_{ij}\le\frac{\sqrt R(\sqrt R-1)^2}{2(i-j)^2}.
\]

A soma em cada linha é majorada por \(\sqrt R(\sqrt R-1)^2\sum_{k\ge1}1/k^2\), pois há no máximo um ponto a cada distância de índice positiva em cada lado. Para matrizes simétricas com entradas não negativas, a norma operador não ultrapassa a maior soma de linha. Isso prova (I1), sem argumento numérico.

## 2. Teorema I.2 — composição de **blocos não uniformes**

Particione os índices em \(K\ge1\) **blocos consecutivos não vazios** \(I_1,\ldots,I_K\). Seja

\[
h_r=\min_{i\in I_r}\delta_i,\quad H_r=\max_{i\in I_r}\delta_i,
\quad R_r=H_r/h_r\le4.
\]

Os blocos podem ter tamanhos arbitrários, espaçamentos internos **irregulares**, escalas diferentes e número arbitrário de pontos. Para \(r<s\), defina \(S_{rs}\) como a distância do **último** ponto de \(I_r\) ao **primeiro** ponto de \(I_s\), e

\[
\kappa_{rs}=\frac12\sqrt{H_rH_s}(\sqrt{H_r}+\sqrt{H_s})^2.
\]

Forme uma matriz simétrica não negativa de dimensão \(K\) por

\[
\begin{aligned}
A_{rr}&=\zeta(2)\sqrt{R_r}(\sqrt{R_r}-1)^2,\\
A_{rs}&=\kappa_{rs}\left[\frac1{S_{rs}^4}+\frac{h_r^{-1}+h_s^{-1}}{3S_{rs}^3}
+\frac1{6h_rh_sS_{rs}^2}\right]^{1/2}\quad(r\ne s).
\end{aligned}\tag{I4}
\]

Então, **para toda configuração real que satisfaça os limites intra-bloco**, com separações entre blocos positivas mas sem exigir progressões aritméticas, vale

\[
\boxed{\|\Delta\|_2\le\|A\|_2\le\max_r\sum_s A_{rs}.}\tag{I5}
\]

O último termo é um **teste condicional por configuração**. Nenhuma constante independente de \(K\) é afirmada sem controle adicional das interações.

### Prova da comparação entre blocos

No bloco diagonal \(r\), a prova do Teorema I.1 vale usando apenas os pontos do bloco e as distâncias \(\delta_i\) calculadas na configuração inteira: cada intervalo interno é \(\ge h_r\). Por isso \(\|\Delta_{rr}\|_2\le A_{rr}\).

Para blocos \(r<s\), enumere pontos do primeiro a partir de sua extremidade direita por \(p=0,1,\ldots\), e pontos do segundo a partir da esquerda por \(q=0,1,\ldots\). As distâncias respectivas satisfazem

\[
|x_i-x_j|\ge S_{rs}+ph_r+qh_s,
\]

pois as lacunas entre pontos adjacentes de um mesmo bloco são \(\ge h_r\) ou \(\ge h_s\). O numerador em \(\Delta_{ij}\) é no máximo \(\kappa_{rs}\), já que \(\delta_i\le H_r\), \(\delta_j\le H_s\) e \((\sqrt a-\sqrt b)^2\le(\sqrt a+\sqrt b)^2\).

Assim, a norma Hilbert–Schmidt do bloco fora da diagonal satisfaz

\[
\|\Delta_{rs}\|_2^2\le\|\Delta_{rs}\|_{HS}^2
\le\kappa_{rs}^2\sum_{p,q\ge0}\frac1{(S_{rs}+ph_r+qh_s)^4}.
\]

A comparação por integrais de funções decrescentes, separando o termo \((p,q)=(0,0)\), os eixos e o interior, dá

\[
\sum_{p,q\ge0}(S+ph+qg)^{-4}
\le S^{-4}+\frac{h^{-1}+g^{-1}}{3S^3}+\frac1{6hgS^2}.
\]

Isso é precisamente \(A_{rs}^2\). Para vetor \(v\) dividido em subvetores \(v_r\), ponha \(t_r=\|v_r\|_2\). A desigualdade triangular fornece \(\|(\Delta v)_r\|_2\le\sum_s A_{rs}t_s\). A soma de quadrados e a simetria de \(A\) completam a prova de (I5). **A P92-E exigia progressões aritméticas internas; aqui não são necessárias.**

## 3. Certificação racional a posteriori

Para qualquer repartição em blocos, o teste rigoroso é calcular um **majorante superior**, não um autovalor flutuante, da maior soma de linha de \(A\). Se ele ficar abaixo de

\[
L_{\mathrm{P88}}:=\frac{2725}{768}
=\frac{(25/8)^2}{3}\frac{109}{100}
<\frac{\pi^2}3\sqrt{6/5},
\]

então está provado que a estimativa residual P88 vale para essa configuração. As direções dos arredondamentos são explícitas: \(\pi<22/7\Rightarrow\zeta(2)<5/3\) para **majorantes**; \(\pi>25/8\) e \(\sqrt{6/5}>109/100\) para o **limite inferior** do alvo P88; raízes de \(A_{rs}^2\) são arredondadas **para cima** com \(\mathrm{isqrt}\) inteiro.

O [certificate.py](certificate.py) faz regressões exatas sobre 6 configurações racionais de 1, 2 ou 3 blocos, com razões internas distintas e contraste global de \(\delta\) até **4225**. Em quatro configurações são executadas também fatorações LDL exatas **para ambos os sinais** de \(L_{\mathrm{P88}}I\pm\Delta\), independentes da lógica de soma de normas. Para o verificador, as lacunas internas foram escolhidas quadrados de racionais; essa restrição **não aparece nos teoremas**, que valem para reais.

Reprodução:

    python -m unittest -v NUMERIA_P92I_SLOW_BANDS.test_p92i
    python -O -m unittest -v NUMERIA_P92I_SLOW_BANDS.test_p92i
    python -m NUMERIA_P92I_SLOW_BANDS.certificate --out p92i-bands.json

O JSON contém os dados completos de cada bloco, norma HS quadrática efetiva, majorantes integrais, certificado Schur em frações e pivôs racionais LDL. Os exemplos são regressões; os teoremas seguem das demonstrações das seções anteriores.

## 4. Fronteira e próximo passo P92-J

- **Provado:** para qualquer conjunto de pontos com \(H/h\le4\), \(\|\Delta\|_2\le\pi^2/3\); para \(H/h\le(2037/1000)^2\), a estimativa residual P88 completa vale. Sem limite em \(N\), número de valores de \(\delta\) ou disposição dos pontos.
- **Provado:** composição de blocos irregulares de faixa interna \(\le4\) pela matriz \(A\), com certificado a posteriori e sem dependência do número de pontos dentro dos blocos.
- **Não provado:** uma cota residual uniforme de \(\Delta\) para **todas** as configurações nem para todo \(K\); controle uniforme dos blocos sem hipótese sobre suas separações; Lemas 5–6 gerais de Preissmann; P93 geométrica de cinco lacunas; RH.
- **Próxima prioridade:** encontrar um algoritmo de partição por escalas **com sobreposição controlada** e obter uma cota de \(\|A\|\) independente de \(K\), sem presumir separações que não existem. Considerar blocos com fronteiras dinâmicas ou espectrais, especialmente nas transições entre regiões quase uniformes e saltos grandes.

Preservar PR como **draft**, sem merge automático. Não atribuir prova de RH aos certificados da matriz \(\Delta\).