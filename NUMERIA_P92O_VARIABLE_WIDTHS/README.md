# NUMERIA–RIEMANN P92-O — Larguras macroscópicas adaptativas e porta LDL exata

**Origem:** continuação de [P92-N / PR #103](https://github.com/jbarumby/scicalc/pull/103). Esta é a trilha espectral da matriz `Delta` da P88, **não** a prova da desigualdade de Preissmann para todas as configurações, nem da Hipótese de Riemann. A análise neste arquivo não depende de autovalores calculados em ponto flutuante.

## Teorema O.1 — larguras variáveis sem limite global comum

Considere qualquer número finito `K>=1` de agrupamentos ordenados e disjuntos de pontos reais. Cada grupo `r` é uma progressão aritmética de `m_r>=2` pontos e passo `h_r>0`; sua largura física é

\[
\ell_r=(m_r-1)h_r>0.
\]

Para `r<s`, seja `D_rs` a distância física entre a extremidade direita do grupo `r` e a extremidade esquerda do grupo `s`.

**Hipótese (O-pair) — todos os pares**:

\[
\boxed{D_{rs}\ge |r-s|^{3/4}(\ell_r+\ell_s)\quad\text{para todos }r<s.}\tag{O-pair}
\]

Então, para as distâncias aos vizinhos mais próximos `delta_i=min_{j!=i}|x_i-x_j|` e a matriz P88

\[
\Delta_{ij}=\frac{\sqrt{\delta_i\delta_j}(\sqrt{\delta_i}-\sqrt{\delta_j})^2}{2(x_i-x_j)^2}
\quad (i\ne j),\qquad \Delta_{ii}=0,
\]

**Teorema P92-O.1 (analítico para todos os parâmetros reais):**

\[
\boxed{\|\Delta\|_2<3<\frac{\pi^2}{3}\sqrt{\frac65}.}\tag{O1}
\]

Não existe limite superior para `K`, `m_r`, `ell_r` nem para os contrastes `ell_r/ell_s`. Diferentemente da P92-N, não é preciso encontrar uma escala física comum `L` que majorize todas as larguras. **Contudo, O-pair é uma condição de TODOS OS PARES, não uma condição meramente local entre grupos consecutivos; não afirmar que substitui a P92-N integralmente.**

### Prova analítica, independente do código

Em particular, para grupos consecutivos, O-pair impõe
`G_r=D_{r,r+1}>=ell_r+ell_{r+1}>=max(h_r,h_{r+1})`. Assim, cada ponto de um grupo tem distância ao vizinho mais próximo `delta_i=h_r`, e os blocos diagonais de `Delta` são nulos.

Fixe dois grupos `r<s` com passos `a=h_r`, `b=h_s`, larguras `ell_r=(m_r-1)a`, `ell_s=(m_s-1)b` e separação `D=D_rs`. Revertendo a ordem dos índices do primeiro grupo, as entradas do bloco cruzado são

\[
(\Delta_{rs})_{jk}=\frac{\beta}{(D+ja+kb)^2},
\quad \beta=\tfrac12\sqrt{ab}(\sqrt a-\sqrt b)^2,
\quad 0\le j<m_r,\quad 0\le k<m_s.
\]

Como `m_r a=ell_r+a<=2ell_r`, analogamente `m_s b<=2ell_s`, e `(sqrt(a)-sqrt(b))²<=a+b<=ell_r+ell_s`, segue da norma de Hilbert–Schmidt:

\[
\begin{aligned}
\|\Delta_{rs}\|_2
&\le \|\Delta_{rs}\|_{\mathrm{HS}}
\le \frac{(\sqrt a-\sqrt b)^2\sqrt{(m_ra)(m_sb)}}{2D^2}\\
&\le \frac{(a+b)\sqrt{\ell_r\ell_s}}{D^2}
\le \frac{(\ell_r+\ell_s)\sqrt{\ell_r\ell_s}}{D^2}
\le \boxed{\frac{(\ell_r+\ell_s)^2}{2D^2}}.
\end{aligned}
\tag{O-HS}
\]

A última desigualdade é somente `2sqrt(ell_r ell_s)<=ell_r+ell_s`. Portanto, definindo uma **matriz racional de comparação**

\[
A_{rr}=0,\qquad A_{rs}=\frac{(\ell_r+\ell_s)^2}{2D_{rs}^2}\quad(r\ne s),
\]

para qualquer vetor particionado `v=(v_1,...,v_K)` e `eta_r=||v_r||_2` vale

\[
\|\Delta v\|_2\le \|A\eta\|_2\le\|A\|_2\|v\|_2,
\qquad\boxed{\|\Delta\|_2\le\|A\|_2}.\tag{O-block}
\]

Pela hipótese O-pair,

\[
0\le A_{rs}\le\frac1{2|r-s|^{3/2}}.
\]

Todas as somas de linha de `A` são majoradas por

\[
\sum_{s\ne r}A_{rs}\le\sum_{k=1}^{\infty}k^{-3/2}
<1+\int_1^\infty t^{-3/2}\,dt=3.
\]

A desigualdade estrita vale porque `t^-3/2` é estritamente decrescente. Como `A` é simétrica e tem entradas não negativas, `||A||_2<=max(row sums)<3`. Isto prova O1 **para todo K finito e todos os parâmetros reais admissíveis**. A comparação `3<pi²/3 sqrt(6/5)` é imediata de `pi>3` e `sqrt(6/5)>1`. QED.

> Uma generalização automática exige `D_rs>=c|r-s|^alpha(ell_r+ell_s)` com `alpha>1/2` e fornece `||Delta||<=zeta(2alpha)/c²`. **Esta é somente uma condição suficiente pelo método das somas; para `alpha<=1/2`, divergência da série não constitui contraexemplo à desigualdade verdadeira.**

## Teorema O.2 — critério espectral racional adaptativo quando O-pair falha

A identidade O-HS e a comparação O-block **não dependem** de O-pair. Basta que todos os grupos tenham passo positivo e que as lacunas de fronteira satisfaçam

\[
G_r\ge\max(h_r,h_{r+1})
\]

para garantir `delta_i=h_r` em cada grupo.

Tome `T=18/5`. Se uma configuração com larguras `ell_r` e separações `D_rs` **racionais** produzir uma matriz racional `A` tal que

\[
\boxed{TI_K-A\succ0},
\]

então a fatoração `LDL^T` com pivôs exatamente positivos certifica

\[
\boxed{\|\Delta\|_2<18/5<\frac{\pi^2}{3}\sqrt{6/5}.}\tag{O2}
\]

Para a última comparação, valem as desigualdades racionais verificáveis `pi>157/50`, `sqrt(6/5)>5477/5000`, e

\[
\frac{(157/50)^2}{3}\frac{5477}{5000}>\frac{18}{5}.
\]

A P92-N documenta a cota inferior de pi por Machin; o script P92-O verifica apenas o último comparativo racional e a desigualdade envolvendo `sqrt(6/5)`. A prova de O2 é direta: o máximo autovalor de `A` é menor que `T`, `A` é simétrica não negativa e `||A||_2=lambda_max(A)`, logo O-block implica a conclusão. A condição positiva definida é um certificado **para a configuração racional fornecida**; não é automaticamente satisfeita por todos os pontos reais.

### Exemplo em que a soma de linha falha, mas a LDL tem sucesso

Uma configuração de `K=13` grupos, `26` pontos, possui uma âncora uniforme com dois pontos e passo `100`, seguida de `12` grupos de dois pontos com passo `1/1024`. Entre o primeiro e o segundo grupo, use `G=100`; entre grupos pequenos, `G=1`. Os passos são quadrados racionais, permitindo auditorias independentes do núcleo `Delta` quando desejado.

Nessa configuração, a matriz de comparação **racional** tem

\[
\max_r\sum_s A_{rs}>5.4>18/5,
\]

portanto o critério ingênuo de soma máxima **não** serve. Porém `TI-A` possui pivôs LDL racionais estritamente positivos, com menor pivô acima de `3.53` (valor decimal só para exibição). Assim O2 certifica `||Delta||<18/5` apesar da falha da soma de linha. A hipótese O-pair também falha.

### Exemplo que não obedece à P92-N, mas obedece à O.1

Os seis grupos com `(m,h)=(2,100),(2,1),(2,1/100),(5,4),(10,1/16),(3,9)` apresentam larguras `100,1,1/100,16,9/16,18`. Escolhendo lacunas entre grupos consecutivos iguais a `7(ell_r+ell_{r+1})`, obtemos O-pair por verificação exata da quarta potência. Entre dois agrupamentos pequenos, a lacuna fica abaixo de `(13/16)*100`, portanto não existe hipótese N-local com `L>=100`, embora O.1 se aplique.

## Certificado computacional e limites

O [certificate.py](certificate.py) usa `fractions.Fraction`, com:

- verificações geométricas da distância mínima e do teste de quarta potência `D_rs^4 >= (s-r)^3(ell_r+ell_s)^4`;
- comparações **exatas** de normas de Hilbert–Schmidt ao quadrado, quando os passos têm raízes racionais;
- construção exata da matriz racional de comparação `A`;
- fatoração `LDL^T` exata de `(18/5)I-A`;
- verificação de `(18/5)I-Delta` por fatoração LDL de uma matriz completa de `13` pontos como regressão independente;
- um exemplo de `K=30` grupos e `93.994` pontos em que **somente** a matriz de comparação `30x30` é fatorada; a matriz total de `93.994x93.994` pontos **não é construída**.

A saída registra explicitamente que **os lemas P88/Preissmann gerais e RH permanecem abertos**. CI e amostragens são verificações da implementação; a afirmação para todos os parâmetros reais em O.1 vem da prova analítica acima.

Comandos na raiz do repositório:

```bash
python -m unittest -v NUMERIA_P92O_VARIABLE_WIDTHS.test_p92o
python -O -m unittest -v NUMERIA_P92O_VARIABLE_WIDTHS.test_p92o
python -m NUMERIA_P92O_VARIABLE_WIDTHS.certificate --out p92o-certificate.json
```

## Próxima investigação P92-P

Explorar hipóteses **somente locais** entre grupos consecutivos, mas com larguras físicas variáveis, que possam implicar uma cota global para `||Delta||`, sem impor O-pair explicitamente. Uma direção é usar a distribuição de comprimentos e a ortogonalidade aproximada de blocos, em vez da soma isolada das cotas `A_rs`.

**Governança:** PR como rascunho; sem merge automático. Não alegar prova geral de Preissmann, dos zeros da função zeta, ou RH.