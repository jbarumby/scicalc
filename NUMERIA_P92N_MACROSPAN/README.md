# NUMERIA–RIEMANN P92-N — Colagem espectral por largura macroscópica e separação apenas local

**Origem:** [P92-M / PR #102](https://github.com/jbarumby/scicalc/pull/102), com correções e auditoria da estrutura de blocos da P92-K e P92-L. Este é um **resultado condicional para uma classe de configurações**, não a demonstração universal do Lema 6 de Preissmann nem da Hipótese de Riemann.

## N.1 — Teorema com larguras físicas limitadas, cardinalidades e escalas arbitrárias

Sejam \(K\ge1\) agrupamentos de pontos reais ordenados e disjuntos. O agrupamento \(r\) é uma progressão aritmética de passo \(h_r>0\) e tem \(m_r\ge2\) pontos. Não existe limite superior imposto a \(K\), a \(m_r\), nem às razões \(h_r/h_s\). Seja \(L>0\) tal que **a largura física** de cada agrupamento satisfaça

\[
\ell_r=(m_r-1)h_r\le L.\tag{N-width}
\]

Denote por \(G_r\) a distância da última posição do grupo \(r\) à primeira do grupo \(r+1\), e suponha, **somente para grupos consecutivos**, que

\[
\boxed{G_r\ge\max\left(\frac{13}{16}L,\,h_r,\,h_{r+1}\right).}\tag{N-local}
\]

Como antes, seja \(\delta_i=\min_{j\ne i}|x_i-x_j|\) e defina \(\Delta\) pela matriz simétrica não negativa

\[
\Delta_{ij}=\frac{\sqrt{\delta_i\delta_j}(\sqrt{\delta_i}-\sqrt{\delta_j})^2}{2(x_i-x_j)^2}\quad(i\ne j),\qquad\Delta_{ii}=0.
\]

**Teorema P92-N.1 (prova analítica abaixo, sem computação de autovalores).** Nas hipóteses (N-width) e (N-local), para todo número finito de pontos,

\[
\boxed{\|\Delta\|_2<\frac{619520}{173901}<\frac{18}{5}
<\frac{\pi^2}{3}\sqrt{\frac65}.}\tag{N1}
\]

**O ganho:** a P92-M exigia uma desigualdade \(D_{rs}\gtrsim(s-r)^{5/4}\max(h_r,h_s)\) para cada par de grupos. Agora usamos apenas lacunas **locais** de ordem \(L\) entre vizinhos e largura física limitada. A condição não controla o número de pontos nem o número de escalas; espaçamentos microscópicos podem implicar grupos de cardinalidade arbitrariamente grande.

## N.2 — Lema de Hilbert–Schmidt com largura finita

Sejam dois grupos uniformes com passos \(a,b>0\), cardinalidades \(m,n\ge2\), larguras \((m-1)a\le L\) e \((n-1)b\le L\), e distância entre extremidades \(D>0\). Supondo que as distâncias ao vizinho mais próximo em cada grupo são seus próprios passos, o bloco cruzado \(W\) da matriz \(\Delta\) tem entradas, após inverter a ordem dos pontos no primeiro grupo,

\[
W_{jk}=\frac{\beta}{(D+ja+kb)^2},\qquad 0\le j<m,\quad0\le k<n,
\quad\beta=\frac{\sqrt{ab}(\sqrt a-\sqrt b)^2}{2}.
\]

Logo, pela norma de Hilbert–Schmidt,

\[
\|W\|_2\le\|W\|_{\rm HS}
\le\frac{\beta\sqrt{mn}}{D^2}
=\frac{(\sqrt a-\sqrt b)^2\sqrt{(ma)(nb)}}{2D^2}
\le\frac{(\sqrt a-\sqrt b)^2\sqrt{(L+a)(L+b)}}{2D^2}.
\tag{N-pair-prelim}
\]

**Lema elementar:** para quaisquer \(0<a,b\le L\),

\[
\boxed{(\sqrt a-\sqrt b)^2\sqrt{(L+a)(L+b)}\le\sqrt2\,L^2.}\tag{N-elementary}
\]

*Prova do lema.* Por simetria assuma \(a\ge b\). Então \((\sqrt a-\sqrt b)^2\le a-b\). Para \(t=b/a\in[0,1]\), como \(a/L\le1\),

\[
\frac{(a-b)^2(L+b)}{a^2L}
=(1-t)^2\left(1+\frac aL t\right)
\le(1-t)^2(1+t)\le1,
\]

pois \(1-(1-t)^2(1+t)=t+t^2-t^3\ge0\). Assim \((a-b)\sqrt{L+b}\le a\sqrt L\); multiplicando por \(\sqrt{L+a}\le\sqrt{2L}\) e usando \(a\le L\), obtemos (N-elementary).

Com isso,

\[
\boxed{\|W\|_2\le\|W\|_{\rm HS}\le\frac{L^2}{\sqrt2\,D^2}
<\frac{5L^2}{7D^2}.}\tag{N-pair}
\]

A última desigualdade é **exatamente racionalizável**, já que \(\sqrt2>7/5\) (porque \(2>49/25\)). A cota depende da distância física como \(D^{-2}\) e **não depende** de \(m,n\). Não contradiz a obstrução da P92-M.2 para \(D^{-1-\varepsilon}\) quando as larguras dos grupos são ilimitadas: aqui impusemos a hipótese adicional \(\ell\le L\).

## N.3 — Prova do teorema com apenas lacunas locais

A hipótese \(G_r\ge\max(h_r,h_{r+1})\) garante que **todos** os pontos do agrupamento \(r\) possuem distância ao vizinho mais próximo exatamente \(\delta_i=h_r\). Consequentemente, os blocos diagonais da matriz \(\Delta\) são nulos.

Para qualquer par \(r<s\), a distância \(D_{rs}\) entre a última posição do primeiro e a primeira do segundo é uma soma de \(s-r\) lacunas de fronteira e das larguras não negativas dos grupos intermediários; portanto,

\[
D_{rs}\ge\sum_{t=r}^{s-1}G_t
\ge\frac{13L}{16}(s-r).
\]

Pelo lema (N-pair), os blocos cruzados satisfazem

\[
\|\Delta_{rs}\|_2<\frac{5L^2}{7D_{rs}^2}
\le\frac{5}{7c^2|r-s|^2},\qquad c=\frac{13}{16}.
\]

Se \(z=(z_1,\ldots,z_K)\) é particionado por grupos, \(\eta_r=\|z_r\|_2\), e \(A\) é a matriz escalar de comparação \(A_{rr}=0\), \(A_{rs}=5L^2/(7D_{rs}^2)\), então

\[
\|\Delta z\|_2\le\|A\eta\|_2
\le\|A\|_2\|z\|_2,
\quad\Rightarrow\quad
\|\Delta\|_2\le\|A\|_2
\le\max_r\sum_{s\ne r}A_{rs}.
\]

Somando a série quadrática dos dois lados de cada linha,

\[
\|\Delta\|_2
\le\frac{5}{7c^2}\,2\sum_{k=1}^\infty\frac1{k^2}
=\frac{5}{7c^2}\frac{\pi^2}{3}.
\]

Pela identidade de Machin e estimativas alternadas, \(\pi<22/7\); substituindo \(c=13/16\):

\[
\|\Delta\|_2
<\frac{5}{7(13/16)^2}\frac{(22/7)^2}{3}
=\boxed{\frac{619520}{173901}}<\frac{18}{5}.
\]

A comparação final é igualmente rigorosa: as mesmas séries de Machin dão \(\pi>157/50\), e \(\sqrt{6/5}>5477/5000\), pois \((5477/5000)^2<6/5\). Logo

\[
\frac{\pi^2}{3}\sqrt{\frac65}
>\frac{(157/50)^2}{3}\frac{5477}{5000}
>\frac{18}{5}.
\]

Isto prova (N1) com **frações e desigualdades analíticas**, sem extrapolar testes finitos. QED.

### Generalização com uma separação local arbitrária \(cL\)

O mesmo argumento fornece, para todo \(c>0\) e lacunas \(G_r\ge\max(cL,h_r,h_{r+1})\),

\[
\boxed{\|\Delta\|_2\le\frac{\pi^2}{3\sqrt2\,c^2}.}\tag{N-c}
\]

É suficiente para a meta residual da P88 tomar \(c^4>5/12\); \(c=13/16\) é apenas uma escolha racional simples, **não ótima**.

## N.4 — Contraexemplo à necessidade da separação P92-M, com 40 escalas

Para \(r=1,\ldots,40\), tome \(h_r=1/r^2\), \(m_r=r^2+1\) e largura física \(\ell_r=1\). Entre grupos vizinhos, use lacuna \(G_r=\max(13/16,h_r,h_{r+1})\). Esta configuração finita é explícita, com

- **40 grupos e 40 valores distintos** de \(\delta_i\);
- **22.180 pontos**;
- tamanho máximo de grupo \(1.601\) pontos;
- distância entre as extremidades dos grupos 1 e 40 igual a \(559/8\) (exata);
- a separação da P92-M exigiria \(D_{1,40}^4\ge39^5\max(h_1,h_{40})^4\), mas \((559/8)^4<39^5\).

Portanto este exemplo está fora das hipóteses da P92-M e **dentro** das da P92-N. Certificamos a matriz de comparação de tamanho \(40\times40\) com aritmética racional, **sem construir nem diagonalizar** uma matriz de \(22.180\times22.180\). A validação universal do teorema é analítica, não depende dessa construção.

## N.5 — Reproduzibilidade independente da prova universal

O script [certificate.py](certificate.py) implementa apenas `fractions.Fraction` e verifica:

1. as cotas inferior e superior de \(\pi\) via Machin e os limiares racionais;
2. seis configurações explícitas, a geometria das fronteiras e a matriz de comparação;
3. **56 interações de dois grupos**, calculando exatamente \(\|W_{rs}\|_{\rm HS}^2\) por somas de racionais e comparando com a cota (N-pair);
4. três fatorações **LDLᵀ completas** de \((18/5)I-\Delta\) e seis fatorações LDL de matrizes de comparação;
5. a família de 40 grupos, 22.180 pontos e 40 escalas por certificados de geometria e comparação, sem reivindicar fatoração LDL da matriz completa.

Executar na raiz do repositório:

```bash
python -m unittest -v NUMERIA_P92N_MACROSPAN.test_p92n
python -O -m unittest -v NUMERIA_P92N_MACROSPAN.test_p92n
python -m NUMERIA_P92N_MACROSPAN.certificate --out p92n-certificate.json
```

**O que não é provado:** a cota universal de \(\Delta\) sobre todas as configurações; os Lemas 5 e 6 gerais de Preissmann; a extensão geométrica de cinco lacunas do programa P93; RH. O teorema P92-N depende da hipótese \(\ell_r\le L\) e de lacunas locais entre agrupamentos; não a remover silenciosamente.

**Próximo passo P92-O:** estudar *larguras físicas não uniformemente limitadas* por meio de partições adaptativas e operadores de comparação \(A_{rs}\) com controle espectral; a obstrução P92-M.2 mostra que é preciso explorar geometria compartilhada, não substituir formalmente \(1/D\) por \(1/D^2\) sem condições.

**Governança:** PR draft, sem merge automático; provas gerais continuam abertas.