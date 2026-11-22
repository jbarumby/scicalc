# NUMERIA–RIEMANN P92-D — Cota espectral exata para âncora + cluster

**Trilha:** aprofundamento espectral da P88, sem substituir o programa P93 de cinco lacunas interiores do Lema 6. Baseada na [PR P92](https://github.com/jbarumby/scicalc/pull/92).

## Teorema novo: família infinita, não somente testes finitos

Sejam L>=h>0 reais arbitrários e M>=2 um inteiro. Tome pontos distintos

\[
x_0=0,\qquad x_j=L+(j-1)h,\quad j=1,\ldots,M.
\]

Cada lacuna de vizinho mais próximo da configuração é
\[
\delta_0=L,\qquad \delta_j=h\quad(1\le j\le M).
\]

Defina a matriz real simétrica de diagonal zero da P88:

\[
\Delta_{ij}=\frac{\sqrt{\delta_i\delta_j}
 (\sqrt{\delta_i}-\sqrt{\delta_j})^2}{2(x_i-x_j)^2},\quad i\ne j.
\]

**Teorema P92-D.1 (prova analítica completa).**
Para TODOS L>=h>0 e TODOS M>=2,
\[
\boxed{\|\Delta\|_2\le\frac1{\sqrt{12}}\approx0.2886751345948129.}
\]
A constante 1/sqrt(12) é **ótima como supremo para toda essa família**, embora não seja atingida por configurações finitas.

Este resultado controla uma família adversária de duas escalas que causa divergência em certos testes Schur–Collatz da P88. **Não implica** uma cota uniforme para Δ sobre configurações de pontos arbitrárias.

## Demonstração independente de cálculos numéricos

Seja λ=L/h>=1. Para 1<=j<k<=M temos δ_j=δ_k=h, logo
\[
\Delta_{jk}=0.
\]
O único bloco não trivial liga a âncora aos pontos do cluster:
\[
\Delta=\begin{pmatrix}0&w^\top\\w&0_{M\times M}\end{pmatrix},
\qquad
w_j=\frac{\sqrt{\lambda}(\sqrt{\lambda}-1)^2}
             {2(\lambda+j-1)^2}.
\]

Consequentemente Δ tem posto no máximo 2 e, para w não nulo,
os seus únicos autovalores não nulos são +||w||₂ e -||w||₂. Assim

\[
\boxed{\|\Delta\|_2^2=
\frac{\lambda(\sqrt{\lambda}-1)^4}{4}
\sum_{k=0}^{M-1}\frac1{(\lambda+k)^4}.}
\tag{1}
\]

A função t↦(λ+t)^-4 é decrescente e positiva. Comparando a soma com sua integral,
\[
\sum_{k=0}^{M-1}(\lambda+k)^{-4}
\le\lambda^{-4}+\int_0^\infty(\lambda+t)^{-4}\,dt
=\lambda^{-4}+\frac1{3\lambda^3}.
\]
Com u=1/sqrt(λ), onde 0<u<=1, (1) implica
\[
\|\Delta\|_2^2\le\frac{(1-u)^4(1+3u^2)}{12}.
\]

Basta observar que f(u)=(1-u)^4(1+3u²) é não crescente em [0,1]:
\[
f'(u)=-(1-u)^3(4-6u+18u^2)\le0
\]
pois o polinômio 4-6u+18u²=18(u-1/6)²+7/2 é estritamente positivo. Como f(0)=1, segue
\[
\|\Delta\|_2^2\le1/12
\]
para todo L,h,M admissível. **Isso é uma prova matemática universal nessa família**, independentemente das regressões.

### Otimalidade no limite de dois parâmetros

Escolhendo λ→∞ e M/λ→∞, temos
\[
\lambda^3\sum_{k=0}^{M-1}(\lambda+k)^{-4}
=\frac1\lambda\sum_{k=0}^{M-1}(1+k/\lambda)^{-4}
\longrightarrow\int_0^\infty\frac{dt}{(1+t)^4}=\frac13.
\]
Além disso, (\sqrtλ-1)^4/λ²→1. Da fórmula (1) resulta
\[
\lim\|\Delta\|_2^2=\frac1{12}.
\]

### Explicação da obstrução ao Schur–Collatz da P88

A soma da linha da âncora é
\[
\sum_{j=1}^M\Delta_{0j}
=\frac{\sqrtλ(\sqrtλ-1)^2}{2}
\sum_{k=0}^{M-1}(\lambda+k)^{-2}.
\]

Para λ inteiro e M>=λ, os primeiros λ termos satisfazem
\[
\sum_j\Delta_{0j}\ge
\frac{(\sqrtλ-1)^2}{8\sqrtλ}\longrightarrow\infty.
\]
Portanto o teste ingênuo de soma máxima de linha pode divergir, **ao mesmo tempo que a norma espectral permanece abaixo de 1/sqrt(12)**. É a estrutura de posto dois, e não o teste de linha, que compensa os termos. Isso fortalece a interpretação da observação L4 da P88 sem afirmar que outros métodos de Schur são impossíveis.

## Regressões racionais auditáveis (não são a prova universal)

O [certificate.py](certificate.py) emprega somente Fraction e reconstrói exatamente a fórmula (1) quando sqrtλ é racional. Foram testadas 8 razões de escala e 7 tamanhos de cluster (56 pares), incluindo δ uniforme e valores de grande contraste. A CI executa sete testes em Python normal e em modo -O, gera um JSON com os valores e margens como frações exatas e registra um SHA-256.

Comandos:
    python -m unittest -v NUMERIA_P92D_DELTA_STAR.test_p92d
    python -O -m unittest -v NUMERIA_P92D_DELTA_STAR.test_p92d
    python -m NUMERIA_P92D_DELTA_STAR.certificate --out p92d-star-certificate.json

O artefato e o código evitam promover testes finitos a afirmações gerais. O teorema universal sobre a família é a demonstração algébrica/analítica acima.

## Fronteiras explícitas e próximos passos

- **Demonstrado:** ||Δ||₂<=1/sqrt12 para cada família âncora+cluster uniforme, com parâmetros L>=h>0 e M>=2 arbitrários; o limite é ótimo para essa classe.
- **Ainda não demonstrado:** ||Δ||₂<=(π²/3)sqrt(6/5) para qualquer conjunto de pontos; controle de interações simultâneas entre vários blocos e escalas; Lema 5/6 geral, limite espectral de Preissmann totalmente independente e RH.
- **Trabalho seguinte:** generalizar a compensação de posto dois a acoplamentos entre **múltiplos clusters**. O controle individual de cada bloco não prova sozinho a soma global; a proposta é tratar uma matriz de blocos e obter limites de ortogonalidade/overlap.
- **P93 geométrica:** cinco lacunas internas continua uma tarefa separada, não resolvida aqui.
- Manter PR como rascunho, sem merge automático.
