# NUMERIA–RIEMANN P92-G — Norma de Delta para progressões geométricas com escalas ilimitadas

**Base:** P92-F [PR #95](https://github.com/jbarumby/scicalc/pull/95). Este é um resultado de trilha espectral independente da P93 geométrica. **Não é** a prova do caso geral de Preissmann ou da Hipótese de Riemann.

## Teorema P92-G.1 — todos os parâmetros reais e dimensões

Para qualquer q>1 REAL e qualquer N>=2, sejam os N pontos
\[
x_j=\sum_{\ell=0}^{j-1}q^\ell=\frac{q^j-1}{q-1},\qquad j=0,\ldots,N-1.
\]
As lacunas consecutivas são 1,q,q²,...,q^(N-2). As distâncias aos vizinhos mais próximos são
\[
\delta_0=\delta_1=1,\qquad \delta_j=q^{j-1}\ (j\ge2).
\]
A matriz P88 é simétrica, de diagonal nula, e tem entradas
\[
\Delta_{ij}=\frac{\sqrt{\delta_i\delta_j}(\sqrt{\delta_i}-\sqrt{\delta_j})^2}{2(x_i-x_j)^2}\ge0\quad(i\ne j).
\]

**Teorema (prova analítica abaixo, sem float):**
\[
\boxed{\|\Delta\|_2 < \frac{15}{16}=0.9375.}\tag{G1}
\]

A cota vale para **TODOS** os q>1 e todos os inteiros N>=2; não depende de N, da razão geométrica nem do contraste das escalas. Portanto, o número de valores distintos das delta_i pode crescer sem limite (m=N-1), e a razão max(delta)/min(delta)=q^(N-2) também pode crescer sem limite.

Esta é uma classe de pontos não uniformes estruturalmente especial; nada aqui demonstra a estimativa universal da P88 para configurações arbitrárias.

## 1. Núcleo exatamente Toeplitz para índices i,j>=1

Escreva r=sqrt(q)>1. Para quaisquer 1<=i<j<=N-1 e k=j-i, temos
\[
\delta_i=q^{i-1},\quad \delta_j=q^{j-1},\quad
x_j-x_i=\frac{q^i(q^{j-i}-1)}{q-1}.
\]

Substituindo e cancelando (r^k-1)^2 contra (r^(2k)-1)^2, obtemos a identidade **exata**:
\[
\boxed{\Delta_{ij}=\kappa_q(k),\quad
\kappa_q(k)=\frac{(q-1)^2}{2q^2}\frac{r^k}{(1+r^k)^2}.}\tag{G2}
\]

O bloco T=(Delta_ij)_{1<=i,j<N} é Toeplitz simétrico, positivo entrada a entrada e de diagonal zero. Todos os termos dependem apenas de |i-j|.

Pelo limite elementar r^k/(1+r^k)^2<=r^(-k):
\[
S_q:=\sum_{k=1}^{\infty}\kappa_q(k)
\le\frac{(q-1)^2}{2q^2}\sum_{k=1}^{\infty}r^{-k}
=\frac{(q-1)^2}{2q^2(r-1)}
=\frac12\,t(1-t)(1+t)^2,\quad t=1/r\in(0,1).
\tag{G3}
\]
Todas as somas de linha de T são <=2S_q, então ||T||_2<=2S_q.

## 2. A primeira linha é dominada por uma linha do bloco Toeplitz

Como delta_0=delta_1=1 e x_0=0<x_1=1<x_j para j>=2, os numeradores de Delta_0j e Delta_1j são idênticos, enquanto (x_j-x_0)^2>(x_j-x_1)^2. Logo:
\[
0\le\Delta_{0j}\le\Delta_{1j}=\kappa_q(j-1),\quad j\ge2,
\qquad \Delta_{01}=0.
\tag{G4}
\]

O vetor da âncora w=(Delta_01,...,Delta_0,N-1) satisfaz ||w||_2<=||w||_1<=S_q.

Decompondo
\[
\Delta=\begin{pmatrix}0&w^\top\\w&T\end{pmatrix}
=\begin{pmatrix}0&w^\top\\w&0\end{pmatrix}
+\begin{pmatrix}0&0\\0&T\end{pmatrix},
\]
a norma da primeira matriz é ||w||_2 (posto<=2), e a norma da segunda é ||T||_2. Da desigualdade triangular:
\[
\boxed{\|\Delta\|_2\le3S_q
\le\frac32\,t(1-t)(1+t)^2.}\tag{G5}
\]

Há igualmente uma prova por máximas somas de linha: as linhas 0 e 1 somam <=S_q, e qualquer linha i>=2 soma no máximo 2S_q+Delta_i0<=3S_q. Nenhuma dessas etapas depende de N ou de limite de amostragem.

## 3. Otimização global exata da envoltória

Defina f(t)=t(1-t)(1+t)^2 para 0<=t<=1. Sua derivada é
\[
f'(t)=(1+t)(1+t-4t^2).
\]
O máximo global ocorre em
\[
t_*=\frac{1+\sqrt{17}}8\quad (0<t_*<1),
\]
e, usando 4t_*²=t_*+1,
\[
f(t_*)=\frac{107+51\sqrt{17}}{512}.
\]

A desigualdade racional exata
\[
\sqrt{17}<\frac{71}{17}
\quad\text{porque}\quad17\cdot17^2=4913<5041=71^2
\]
fornece
\[
f(t_*)<\frac{107+51(71/17)}{512}=\frac58.
\]
Portanto, por (G5),
\[
\boxed{\|\Delta\|_2
\le\frac32 f(t)<\frac{3}{2}\frac58
=\frac{15}{16}.}
\]
A demonstração completa é em aritmética real elementar. A cota 15/16 não é declarada ótima; apenas válida e uniforme.

A meta residual da P88 é (pi²/3)sqrt(6/5)>3, logo **G1 implica esse limite com ampla margem para toda a família geométrica**, sem condições sobre N ou q.

## 4. Relação com P92-F.2

O critério espectral P92-F.2 pode se tornar inconclusivo para muitos valores de delta com contrastes elevados. A estimativa P92-G emprega a **localização dos pontos e a separação geométrica**, não somente os valores possíveis de delta_i.

Isso demonstra por que a soma das cotas independentes por escala não pode ser usada como critério de impossibilidade para a estimativa universal.

Não basta, porém, decompor uma configuração arbitrária em progressões geométricas e somar as respectivas normas: os **acoplamentos entre progressões** não estão controlados aqui. A ponte necessária é uma estimativa de blocos ou um argumento de interpolação que preserve a norma.

## 5. Verificador racional e CI

O arquivo [certificate.py](certificate.py) usa apenas Fraction e compara todas as entradas da matriz com a fórmula Toeplitz G2, a dominação G4 e a cota por somas de linha G5. Para sqrt(q) racional, todas as contas são exatas.

- 7 razões r=sqrt(q) e 8 dimensões N=2,3,4,6,9,13,17,25: 56 configurações;
- uma configuração adicional de 60 pontos: **57 exemplos racionais**;
- **42 fatorações LDL** de (15/16)I-Delta feitas com Fraction, independentes da fórmula de Schur, como testes adicionais;
- a fórmula derivativa e seu máximo são a **prova** válida para todo r REAL, não os 57 exemplos.

Reproduzir:
    python -m unittest -v NUMERIA_P92G_GEOMETRIC_GAPS.test_p92g
    python -O -m unittest -v NUMERIA_P92G_GEOMETRIC_GAPS.test_p92g
    python -m NUMERIA_P92G_GEOMETRIC_GAPS.certificate --out p92g-geometric.json

O JSON separa explicitamente os resultados verdadeiros para a classe geométrica dos problemas ainda abertos.

## 6. Fronteira científica

- **Provado sem importações gerais:** ||Delta||_2<15/16 para TODOS os N e q>1 na família geométrica, apesar de escala e contraste ilimitados.
- **Não provado:** mesmo limite para todas as configurações, nem a cota P88 original para todas as configurações; Preissmann Lema 5/6 gerais; P93 geométrica de cinco lacunas; Hipótese de Riemann.
- **Próxima investigação:** interações entre **várias progressões geométricas intercaladas**, usando somas de linha entre sequências, geometria e controle da norma de blocos. Não extrapolar G1 por mera soma de normas.
