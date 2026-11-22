# NUMERIA P88 — Lema 6: identidade mista, malha uniforme exata e obstáculo residual

**Objetivo:** ir além dos 223 testes finitos da P87, sem repetir DOC1/DOC2, os patches S3 ou a P84. O alvo é a demonstração independente do **Lema 6 uniforme** de Preissmann para quaisquer espaçamentos.

**Estado científico:** obtivemos duas afirmações matemáticas gerais novas e uma igualdade exata no caso uniforme; **não** provamos ainda a desigualdade afiada para toda configuração arbitrária. O Lema 6 permanece como teorema externo até fechar o termo misto identificado abaixo.

## 1. Enunciado e fonte exatos

Sejam \(x_1<\cdots<x_N\) pontos reais distintos, \(\delta_j=\min_{k\ne j}|x_j-x_k|>0\), e \(a\ne b\) índices. Escreva \(h=x_b-x_a\ne0\). A desigualdade afiada é

\[
\boxed{\quad
\sum_{k\notin\{a,b\}}
\frac{\delta_k}{(x_k-x_a)^2(x_k-x_b)^2}
\le
\frac{\pi^2}{3h^2}\left(\frac1{\delta_a}+\frac1{\delta_b}\right)
-\frac{3(\delta_a+\delta_b)}{h^4}.
\quad}\tag{L6}
\]

- [Yangjit, *On the Montgomery–Vaughan weighted generalization of Hilbert's inequality*, arXiv:2203.14950v1, Lema 3.2, equação (3.6), p. 9](https://arxiv.org/pdf/2203.14950) reproduz o enunciado e atribui a prova a Preissmann.
- [E. Preissmann, *Sur une inégalité de Montgomery–Vaughan*, Enseign. Math. 30 (1984), pp. 108–110, Lemme 6](https://www.e-periodica.ch/cntmng?pid=ens-001:1984:30::31). O PDF indexado pela E-Periodica confirma a estrutura da prova: uniformização dos pontos exteriores, decomposição por frações parciais, aplicação do Lema 5 às somas de reciprocais e limite das caudas. A digitalização tem fórmulas de extração textual degradadas, portanto NÃO reivindicamos ter rederivado o Lema 5 ou o passo final de pp. 109–110 apenas a partir desse texto.
- P87 [PR #87](https://github.com/jbarumby/scicalc/pull/87) demonstrou \(\sum_{k\ne a}\delta_k/(x_k-x_a)^2\le\pi^2/(3\delta_a)\) **para todo N finito**, reconstruindo o método convexo de Shan.

## 2. Redução exata a um compensador de sinal misto

Defina \(\mathcal A_a=\sum_{k\ne a}\delta_k/(x_k-x_a)^2\), \(\mathcal A_b\) analogamente, e

\[
C_{a,b}=\sum_{k\notin\{a,b\}}
\frac{\delta_k}{(x_k-x_a)(x_k-x_b)}.
\]

Como \(h=(x_k-x_a)-(x_k-x_b)\), temos a **identidade algébrica exata**
\[
\frac{1}{(x_k-x_a)^2(x_k-x_b)^2}
=\frac{1}{h^2}\left(
\frac{1}{(x_k-x_a)^2}
+\frac{1}{(x_k-x_b)^2}
-\frac{2}{(x_k-x_a)(x_k-x_b)}
\right).
\]
Ao excluir os termos \(k=a,b\) das duas somas diagonais, aparece exatamente \((\delta_a+\delta_b)/h^2\). Logo,
\[
h^2 S_{a,b}=\mathcal A_a+\mathcal A_b
-\frac{\delta_a+\delta_b}{h^2}-2C_{a,b}.
\]

Escreva os **déficits diagonais não negativos**, provados na P87:
\[
\Delta_a=\frac{\pi^2}{3\delta_a}-\mathcal A_a\ge0,\qquad
\Delta_b=\frac{\pi^2}{3\delta_b}-\mathcal A_b\ge0.
\]

**Proposição (equivalência rigorosa):** a desigualdade afiada (L6) é **exatamente equivalente** a
\[
\boxed{\quad
E_{a,b}:=\Delta_a+\Delta_b+2C_{a,b}
-\frac{2(\delta_a+\delta_b)}{h^2}\ \ge0.
\quad}\tag{P88-GATE}
\]
De fato,
\[
h^2\bigl(\mathrm{RHS}_{L6}-S_{a,b}\bigr)=E_{a,b}
\]
é identidade, antes de qualquer aproximação ou teste.

**Obstáculo real:** o termo \(C_{a,b}\) **não tem sinal positivo em geral**: na configuração \((x_0,x_1,x_2,x_3)=(0,1,2,3)\), escolhendo os extremos \(a=0\), \(b=3\), temos \(C_{a,b}=-1\). Portanto, não basta provar separadamente \(\Delta_a,\Delta_b\ge0\), e seria *incorreto* usar \(C_{a,b}\ge0\) como atalho. A prova completa exige compensação quantitativa entre os déficits diagonais e o termo misto negativo.

## 3. Um limite universal provado, mas mais fraco

A identidade real \((u-v)^2\le2(u^2+v^2)\), com \(u=x_k-x_a\), \(v=x_k-x_b\), fornece
\[
\frac1{u^2v^2}\le
\frac2{h^2}\left(\frac1{u^2}+\frac1{v^2}\right).
\]
Somando e aplicando o Lema 1 já provado na P87, com as exclusões explícitas:
\[
\boxed{\quad
S_{a,b}\le\frac2{h^2}\left[
\frac{\pi^2}{3}\left(\frac1{\delta_a}+\frac1{\delta_b}\right)
-\frac{\delta_a+\delta_b}{h^2}
\right].
\quad}\tag{P88-WEAK}
\]
Esta desigualdade vale **para quaisquer N finitos e pontos reais distintos**, sem importar o Lema 6. É comprovada por uma demonstração analítica geral, não extrapolação de testes. Sua constante é insuficiente para manter a mesma \(m_0\) na cadeia de Preissmann. Por isso a prova afiada continua pendente.

## 4. Igualdade exata na malha uniforme e necessidade do termo −3

Tome a malha infinita \(x_k=kd\), \(k\in\mathbb Z\), \(d>0\). Para índices \(0,m\), \(m\ge1\), todas as \(\delta_k=d\). No caso \(d=1\), a fração parcial do §2 implica:
\[
\sum_{k\ne0,m}\frac1{k^2(k-m)^2}
=\frac1{m^2}\left[
\sum_{k\ne0,m}\frac1{k^2}
+\sum_{k\ne0,m}\frac1{(k-m)^2}
-2\sum_{k\ne0,m}\frac1{k(k-m)}
\right].
\]

A soma cruzada é **exatamente** \(2/m^2\). Para demonstrar, dividimos em três partes absolutamente convergentes:
\[
\sum_{k=m+1}^\infty\frac1{k(k-m)}=\frac{H_m}m,\quad
\sum_{k=-\infty}^{-1}\frac1{k(k-m)}=\frac{H_m}m,\quad
\sum_{k=1}^{m-1}\frac1{k(k-m)}=-\frac{2H_{m-1}}m.
\]
A soma é \(2(H_m-H_{m-1})/m=2/m^2\). Como
\[
\sum_{k\ne0,m}k^{-2}=2\zeta(2)-m^{-2}
\]
e \(\zeta(2)=\pi^2/6\), segue a **identidade exata**
\[
\boxed{\quad
\sum_{k\in\mathbb Z\setminus\{0,m\}}
\frac{d}{(kd)^2((k-m)d)^2}
=\frac1{d^3}\left(\frac{4\zeta(2)}{m^2}-\frac6{m^4}\right)
=\frac{2\pi^2}{3m^2d^3}-\frac6{m^4d^3}.
\quad}\tag{P88-UNIFORM}
\]

Isso é exatamente a RHS do Lema 6, pois \(\delta_a=\delta_b=d\), \(h=md\).
Portanto:
- **A igualdade em malha uniforme infinita é provada analiticamente.**
- Para o segmento finito \(\{-M,\dots,m+M\}d\), a soma é estritamente menor.
- A diferença positiva entre o valor infinito e a soma finita é majorada uniformemente por \(2/(3M^3d^3)\), pois as duas caudas omissas são no máximo \(\sum_{k>M} k^{-4}\le1/(3M^3)\).
- O coeficiente \(-3\) é **ótimo no limite uniforme** para o termo corretivo fixando a parte positiva: substituí-lo por \(-(3+\varepsilon)\) falharia quando \(M\to\infty\). Não basta ignorar o termo negativo na continuação da prova.

`compensation.py` verifica exatamente as três parcelas harmônicas, a decomposição racional, a fórmula da soma finita e os limites da cauda por **frações exatas**. Os limites inferiores e superiores de \(\pi\) são obtidos pela identidade de Machin e pelas séries alternadas de arctan; nenhum `float` decide um gate matemático.

## 5. Testes independentes e escopo

O código percorre configurações fornecidas em P86, grades de até 31 pontos, uma progressão exponencial e configurações com forte aglomeração. Para cada par \((a,b)\), registra:
- \(\mathcal A_a\), \(\mathcal A_b\), \(C_{a,b}\) e \(S_{a,b}\) exatamente como frações;
- o compensador \(E_{a,b}\) com o limite **inferior racional** para \(\pi\);
- a identidade \(\,E_{a,b}=h^2(\mathrm{RHS}_{L6}-S_{a,b})\);
- a validade do limite relaxado P88-WEAK;
- exemplos com \(C_{a,b}<0\) que impedem provas por simples positividade do termo cruzado.

**Advertência:** verificar um número finito de pares com folga positiva, mesmo em aritmética exata, **não prova** que \(E_{a,b}\ge0\) para qualquer sequência real. A prova analítica geral de \(\mathrm{P88-GATE}\) continua aberta nesta implementação. O Lema 6 como teorema publicado é válido, mas ainda **importado** na cadeia independente.

Os arquivos S3 reenviados anteriormente foram deduplicados na P84; não se alteram nem se executam novamente os certificados B1. Os certificados P84/P85 e o intervalo Suzuki \(18\le t\le79\) continuam preservados e **condicionais**.

## 6. Próximo gate matemático P89

A prova primária de Preissmann utiliza seu **Lema 5** para comparar somas de recíprocos após uniformização dos pontos exteriores (Lema 6, pp. 108–110). A prioridade seguinte é reconstruir e verificar essa comparação de somas, resolvendo \(E_{a,b}\ge0\) para espaçamentos arbitrários. Somente depois a P87 e a P88 poderiam retirar o status externo do Lema 6 e da constante \(m_0\) para \(F_X\).

**Outras lacunas:** resultados de Ramaré–Zúñiga para \(S_\varepsilon(X)\) e a faixa \(3\times10^{10}\le X<3\times10^{12}\). Nenhuma parte desta etapa demonstra RH global.
