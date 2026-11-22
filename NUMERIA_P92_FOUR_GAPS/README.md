# NUMERIA–RIEMANN P92 — Preissmann Lema 6 afiado com quatro lacunas internas

**Base:** [PR #91](https://github.com/jbarumby/scicalc/pull/91); P90 demonstrou pares de separação 1 e 2; P91 acrescentou separação 3. **P92 acrescenta separação 4**, ou três pontos interiores entre os índices extremos, com número e espaçamentos arbitrários de pontos exteriores.

## Teorema P92 (escopo matemático exato)

Para toda configuração finita de pontos reais distintos x1<...<xn, delta_i = distância ao vizinho mais próximo, e todo par a<b com b-a<=4, vale a desigualdade afiada:

\[
\sum_{k\notin\{a,b\}}\frac{\delta_k}{(x_k-x_a)^2(x_k-x_b)^2}
\le\frac{\pi^2}{3h^2}\left(\frac1{\delta_a}+\frac1{\delta_b}\right)
-\frac{3(\delta_a+\delta_b)}{h^4},\qquad h=x_b-x_a.
\]

Este resultado, em sua classe de índices, dispensa a importação do Lema 5 ou Lema 6 geral de Preissmann: a prova é uma redução exata a 133.363 sinais racionais, com código, geometria e reconstrução explícitos abaixo. **Não prova separação >=5 nem RH.**

## 1. Redução interior P89/P90

Normalize h=1. Para 4 lacunas interiores positivas, escreva a+b+c+d=1. Os três pontos interiores estão em a, a+b, 1-d. Seu núcleo é:

\[
M_4=
\frac{\min(a,b)}{a^2(1-a)^2}
+\frac{\min(b,c)}{(a+b)^2(c+d)^2}
+\frac{\min(c,d)}{d^2(1-d)^2}.
\]

Das caudas exteriores de P89 e da monotonicidade do orçamento D em P90, basta provar M_4<=D(a)+D(d), com

\[
Z=\zeta(2)=\pi^2/6,\quad
S(s)=\sum_{j\ge1}\frac1{j^2(1+js)^2},\quad
D(s)=\frac{2Z-S(s)}s-3s.
\]

As distâncias dos extremos normalizadas são <=a e <=d. Como D decresce, a substituição por a,d é a pior situação para o orçamento. A contribuição dos pontos exteriores já foi coberta em P89.

## 2. Minorante de D especialmente afiado para quatro lacunas

A interpolação de P91 sem s=1/4 não preserva a igualdade uniforme em m=4. Inserimos o novo nó usando a identidade exata P90:

\[
S(0)=Z,\quad S(1/4)=2Z-355/144,\quad
S(1/3)=2Z-31/12,\quad S(1/2)=2Z-11/4,\quad S(1)=2Z-3.
\]

Convexidade de S: cada corda entre nós sucessivos majoriza S. No primeiro intervalo o coeficiente de Z é 1/s-4>=0, permitindo o minorante rigoroso Z>3/2. Nos demais intervalos Z se cancela exatamente. Resulta:

\[
\boxed{D(s)\ge d_4(s)=
\begin{cases}
3/(2s)+139/36-3s,&0<s\le1/4,\\
19/(9s)+17/12-3s,&1/4\le s\le1/3,\\
9/(4s)+1-3s,&1/3\le s\le1/2,\\
5/(2s)+1/2-3s,&1/2\le s\le1.
\end{cases}}
\]

As quatro expressões concordam nas fronteiras. Para a malha uniforme, M_4(1/4,1/4,1/4,1/4)=164/9=2d_4(1/4)=2D(1/4), preservando exatamente o termo afiado -3.

O alvo P92 é, portanto, a **desigualdade puramente racional para todos os reais admissíveis**:

\[
\boxed{M_4(a,b,c,d)\le d_4(a)+d_4(d)},\quad a,b,c,d>0,\ a+b+c+d=1.
\]

## 3. Domínio 3D e triangulação racional completa

Por reflexão (a,b,c,d) para (d,c,b,a) basta considerar a<=d. Use coordenadas (a,b,d), c=1-a-b-d. O domínio fechado Omega é a metade do simplex: a,b,d>=0, a+b+d<=1, a<=d, volume exato **1/12**.

Partimos pelas três comparações adjacentes a<=b ou a>=b; b<=c ou b>=c; c<=d ou c>=d; e pelas faixas a em [0,1/4],[1/4,1/3],[1/3,1/2]; d em [0,1/4],[1/4,1/3],[1/3,1/2],[1/2,1]. São 96 células formais; **35 têm interior tridimensional não vazio**. Elas são interseções de doze semiespaços definidos por formas lineares racionais, com interiores disjuntos.

No [certificate.py](certificate.py), todas as triplas de planos são intersectadas **com determinantes de Fraction**, e somente as soluções que satisfazem cada semiespaço são retidas. Cada poliedro convexo é triangulado por uma decomposição estrela: ordenar as faces por comparações polares **racionais sem atan2/float**, triangular cada face com um leque, unir ao centroide dos vértices. O centroide de uma região tridimensional convexa fica em seu interior. Os tetraedros de interiores disjuntos cobrem o poliedro. O volume de cada tetraedro é o determinante racional absoluto dividido por 6.

O algoritmo encontra **246 tetraedros**, verifica que seus vértices pertencem às células e que a soma exata dos volumes é **1/12**, o volume total exigido. Essas condições, combinadas à partição por semiespaços, certificam a cobertura integral.

## 4. Polinômios exatos homogêneos de grau 13

Em cada tetraedro, sejam u1,u2,u3,u4>=0, H=sum(ui)>0, e

\[
A=\sum_i a_i u_i,\quad B=\sum_i b_i u_i,\quad
D_0=\sum_i d_i u_i,\quad C=H-A-B-D_0.
\]

Aqui D_0 é a forma homogênea da última lacuna, distinta da função D(s). Escreva L=H-A, R=H-D_0, T=A+B, V=C+D_0. Dentro de cada região, min(A,B)=m1, min(B,C)=m2, min(C,D_0)=m3 são formas **lineares fixas** por suas comparações. A forma do minorante em cada faixa é alpha/s+beta-3s para a primeira extremidade e gamma/s+eta-3s para a segunda.

Multiplicando [d_4(a)+d_4(d)-M_4] pelo fator **estritamente positivo**

\[
H^{13}a^2d^2(1-a)^2(1-d)^2(a+b)^2(c+d)^2
\]

e substituindo a=A/H, b=B/H, c=C/H, d=D_0/H, resulta o polinômio homogêneo de grau 13:

\[
\begin{aligned}
P={}&\alpha H^2 AD_0^2L^2R^2T^2V^2+
\gamma H^2A^2D_0L^2R^2T^2V^2\\
&+(\beta+\eta)HA^2D_0^2L^2R^2T^2V^2
-3(A+D_0)A^2D_0^2L^2R^2T^2V^2\\
&-H^4m_1D_0^2R^2T^2V^2
-H^4m_2A^2D_0^2L^2R^2
-H^4m_3A^2L^2T^2V^2.
\end{aligned}
\]

O verificador expande todos os monômios nas quatro coordenadas baricêntricas **usando somente Fraction**: em 246 tetraedros produz **133.363 coeficientes não nulos e todos são estritamente positivos**, sem coeficientes negativos. Ausências equivalem a coeficientes zero. Há controles racionais independentes que comparam P expandido ao denominador positivo vezes a expressão original num ponto interior de cada tetraedro.

Para u_i>=0 todo monômio é não negativo. Assim P>=0 em cada tetraedro e a cobertura certifica d_4(a)+d_4(d)-M_4>=0 para **todos os reais positivos**, não somente exemplos racionais. Pela redução P89/P90, segue o Lema 6 para b-a=4.

**Natureza científica:** certificação algébrica finita executável e auditável; não afirmar que foi verificada por um assistente formal de provas, nem extrapolar o resultado para todos os m.

## 5. Reproduzir

Execute a partir da raiz do repositório:

    python -m unittest -v NUMERIA_P92_FOUR_GAPS.test_p92
    python -O -m unittest -v NUMERIA_P92_FOUR_GAPS.test_p92
    python -m NUMERIA_P92_FOUR_GAPS.certificate --out p92-summary.json
    python -m NUMERIA_P92_FOUR_GAPS.certificate --full --out p92-full.json

O JSON resumido contém vértices de cada tetraedro, volume, quantidade de coeficientes e SHA-256 canônico dos coeficientes. A opção --full inclui também todos os expoentes e valores exatos dos 133.363 coeficientes. O GitHub Actions executa testes em Python normal e com -O, verifica os marcadores de limitação e exporta manifesto SHA-256 + JSON completo.

**Fronteira restante:** b-a>=5, Lema 5 e Lema 6 gerais, constante espectral global sem teorema externo, Ramaré–Zúñiga pendente e RH. **Próxima P93:** cinco lacunas interiores ou uma prova indutiva por compensação variacional; a P92 NÃO constitui indução geral. Manter PR draft, sem merge automático.
