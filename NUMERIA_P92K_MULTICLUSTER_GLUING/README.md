# NUMERIA–RIEMANN P92-K — Colagem espectral de agrupamentos com escalas diferentes

**Base:** P92-J, [PR #99](https://github.com/jbarumby/scicalc/pull/99). Esta fase estende a trilha da matriz \(\Delta\) da P88 para **um número arbitrário de agrupamentos, com comprimentos e espaçamentos internos arbitrários**. Não utiliza nem demonstra o Lema 5/6 geral de Preissmann ou RH.

## K.1 — Teorema analítico uniforme de separação quadrática

Sejam \(K\ge1\) e \(m_1,\dots,m_K\ge2\) inteiros. Forme \(K\) agrupamentos ordenados de pontos reais distintos; o grupo \(r\) (\(1\le r\le K\)) é uma progressão aritmética de \(m_r\) pontos com passo \(h_r>0\). Para \(1\le r<s\le K\), seja

\[
D_{rs}=\min\{y-x:x\text{ no grupo }r,\ y\text{ no grupo }s\},
\]

a distância entre a extremidade direita do primeiro grupo e a extremidade esquerda do segundo. Suponha a condição **quantitativa e verificável**

\[
\boxed{D_{rs}\ge(s-r)^2\max(h_r,h_s)\qquad(1\le r<s\le K).}\tag{K-sep}
\]

Em particular, para grupos adjacentes \(D_{r,r+1}\ge\max(h_r,h_{r+1})\). Como todo agrupamento tem ao menos dois pontos, a distância do ponto \(x_i\) ao seu vizinho mais próximo é exatamente \(\delta_i=h_r\) para todos os pontos do grupo \(r\). Considere a matriz simétrica, de diagonal nula, da P88:

\[
\Delta_{ij}=\frac{\sqrt{\delta_i\delta_j}(\sqrt{\delta_i}-\sqrt{\delta_j})^2}{2(x_i-x_j)^2} \quad (i\neq j).
\]

**Teorema P92-K (prova completa abaixo):** para todo \(K\), todas as cardinalidades e quaisquer passos **reais** positivos sob (K-sep),

\[
\boxed{\|\Delta\|_2\le\frac{\zeta(2)}{\sqrt2}=\frac{\pi^2}{6\sqrt2}<\frac{10}{7}.}\tag{K-global}
\]

Não há dependência de \(K\), do número de pontos dentro de cada grupo, do maior contraste entre \(h_r\) ou da dimensão total da matriz. **A hipótese de separação quadrática é indispensável ao argumento; (K-global) não foi demonstrado para pontos arbitrários.**

## K.2 — Lema Hilbert–Schmidt para *dois* grupos arbitrários

Considere dois agrupamentos aritméticos, de cardinalidades \(m,n\ge2\), passos \(a,b>0\), com distância entre suas extremidades enfrentadas \(D\ge\max(a,b)\). Para os índices locais contados *a partir da fronteira* (\(j,k\ge0\)), a matriz de interação tem entradas exatas

\[
W_{jk}=\frac{\beta}{(D+ja+kb)^2},\qquad
\beta=\frac{\sqrt{ab}(\sqrt a-\sqrt b)^2}{2}.
\]

Seja \(H=\max(a,b)\), \(\rho=D/H\ge1\). Então

\[
\boxed{\|W\|_2^2\le\|W\|_{HS}^2\le\frac1{8\rho^2},\qquad
\|W\|_2\le\frac1{2\sqrt2\,\rho}.}\tag{K-pair}
\]

**Prova com integrais elementares.** A função \((D+ax+by)^{-4}\) é não crescente nas duas variáveis não negativas. Separando os termos de borda do reticulado positivo e comparando os termos interiores com integrais sobre quadrados unitários, temos

\[
\sum_{j,k\ge0}\frac1{(D+ja+kb)^4}\le\frac1{D^4}+\frac1{3aD^3}+\frac1{3bD^3}+\frac1{6abD^2}.\tag{K-sum}
\]

A soma finita para \(m,n\) é menor que a soma infinita. Sem perda, tome \(a=H\ge b\), e ponha \(u=\sqrt{b/a}\in(0,1]\). Então \(\beta=a^2u(1-u)^2/2\) e \(D=\rho a\). Multiplicando (K-sum) por \(\beta^2\), obtemos a desigualdade **exata para todos os parâmetros reais**:

\[
\begin{aligned}
\|W\|_{HS}^2
&\le(1-u)^4\left(\frac{u^2}{4\rho^4}
+\frac{1+u^2}{12\rho^3}+\frac1{24\rho^2}\right)\\
&\le\frac{(1-u)^4(3+8u^2)}{24\rho^2}
\le\frac1{8\rho^2}.
\end{aligned}
\]

Na última passagem, \(F(u)=(1-u)^4(3+8u^2)\) é decrescente em \([0,1]\), porque

\[
F'(u)=-(1-u)^3(12-16u+48u^2),\qquad
12-16u+48u^2=48(u-1/6)^2+\frac{32}{3}>0.
\]

Portanto \(F(u)\le F(0)=3\), e (K-pair) está provado. A constante \(1/(2\sqrt2)\) não é reivindicada como ótima.

## K.3 — Comparação de blocos e fechamento sem dependência de K

Separe \(\Delta\) na partição por agrupamentos: os blocos diagonais são identicamente zero porque \(\delta_i=h_r\) em cada grupo. Denote os blocos cruzados por \(W_{rs}\), e escreva o vetor original como \(z=(z_1,\dots,z_K)\), com \(z_r\) o vetor das coordenadas do grupo \(r\).

Para \(r\neq s\), pelo Lema K.2 e (K-sep):

\[
\|W_{rs}\|_2\le\frac{\max(h_r,h_s)}{2\sqrt2\,D_{rs}}
\le\frac1{2\sqrt2\,(s-r)^2}.
\]

Considere a matriz real não negativa \(A\) de dimensão \(K\) com \(A_{rr}=0\) e \(A_{rs}=1/[2\sqrt2\,(r-s)^2]\). A comparação dos blocos, por Cauchy–Schwarz e o vetor das normas \(\eta_r=\|z_r\|_2\), mostra que \(\|\Delta z\|_2\le\|A\eta\|_2\). Consequentemente

\[
\|\Delta\|_2\le\|A\|_2\le\max_r\sum_{s\ne r}A_{rs}
\le\frac1{2\sqrt2}\,2\sum_{k\ge1}\frac1{k^2}
=\frac{\zeta(2)}{\sqrt2}=\frac{\pi^2}{6\sqrt2}.
\]

A estimativa racional complementar \(\pi^2/(6\sqrt2)<10/7\) pode ser justificada sem números decimais: \(\zeta(2)<2\) e \(\sqrt2>7/5\); logo \(\zeta(2)/\sqrt2<2/(7/5)=10/7\). A cota \(\zeta(2)<2\) segue de \(1/n^2<1/[n(n-1)]\) para \(n\ge2\) e telescopagem. Assim \(\|\Delta\|_2<10/7\), muito abaixo do residual \((\pi^2/3)\sqrt{6/5}\) da P88. Isso conclui K-global.

## K.4 — Certificado racional independente

[certificate.py](certificate.py) reconstrói diretamente as distâncias aos vizinhos e todas as entradas de \(\Delta\) a partir das coordenadas com `Fraction`. Para obter raízes exatas das \(\delta_i\), os testes escolhem \(h_r=s_r^2\) com \(s_r\in\mathbb Q_{>0}\). Para cada par de grupos, recalcula o **quadrado exato da norma de Hilbert–Schmidt** e confere a cota racional \(\|W_{rs}\|_{HS}^2\le1/(8\rho_{rs}^2)\), com \(\rho_{rs}=D_{rs}/\max(h_r,h_s)\).

Uma matriz auxiliar **racional** usa entradas \(A'_{rs}=5/(14\rho_{rs})\), majorantes estritos de \(1/(2\sqrt2\rho_{rs})\), porque \((5/14)^2>1/8\). Os gates conferem o limite de soma de linha \(<10/7\), uma fatoração LDL exata de \((3/2)I_K-A'\) e, nos casos menores, uma fatoração LDL completa de \((3/2)I_N-\Delta\).

A configuração gerada para regressão usa separações adjacentes \(D_{r,r+1}=H r^2\), com \(H=\max h_r\). Para qualquer \(r<s\), \(D_{rs}\ge H\sum_{j=r}^{s-1}j^2\ge(s-r)^2\max(h_r,h_s)\), já que \(\sum_{j=1}^k j^2\ge k^2\) para todo \(k\ge1\).

**Reproduzir na raiz do repositório:**

    python -m unittest -v NUMERIA_P92K_MULTICLUSTER_GLUING.test_p92k
    python -O -m unittest -v NUMERIA_P92K_MULTICLUSTER_GLUING.test_p92k
    python -m NUMERIA_P92K_MULTICLUSTER_GLUING.certificate --out p92k-multicluster.json

Os testes finitos são regressões; **a prova universal sob a condição (K-sep)** é o argumento analítico das seções K.2–K.3. O JSON de saída documenta explicitamente o escopo e as quantidades efetivamente certificadas.

## K.5 — Status científico e próxima prioridade

- **Provado:** uma cota independente de dimensão, número de grupos e contrastes de escala sob a separação quadrática entre grupos; todas as cardinalidades e passos reais positivos.
- **Não provado:** uma partição canônica de toda configuração arbitrária satisfazendo (K-sep), nem a mesma cota para grupos arbitrariamente próximos; cota universal de \(\|\Delta\|_2\), Lemas 5/6 gerais de Preissmann, P93 geométrica de cinco lacunas e Hipótese de Riemann.
- **Próxima fase P92-L:** reduzir a exigência de separação quadrática, combinando distâncias reais entre grupos com uma matriz de comparação por norma, usando certificação de blocos e um núcleo descontínuo/sparsificação. Uma extensão a toda configuração requer controlar as interações de muitos blocos **sem presumir** (K-sep).
- **Política:** PR *draft*, sem merge automático, aguardando auditoria independente.