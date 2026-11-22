# NUMERIA–RIEMANN P92-E — Lema espectral multicluster para a matriz Δ

**Base:** [NUMERIA P92-D, PR #93](https://github.com/jbarumby/scicalc/pull/93). A P92-D provou o caso especial *âncora + um cluster uniforme*, de posto no máximo dois. Esta fase controla **vários clusters uniformes diferentes**, com comprimentos arbitrários, por estimativas de Hilbert–Schmidt e uma matriz pequena de normas de blocos.

## 1. Teorema P92-E (resultado analítico, classes de parâmetros reais)

Sejam \(K\ge1\) agrupamentos de pontos reais disjuntos, ordenados, cada agrupamento \(r\) contendo \(n_r\ge2\) pontos em progressão aritmética de passo \(h_r>0\). Denote por \(G_r\) a distância do último ponto do agrupamento \(r\) ao primeiro de \(r+1\), e suponha

\[
G_r\ge\max(h_r,h_{r+1}),\qquad r=1,\ldots,K-1.
\tag{separação}
\]

Em toda a configuração, a distância de qualquer ponto do grupo \(r\) ao vizinho mais próximo é exatamente \(\delta_i=h_r\). Defina a matriz simétrica de diagonal zero da P88:

\[
\Delta_{ij}=\frac{\sqrt{\delta_i\delta_j}(\sqrt{\delta_i}-\sqrt{\delta_j})^2}{2(x_i-x_j)^2},\quad i\ne j.
\]

**Teorema P92-E.1 (demonstrado):** para TODOS os valores reais que satisfaçam as hipóteses, sem restrições sobre \(n_r\) nem sobre a razão entre os \(h_r\),

\[
\boxed{\|\Delta\|_2\le\frac{K-1}{\sqrt8}.}
\tag{P92-E1}
\]

Além disso, para qualquer quantidade \(K\) de agrupamentos, a norma tem uma estimativa **a posteriori**, por uma matriz \(K\times K\) de coeficientes explícitos descrita abaixo. A P92-E.1, em particular, demonstra o **limite residual desejado pela P88** para todos os conjuntos da classe com \(K\le11\):

\[
\boxed{\|\Delta\|_2\le\frac{10}{\sqrt8}
<\frac{\pi^2}{3}\sqrt{\frac65}.}
\tag{P92-E2}
\]

Com a energia universal \(\max_i Q_i\le\pi^2/3\) da P87 (com errata analítica) e o argumento corrigido de Rayleigh da P88, também se obtém a **constante de Preissmann** para essa classe restrita: \(\|B\|_2\le\pi\sqrt{1+\frac23\sqrt{6/5}}\). Isso NÃO é uma demonstração do Lema 6 geral de Preissmann nem de RH.

## 2. Cota exata para o bloco entre dois agrupamentos

Fixe dois grupos \(r<s\). Seja \(S_{rs}\) a distância do último ponto de \(r\) ao primeiro de \(s\), e escreva \(h=h_r\), \(g=h_s\). A hipótese de separação implica \(S=S_{rs}\ge\max(h,g)\), inclusive para grupos não adjacentes.

Se a posição no primeiro grupo é indexada de **direita para esquerda** por \(p\) e, no segundo, de **esquerda para direita** por \(q\), então a diferença entre os pontos vale exatamente \(S+ph+qg\), com \(p=0,\ldots,n_r-1\) e \(q=0,\ldots,n_s-1\). Portanto o bloco de \(\Delta\) possui entradas positivas

\[
T_{pq}=\frac{\kappa}{(S+ph+qg)^2},\qquad
\kappa=\frac{\sqrt{hg}(\sqrt h-\sqrt g)^2}{2}.
\]

Sua norma de operador é limitada pela norma Hilbert–Schmidt:

\[
\|T\|_2^2\le\|T\|_{HS}^2
=\kappa^2\sum_{p=0}^{n_r-1}\sum_{q=0}^{n_s-1}
\frac1{(S+ph+qg)^4}.
\]

**Comparação integral bidimensional exata:** como \(t\mapsto(S+t)^{-4}\) é decrescente, cada soma com índice positivo é majorada por uma integral no retângulo que precede o ponto da malha. Separando \((0,0)\), os dois eixos positivos e o interior positivo, temos:

\[
\begin{aligned}
\sum_{p,q\ge0}(S+ph+qg)^{-4}
&\le S^{-4}+\frac1h\int_0^\infty(S+t)^{-4}dt
+\frac1g\int_0^\infty(S+t)^{-4}dt\\
&\quad+\frac1{hg}\int_0^\infty\!\int_0^\infty(S+t+u)^{-4}dt\,du\\
&=\frac1{S^4}+\frac{h^{-1}+g^{-1}}{3S^3}+\frac1{6hgS^2}.
\end{aligned}
\]

Essa expressão é o **quadrado da cota de norma entre os dois grupos**, denominada \(W_{rs}^2\):

\[
\boxed{W_{rs}^2=\kappa^2\left(\frac1{S^4}
+\frac{h^{-1}+g^{-1}}{3S^3}+\frac1{6hgS^2}\right).}
\tag{bloco}
\]

## 3. A constante universal por par \(1/\sqrt8\)

Normalize \(x=h/S\), \(y=g/S\). Assim \(0<x,y\le1\) e

\[
W_{rs}^2=\frac{(\sqrt x-\sqrt y)^4}{4}
\left[xy+\frac{x+y}{3}+\frac16\right].
\]

Escreva \(p=\max(x,y)\le1\) e \(q=\sqrt{\min(x,y)/p}\in[0,1]\). Então

\[
\begin{aligned}
W_{rs}^2
&=\frac{p^2(1-q)^4}{4}
\left[p^2q^2+\frac{p(1+q^2)}3+\frac16\right]\\
&\le\frac{(1-q)^4(3+8q^2)}{24}\le\frac18.
\end{aligned}
\]

A segunda desigualdade segue de que \(f(q)=(1-q)^4(3+8q^2)\) é decrescente em \([0,1]\):

\[
f'(q)=-4(1-q)^3(3-4q+12q^2)\le0,
\]

pois \(3-4q+12q^2=12(q-1/6)^2+8/3>0\), e \(f(0)=3\). Portanto, \(\|T\|_2\le W_{rs}\le1/\sqrt8\) **independentemente de ambos os tamanhos dos grupos e de todos os parâmetros reais**.

## 4. Matriz comparadora e teorema para K agrupamentos

Tome a matriz simétrica \(W\in\mathbb R^{K\times K}\), com \(W_{rr}=0\) e \(W_{rs}\) dado pelo limite quadrático **(bloco)**. Para um vetor dividido em subvetores \(v=(v_1,\ldots,v_K)\), ponha \(t_r=\|v_r\|_2\). A desigualdade triangular e o limite por bloco dão

\[
\|(\Delta v)_r\|_2\le\sum_{s\ne r}W_{rs}t_s.
\]

Somando os quadrados e usando \(\sum_rt_r^2=\|v\|_2^2\):

\[
\boxed{\|\Delta\|_2\le\|W\|_2\le
\max_r\sum_{s\ne r}W_{rs}\le\frac{K-1}{\sqrt8}.}
\]

O limite intermediário \(\max_r\sum_{s\ne r}W_{rs}\) pode ser usado **mesmo para K>11**, contanto que seja suficientemente pequeno; a expressão para cada \(W_{rs}\) é explícita e rigorosa. Isso fornece uma certificação específica por configuração sem restrição no número de pontos de cada agrupamento.

Para \(K\le11\), \((K-1)^2/8\le100/8=25/2\). Usando a desigualdade elementar conhecida \(\pi>25/8\), e a verificação racional

\[
(25/8)^4-375/4=6625/4096>0,
\]

obtemos

\[
\frac{25}{2}<\frac{2}{15}(25/8)^4
<\frac{2\pi^4}{15}
=\left(\frac{\pi^2}{3}\sqrt{\frac65}\right)^2,
\]

concluindo rigorosamente **(P92-E2)**. A derivação usa a **cota inferior** de \(\pi\), na direção adequada, corrigindo a deficiência dos testes B3 da P88.

## 5. Consequência espectral restrita

Para \(B\) da P87 e a matriz \(C=D+\Delta\) da P88, o argumento de Rayleigh com o comutador corrigido diz

\[
\|B\|_2^2\le\max_i Q_i+2\|C\|_2
\le\frac{\pi^2}{3}+2\left(\frac{\pi^2}{3}+\|\Delta\|_2\right)
=\pi^2+2\|\Delta\|_2.
\]

Assim, na classe de \(K\) agrupamentos acima:

\[
\boxed{\|B\|_2^2\le\pi^2+\frac{K-1}{\sqrt2}.}
\]

Se \(K\le11\), a cota P92-E2 também conduz à constante geral de Preissmann \(\pi\sqrt{1+(2/3)\sqrt{6/5}}\), **mas somente para esta classe estruturada**. Para \(K>11\), qualquer certificado por configuração que prove \(\max\sum W_{rs}\le(\pi^2/3)\sqrt{6/5}\) produz a mesma conclusão.

## 6. Certificados racionais e CI

O [certificate.py](certificate.py) constrói oito famílias de exemplo, contendo de 1 a 14 grupos, e verifica **162 pares de blocos** com `fractions.Fraction` puro. Para manter os coeficientes de \(\Delta\) racionais, os passos nos exemplos são quadrados de racionais: \(h_r=s_r^2\), embora o teorema analítico seja válido para qualquer passo real. Os limites de raízes \(W_{rs}\) são *superiores* racionais gerados por \(\mathrm{isqrt}\), sem flutuação. O JSON inclui somas exatas de Hilbert–Schmidt, limites integrais e os limites das somas de linhas.

O arquivo [test_p92e.py](test_p92e.py) contém **8 testes** repetidos em Python normal e com `-O`, incluindo uma regressão espectral independente por **LDLᵀ exata em frações** para matrizes \(\Delta\) pequenas. Essas regressões validam a implementação, mas **não são a prova universal**, que está nas Seções 2–5.

Reproduzir na raiz do repositório:

    python -m unittest -v NUMERIA_P92E_MULTICLUSTER.test_p92e
    python -O -m unittest -v NUMERIA_P92E_MULTICLUSTER.test_p92e
    python -m NUMERIA_P92E_MULTICLUSTER.certificate --out p92e-certificate.json

## 7. Escopo e próximos gates

- **Demonstrado:** cota de Δ para qualquer K de agrupamentos uniformes separados, \(\|\Delta\|\le(K-1)/\sqrt8\); para \(K\le11\), a estimativa da P88 e a constante espectral de Preissmann seguem nessa **classe restrita**; para K arbitrário há um critério explícito por matriz W.
- **Não demonstrado:** controle universal de Δ para *todo* conjunto de pontos sem estrutura de blocos; Lema 5 e Lema 6 completos de Preissmann; P93 geométrica de cinco lacunas; qualquer ampliação da cobertura Ramaré–Zúñiga; RH.
- **Próximo gate:** procurar uma decomposição em clusters para uma configuração geral com erro de aproximação controlado e uma estimativa de \(\|W\|\) independente de K. **Não** presumir que a decomposição de qualquer configuração finita satisfaça a condição de separação nem que o limite \((K-1)/\sqrt8\) permaneça controlado para K arbitrário.
- Manter PR como **draft**, sem merge automático.
