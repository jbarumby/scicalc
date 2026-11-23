# NUMERIA–RIEMANN P92-M — Separação 5/4 sem cardinalidade máxima e barreira para o expoente de dois blocos

**Continuação rigorosa da P92-L, PR #101.** A P92-M trata somente a matriz simétrica de entradas não negativas 

\[
\Delta_{ij}=\frac{\sqrt{\delta_i\delta_j}(\sqrt{\delta_i}-\sqrt{\delta_j})^2}{2(x_i-x_j)^2}\ (i\ne j),\quad \Delta_{ii}=0,
\qquad \delta_i=\min_{k\ne i}|x_i-x_k|.
\]

Os resultados seguintes não demonstram a cota universal de \(\|\Delta\|_2\), os Lemas 5–6 completos de Preissmann nem a Hipótese de Riemann. **A prova é analítica para todos os parâmetros reais da classe especificada**; testes computacionais servem apenas para auditoria.

## M.1 — Separação subquadrática com cardinalidades inteiramente livres

Sejam \(K\ge1\) agrupamentos ordenados de pontos distintos. O agrupamento \(r\) possui \(m_r\ge2\) pontos em progressão aritmética de passo \(h_r>0\), com **nenhuma cota sobre** \(K\), \(m_r\) ou \(h_r/h_s\). Para \(r<s\), denote por \(D_{rs}\) a distância da última posição do agrupamento \(r\) à primeira do agrupamento \(s\).

**Teorema P92-M.1.** Se, para todo \(r<s\),

\[
\boxed{D_{rs}\ge (s-r)^{5/4}\max(h_r,h_s),}\tag{M-sep}
\]

então, para *todos* os passos reais e cardinalidades:

\[
\boxed{\|\Delta\|_2<\frac5{\sqrt2}<\frac{\pi^2}{3}\sqrt{\frac65}.}\tag{M1}
\]

**Demonstração.** A condição com \(s=r+1\) implica isolamento de todos os agrupamentos: cada ponto do grupo \(r\) tem exatamente \(\delta_i=h_r\). Assim, os blocos diagonais de \(\Delta\) são nulos. A estimativa **analítica de dois blocos da P92-K** (revalidada na P92-L) estabelece, sem restrições de cardinalidade,

\[
\|\Delta_{rs}\|_2\le\frac{\max(h_r,h_s)}{2\sqrt2\,D_{rs}}
\le\frac{1}{2\sqrt2\,|r-s|^{5/4}}.\tag{M-pair}
\]

Para um vetor particionado \(z=(z_1,\ldots,z_K)\), escreva \(\eta_r=\|z_r\|_2\). A desigualdade triangular fornece \(\|\Delta z\|_2\le\|A\eta\|_2\), onde \(A\) é simétrica, com entradas nulas na diagonal e \(A_{rs}=1/[2\sqrt2|r-s|^{5/4}]\) fora dela. Como \(A\) é simétrica de entradas não negativas,

\[
\|\Delta\|_2\le\|A\|_2\le\max_r\sum_{s\ne r}A_{rs}
\le\frac1{\sqrt2}\sum_{k\ge1}k^{-5/4}
=\frac{\zeta(5/4)}{\sqrt2}.
\]

A função \(x^{-5/4}\) é estritamente decrescente; logo a comparação soma–integral dá

\[
\zeta(5/4)=1+\sum_{k=2}^\infty k^{-5/4}
<1+\int_1^\infty x^{-5/4}\,dx=5.
\]

Assim \(\|\Delta\|_2<5/\sqrt2\), sem depender de \(K\), do número total de pontos, da maior cardinalidade ou do contraste entre passos.

Para verificar **rigorosamente** a comparação com a meta da P88, elevar ao quadrado mostra que

\[
\frac5{\sqrt2}<\frac{\pi^2}{3}\sqrt{\frac65}
\quad\Longleftrightarrow\quad \pi^4>\frac{375}{4}.
\]

Basta \(\pi>25/8\), pois \((25/8)^4=390625/4096>375/4\). Esta cota inferior elementar segue da identidade de Machin \(\pi=16\arctan(1/5)-4\arctan(1/239)\) e das séries alternadas:

\[
\pi>16\Big(\frac15-\frac1{375}\Big)-\frac4{239}>\frac{25}{8}.
\]

A prova não usa amostragem numérica nem supõe \(K\) finito fixado a priori; ela vale para cada dimensão finita arbitrária. QED.

### M.1.1 — Princípio geral de decaimento

Com a hipótese mais geral \(D_{rs}\ge c|r-s|^p\max(h_r,h_s)\), para \(c\ge1\) e \(p>1\), a mesma prova produz

\[
\boxed{\|\Delta\|_2\le\frac{\zeta(p)}{\sqrt2\,c}
<\frac{1+1/(p-1)}{\sqrt2\,c}.}\tag{M-general}
\]

A escolha \(p=5/4,c=1\) é a primeira condição simples aqui certificada abaixo da meta residual. Não se afirma que \(5/4\) seja o melhor expoente possível.

## M.2 — Obstrução exata: a potência \(1/\rho\) por par não pode melhorar uniformemente

**Teorema P92-M.2 (negativo).** Não existem constantes \(C<\infty\) e \(\varepsilon>0\) tais que, **para todos** os pares de agrupamentos aritméticos de tamanhos e passos arbitrários e separação \(D\ge H=\max(a,b)\),

\[
\|W\|_2\le C\left(\frac{H}{D}\right)^{1+\varepsilon}.
\]

**Prova por família racional explícita.** Para cada inteiro \(n\ge2\), tome um agrupamento de \(n\) pontos com passo \(a=1\) e outro de \(n^3\) pontos com passo \(b=1/n^2\), separados entre as extremidades por \(D=n\). Assim \(H=1\), \(\rho=D/H=n\). A distância ao vizinho mais próximo é realmente \(a\) no primeiro bloco e \(b\) no segundo, porque a separação supera os passos internos.

Para índices a partir das extremidades, \(0\le j<n\), \(0\le k<n^3\), vale

\[
W_{jk}=\frac{\beta}{(n+j+k/n^2)^2},\qquad
\beta=\frac1{2n}\left(1-\frac1n\right)^2.
\]

Como \(n+j+k/n^2<3n\), **cada entrada** de \(W\) é maior que \(\beta/(9n^2)\). Aplicando \(W\) aos vetores constantes unitários de dimensões \(n\) e \(n^3\) (cujos produtos de normalização são \(1/n^2\)), obtemos

\[
\boxed{\|W\|_2\ge \frac{1}{n^2}\sum_{j,k}W_{jk}
>\frac{n^2\beta}{9n^2}
=\frac{(1-1/n)^2}{18n}
\ge\frac1{72\rho}.}\tag{M-barrier}
\]

Se existisse \(C\rho^{-1-\varepsilon}\) para todos os tamanhos, então \(1/(72n)\le Cn^{-1-\varepsilon}\), isto é, \(n^\varepsilon\le72C\) para todos os \(n\), uma contradição. QED.

**Interpretação crucial:** esse resultado descarta a melhoria do *expoente de uma estimativa universal isolada de dois blocos*. **Não** exclui uma prova de separação linear ou de uma cota espectral geral, pois a posição relativa de vários blocos pode proporcionar compensações e compatibilidades geométricas não exploradas pelo majorante par-a-par.

## M.3 — O que a P92-M efetivamente avançou

- **P92-K:** exigia \(D_{rs}\ge|r-s|^{2}H_{rs}\) para cardinalidades sem limite.
- **P92-L:** permitia separação sublinear \(|r-s|^{3/4}\), mas com \(m_r\le4\); com separação linear \(|r-s|\), permitia \(m_r\le14\).
- **P92-M:** exige apenas separação \(|r-s|^{5/4}\) para **cardinalidades sem limite**, além de demonstrar que melhorar o expoente de dois blocos para \(1+\varepsilon\) é impossível sem hipóteses extras.

A P92-M **não** demonstra uma classe que contenha integralmente as classes da P92-L, mas estabelece um novo ponto rigoroso no equilíbrio entre separação e cardinalidade.

## M.4 — Verificação e limites dos certificados

O pacote [certificate.py](certificate.py) reconstrói distâncias com `Fraction`, verifica cada desigualdade de separação \(D_{rs}^4\ge|r-s|^5H_{rs}^4\) exatamente (sem raízes irracionais), confere o Hilbert–Schmidt real dos blocos contra a cota da P92-K e executa LDL racional em matrizes completas selecionadas. Os exemplos com passos \(h_r=s_r^2\) para \(s_r\) racional permitem calcular as entradas de \(\Delta\) sem arredondamento.

Também testa o contraexemplo M-barrier por frações exatas e confirma uma cota inferior para vetores constantes. O código **não usa cálculos finitos para demonstrar os teoremas universais**. O JSON contém flags falsas para configurações arbitrárias, Lema 6 completo e RH.

Comandos na raiz do repositório:

```bash
python -m unittest -v NUMERIA_P92M.test_p92m
python -O -m unittest -v NUMERIA_P92M.test_p92m
python -m NUMERIA_P92M.certificate --out p92m-certificate.json
```

## M.5 — Próxima investigação P92-N

O limite de expoente em M.2 mostra que o caminho mais promissor não é substituir \(1/\rho\) por \(1/\rho^{1+\varepsilon}\) sem condições. Procurar **somas de operadores com geometria compartilhada**, regimes de blocos esparsos e **estimativas de quase ortogonalidade**: duas matrizes \(W_{rs}\) podem compartilhar vetores principais e exigir controle de correlações, não apenas das suas normas individuais.

Em particular, investigar se o arranjo linear ordenado e a aditividade das distâncias entre agrupamentos fornecem estimativas do operador global mais fortes que o máximo de soma das normas par-a-par.

**Status científico:** prova analítica apenas das classes e da obstrução descritas. A cota de \(\Delta\) em configurações inteiramente arbitrárias, os Lemas 5–6 gerais e a RH permanecem em aberto. Manter PR draft, sem merge automático.
