# NUMERIA–RIEMANN P92-F — Controle espectral de Δ por número de escalas

**Base:** P92-E, [PR #94](https://github.com/jbarumby/scicalc/pull/94). A P92-E usa agrupamentos uniformes separados e obtém uma cota em função do **número de agrupamentos** K. Esta nova etapa elimina completamente a necessidade de blocos, desde que o número de **valores distintos das folgas ao vizinho mais próximo** seja pequeno.

## 1. Teorema exato e novo (todos os pontos e parâmetros reais)

Sejam N>=2, x1<...<xN reais distintos, e `delta_i = min_{j != i}|x_i-x_j|`. Seja **m** o número de valores *distintos* entre delta1,...,deltaN. Defina a matriz real simétrica de diagonal nula da P88:

\[
\Delta_{ij} = \frac{\sqrt{\delta_i\delta_j}(\sqrt{\delta_i}-\sqrt{\delta_j})^2}{2(x_i-x_j)^2}\quad(i\ne j).
\]

**Teorema P92-F.1:** sem supor uniformidade, agrupamento ou separação especial,

\[
\boxed{\|\Delta\|_2 \le (m-1)\zeta(2)=(m-1)\frac{\pi^2}{6}.}
\tag{F1}
\]

Para **m<=3**, em particular,
\[
\boxed{\|\Delta\|_2\le\pi^2/3
< (\pi^2/3)\sqrt{6/5}.}
\tag{F2}
\]

O número total N de pontos e o número de agrupamentos visuais K podem ser arbitrariamente grandes. As folgas podem ser distribuídas irregularmente e as três razões entre escalas podem ser arbitrárias. O teorema **não afirma** que m<=3 para configurações gerais nem controla m arbitrário com a constante almejada.

### Resultado refinado por linha

Se o ponto i pertence à escala `a=delta_i`, sejam `m_<(a)` e `m_>(a)` os números de *valores de escala* estritamente menores e maiores que a, respectivamente. Então a prova oferece:

\[
\rho_i:=\frac{(\Delta p)_i}{p_i}
\le m_<(a)+\zeta(2)m_>(a),\qquad p_i=\sqrt{\delta_i}.
\tag{F3}
\]

Como 1<=zeta(2), F1 decorre imediatamente.

## 2. Demonstração do teorema P92-F.1

**Lema de separação:** para qualquer valor b que apareça entre as `delta_j`, dois índices distintos da classe `J_b={j: delta_j=b}` têm `|x_j-x_k|>=b` (pois b é a distância ao vizinho mais próximo de cada um). Fixado i de escala a, todo ponto `j in J_b`, `j != i`, obedece `|x_i-x_j|>=max(a,b)`. Em cada lado de x_i, os pontos de J_b ficam ordenados a distâncias pelo menos `max(a,b)+kb`, k=0,1,... . Estas afirmações independem da disposição dos outros pontos.

**Schur com peso p_i=sqrt(delta_i):** como Δ é simétrica, não negativa entrada a entrada, e p é positivo,

\[
\|\Delta\|_2\le\max_i\rho_i,\quad
\rho_i=\sum_{j\ne i}\frac{\delta_j(\sqrt{\delta_i}-\sqrt{\delta_j})^2}{2(x_i-x_j)^2}.
\tag{F4}
\]

A classe b=a não contribui. Para uma classe b!=a, a contribuição de J_b a rho_i é

\[
R_{a\to b}=\frac{b(\sqrt a-\sqrt b)^2}{2}
\sum_{j\in J_b}\frac1{(x_i-x_j)^2}.
\]

**Caso 1: b>a.** As distâncias são >=b e os pontos de J_b ficam b-separados. Há no máximo um por lado em cada posição da malha de referência; assim,

\[
\sum_{j\in J_b}|x_i-x_j|^{-2}
\le 2\sum_{k=1}^{\infty}(kb)^{-2}=\frac{2\zeta(2)}{b^2}.
\]

Com u=sqrt(a/b) em [0,1), obtemos

\[
R_{a\to b}\le\zeta(2)(1-u)^2\le\zeta(2).
\tag{F5}
\]

**Caso 2: 0<b<a.** Agora as distâncias iniciais são >=a e a separação dos pontos de J_b é >=b. A soma decrescente é majorada por

\[
\sum_{j\in J_b}|x_i-x_j|^{-2}
\le 2\sum_{k\ge0}(a+kb)^{-2}
\le 2\left(a^{-2}+\int_0^\infty(a+bt)^{-2}dt\right)
=2(a^{-2}+(ab)^{-1}).
\]

Com u=sqrt(b/a) em (0,1), resulta

\[
R_{a\to b}\le (1-u)^2(1+u^2)\le1.
\tag{F6}
\]

A última desigualdade segue de `f(u)=(1-u)^2(1+u^2)`, pois

\[
f'(u)=-2(1-u)(1-u+2u^2)\le0\quad(0\le u\le1),\qquad f(0)=1.
\]

Somando (F5) e (F6) sobre **valores distintos** b, obtemos F3, e o teste de Schur implica F1. Para m<=3, F2 é imediata. **Essa é a prova universal da classe de até três escalas: nenhum teste numérico substitui esses lemas.**

## 3. Consequência da P88/P87 para a desigualdade de Hilbert ponderada

Na P88 corrigida, o argumento de Rayleigh com o sinal apropriado do comutador fornece

\[
\|B\|_2^2\le\max_i Q_i+2\|C\|_2,\qquad C=D+\Delta,\quad \|D\|_2\le\max_i Q_i.
\]

Com o **lema de energia** universal `max_i Q_i <= pi^2/3` da P87 (cuja demonstração Laplace exige a errata registrada nas P89/P91), e F2, segue para configurações com m<=3:

\[
\boxed{\|B\|_2^2\le 3\,(\pi^2/3)+2\,(\pi^2/3)=\frac53\pi^2,\qquad
\|B\|_2\le\pi\sqrt{5/3}.}
\tag{F7}
\]

Como `5/3 < 1 + (2/3)*sqrt(6/5)`, trata-se de uma constante ainda mais forte que o limite de Preissmann **somente nessa classe de até três valores distintos das folgas delta_i**. Essa dedução é condicional à validade da prova analítica reparada do elo de energia P87; os testes de π da P87/P88 que usavam a cota SUPERIOR não devem ser citados como prova desse elo.

## 4. Certificação mecânica independente e limitações

[certificate.py](certificate.py) usa somente `fractions.Fraction`, `math.isqrt`, funções combinatórias inteiras e aritmética exata. Como regressão, toma folgas consecutivas em `{1,4,9}`: raízes quadradas das folgas e delta_i são então racionais, tornando cada entrada de Δ, cada soma de Schur e cada pivô LDL uma fração exata. Há 360 sequências de 2 a 5 folgas e 8 casos adicionais maiores; **368 configurações**, com **40 fatorações LDL exatas** e comparação individual contra (F5)/(F6). Os testes com folgas quadradas não restringem a prova analítica, válida para qualquer real positivo.

Para auditoria sem invocar π como float, o verificador usa os limites racionais estritos `zeta(2)<5/3` (pois π²<10), `π>25/8` e `sqrt(6/5)>109/100`, garantindo a folga da cota residual da P88 também sob um teto racional. Esses limites são auxiliares de regressão, não alteram a constante exata em F1.

A CI executa testes Python normal e `-O`, gera JSON com as margens e um SHA-256. Reproduzir:

    python -m unittest -v NUMERIA_P92F_THREE_SCALES.test_p92f
    python -O -m unittest -v NUMERIA_P92F_THREE_SCALES.test_p92f
    python -m NUMERIA_P92F_THREE_SCALES.certificate --out p92f-certificate.json

**Não demonstrados:** estimativa residual global de Δ para arbitrariamente muitos valores distintos de delta_i; generalização universal de Preissmann sem hipóteses sobre as escalas; Lema 5/6 completo; P93 geométrica de cinco lacunas; RH. A P92-F é uma extensão espectral independente de K do P92-E, e a PR deve permanecer draft.

**Próximo gate:** estudar como combinar muitos valores distintos de delta_i sem fazer a cota crescer proporcionalmente a m, aproveitando a distribuição real dos seus valores e a localização dos pontos. Não presumir que qualquer configuração pode ser aproximada por três escalas com erro de norma uniformemente pequeno.
