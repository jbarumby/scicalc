# NUMERIA–RIEMANN P86 — Preissmann: decomposição espectral e auditoria independente

**Objetivo:** continuar a prioridade aberta na P85: ler a prova original de Emmanuel Preissmann, isolar exatamente seus lemas analíticos uniformes e auditar separadamente as identidades de matriz. **Não repetir** os 42 limiares b4* já certificados no PR #85.

## Fontes examinadas

1. [E. Preissmann, *Sur une inégalité de Montgomery–Vaughan*, L'Enseignement Mathématique 30 (1984), 95–113](https://www.e-periodica.ch/cntmng?pid=ens-001:1984:30::29), fonte **primária**: Lema 1, pp. 96–99; Lemas 2–5, pp. 99–108; Lema 6, pp. 108–110; prova espectral, pp. 111–113. O texto extraído da digitalização contém degradações em fórmulas; não inferir igualdades com índices ilegíveis.
2. [E. Carneiro e F. Littmann, *Monotone extremal functions and the weighted Hilbert's inequality*, Pacific J. Math. (2023), p. 3](https://ems.press/content/serial-article-files/44548): transcreve explicitamente a constante de Preissmann `C=pi*sqrt(1+(2/3)*sqrt(6/5))` e demonstra uma versão **mais fraca**, de constante `2*pi`, por Fourier. Não substituir a constante necessária da P84 por essa versão.
3. [Ramaré, Lema 6.5](https://ramare-olivier.github.io/Maths/DensityEstimate-34.pdf): usa exatamente `m0=sqrt(1+2/3 sqrt(6/5))` e a melhoria para coeficientes reais já auditada na P85.

## 1. Isolamento da matriz

Para números reais distintos `x_1,...,x_N`, com `N>=2`, defina

```
delta_i = min_(j!=i) |x_i-x_j| > 0,
D = diag(delta_1,...,delta_N),
K_ij = 1/(x_i-x_j) se i!=j, K_ii=0,
B = D^(1/2)*K*D^(1/2).
```

A matriz `K` é **antissimétrica**; `B` também. A desigualdade ponderada de Preissmann se reduz a provar `||B||_2<=pi*m0`, pois `w_i=a_i/sqrt(delta_i)` dá

```
|sum_(i!=j) a_i*conj(a_j)/(x_i-x_j)|
   <= ||B||_2 * sum_i |a_i|^2/delta_i.
```

Esta redução, usando a desigualdade de Cauchy e o operador normal `iB`, é algébrica e **não depende** dos lemas de espaçamento de Preissmann.

## 2. Identidade racional para a verificação espectral

Para `x_i` **racionais**, temos `delta_i` e `K` racionais. A matriz

```
H(rho) := rho*D^(-1) + K*D*K
```

é racional e simétrica. A identidade

```
D^(-1/2) * (rho I + B^2) * D^(-1/2) = H(rho)
```

demonstra

```
||B||_2^2 < rho   <=>   H(rho) é definida positiva.
```

O cálculo dirigido não precisa de raízes quadradas ou autovalores aproximados. Escrevendo `A=K*D` (similar a B), confirma-se exatamente a identidade `H(rho)*D=rho*I+A^2`.

O programa `spectral_rational.py` executa **LDL sem pivotamento sobre frações exatas**. Se todos os pivôs forem positivos, ele reconstrói `H=L diag(pivots)L^T` **entrada por entrada**, sem `float`, tolerâncias ou limites numéricos supostos.

Para os exemplos finitos, adota-se o **limite racional mais rigoroso**

```
rho = ((157/50)*(263/200))^2 < (pi*m0)^2,
```

pois `pi>157/50` (série alternada de arctan e identidade de Machin, ver P83) e `m0>263/200` por comparação de quadrados (ver P85).

O certificado de uma configuração concreta com `H(rho)>0` prova rigorosamente que `||B||_2<pi*m0` para aquela configuração — **não** para todos os `N`.

## 3. Um resultado geral que não usa Preissmann

Se `2<=N<=5`, então para quaisquer `x_i` reais distintos:

```
delta_i<=|x_i-x_j|, delta_j<=|x_i-x_j| (i!=j)
=> |B_ij|=sqrt(delta_i delta_j)/|x_i-x_j|<=1.
```

Assim, `||B||_infinity<=N-1`, `||B||_1<=N-1` porque B é antissimétrica, e

```
||B||_2 <= sqrt(||B||_1||B||_infinity)
        <= N-1 <= 4
        < (157/50)*(263/200) < pi*m0.
```

**Esta prova é válida uniformemente para todo conjunto real de 2 a 5 pontos.** É uma contribuição independente verificável; não fecha o caso de dimensão arbitrária necessário ao teorema de valor médio.

## 4. Cadeia original de Preissmann: registro das lacunas

| Elo na fonte original | Papel matemático | P86 |
|---|---|---|
| Lema 1 (pp. 96–99) | maximização de `sum delta_k f(x_k)` por espaçamento uniforme, para f decrescente/convexa sob hipóteses de derivadas; gera cotas de energia com `pi^2/3` | **não demonstrado independentemente**; identificado como gate |
| Lemas 2–5 (pp. 99–108) | análise variacional de padrões de espaçamento e extremizadores | **não demonstrados independentemente** |
| Lema 6 (pp. 108–110) | controle da soma de três denominadores e de termos cruzados necessários à cota final | **não demonstrado independentemente** |
| Prova espectral (pp. 111–113) | antissimetria, diagonalização de `iB` e estimativas dos termos de `B^2` | fatoração e congruência **auditadas exatamente** para entradas racionais |
| Constante `m0` | consequência dos limites uniformes anteriores, não apenas da positividade finita | `pi*m0` e seus limites racionais **verificados**; teorema geral **permanece importado** |

Uma observação útil: o *post* de [MathOverflow sobre a desigualdade ponderada](https://mathoverflow.net/questions/444337/on-the-weighted-hilbert-inequality-a-la-montgomery-and-vaughan) destaca a cadeia `||B||^2 <= pi^2/3 + 2||C||` e a dificuldade em controlar `C_ij=delta_i^(3/2) delta_j^(1/2)/(x_i-x_j)^2`. É uma explicação informal, **não** substitui os Lemas 1–6 publicados.

## 5. O que foi testado

A CI executa o auditor racional para **oito configurações finitas**, incluindo grades uniformes, triangulares e quadráticas e pontos fortemente agrupados. Cada configuração inclui seus `delta_i`, pivôs LDL exatos e SHA-256. A CI não aceita `status=THEOREM_PROVED` para o caso arbitrário: exige as bandeiras analíticas `false`.

Verificações adicionais: simetria de `H`, antissimetria de `K`, identidade de similaridade, reconstrução do LDL, entrada adversarial com limiar `rho` deliberadamente pequeno, ausência de nós repetidos, e limite geral para `N<=5`.

Os certificados não empregam `python-flint` porque essa etapa é **exatamente racional**, mais forte do que intervalos para essas identidades finitas; a computação Arb de b4* em P85 permanece preservada.

## 6. Implicação para NUMERIA–RIEMANN

A entrada de Preissmann usada no limite `F_X` é rigorosamente um teorema publicado, mas esta auditoria **não a redemonstra completamente**. A P86 identifica as partes matemáticas essenciais ainda abertas: sobretudo o **Lema 1** sobre energia dos espaçamentos e o **Lema 6** sobre três denominadores.

A próxima fase útil é **P87**, priorizando a reconstrução independente e mecanicamente verificável do Lema 1 de Preissmann (cota uniforme `sum_(j!=i) delta_i delta_j/(x_i-x_j)^2 <= pi^2/3`). Só depois devemos atacar o Lema 6 e a combinação espectral em dimensões arbitrárias.

Os resultados P78–P85, a tabela S2 robusta, a cobertura de Suzuki até `t=79` (finita e condicional), e a restrição `X>=3e12` continuam intactos.

**Nenhum resultado estabelece a Hipótese de Riemann global.**
