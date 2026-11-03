# NUMERIA–RIEMANN — ES8-P: projeção de polos, núcleo dos primos e cotejo V17/V24

**Data:** 03/11/2026  
**Status:** auditoria estrutural parcial, provas algébricas exatas e envelopes racionais de erro **condicionais ao modelo V24**; verificações numéricas independentes do termo primo; **sem certificação da matriz global de Weil nem de RH**.

## 1. Escopo e o que foi efetivamente executado

1. Prova simbólica de cancelamento das **duas avaliações de polo** nas colunas da projeção construída a partir das funções seno de V17, para `N=6,8,10`.
2. Construção analítica, pela integral de duas funções seno em intervalos sobrepostos, do núcleo simétrico de convolução que seleciona **somente** a potência de primo `2` para `L=9/20`.
3. Verificação numérica independente das fórmulas de convolução por quadratura de alta precisão, sem afirmar que sejam certificados intervalares.
4. Propagação **em frações exatas com arredondamento para fora** de limites das entradas do **modelo de erro V24** pela projeção `B* E B`, para `N=6,8,10` e `delta=10^{-4},10^{-5},10^{-6}`.
5. Detecção de **uma diferença real na fórmula arquimediana de V24**, quando confrontada com a simetrização de V17.
6. Triagem dos **12 registros do CSV (11 DOIs distintos, uma duplicação)** e identificação de referências espectrais e adélicas pertinentes.

**Não foram reconstruídas** todas as entradas `G+P` originais do arquivo E7 nem calculadas rigorosamente todas as integrais cruzadas arquimedianas para os modos V24. O ZIP E7 original e as matrizes E7/V16 completas não foram importados para este experimento. Por essa razão, o gate de identidade entre matrizes continua fechado.

## 2. Projeção dos dois polos: igualdade algébrica exata

Seja `L=9/20`, `f_k(x)=sin(k*pi*(x+L)/(2L))/sqrt(L)` em `[-L,L]`, `k>=1`. Em V17, os momentos dos dois polos, para o mesmo grupo de paridade, têm dependência no índice proporcional a

\[
 C_k(X)=\frac{k}{k^2X+L^2},\qquad X=\pi^2.
\]

Para `k>=3` escolhemos `p=1` se `k` é ímpar e `p=2` se `k` é par. A coluna projetada é `e_k+b_{p,k}e_p`, com

\[
 b_{p,k}(X)=-\frac{k}{p}\frac{p^2X+L^2}{k^2X+L^2}.
\]

A identidade racional

\[
 C_k(X)+b_{p,k}(X)C_p(X)=0
\]

é exata para qualquer `X>0`. Como os momentos em `s=0` e `s=1` diferem somente pelo fator de paridade comum em cada grupo, ambos se anulam em cada coluna. Resultado: dimensões projetadas `N-2=4,6,8` para `N=6,8,10`.

Isso **não significa** que a soma inteira de polos e primos foi transportada do E7 original: a prova se refere à construção analítica de V17, independente de arquivos E7 ausentes.

## 3. Primo `2`: núcleo real e simetrização

A desigualdade `log 2 < 2L=9/10 < log 3` foi demonstrada novamente com cotas racionais, sem números de ponto flutuante. Só o termo `n=2` pode entrar na soma de potências de primos com `|u|<2L`.

Defina a função de sobreposição para `0<=u<=2L` por

\[
 C_{ij}(u)=\int_{-L+u}^L f_i(x)f_j(x-u)\,dx.
\]

Expandir `sin A sin B = [cos(A-B)-cos(A+B)]/2` fornece uma fórmula fechada implementada, inclusive o limite contínuo para `i=j`. A parte simétrica é

\[
 C^{\rm sym}_{ij}(u)=\frac12[C_{ij}(u)+C_{ji}(u)].
\]

Na convenção de V17, a contribuição matricial do primo é

\[
 \mathrm{Prime}_{ij}=-2\frac{\log2}{\sqrt2}\,C^{\rm sym}_{ij}(\log2).
\]

A simetria `f_k(-x)=(-1)^{k+1}f_k(x)` torna a parte simétrica identicamente zero quando `i,j` têm paridades opostas. Para `i=j=1`, obtemos numericamente `Prime_11≈−0,0371502796970`; para `i=1,j=3`, `Prime_13≈−0,0900844624602`. Essas **entradas são diagnósticas mpmath**, não intervalos certificados; os testes de comparação direta via quadratura concordam a ~`10^-66` na precisão de trabalho.

## 4. Achado crítico: produto arquimediano na semirreta não foi simetrizado na V24

O código V17 fornece o integrando simétrico real

\[
 H_{ij}(r)=\tfrac12\{F_i(r)F_j(-r)+F_j(r)F_i(-r)\}
         =\operatorname{Re}(F_i(r)\overline{F_j(r)})
\]

para funções reais. Na documentação da V24 aparece, literalmente,

\[
 A^{V24}_{ij}=\frac1\pi\int_0^\infty\Omega(r)F_i(r)\overline{F_j(r)}\,dr,
\]

**sem** o operador `Re` ou simetrização. Para índices de paridades opostas, `F_i` é real e `F_j` é puramente imaginária (ou vice-versa), de modo que o integrando correto de V17 é identicamente zero. O produto sem simetrizar é, em geral, imaginário e não nulo.

Testemunho pontual `i=1,j=2,r=1`: `F_1(1) conjugate(F_2(1)) ≈ +0,1588947886567895863 i`, enquanto `H_12(1)=0` **exatamente por paridade**. A integral imaginária não simetrizada em `r∈[0,4]` é numericamente ~`−0,3511398118087`; a integral global, com cauda, não foi certificada. Isso demonstra uma **falha de identificação pontual dos kernels**, não, isoladamente, um erro numérico certificado para a entrada final. É necessário revisar a fórmula e o seu domínio (matriz real versus complexa), não apenas assumir equivalência entre V17 e V24.

Os envelopes **absolutos** de perturbação da V24 podem continuar úteis, desde que sejam reatribuídos à **mesma** forma bilinear corrigida. Os certificados prévios de positividade não podem ser transportados para a V24 sem essa revisão.

## 5. Propagação condicional do erro V24 ao espaço projetado

Se `|Delta K_ij|<=E_ij` for válido para **a mesma matriz modelada pela V24 e pela V17 corrigida**, com `E_ij>=0`, então para a projeção `B` definida acima,

\[
 |(B^*\Delta K B)_{ab}|\le (|B|^T E |B|)_{ab},\quad
 \|B^*\Delta K B\|_2\le\max_a\sum_b (|B|^T E |B|)_{ab}.
\]

O programa reconstrói a cota modelo V24

\[
E_{ij}=\underbrace{\frac{12B_0(e_i+e_j)+16(e_iD'_j+e_jD'_i)}{\pi}}_{\rm arquimediano}
+\underbrace{2B_0q^2(e_i+e_j)}_{\rm polos}
+\underbrace{\frac{2\log2}{\sqrt L}(e_i+e_j)}_{\rm primos}
\]

com `e_k=2δ² kπ/(2L√L)`, `D'_k=2L(kπ/(2L))²/√L+288(kπ/(2L))/√L`, `B_0=2√L` e `q=(1-L/2)^−1`. Limites de `π`, `log 2`, `log 3` e `sqrt L` vêm de séries racionais e raízes inteiras. Os valores da tabela são **limites superiores arredondados para fora**, sob as premissas do modelo V24.

| N | dimensão projetada | `δ=10^-4` | `δ=10^-5` | `δ=10^-6` |
|---:|---:|---:|---:|---:|
| 6 | 4 | 0,226025 | 0,00226025 | 0,0000226025 |
| 8 | 6 | 0,519805 | 0,00519805 | 0,0000519805 |
| 10 | 8 | 0,998524 | 0,00998524 | 0,0000998524 |

**Advertência:** os valores não estabelecem que a matriz projetada completa seja positiva definida. Isso exigiria **uma margem rigorosa do menor autovalor da mesma matriz** maior que os erros, além da verificação de todos os termos, sinais, fases, bases e coeficientes.

## 6. Triagem dos artigos recebidos (References_LeapSpace (9).csv)

São **12 registros e 11 títulos/DOIs distintos**. A referência `384` repete o artigo de Khuri `380`; a referência `385` não traz título no CSV, embora seu resumo identifique estudo das S-matrizes p-ádicas e adélicas. A tabela detalhada com DOI, utilidade e ressalvas está em `REFERENCIAS_TRIADAS_ES8P.csv`. Síntese:

- **Sierra 2007 (ref. 383)**, [DOI 10.1016/j.nuclphysb.2007.03.049](https://doi.org/10.1016/j.nuclphysb.2007.03.049): prioridade alta para o eixo espectral. Função de Jost e interação não local no `H=xp`, mas a identificação exata com todos os zeros **não** é estabelecida.
- **Sierra 2008 (ref. 386)**, [arXiv:0711.1063](https://arxiv.org/abs/0711.1063): prioridade alta para a reconstrução espectral; o modelo em `σ>1` não demonstra a construção desejada em `σ=1/2`.
- **S-matrizes adélicas (ref. 385)**, [DOI 10.1016/0370-2693(91)90868-Q](https://doi.org/10.1016/0370-2693(91)90868-Q): prioridade contextual para estruturas p-ádicas e fator local-global, **não** uma identidade pronta para `G+P` E7/V24.
- **LeClair 2025 (ref. 388)**, [arXiv:2406.01828](https://arxiv.org/abs/2406.01828): fluxo espectral e modelos de espalhamento com produto de Euler; verificar rigorosamente condições de limites e unitariedade.
- **Bender–Brody–Müller 2017 (ref. 389)**, [DOI 10.1103/PhysRevLett.118.130201](https://doi.org/10.1103/PhysRevLett.118.130201): operador com condição de contorno relacionada aos zeros; a autoadjunção é heurística e permanece obstáculo.
- **Khuri 1990 (refs. 380/384)**, [DOI 10.1016/0003-4916(90)90342-L](https://doi.org/10.1016/0003-4916(90)90342-L): potencial que **codifica** zeros em polos redundantes de espalhamento; não transforma essa codificação em prova de RH.
- **Khare 1997 (ref. 381)** e **Wardlaw–Jaworski 1989 (ref. 382)**: analogias de ressonância, fase e fórmulas de traço; úteis para controles espectrais, não para as entradas rigorosas E7/V24.
- **Tamburini–Licata 2021 (ref. 379)**: hipótese física Majorana/string; exploratória.
- **Savvidy 2020 (ref. 378)** e **Joffily 2007 (ref. 387)**: prioridade menor para o gate matricial atual.

**Nenhum dos artigos analisados fornece uma demonstração aceita da Hipótese de Riemann ou o certificado de normalização/primos/arquimediano das matrizes NUMERIA.** As classificações baseiam-se nos metadados e resumos do CSV, com verificação externa dos resumos e conclusões principais de artigos selecionados; não se afirma leitura integral de todos os artigos.

## 7. Gate de transição

**Aprovados:** prova algébrica das projeções de polos; suporte racional de primos; checagens cruzadas do núcleo de convolução; propagação condicional rigorosa de envelopes de erro V24 pela projeção; identificação de incompatibilidade pontual de fases.

**Não aprovados / próximos passos (ES8-P2):** correção da convenção arquimediana e testes de paridade em V24; certificação intervalar da própria matriz `G+P` com suas integrais cruzadas; confronto *por SHA* das fontes E7/V16; majorante de cauda de Fourier e da quadratura para todas as entradas; certificação LDL/eigenvalue com margem excedendo a perturbação projetada. Até lá, `NOT_ESTABLISHED` para positividade da matriz E7/V24, para positividade universal de Weil e para RH.

## Fontes principais

- `jbarumby/scicalc/v17/AUDITAR_E8K_V17_NORMALIZACAO.py`, fórmula de `H_ij`, projeções e termo primo.
- `jbarumby/scicalc/v24/README_V24.md` e `v24/AUDITAR_V24_CONTINUIDADE.py`, envelopes globais de perturbação **condicionais ao modelo**.
- J. C. Lagarias, *The Riemann Hypothesis: Arithmetic and Geometry*, Teorema 3.2, [texto](https://www.math.lsa.umich.edu/~lagarias/doc/mt-holyoke-rev.pdf): critério de Weil como equivalência incondicional conhecida.
- [Fórmula explícita e funções-teste](https://openmaths.org/courses/NT-EXPF/NT-EXPF-01.html), referência didática para sinais e polos.
- Fontes específicas nos DOIs das referências triadas.

**Conclusão:** a ES8-P produz novo diagnóstico matemático verificável e um limite projetado condicional para `N≤10`; ela não certifica a forma projetada de Weil E7/V24 nem a Hipótese de Riemann.
