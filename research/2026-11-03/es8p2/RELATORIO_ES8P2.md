# NUMERIA–RIEMANN — ES8-P2: paridade e certificação racional do termo primo projetado

**Data:** 03/11/2026  
**Repositório-base:** `jbarumby/scicalc`, etapa ES8-P, commit `2bc6a293a67213bc6cdf8321bab2ab2f0bf4cbdc`.  
**Status científico:** resultado **finito e parcial**, com intervalos racionais para **somente a parcela do primo 2** na convenção V17. Identificação com a matriz integral original E7/V16 e positividade global de Weil: **NÃO ESTABELECIDAS**. Nenhuma prova de RH.

## 1. Entradas, fórmulas e convenção

Com `L=9/20`, para `k=1,...,N` a base real é

\[
 f_k(x)=\frac{\sin\big(k\pi(x+L)/(2L)\big)}{\sqrt L}\,\mathbf1_{[-L,L]}(x),\quad
 F_k(r)=\int_{-L}^{L}f_k(x)e^{irx}\,dx.
\]

Para funções reais, `F_k(-r)=\overline{F_k(r)}`. Como `Ω(-r)=Ω(r)` e as transformadas têm decaimento suficiente, segue para cada par de modos:

\[
\frac1{2\pi}\int_{-\infty}^\infty \Omega(r)F_i(r)\overline{F_j(r)}\,dr
=\frac1\pi\int_0^\infty \Omega(r)\operatorname{Re}\!\left(F_i(r)\overline{F_j(r)}\right)dr.
\]

**Correção da ES8-P:** na metade positiva da reta, a expressão real correspondente à convenção V17 tem de usar `Re(F_i conjugate(F_j))`. A V24 escreve o produto complexo não simetrizado. Não se presume a igualdade desses *integrandos* no semieixo; a igualdade é da integral bilateral com a integral unilateral **da parte real**, sob convergência adequada.

Sob reflexão, `f_k(-x)=(-1)^(k+1) f_k(x)`. Assim, se `i` e `j` têm paridades opostas, `F_i conj(F_j)` é puramente imaginário no eixo real e a contribuição **real** arquimediana entre essas paridades é zero. A parte **simétrica** do termo primo também é exatamente zero entre as duas paridades. Como a projeção de polo preserva paridade, a matriz real projetada deve se decompor em blocos ímpar e par. Isto é uma identidade estrutural, **não** uma garantia de positividade.

## 2. Cálculo analítico do termo aritmético

O suporte duplo é `2L=9/10`, e os logaritmos foram confinados por frações racionais para verificar

\[
\log 2 <\frac9{10}<\log 3.
\]

Na convenção e coeficientes herdados **explicitamente** de V17, só a potência de primo `n=2` entra na janela. Defina

\[
 C_{ij}(u)=\int_{-L+u}^L f_i(x)f_j(x-u)dx,
 \quad C^{\rm sym}_{ij}(u)=\tfrac12(C_{ij}(u)+C_{ji}(u)),\quad u=\log2.
\]

Sua forma fechada, para `ω_k=kπ/(2L)` e `0≤u≤2L`, é

\[
 C^{\rm sym}_{ij}(u)=
\begin{cases}
0,& i-j \text{ ímpar},\\[1mm]
\left(1-\frac{u}{2L}\right)\cos(\omega_i u)+\frac{\sin(\omega_i u)}{i\pi}, &i=j,\\[1mm]
\frac{2}{\pi(i^2-j^2)}\left[i\sin(\omega_j u)-j\sin(\omega_i u)\right],&i\ne j,\ i-j\text{ par}.
\end{cases}
\]

A parcela que auditamos é

\[
 \mathrm{Prime}_{ij}=-2\frac{\log2}{\sqrt2}\,C^{\rm sym}_{ij}(\log2)
 =-\sqrt2\log2\,C^{\rm sym}_{ij}(\log2).
\]

Os momentos dos dois polos são cancelados algebricamente nas colunas da projeção `B`, com índices `k=3,...,N` e pivô `p=1` se `k` ímpar, senão `p=2`:

\[
 B_{kk}=1,\qquad
 B_{pk}=-\frac{k}{p}\frac{p^2\pi^2+L^2}{k^2\pi^2+L^2}.
\]

Essa é a construção analítica V17; **a equivalência literal com a matriz original E7/V16 permanece um gate externo**. Certificamos `B^T Prime B` para `N=6,8,10`, não `B^T(G+P)B`.

## 3. Mecanismo de certificação

1. `π`: Machin `π=16 atan(1/5)-4 atan(1/239)` com resto de série alternada dirigido corretamente. `log 2` e `log 3`: série positiva `atanh` com majorante de cauda racional. `sqrt 2`: raiz inteira sobre escala `10^40`.
2. `sin`, `cos`: polinômios de Taylor exatos de ordens 161/160 no argumento racional central, com resto de Taylor majorado por fatoriais. O afastamento entre o argumento real e o centro é propagado pela desigualdade global de Lipschitz `|sin'|,|cos'|≤1`.
3. Combinação das fórmulas fechadas e das colunas de projeção usando `fractions.Fraction` e arredondamento externo em cada operação. Todos os extremos da **matriz inteira** ficam no JSON como inteiros multiplicados por `10^40`; há também apresentações decimais arredondadas para fora.
4. Os **18 testes** cobrem direção do resto de arctangente (par e ímpar), enclosures de constantes, avaliação independente `mpmath` (somente diagnóstico), sinais de entradas, simetria, anulamento exato de paridades cruzadas, projeção e controles adversariais.

**Achado de auditoria:** a implementação ES8-P anterior tinha o sentido do intervalo da série alternada `atan` invertido antes do arredondamento. O teste adicional ES8-P2 prova em frações que, **para os parâmetros específicos da ES8-P (45 termos e arredondamento para 30 casas)**, o intervalo final arredondado de `π` ainda engloba o intervalo rigoroso corrigido. Isso mitiga o defeito naquele número, mas não dispensa auditar outros lemas ou transferir a equivalência de matrizes.

## 4. Resultados finitos certificados

Para cada `N`, a dimensão projetada é `N-2`. As entradas de paridade oposta são **literalmente** `[0,0]` no JSON, não valores pequenos comparados a uma tolerância.

| N | dimensão | pares cruzados anulados (entradas ordenadas) | largura máxima interna dos intervalos |
|---:|---:|---:|---:|
| 6 | 4 | 8 | < 2 × 10^-24 |
| 8 | 6 | 18 | < 2 × 10^-24 |
| 10 | 8 | 32 | < 2 × 10^-24 |

Duas entradas do **primo projetado**, independentes de aumentar `N` dentro desta família, são:

- `Prime_proj[3,3]` entre **−0,155593305650730642** e **−0,155593305650730641**.
- `Prime_proj[4,4]` entre **+0,088507025352479562** e **+0,088507025352479563**.
- `Prime_proj[3,4] = 0` **exatamente**, por paridade.

As extremidades decimais são intencionalmente mais largas que os intervalos internos, que ficam preservados como inteiros escalonados no JSON. A entrada negativa **não** refuta Weil: ela é da parcela prima **isolada**, sem a contribuição arquimediana e demais normalizações. A entrada positiva também não prova positividade da matriz completa.

As matrizes real-simetrizadas projetadas admitem blocos de paridade `2+2` (N=6), `3+3` (N=8) e `4+4` (N=10), permitindo em princípio certificar LDL por blocos **após** calcular o componente arquimediano completo.

## 5. Literatura recebida: três CSVs

`References_LeapSpace (10).csv`, `(11).csv`, `(12).csv` totalizam **51 registros**, **44 DOIs distintos**, **7 duplicações excedentes**; `REFERENCIAS_TRIADAS_ES8P2.csv` preserva os autores/DOIs de origem e corrige o nome de **Alain Connes** em uma cópia normalizada, sem alterar a planilha original. Classificação: 4 prioritários, 12 relevância média, 23 contextuais, 5 baixa prioridade. Não foi realizada leitura integral dos 44 artigos.

- **Connes (1999), DOI 10.1007/s000290050042**: primeira prioridade de fundamentação da fórmula de traço adélica e das parcelas locais. A proposição dele sobre RH depende de estabelecer uma identidade global apropriada; não resolve o gate das matrizes numéricas.
- **Carneiro–Chirre–Milinovich (2019), DOI 10.5565/PUBLMAT6321906**: importante para aproximação com Fourier limitado e desenho de funções-teste; as cotas sobre zeta descritas no resumo são **sob RH**, o que exige precaução contra circularidade.
- **Chirre–Gonçalves–de Laat (2020), DOI 10.1016/j.aim.2019.106926**: métodos de programação semidefinida e correlação de pares; ferramenta conceitual para otimização certificada. O uso de RH como hipótese em resultados derivados deve ser explicitamente separado de uma tentativa de prová-la.
- **Newman–Ryavec–Shure (1979), DOI 10.1016/0022-1236(79)90049-1**: candidatos a operadores integrais na teoria dos números; é preciso localizar a construção específica antes de relacioná-la ao NUMERIA.
- **Katz–Sarnak (1999), DOI 10.1090/S0273-0979-99-00766-1**: contexto de simetrias espectrais; não é uma fórmula para as entradas V17/V24.
- **Yakaboylu (2024), DOI 10.1088/1751-8121/ad4c2d**: candidato de Hamiltoniano; a hipótese relativa ao domínio/autoadjunção ainda requer auditoria rigorosa e não foi incorporada como prova do NUMERIA.

## 6. Gate matemático ES8-P3

O próximo passo NÃO é aumentar `N` sem controle, mas computar rigorosamente, para o **mesmo espaço projetado real**, cada

\[
 A_{ij}=\frac1\pi\int_0^\infty\Omega(r)\operatorname{Re}\big(F_i(r)\overline{F_j(r)}\big)\,dr
\]

com majorantes de cauda, erros intervalares e, para as versões suavizadas, a mesma função de corte de V24. Após conferir a base e a identidade da matriz original E7/V16, combinar `A` + `Prime` + polos (que devem cancelar na projeção), testar autovalores/LDL com intervalos e **só então** confrontar as margens com o erro de suavização V24. Sem esses passos, a porta `FULL_WEIL_MATRIX_POSITIVE` fica fechada.

## Conclusão

**Certificado:** entradas da parcela do primo `2` projetada, com paridade e intervalos racionais, para três dimensões finitas da base V17.  
**Não certificado:** integrais arquimedianas matriciais, identidade com os dados originais E7/V16, positividade de todas as matrizes, classe infinita de funções-teste e a Hipótese de Riemann.
