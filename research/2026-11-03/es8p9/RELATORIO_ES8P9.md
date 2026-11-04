# NUMERIA–ES8-P9 — Integração independente, cauda e projeção finita

**Data:** 03/11/2026  
**Escopo:** 30 entradas arquimedianas de modos de mesma paridade; malha racional [0,160]; `N=6,8,10`.  
**Status:** reprocessamento de 768.000 valores em somas exatas de intervalos, sob lemas analíticos identificados; matrizes finitas positivas condicionais. **Nenhuma afirmação sobre RH.**

## 1. Arquivos originais e independência

A versão E7/V16 histórica foi disponibilizada pelo usuário como ZIP e BIN. Da subpasta Fase 4 extraiu-se o conjunto de 30 registros, validando cada SHA-256 contra o resumo original. A fonte Fase 4 tem SHA-256:

`b0d7d335534239a6863e655295a658282551d334b0808ee801710a4b5699d8fe`.

A ES8-P8 comparara dois algoritmos que produzem intervalos nodais. Nesta etapa a soma foi reconstruída usando **apenas o avaliador independente** `IndependentFixed` de 35 casas, com 128 parcelas de digama, sem importar nenhum dos valores das 30 somas antigas. O arquivo histórico só entra na verificação de interseção posterior.

A independência é **numérica e algébrica** no arquimediano; as identidades de funções especiais e os parâmetros de Fourier são os mesmos. Os termos de primo e de projeção da etapa final reutilizam a convenção histórica V16/V17 e **não** são uma rederivação independente da fórmula global de Weil.

## 2. Quadratura racional com pesos assinados

Se `w_j` são os 16 pesos que integram exatamente os polinômios interpoladores nos pontos `x_j=j/15`, então

\[
\int_{0}^{160}q_{ik}(r)\,dr \in
  \frac1{10}\sum_{c=0}^{1599}\sum_{j=0}^{15}w_j\,Q_{ik}\left(\frac{c}{10}+\frac{j}{150}\right)
  +[-E,E],
\]

onde `Q_ik(r)` é o intervalo de valor nodal obtido com arredondamento dirigido e

\[
E=\frac{4154891126177792}{209922127267236328125}
\approx 1{,}979253535711606\times10^{-5}
\]

é a cota de interpolação de Cauchy deduzida na ES8-P6.

Os pesos são determinados *de novo* resolvendo em frações exatas o sistema

\[
\sum_{j=0}^{15}w_j\Bigl(\frac{j}{15}\Bigr)^k=\frac1{k+1}
\qquad(0\le k\le15).
\]

Uma implementação independente das integrais das bases de Lagrange fornece os mesmos 16 pesos. Dez desses pesos são positivos e **seis negativos**; por isso, se `w_j < 0`, a contribuição ao extremo inferior usa o extremo **superior** do valor nodal. As operações inteiras preservam todas as inclusões.

Cada lote de 50 células grava SHA-256 da sequência de 30 somas inteiras do lote. O resultado JSON permite reproduzir e conferir as 32 sequências sem arredondamentos decimais.

## 3. Comparação das 30 integrais

- `768.000/768.000` parcelas da quadratura somadas independentemente em 32 lotes.
- `30/30` intervalos de integrais nodais têm interseção com os históricos.
- `30/30` intervalos com erro de interpolação de Cauchy têm interseção com os históricos.
- Maior largura nodal independente (entre as 30 entradas): aproximadamente `2,24809 × 10^-6`, observada na entrada `G_10,10`; a largura histórica é menor por usar avaliador distinto.

Duas entradas ilustrativas, limites decimais arredondados para fora:

| Entrada | Intervalo nodal independente em [0,160] | Intervalo com erro de Cauchy |
|---|---|---|
| G11 | [-1,428208480997293; -1,428208457813065] | [-1,428228273532650; -1,428188665277708] |
| G22 | [-0,102578990877565; -0,102578898271258] | [-0,102598783412922; -0,102559105735901] |

O arquivo `RESULTADOS_ES8P9_INTEGRAIS.json` guarda os limites dirigidos com 24 casas decimais. A cauda `r>=160` é incluída separadamente via fórmula da Proposição 6, rededuzida na ES8-P6. Nenhum número proveniente de `mpmath` foi usado como extremidade de intervalo.

## 4. Propagação à matriz projetada

Para a matriz real simetrizada, os pares de paridades opostas são exatamente zero. Em cada entrada restante:

\[
G_{ij}\in [I^{\rm novo}_{ij,-}-T_{ij},\ I^{\rm novo}_{ij,+}+T_{ij}],
\quad T_{ij}\text{ é cota racional da cauda}. 
\]

A contribuição aritmética é a do primo 2, com normalização V16/V17:

\[
P_{ij}=-\frac{2\log2}{\sqrt2}C_{ij}(\log2),
\qquad M_N=B^{\mathsf T}(G+P)B.
\]

As colunas de `B` cancelam os dois momentos de polos, por identidades algébricas verificadas nas etapas anteriores. A projeção utiliza a base das funções-teste originais, não uma base ortonormal.

Para uma matriz intervalar simétrica `M` e ponto médio `C`, o cálculo usa

\[
\|M-C\|_2\le \delta:=\max_i\sum_j\operatorname{rad}(M_{ij}).
\]

A decomposição LDL racional verifica que `C-(δ+0,05)I` é positiva definida:

| N | Dimensão projetada | Cota superior de δ | Menor pivô LDL deslocado (aprox.) | Status |
|---|---:|---:|---:|---|
| 6 | 4×4 | 0,001995097 | 0,179596922 | Positiva condicionalmente |
| 8 | 6×6 | 0,004420901 | 0,177171118 | Positiva condicionalmente |
| 10 | 8×8 | 0,008356446 | 0,173235573 | Positiva condicionalmente |

**Importante:** `δ` é o limite da norma do erro propagado no *modelo matricial*; não é erro estatístico, e a margem `0,05` é medida nas coordenadas dessa base. A passagem dessa margem à norma funcional exigiria a matriz de Gram (tratada à parte na ES8-P4).

## 5. O que foi estabelecido, e o que não foi

**Estabelecido computacionalmente, sob os lemas externos explicitados:**

- Integridade dos dados de origem e procedimento reproduzível de soma dos 30 intervalos.
- Verificação exata dos 16 pesos e tratamento adequado de seus sinais.
- Verificação independente das 30 integrais com interseção entre dois métodos.
- Certificados racionais condicionais para as três matrizes finitas projetadas, com termo primo e cauda separados.

**Não estabelecido:** formalização integral em Lean dos lemas de digama, restos e interpolação; equivalência desta família finita ao conjunto de todas as funções-teste de Weil; controle uniforme em `N→∞`; positividade universal de Weil; Hipótese de Riemann.

## 6. Próxima etapa ES8-P10

Separar a prova analítica dos programas de soma: formalizar ou revisar independentemente as cotas de digama e da interpolação em uma segunda ferramenta, além de tentar expandir as dimensões da família com limites uniformes, **sem** extrapolar a partir dos casos `N≤10`.
