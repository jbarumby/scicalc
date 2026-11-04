# NUMERIA–RIEMANN — Relatório matemático ES8-P13

**Data:** 04/11/2026. **Escopo:** matrizes projetadas concretas, \(N=6,8,10\).

## 1. Questão verificada

A ES8-P12 provou no Lean teoremas gerais de LDL. A ES8-P13 produz dados racionais *concretos* a partir das inclusões ES8-P10 e exige que a máquina confira a igualdade da fatoração e a positividade dos pivôs, sem substituir isso por simples impressão de números positivos.

## 2. Fonte e reconstrução

São lidas as 30 entradas originais da ES8-P9; a ES8-P10 reduz a cota de interpolação sem modificar os valores nodais. A aritmética projetada é reconstruída por intervalos de escala \(10^{90}\), incluindo cauda, primo \(2\) e projeção \(B^\mathsf{T}(G+P)B\). A identidade das matrizes de dimensões 4 e 6 como submatrizes principais do caso de dimensão 8 foi verificada **exatamente, nos extremos brutos da aritmética de ponto fixo**.

Definindo \(A_N\) como intervalo simétrico e \(C_N\) como aproximação de ponto médio arredondada a múltiplos de \(10^{-4}\), usamos as entradas racionais explícitas de \(C_N-3I/50\).

| Parâmetro N | Dimensão | Cota máxima exata de erro por linha (aprox. de leitura) | Menor pivô LDL (aprox. de leitura) |
|---|---:|---:|---:|
| 6 | 4 | 0,001936050 | 0,171600000 |
| 8 | 6 | 0,004378245 | 0,171600000 |
| 10 | 8 | 0,008330007 | 0,171600000 |

Todos os **valores decisivos** estão armazenados como frações exatas em `RESULTADOS_ES8P13.json`. Os números da tabela são apresentados apenas para leitura.

## 3. Certificado Lean

Para cada tamanho, o script deriva \(D=\operatorname{diag}(d_i)\), \(U\) triangular superior com diagonal unitária e verifica exatamente

\[
C_N-\frac{3}{50}I=U^\mathsf{T}DU,\qquad \forall i\; d_i>0.
\]

O arquivo Lean contém também o argumento de injetividade da matriz triangular e a igualdade entre a forma quadrática e a soma de quadrados. Isto prova \(x^\mathsf T(C_N-3I/50)x>0\) para todo vetor racional não nulo \(x\). Como os coeficientes são racionais e a forma é contínua, a mesma positividade estrita na decomposição vale para vetores reais não nulos (os pivôs não se alteram). A prova formal exibida foi escrita para vetores racionais; o argumento de extensão a reais não foi formalizado nesse arquivo.

## 4. Ponte de perturbação e limites

Em Python com `fractions.Fraction`, verifica-se, linha por linha, que

\[
\delta_N:=\max_i \sum_j \sup_{a_{ij}\in A_{ij}}|a_{ij}-(C_N)_{ij}|<1/100.
\]

Para qualquer matriz real simétrica \(M\) com \(M_{ij}\in A_{ij}\), segue \(\|M-C_N\|_2\le \delta_N\). Como \(C_N-0,06I\succ0\) e \(\delta_N<0,01\), concluímos \(M-0,05I\succ0\). Isto depende da **validade externa dos intervalos** da fase numérica. A identidade do verdadeiro operador de Weil com essa forma matricial também depende da normalização arquimediana, primo 2 e projeção herdadas.

## 5. Pendências científicas

- Traduzir para Lean o elo de inclusão analítica entre \(G_{ij}\), os restos da digama, a quadratura e os extremos intervalares.
- Traduzir a verificação dos extremos de aritmética fixa e a robustez espectral ao próprio Lean.
- Investigar a positividade **uniforme em N** e a densidade/admissibilidade de funções com suportes arbitrários. A prova dessas propriedades não decorre de três matrizes positivas, mesmo que o Lean as compile.

**Nenhuma afirmação de demonstração da Hipótese de Riemann.**
