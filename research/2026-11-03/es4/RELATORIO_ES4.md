# NUMERIA–RIEMANN ES4 — Operador não quadrático com lei espectral de Riemann–von Mangoldt

**Data:** 3 de novembro de 2026  
**Escopo:** continuação da auditoria ES3. Modelo de controle explícito; nenhuma afirmação de demonstração de RH, identificação com a zeta ou certificação intervalar de seus zeros.

## 1. Contexto e critério de aceitação

A ES3 descartou a família de Hamiltonianos quadráticos confinantes proposta no documento Spin–Zeta: espectro afim em `n`, enquanto o número de zeros não triviais da zeta até a altura `T` tem termo dominante `(T/2π) log(T/2π)-T/(2π)`. O objetivo de ES4 é construir **um controle espectral autoadjunto não quadrático** que supere somente essa barreira de crescimento, para depois testar a barreira mais exigente da fórmula explícita de Weil. O modelo foi construído **sem inserir zeros conhecidos nos seus elementos matriciais**.

## 2. Construção analítica (lema elementar)

Defina para `T>0` o termo suave da contagem

\[
N_0(T)=\frac{T}{2\pi}\left(\log\frac{T}{2\pi}-1\right)+\frac78.
\]

A derivada é `N_0'(T)=(2π)^{-1}log(T/(2π))`, positiva para `T>2π`. Para inteiros `n≥1`, a equação `N_0(E_n)=n` tem a solução única com `E_n>2π`, expressa por

\[
E_n=\frac{2\pi (n-7/8)}{W_0((n-7/8)/e)},
\]

onde `W_0` é o ramo real principal da função de Lambert (`W e^W=z`). O valor de `E_1` é aproximadamente `17,8478365128496`. A interpolação contínua obedece

\[
E'(x)=\frac{2\pi}{1+W_0((x-7/8)/e)}>0,
\]

que decresce para zero. Assim, a separação espectral diminui; esse comportamento difere dos níveis igualmente espaçados do modelo quadrático estável da ES3.

Seja `\mathcal H=\ell^2(\mathbb{N}_{\ge 1})` com base canônica `e_n`. Defina o operador diagonal

\[
H e_n=E_n e_n,\qquad
D(H)=\left\{c\in\ell^2:\sum_{n\ge1}E_n^2|c_n|^2<\infty\right\}.
\]

Esse domínio natural faz `H` autoadjunto (teorema elementar dos operadores diagonais reais). Como `E_n→∞`, o resolvente `(H-i)^{-1}` é compacto: ele é diagonal com entradas `(E_n-i)^{-1}→0`. O espectro é puramente discreto, de multiplicidade um, e para `T≥E_1` a contagem de autovalores satisfaz **exatamente**

\[
N_H(T)=\lfloor N_0(T)\rfloor;
\quad\text{para }T<E_1,\ N_H(T)=0.
\]

Consequentemente, `N_H(T)=N_0(T)+O(1)` para `T→∞`, que reproduz o **termo principal** de Riemann–von Mangoldt. Não há dependência de `\zeta`, primos, convoluções adélicas ou funções de teste de Weil nesta construção; portanto ela não estabelece identidade espectral com os zeros da zeta.

## 3. Teste adversarial: a contagem assintótica não identifica os zeros

Considere `H_1=H+I`, também autoadjunto com resolvente compacto. Seu espectro é `{E_n+1}`, estritamente distinto do de `H`. Para `T` suficientemente grande,

\[
N_{H_1}(T)=\lfloor N_0(T-1)\rfloor=N_0(T)+O(\log T).
\]

A segunda igualdade segue do teorema do valor médio e de `N_0'(T)=O(log T)`. Logo, **dois operadores diferentes**, com espectros distintos, compartilham a lei de contagem `N(T)=(T/2π)log(T/2π)-T/(2π)+O(log T)` válida para a zeta. Ter a densidade correta não seleciona os zeros individualmente.

## 4. Experimento computacional e limites

Executamos 12 testes automatizados e construímos os primeiros `J=25` níveis independentes dos zeros. Só depois calculamos como referência numérica as ordenadas `Im(zetazero(j))` com `mpmath 1.3.0`, usando 60 dígitos de trabalho; essas ordenadas **não estão certificadas por Arb/intervalos nesta etapa**.

| Índice | Nível do modelo `E_n` | Ordenada do zero conhecido | `E_n−γ_n` |
|---|---:|---:|---:|
| 1 | 17,847837 | 14,134725 | 3,713111 |
| 2 | 23,171661 | 21,022040 | 2,149621 |
| 3 | 27,671198 | 25,010858 | 2,660340 |

O erro quadrático médio das primeiras 25 ordenadas é aproximadamente `1,7563516051`. Não houve ajuste de parâmetros a esses zeros. O modelo captura qualitativamente a densidade média, mas **falha em reproduzir as oscilações aritméticas que determinam os zeros individuais**.

Testamos ainda somas gaussianas finitas (ambas truncadas em 25 termos), com a função `g_T(x)=exp(-(x/T)^2)`:

| T | Soma para níveis do modelo | Soma para zeros conhecidos | Diferença |
|---|---:|---:|---:|
| 20 | 1,027592 | 1,373060 | −0,345468 |
| 40 | 5,305121 | 5,759337 | −0,454216 |
| 80 | 14,610366 | 14,966904 | −0,356538 |

**Aviso:** essas somas parciais **não são** uma fórmula do traço global, e os números não estão certificados por intervalos. A diferença finita registrada não fornece por si só limite da diferença entre traços infinitos.

## 5. Armadilha dos coeficientes de Li

Para quaisquer energias reais positivas `E_n`, defina `θ(E)=2 arctan(1/(2E))`. A soma finita `2 Σ_{n≤J}(1−cos(k θ(E_n)))` é sempre não negativa, **mesmo quando os níveis não são zeros da zeta**. Portanto, repetir essa positividade com `H` não testa o critério completo de Li. A identidade exigiria correspondência exata de todos os zeros com multiplicidades, inclusive zeros fora da reta se existirem, e justificativa da soma completa.

## 6. Ponte proposta com Weil e E7/V24 — pendências

Nenhuma ponte foi provada. Uma eventual ligação exige **primeiro**:

1. Definir em domínio denso um operador com informação aritmética não tautológica (primos e polo arquimediano), que não seja simplesmente uma sequência escolhida para reproduzir a contagem suave.
2. Provar uma fórmula de traço para uma classe explícita de funções de teste, distinguindo lado espectral, termos de primos, arquimediano e polos com suas normalizações.
3. Comparar cada termo com o implementado em E7/V24 (inclusive sinais, projeções e normalizações); estabelecer majorantes de erro e controle de caudas independentes de RH.
4. Somente depois investigar a positividade para todas as funções admissíveis. Positividade de um operador autoadjunto/forma construída por quadrados não é a positividade de Weil.

**Gate ES5 proposto:** reproduzir de maneira independente a fórmula explícita para *funções de teste conhecidas e admissíveis*, antes de acoplar espins ou otimizar parâmetros. Caso não exista identificação rigorosa, o operador ES4 permanece um **controle negativo da suficiência da contagem**.

## 7. Proveniência e referências

- Código e resultados: `es4_spectral_count.py`, `test_es4_spectral_count.py`, `RESULTADOS_ES4.json`.
- Referência da etapa anterior: NUMERIA–ES3, commit `83c58fd41df3748f3ba0d27be21266b0813fda55` no repositório `jbarumby/scicalc`.
- NIST DLMF, §25.10, sobre zeros e contagem: https://dlmf.nist.gov/25.10
- Berry, M. V. e Keating, J. P. (1999), *The Riemann Zeros and Eigenvalue Asymptotics*, SIAM Review **41**, DOI https://doi.org/10.1137/S0036144598347497
- Hall, R. R. (1999), *On the Zeros of the Riemann Zeta-Function*, J. London Math. Soc. **59**, DOI https://doi.org/10.1112/S0024610798006942

**Status final ES4:** um modelo de controle autoadjunto, com resolvente compacto e lei de contagem correta foi construído; ele **não** realiza o espectro real da zeta, não reproduz a fórmula explícita de Weil e não constitui prova da Hipótese de Riemann.
