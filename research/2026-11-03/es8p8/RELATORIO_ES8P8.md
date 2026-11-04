# NUMERIA–ES8-P8 — Auditoria integral dos 768.000 nós históricos

**03/11/2026 — Relatório técnico de escopo finito.**

## 1. Pergunta resolvida

A Fase 4 histórico-E7/V16 avaliou 30 integrais do integrando real associado aos produtos de Fourier de modos em `[-L,L]`, com `L = 9/20`, usando células na faixa `0 ≤ r ≤ 160`. A ES8-P7 verificou 450 interseções com um avaliador racional separado. A ES8-P8 amplia isso **à totalidade dos 1.600 × 16 × 30 nós** e registra hashes por lote.

Os 30 pares são `1 <= i <= j <= 10`, de mesma paridade. Para paridades diferentes, a combinação simetrizada é identicamente nula, como documentado na ES8-P2.

## 2. Duas implementações e escopo de independência

- A implementação histórica (`fase4_fixedpoint_sem_arb.py`) usa inteiros em escala `10^40`, 64 termos explícitos da digama e uma cauda Euler–Maclaurin. Seus dois módulos possuem SHA-256 fixado, confrontado com o ZIP Fase 4 original.
- A implementação nova (`IndependentFixed`) utiliza escala `10^35`, 128 termos de uma série real positiva e uma cauda derivada do teste da integral mais correção da derivada. É código distinto, com outro esquema de aproximação do seno e do sinc, outro número de parcelas e outra precisão de arredondamento. As constantes independentes são calculadas a partir do módulo ES8-P7, não do módulo histórico de digama.
- **Independência parcial:** os métodos partilham identidades analíticas usuais e a convenção da transformada de Fourier. Concordância entre eles não constitui prova dessas identidades nem auditoria externa de sua normalização em dimensão infinita.

A relação usada para a digama real é, com `b=r/2` e `y_n=n+1/4`,

\[
\Omega(r)=\Omega(0)+\sum_{n=0}^\infty \frac{b^2}{y_n(y_n^2+b^2)},\quad
\Omega(r)=\Re\psi\left(\frac14+\frac{ir}{2}\right)-\log\pi.
\]

O resto a partir de `u=N+1/4`, com `g(u)=b²/[u(u²+b²)]` e `I(u)=0.5 log(1+b²/u²)`, obedece ao controle adotado

\[
I(u)+\frac{g(u)}{2}\le\sum_{n=N}^\infty g(n+1/4)
\le I(u)+\frac{g(u)}2-\frac{g'(u)}{12}.
\]

A justificação utiliza a representação de Laplace de `g` como transformada de `(1-cos bt) >= 0` e os limites clássicos do kernel de somação. Esse lema continua como dependência matemática explícita (não formalizado em Lean).

Para seno usa-se Taylor de grau 63 com resto majorado por `4^65/65!`, após redução inteira verificável ao intervalo `[-4,4]`. Para `sinc` próximo a zero usa-se Taylor com resto `1/57!` e a estimativa `|sinc'(x)| ≤ 1/2`.

## 3. Teste de compatibilidade

Para cada nó `r = cell/10 + j/150` e para cada par `(i,k)`, com intervalos histórico `[a,b]` em escala `10^-40` e independente `[c,d]` em escala `10^-35`, o teste usa somente multiplicação inteira:

\[
 a\,10^{35}\le d\,10^{40}\quad\text{e}\quad
 c\,10^{40}\le b\,10^{35}.
\]

Essa condição é exatamente equivalente à interseção dos dois intervalos, sem conversão aproximada de `float`. O resultado de cada avaliação, inclusive seus quatro extremos inteiros, entra no hash SHA-256 do lote correspondente.

## 4. Resultados reproduzidos

| Indicador | Valor |
|---|---:|
| Células | 1.600 |
| Nós por célula | 16 |
| Pares de modos de mesma paridade | 30 |
| Posicionamentos nodais | 25.600 |
| Total de comparações de intervalos | **768.000** |
| Intervalos incompatíveis | **0** |
| Lotes com digest SHA-256 | **32** |
| Largura máxima do intervalo histórico | cerca de `4.07108×10^-11` |
| Largura máxima do intervalo novo | cerca de `3.92120×10^-8` |

Todos os intervalos comparados tiveram interseção. A largura independente é deliberadamente maior, pois utiliza outra cota de cauda, que privilegia simplicidade analítica. O registro `RESULTADOS_ES8P8.json` guarda os valores exatos das larguras e todos os hashes, sem depender dos arredondamentos desta tabela.

## 5. Proveniência e referências

O objeto histórico `E7_original` e o ZIP E7 recuperado pelo usuário possuem SHA-256 `4acf99a2...02a7`; o ZIP que contém os fontes Fase 4 tem SHA-256 `b0d7d3...699d8fe`. O resumo `RESUMO_FASE4_1600CELULAS_D40.json` foi extraído do ZIP autenticado. Hashes completos dos conteúdos e arquivos recuperados estão em `PROVENIENCIA_ES8P8.json`.

**Bibliografia selecionada:** Boyadzhiev (2007), DOI 10.1016/j.amc.2006.08.061, trata de representações por contorno de Hankel para funções especiais, incluindo digama. Kanemitsu–Katsurada–Yoshimoto (2000), DOI 10.1007/PL00000117, trata de continuação analítica da zeta de Hurwitz–Lerch. Essas fontes fundamentam possíveis extensões e auditorias, mas **suas técnicas não foram executadas como terceiro avaliador** nesta etapa. A identidade de série utilizada é clássica e também pode ser consultada no DLMF §5.7.

## 6. Limites e próximos gates

1. **Abertura de código, integridade e compatibilidade:** verificadas para todas as 768.000 avaliações de integrando nesta malha.
2. **Reprodução do valor integrado das 30 entradas:** não é automaticamente implicada pela compatibilidade dos intervalos em cada nó; para isso é necessário reproduzir os pesos racionais e a soma, a cota de interpolação, a cauda `r>160` e a projeção, com verificações separadas.
3. **Verificação externa e formalização dos lemas analíticos:** pendentes. A comparação de dois programas não valida automaticamente uma hipótese matemática comum.
4. **Positividade universal de Weil / RH:** não demonstradas. Os resultados dizem respeito a 30 funções-teste finitas da família V17, não a todas as funções-teste admissíveis.

**Próxima etapa ES8-P9 sugerida:** auditar a soma dos pesos racionais das 30 integrais e o transporte de incerteza através da projeção, com revisão formal ou independente dos majorantes residuais.
