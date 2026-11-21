# ES8-P262 — Confronto direto de fontes: KPS (2022 e 2024)

## Fonte primária, páginas efetivamente conferidas

- Takis Konstantopoulos, Pierre Patie e Rohan Sarkar, *A new class of solutions to the van Dantzig problem, the Lee–Yang property, and the Riemann hypothesis*, **Annales de l'Institut Fourier 74 (2024)**, p. 395 (Eq. 4.8), p. 397 (Eqs. 4.15–4.18): https://www.numdam.org/item/10.5802/aif.3600.pdf
- Mesmos autores, **arXiv:2211.16680v1 (30 nov. 2022)**, p. 16–17, Eqs. 4.17–4.18: https://arxiv.org/pdf/2211.16680

A expressão impressa de (4.18) com sinal de menos está presente **nas duas versões**. O presente relatório identifica uma obrigação de reconciliação, não uma errata publicada nem uma refutação da pesquisa original.

## O que foi certificado com duas cadeias independentes

### Cadeia A: integrar theta

Da Eq. (4.15) dos autores, para inteiro par k≥0,

\[
F(k)=\int_1^\infty(\log x)^k x^{-3/4}
  \sum_{m\ge1}\exp(-\pi m^2 x)\,dx.
\]

O auditor P262 usa acb.integral em u=log x (u∈[0,5], m=1,2,3), com erros garantidos pelo Arb. As caudas m≥4 e u≥5 são majoradas por frações exatas conforme PROVA_P262.md. Os extremos dos intervalos dos três primeiros momentos e os seus SHA-256 são parte do artefato GitHub Actions.

| Momentos | Centro aproximado; sinal rigoroso |
| --- | ---: |
| F(0) | 0.01151688724674356 >0 |
| F(2) | 0.0009892828534669306 >0 |
| F(4) | 0.0003200081096422168 >0 |

Pela definição de G na Eq. (4.18), G(2)≈744.0657627477>0, G(4)≈1186.1093396853>0, ambos com extremos diádicos estritamente positivos. A **expressão impressa** da Eq. (4.18) fornece φ(2)≈−0.0392072708778<0.

### Cadeia B: somente frações exatas, sem nova quadratura

Pela Eq. (4.17) *como impressa*, defina

\[
A_n=\frac{n!}{(2n)!}\,
\frac{32\binom{2n}{2}F(2n-2)-F(2n)}{2^{2n-1}}
=\frac1{W_\phi(n+1)},\quad n\ge1.
\]

A Eq. (4.8), página 395, define W pela recorrência multiplicativa: Wφ(n+2)=φ(n+1)Wφ(n+1). Portanto

\[
\boxed{\phi(n+1)=A_n/A_{n+1}}
\]

se os denominadores forem não nulos. Usando **os mesmos extremos binários exatos** de F(0), F(2), F(4), o replay exclusivamente racional prova:

| Expressão | Intervalo verificado por frações exatas |
| --- | ---: |
| A₁ da Eq. (4.17) | 0.09188777726058054… a 0.09188777726058296… |
| A₂ da Eq. (4.17) | 0.0019752322891250878… a 0.0019752322891250884… |
| φ(2) imposto pela recorrência de W | **46.51998540449205… a 46.51998540449330…** |
| φ(2) da Eq. (4.18) literalmente impressa | **−0.03920727087780425… a −0.03920727087780319…** |

Os dois intervalos para φ(2) são rigorosamente **disjuntos**. Esse é um resultado de consistência *entre expressões impressas*, sem usar RH. O algoritmo rejeita arquivos com hash inválido quando se solicita o pino de SHA-256.

## Terceiro teste de sinais: Xi como função característica

A Eq. (2.7) do mesmo artigo fornece a transformada de Fourier
\[
\Xi(t)=\int_{\mathbb R}e^{itx}\Phi(x)\,dx
\]
com densidade \(\Phi(x)>0\) e decaimento suficientemente rápido para ter segundo momento finito não nulo. Portanto
\[
\Xi''(0)=-\int_{\mathbb R}x^2\Phi(x)\,dx<0.
\]
Se a Eq. (4.15) \(\Theta(z)=\Xi(\sqrt z)\) for interpretada literalmente, o coeficiente de z deve ser \(\Xi''(0)/2<0\). Entretanto a expressão literal para o primeiro coeficiente de (4.16)/(4.17), quando aplicada aos momentos certificados, tem sinal positivo.

Este é um **indício adicional de problema de sinal, convenção ou normalização entre as fórmulas impressas**. Ele não permite concluir onde está a correção sem rederivar a transformação theta a partir da representação original de Xi, e não é prova nem refutação de RH.

## O que NÃO resulta desta auditoria

1. Não provamos que a fórmula **pretendida** pelos autores seja falsa; constatamos incompatibilidade na interpretação literal das expressões publicadas. Poderia haver convenção elidida ou erro de sinal/índice. Verificar errata, corrigenda e rederivação independente é necessário.
2. Não excluímos a existência de representação de Xi na classe de Bernstein, muito menos sua pertença a Laguerre–Pólya.
3. Não concluímos nada sobre zeros de Xi fora da reta nem sobre a positividade global de Weil.
4. Não inferimos que uma análise numérica finita de momentos resolva quantificadores de todos os índices.

## Próxima prioridade científica

Derivar diretamente a série de Taylor de Xi pela representação positiva da densidade theta da Eq. (2.7), com normalização Xi/Xi(0) explícita, conferir os fatores 2, n!, binomial e o sinal (-1)^n, e então reconstruir os candidatos φ(n) **corrigidos** pela recorrência de W; só depois estudar não negatividade, monotonicidade, concavidade e derivadas superiores necessárias para uma função de Bernstein. 
