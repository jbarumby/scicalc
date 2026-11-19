# P246 — escala do corte cofinal e limite da estratégia de força bruta

## Resultado

A sequência cofinal proposta em P244 fecha o sinal da cauda de alta frequência para cada meia-largura finita `L`, mas o custo cresce muito rapidamente. Para

\[
A_L=2\sum_{n<e^{2L}}\frac{\Lambda(n)}{\sqrt n},\qquad
T_L=\max\{504,4\pi e^{A_L}\},
\]

o Teorema dos Números Primos implica

\[
A_L\sim4e^L,\qquad \log T_L=4e^L+o(e^L).
\]

Assim, `T_L` — e o corte de modos `N` suficiente para a estimativa de cauda de Bessel, da ordem de `L T_L` — cresce essencialmente como \(\exp(4e^L)\). Esta é uma obstrução concreta à **exaustão por força bruta** dentro da escolha de corte de P244; não é uma obstrução à positividade matemática nem prova de que outra redução não possa ser mais eficiente.

## Derivação

Defina a função de Chebyshev

\[
\psi(x)=\sum_{n\le x}\Lambda(n).
\]

O PNT equivale a `ψ(x)=x+o(x)` [1]. Aplicando somação por partes a `f(t)=t^{-1/2}` [2],

\[
\sum_{n\le X}\frac{\Lambda(n)}{\sqrt n}
=\frac{\psi(X)}{\sqrt X}
+\frac12\int_1^X\frac{\psi(t)}{t^{3/2}}\,dt.
\]

O primeiro termo é assintótico a `√X`; no segundo, substituir `ψ(t)` por `t+o(t)` dá novamente `√X`. Portanto a soma sem o fator externo 2 é assintótica a `2√X`. Tomando `X=e^{2L}` e recolocando esse fator,

\[
A_L=2\sum_{n<e^{2L}}\frac{\Lambda(n)}{\sqrt n}
\sim 4\sqrt{e^{2L}}=4e^L.
\]

A convenção estrita `n<e^{2L}` não altera a assíntota. Se o extremo atravessar uma potência prima, a diferença é no máximo um termo `2 log(p)/sqrt(p^k)`, que é desprezível frente a `e^L`.

Para `L` suficientemente grande, o máximo com 504 deixa de atuar. Da definição de P244,

\[
\log T_L=\log(4\pi)+A_L=4e^L+o(e^L),
\quad
T_L=\exp\!\bigl(4e^L+o(e^L)\bigr).
\]

O controle de razão do P239 é feito com `x=L T_L` e

\[
\left(\frac{b_{n+2}}{b_n}\right)^2
=\frac{x^4}{(2n+1)(2n+3)^2(2n+5)}.
\]

O programa P246 escolhe, conservadoramente, `N=ceil(L T)` por setor, com `T` racional e maior ou igual ao corte P244. Verifica em frações exatas que a razão no primeiro grau descartado é menor que `(3/40)^2` nos dois setores. Logo essa escolha suficiente tem `N` da mesma escala `L T_L`; ela não testa a positividade do bloco de baixa frequência.

## Recomputação reproduzível

O script calcula a soma finita `A_L` com bolas Arb, isola o inteiro `floor(exp(2L))`, escolhe um corte racional superior e verifica exatamente as razões Bessel. Os resultados incluem:

- `L=1`: `A_L≈5,8525`, `log10(T)≈3,6409` e `N=4.375` modos por setor.
- `L=1,5`: `A_L≈13,0157`, `log10(T)≈6,7518` e `N=8.470.867`.
- `L=2`: `A_L≈24,3833`, `log10(T)≈11,6887` e `N=976.717.466.762`.
- `L=2,5`: `A_L≈42,5516`, `log10(T)≈19,5791` e `N≈9,49×10^{19}`.
- `L=3`: `A_L≈75,1905`, `log10(T)≈33,7540` e `N≈1,70×10^{34}`.
- `L=4`: `A_L≈213,5475`, `log10(T)≈93,8417` e `N≈2,78×10^{94}`.

Os dígitos exibidos são arredondamentos de bolas Arb armazenadas no JSON, não entradas de ponto flutuante usadas na decisão do corte. Para `L=161/200`, o algoritmo genérico P244 resulta em `N=809`; isto é mais conservador que o certificado P239 específico (`T♯=504`, `N=377`) e não altera nem substitui esse certificado.

Reprodução a partir da raiz do repositório:

```bash
python3 numeria/p246_cofinal_scaling.py --bits 512 \
  --out numeria/P246_COFINAL_SCALING.json
python3 -m unittest discover -s numeria -p 'test_p246_cofinal_scaling.py' -v
python3 -O -m unittest discover -s numeria -p 'test_p246_cofinal_scaling.py' -v
```

## Consequência para a próxima etapa

P239 certifica uma janela fixa, `L=0,805`. P244 formula corretamente um critério cofinal de Schur, mas não verifica as suas hipóteses para uma sequência `L_j→∞`. P246 mostra que apenas ampliar as matrizes com o corte global `4πe^{A_L}` torna a dimensão rapidamente impraticável. Para avançar rumo à positividade global será necessário melhorar a redução ou provar, para cada janela da exaustão, a positividade do bloco de baixa frequência e o controle do acoplamento com margens adequadas. Uma possível linha de investigação é substituir a majorização absoluta da soma de primos por um controle rigoroso mais fino do termo oscilatório; tal cancelamento uniforme ainda não foi provado aqui e não deve ser presumido.

O resultado P246 é um diagnóstico assintótico e computacional sobre a escalabilidade do método de corte atual, não um novo critério de RH. Não estabelece positividade global da forma de Weil, nem a Hipótese de Riemann.

## Referências

[1]: https://dlmf.nist.gov/25.16.i "NIST Digital Library of Mathematical Functions, §25.16(i), Distribution of Primes"
[2]: https://www.math.lmu.de/~forster/v/ann/annth_11.pdf "Otto Forster, Analytic Number Theory, Chapter 11, Chebyshev Functions and Abel Summation"

Elaborado por Manus AI para o projeto NUMERIA.
