# NUMERIA–RIEMANN — Etapa ES5: derivadas logarítmicas e traço aritmético

**Data:** 03 de novembro de 2026  
**Status:** duas avaliações finitas intervalares certificadas (mesma grandeza clássica), um teste de traço gaussiano numericamente reproduzido, derivada da zeta completada verificada numericamente. **Sem identificação da forma global de Weil, sem ligação espectral provada com ES4 e sem prova de RH.**

## 1. O problema deixado pela ES4

O operador diagonal ES4 foi construído de modo que sua contagem de níveis reproduza o termo principal de Riemann–von Mangoldt, mas não reproduz os zeros individuais da zeta. Para avançar, a construção precisa de informação **aritmética**, não somente de autoadjunção e de uma densidade espectral imposta. A ES5 testa a fonte clássica mais direta dessa informação, a derivada logarítmica da zeta.

## 2. Uma identidade aritmética incondicional por diferenciação

Para `Re(s)>1`, o produto de Euler converge absolutamente:

\[
\zeta(s)=\prod_p(1-p^{-s})^{-1},\qquad
\log\zeta(s)=-\sum_p\log(1-p^{-s}).
\]

Diferenciar na região de convergência dá

\[
\boxed{-\frac{\zeta'(s)}{\zeta(s)}
=\sum_p\sum_{k\ge1}\frac{\log p}{p^{ks}}
=\sum_{n\ge2}\frac{\Lambda(n)}{n^s}.}
\]

Ela não envolve hipótese alguma sobre zeros na faixa crítica. Os dois métodos de avaliação certificados na ES5 são: (A) somar \(\Lambda(n)n^{-\sigma}\) e limitar sua cauda; (B) somar \(\log n/n^\sigma\) e \(1/n^\sigma\) separadamente, limitar as duas caudas e enquadrar o quociente. Essas duas avaliações usam somas e estratégias de arredondamento diferentes e por isso funcionam como auditoria cruzada, embora dependam da identidade analítica clássica para interpretá-las como a mesma grandeza.

## 3. Cotas de cauda e cálculo com frações exatas

Para σ inteiro ≥2 e N≥3, de \(0\le\Lambda(n)\le\log n\) e da monotonia de \(\log x/x^\sigma\) para x≥N, segue:

\[
0\le \sum_{n>N}\frac{\Lambda(n)}{n^\sigma}
\le \sum_{n>N}\frac{\log n}{n^\sigma}
\le \frac{N^{1-\sigma}}{\sigma-1}
\left(\log N+\frac{1}{\sigma-1}\right),
\]

\[
0\le\sum_{n>N}n^{-\sigma}\le\frac{N^{1-\sigma}}{\sigma-1}.
\]

Todos os logaritmos dos inteiros necessários são confinados em frações racionais com a série

\[
\log x=2\sum_{k=0}^{M-1}\frac{y^{2k+1}}{2k+1}+R_M,\quad
 y=\frac{x-1}{x+1},\quad
0\le R_M\le\frac{2y^{2M+1}}{(2M+1)(1-y^2)},
\]

após decompor `n=2^k u` com `1≤u<2`. A escala fixa é `10^40`; cada divisão e quociente usa arredondamento dirigido por inteiros, sem depender de arredondamento de ponto flutuante. O código preserva os extremos exatos em unidades da escala no JSON. Para `σ=4`, `N=4096`:

| Quantidade | Limite inferior | Limite superior |
|---|---:|---:|
| Série de von Mangoldt | 0,06366976495056302494 | 0,06366976499252638046 |
| Quociente de séries de Dirichlet | 0,06366976491661311567 | 0,06366976495567001957 |
| Interseção | **0,06366976495056302494** | **0,06366976495567001957** |

Cota da cauda da soma com logaritmos: `4.196335552×10^−11`; cota da cauda de \(\zeta(4)\): `4.85063841×10^−12`. Ambas são maiores que os erros efetivos e foram obtidas sem suposições sobre RH. A interseção não é usada para comprovar a identidade conhecida: fornece uma verificação independente e uma estimativa certificada da mesma grandeza quando se assume essa identidade estabelecida.

## 4. Uma identidade gaussiana de traço aritmético

Seja \(\sigma>1\), \(a>0\), \(u\in\mathbb R\). Defina

\[
\mathcal T_{\sigma,a}(u)=\sqrt{\frac a\pi}
\int_{-\infty}^{\infty}e^{-at^2}e^{itu}
\left[-\frac{\zeta'(\sigma+it)}{\zeta(\sigma+it)}\right]dt.
\]

A série para a derivada logarítmica converge absolutamente na semirreta `Re(s)=σ>1`. A identidade da transformada de Fourier da gaussiana e a convergência absoluta permitem trocar integral e soma, obtendo **exatamente**:

\[
\boxed{\mathcal T_{\sigma,a}(u)=
\sum_{n\ge 2}\frac{\Lambda(n)}{n^\sigma}
\exp\left[-\frac{(u-\log n)^2}{4a}\right].}
\]

Essa é uma **identidade de traço aritmético suavizado** (uma filtragem dos primos por escala logarítmica), ainda **não** a fórmula explícita global de Weil que relaciona os zeros, a contribuição arquimediana e os polos. O operador ES4 não foi introduzido na dedução.

Teste com `σ=4`, `a=1/4`, `u=log 2`, integração truncada em `|t|≤18` e soma com `n≤4096`:

| Quantidade | Resultado exploratório (`mpmath`, 45 dps) |
|---|---|
| Integral numérica truncada | 0,0578387187301794498988431637232891869 |
| Soma ponderada de potências de primos | 0,0578387187301794498988431637232891868 |
| Diferença de ponto flutuante | aproximadamente `2,2585×10^-38` |

As **cotas analíticas** para a cauda aritmética (antes do fator gaussiano) e para a integral além de `|t|>18` são, respectivamente, `4,20×10^-11` e `2,09×10^-37`. Para a cota da integral, usamos `|-ζ'/ζ(4+it)|≤Σ_{n≥2} log(n)/n^4 ≤ Σ_{n≥2}n^-3 ≤ 1/2` e a desigualdade gaussiana usual. **A precisão de quadratura e de arredondamento da parte finita não está certificada por intervalos**, portanto a pequena diferença não é uma certificação numérica completa dessa igualdade.

Uma reprodução com `dps=65` confirmou a estabilidade dos dígitos mostrados. Isso também não é certificação intervalar.

## 5. Diferenciais e termo arquimediano da função zeta completada

Para a função inteira completada

\[
\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s),
\]

em pontos que não são zeros ou polos dos fatores, a derivada logarítmica satisfaz

\[
\frac{\xi'(s)}{\xi(s)}=
\frac1s+\frac1{s-1}-\frac12\log\pi
+\frac12\psi(s/2)+\frac{\zeta'(s)}{\zeta(s)}.
\]

Da equação funcional \(\xi(s)=\xi(1-s)\) deduz-se

\[
\frac{\xi'(s)}{\xi(s)}=-\frac{\xi'(1-s)}{\xi(1-s)}.
\]

Com `s=2.25+0.75i`, o resíduo numérico da segunda identidade foi aproximadamente `2.19×10^-47` (45 dps; precisão de ponto flutuante). Este cálculo verifica a normalização do fator \(\Gamma\) no teste, **não demonstra** a identidade de Weil projetada da V24.

## 6. Ponte local com a V24: suporte de primos

No projeto V24, para funções de teste sobre `[-L,L]` com `L=9/20`, os pares de variáveis logarítmicas apresentam suporte em `[-2L,2L]`, portanto `2L=9/10`. Nossa aritmética racional **certifica**:

\[
\log2<\frac9{10}<\log3.
\]

Logo, a soma de potências de primos na janela positiva pode conter apenas `n=2`, e sua amplitude local é proporcional a `Λ(2)/sqrt(2)=log2/sqrt2` na convenção simétrica que aparece na V24. Isso **não** certifica que os fatores de sinal, paridade, projeção e normalização de todas as contribuições `G+P` na V24 foram reproduzidos corretamente. Não se pode usar diretamente as cotas intervalares acima para certificar positividade LDL de E7/V24.

## 7. Critérios de aceitação da próxima fase ES6

1. Escolher uma família explicitamente admissível de funções suaves e compactamente suportadas no domínio de Weil, com transformadas e simetrias conferidas.
2. Reproduzir a **fórmula explícita global**, com parcelas primos, fator arquimediano, polos e zeros e sinais fixados por fonte primária independente, sem pressupor RH.
3. Certificar a quadratura gaussiana/arquimediana com Arb ou prova de erro de quadratura e estimativas de cauda.
4. Confrontar os termos com E7/V24 por transformação de base e projeção documentadas; não reutilizar limites de matriz-modelo como se fossem os da matriz projetada.
5. Testar operadores ES4 como controles adversariais; só declarar identidade espectral se houver um teorema, não um ajuste de níveis.

## 8. Referências e status

- *Encyclopedia of Mathematics*, "Mangoldt function": https://encyclopediaofmath.org/wiki/Mangoldt_function — identidade de Dirichlet clássica.
- DLMF, NIST, *Riemann Zeta Function*: https://dlmf.nist.gov/25 — funções especiais e equação funcional.
- V24, `jbarumby/scicalc/v24/README_V24.md` — envelope de erro condicionado à equivalência das matrizes e fórmulas de E7/V16.
- ES4, commit `5bc7fd34384d37733a0b4c7d121103171aeacbdb` — controle espectral com contagem correta, sem identidade aritmética.

**ES5 aceita apenas como porta aritmética e auditoria independente de normalizações. A Hipótese de Riemann permanece em aberto.**
