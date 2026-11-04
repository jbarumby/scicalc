# NUMERIA–RIEMANN — ES8-P10: lemas analíticos e interpolação otimizada

**Data:** 03/11/2026.  
**Continuação:** ES8-P9 (`research/numeria-es8p9-integrated-2026-11-03`).  
**Estado:** verificação matemática autônoma de um resto da digama, cálculo racional exato de uma cota de interpolação mais forte, e novo transporte dessas cotas às matrizes finitas. **Não é prova da hipótese de Riemann nem de positividade universal de Weil.**

## 1. Resultados quantitativos

Para `L=9/20`, 1600 células de largura `h=1/10`, 16 nós por célula, raio de Cauchy `ρ=3/10`, dados da ES8-P9:

| Grandeza | Resultado |
|---|---|
| Erro global antigo de interpolação | `0.000019792535357116058218` (superior dirigido) |
| **Novo erro global de interpolação** | **`0.000000000003210190573310`** (superior dirigido) |
| Fator de redução estrito | **maior que 6.165.532** |
| Entradas arquimedianas refinadas | 30 de 30, todas contidas nos intervalos antigos |
| Matrizes projetadas | `N=6,8,10`: `M_N - (1/20)I` positiva nas condições descritas |

O fator de redução é **exato**: `16416784384/1601806640625` dividido por `250582204083791/150744971763610839843750`.

Os novos valores do limite da norma da matriz intervalar, `δ=max_i∑_j rad(M_ij)`, são aproximadamente:

| N | Dimensão | δ novo | Menor pivô LDL deslocado |
|---|---:|---:|---:|
| 6 | 4×4 | 0,001920149 | 0,179671870 |
| 8 | 6×6 | 0,004319664 | 0,177272355 |
| 10 | 8×8 | 0,008230725 | 0,173361295 |

O teste LDL é feito sobre `mid(M) − (δ + 1/20)I` em frações racionais, com entradas originais geradas por aritmética intervalar dirigida de 90 casas. Estes números **não** são limites de autovalores do operador infinito nem automaticamente margens na norma funcional, pois a base projetada não é ortonormal.

## 2. Lema demonstrado — resto da série da digama

Pela expansão de frações parciais da digama (§5.7.6 do NIST DLMF), para `r` real,

\[
\Omega(r)-\Omega(0)=\sum_{n=0}^{\infty} f(n+\tfrac14),\quad
f(u)=\frac{b^2}{u(u^2+b^2)},\quad b=r/2.
\]

Para qualquer `u>0`, `b≥0`, a identidade de Laplace é

\[
f(u)=\int_0^\infty e^{-ut}(1-\cos bt)\,dt.
\]

O núcleo é não negativo; a troca entre integral e soma é justificada por Tonelli. Assim,

\[
\sum_{n=0}^\infty f(u+n)-\int_u^\infty f(s)ds-\tfrac12f(u)
=\int_0^\infty e^{-ut}(1-\cos bt)
\left(\frac1{1-e^{-t}}-\frac1t-\frac12\right)dt.
\]

Para `t>0` vale, **sem aproximação assintótica**,

\[
\frac t{12}-\frac{t^3}{720}\le
\frac1{1-e^{-t}}-\frac1t-\frac12
\le\frac t{12}.
\]

**Demonstração dessas desigualdades:** escreva `x=t/2`. As séries de Taylor de `sinh(x)` e `cosh(x)` fornecem

\[
(1+x^2/3)\sinh x-x\cosh x
=\sum_{n=2}^\infty\frac{4n(n-1)}{3(2n+1)!}x^{2n+1}\ge0
\]

 e

\[
x\cosh x-(1+x^2/3-x^4/45)\sinh x
=\sum_{n=3}^\infty\frac{16n(n-1)(n-2)(n+2)}{45(2n+1)!}x^{2n+1}\ge0.
\]

Dividindo por `sinh(x)>0`, obtemos `x/3−x³/45 ≤ coth(x)−1/x ≤ x/3`, equivalente à desigualdade anterior. Integrando contra o núcleo positivo e usando a diferenciação da transformada de Laplace:

\[
I(u)+\tfrac12f(u)-\tfrac1{12}f'(u)+\tfrac1{720}f'''(u)
\le\sum_{n=0}^\infty f(u+n)\le
I(u)+\tfrac12f(u)-\tfrac1{12}f'(u),
\]

onde `I(u)=1/2 log(1+b²/u²)`. As correções implementadas são rigorosamente não negativas:

\[
B_2=-\frac{f'(u)}{12}
=\frac{b^2(3u^2+b^2)}{12u^2(u^2+b^2)^2},
\]

\[
B_4=-\frac{f'''(u)}{720}
=\frac{b^2(10u^6+5b^2u^4+4b^4u^2+b^6)}{120u^4(u^2+b^2)^4}.
\]

**Resultado:** `I+f/2+B₂−B₄ ≤ cauda ≤ I+f/2+B₂`. Isso justifica a direção e os sinais das correções nos avaliadores nodais históricos e independentes. Não demonstra, por si, que cada operação computacional esteja implementada sem erro; os testes de mutação e a auditoria de arredondamentos devem permanecer.

## 3. Lema demonstrado — cota integrada de interpolação

Na célula `[c,c+h]` da Fase 4, defina `x=(r−c)/h` e `x_j=j/15`, `j=0,…,15`. Para a função real analítica `q_ij(r)` e seu polinômio interpolador `P_15`, a fórmula clássica do resto de Lagrange fornece, para um `ξ` real na célula,

\[
q_{ij}(c+hx)-P_{15}(c+hx)
=\frac{q_{ij}^{(16)}(\xi)}{16!}\,h^{16}\prod_{j=0}^{15}(x-x_j).
\]

A cota de Cauchy da ES8-P6 continua válida com um disco de raio `ρ=3/10` ao redor de qualquer ponto da célula: o disco maior centrado no ponto médio tem raio `D=ρ+h/2=7/20` e mantém `Re(1/4±iz/2)≥3/40`. Portanto `|q^(16)|/16! ≤ M_c/ρ^16`, para o mesmo majorante racional `M_c` já auditado. O fator de Fourier é limitado usando `1−2LD=137/200>0`.

A ES8-P6 substituía `|∏(x−x_j)|` pela cota **muito conservadora** `∏max(x_j,1−x_j)`; aqui integramos **exatamente** o valor absoluto:

\[
K_{16}=\int_0^1\left|\prod_{j=0}^{15}(x-j/15)\right|dx
=\frac{250582204083791}{150744971763610839843750}.
\]

O cálculo é puramente racional: constrói-se o polinômio coeficiente por coeficiente, integra-se exatamente em cada um dos 15 intervalos entre raízes, e soma-se o módulo de cada integral, aproveitando a alternância do sinal. O majorante antigo era

\[
K_{16}^{\rm antigo}=\prod_{j=0}^{15}\max(j/15,1-j/15)
=\frac{16416784384}{1601806640625}.
\]

Assim,

\[
|\int_c^{c+h}(q_{ij}-P_{15})dr|
\le h M_c(h/\rho)^{16}K_{16},
\]

 e a soma sobre as 1600 células produz a nova cota global

\[
E_{\rm novo}=E_{\rm antigo}\frac{K_{16}}{K_{16}^{\rm antigo}}
\approx 3,21019057330963\times 10^{-12}.
\]

**Esse é um aperfeiçoamento de uma desigualdade válida para as mesmas funções-teste; não pressupõe RH.** É seguro porque as entradas simetrizadas são reais na reta real, permitindo aplicar o resto real de Lagrange. A referência para o teorema clássico é NIST DLMF §3.3.5.

## 4. Transporte dos certificados

A ES8-P10 não modifica o valor computado em nenhum dos 768 mil nós. Ela reutiliza **apenas** os intervalos nodais auditados da ES8-P9 e substitui `±E_antigo` por `±E_novo`. A precisão de apresentação é de 24 casas, com arredondamento dirigido. As trinta entradas ficam contidas nos intervalos antigos. A cauda `r≥160`, os coeficientes do primo 2 e a base de projeção são os mesmos das ES8-P6/P9.

O arquivo `RESULTADOS_ES8P10.json` registra os erros como frações exatas, a prova algébrica das correções, os parâmetros do disco e os resultados LDL para `N=6,8,10`.

## 5. Fontes e controles

- NIST DLMF [5.7.6](https://dlmf.nist.gov/5.7.E6), frações parciais da digama.
- NIST DLMF [3.3.5](https://dlmf.nist.gov/3.3.E5), resto da interpolação de Lagrange para funções reais.
- `es8p6_independent_audit.py` — cálculo dos majorantes `M_c` e da cota de cauda da Proposição 6.
- `es8p9_projected.py`, `AUDITAR_FASE3_PROJECAO.py` — transporte racional aos casos finitos, sem rederivar fórmula de Weil.
- `RESULTADOS_ES8P9_INTEGRAIS.json` — 30 intervalos nodais vindos da execução integral.
- `test_es8p10.py` — 17 testes, incluindo mutações de dados parciais, extremos forjados, sinais de correção e tratamento de arredondamento negativo.

## 6. Pendências

1. Formalização mecânica em Lean/Coq do teorema da função especial, dos lemas de Cauchy e dos ramos de arredondamento fixo.
2. Demonstração uniforme da positividade de Weil para uma família densa de funções-teste, com limiares independentes de `N`.
3. Revisão externa dos teoremas analíticos herdados, em particular do majorante arquimediano no infinito e da normalização global da fórmula explícita.

**Conclusão:** duas justificativas matemáticas foram fortalecidas, e uma cota de erro de interpolação foi reduzida mais de 6 milhões de vezes, sem perder a positividade das matrizes finitas testadas. **A Hipótese de Riemann permanece em aberto.**
