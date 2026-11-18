# ES8-P239 — positividade da forma de Weil na janela compacta de meia-largura 0,805

## Resultado para os parâmetros P237

Aplicando a redução one-stroke de Zhu à matriz finita P237 e recalculando a cauda em Arb, obtive um certificado de positividade para toda função real, e também complexa, com suporte em `[-0,805, 0,805]`:

```text
Q(f) >= 5 × 10^-19 ||f||₂².
```

O limite é deliberadamente conservador. A margem comum certificada dos dois blocos finitos de 377 modos é `m₀ = 10^-18`; o orçamento de acoplamento infinito fica abaixo de `10^-90`, e a perda diagonal da cauda fica abaixo de `10^-190`. O resultado vale para a janela indicada, não para todos os suportes. Portanto, não prova positividade global da forma de Weil nem a Hipótese de Riemann.

## Identificação com a redução de Zhu

O preprint de Zhu escreve a forma geométrica como um termo de polo de posto um mais uma integral do símbolo de Weil `Ψ_L`. Para um corte `T♯` com

```text
β* = log(T♯/(2π)) − 1/T♯ − A_L > 0,
A_L = Σ_{log n<2L} 2Λ(n)/sqrt(n),
```

O Teorema 1.1 dá `Q(f) >= R(f)` para todo `f` real e par de suporte em `[-L,L]`. A Seção 6, pelo lema de decomposição par/ímpar, estende a redução ao setor ímpar, com sinal oposto no termo de polo. Na base de Legendre, os blocos são `β* I + C + 2ppᵀ` no setor par e `β* I + C − 2ppᵀ` no ímpar; a prova controla a cauda por Gershgorin e o acoplamento por um teste de Schur. A mesma seção prova `Q(f)=Q(Re f)+Q(Im f)` para funções complexas. Portanto, certificar os dois setores reais fornece a mesma margem para funções complexas. [1]

A implementação de montagem P237 usa os mesmos parâmetros e as mesmas peças: a quadratura de `C_nm = (1/π)∫₀^{T♯}(Ψ_L(t)−β*) T̂_n(t)T̂_m(t)dt`, o mesmo `β*`, e os coeficientes de polo obtidos de `2 sqrt(L(n+1/2)) i_n(L/2)`, com o sinal correspondente a cada paridade. A recomputação P239 confere que os `β*` e `A_L` intervalares se sobrepõem aos metadados P237 e fixa por SHA-256 a matriz, o certificado Weyl, a auditoria-sombra e os dois dados `C`. Não se regeneraram neste passo os 1008 painéis de quadratura; o certificado finito P237 continua sendo a fonte desses painéis.

Aqui `L=161/200`, `T♯=504` e `N=377` por paridade. A verificação Arb dá `5 < exp(2L) < 6`, portanto os termos da comb de primos são exatamente `n=2,3,4,5`. Ela também dá `A_L < 5` e `10^-3 < β* < 1`, com `β* ≈ 0,0012165183470747105`. Os graus descartados começam em `n₀=754` no setor par e `n₀=755` no ímpar.

## Cota uniforme do símbolo e das caudas

Para usar um majorante simples no intervalo `[0,504]`, aplico a fórmula de Binet a `ψ(5/4+iy)` com `y=t/2≤252` e a recorrência `ψ(1/4+iy)=ψ(5/4+iy)−1/(1/4+iy)`. Com `a=5/4`, vale `|u²+(a+iy)²|≥a²`; assim, o resto integral de Binet tem módulo no máximo `4/75`, pois `∫₀∞ u/(e^{2πu}−1)du=1/24`. Além disso, `log|5/4+iy|<6`, `|1/(2(5/4+iy))|≤2/5` e o termo da recorrência é menor que `4`. Logo

```text
|Re ψ(1/4+iy)| < 6 + 2/5 + 4/75 + 4.
```

Com `|log π|<2`, `A_L<5` e `β*<1`, resulta

```text
max_[0,T♯] |Ψ_L(t)−β*| < 1384/75 < 20.
```

A cota é verificada por constantes racionais e limites Arb no script P239; não depende de amostragem numérica do símbolo.

A desigualdade de Bessel usada na prova de cauda é `|j_n(x)|≤x^n/(2n+1)!!` para `x≥0`. Assim, para `x=LT♯=10143/25`, defini

```text
b_n = 2 sqrt(L(n+1/2)) x^n/(2n+1)!!.
```

O quadrado da razão `b_(n+2)/b_n` é `x^4/((2n+1)(2n+3)^2(2n+5))`. Nos dois primeiros graus descartados essa razão é menor que `3/40`; o denominador cresce com `n`, então o mesmo limite vale em toda a cauda de cada paridade. Portanto, `Σ_{n≥n₀} b_n ≤ (40/37)b_{n₀}`. Para os modos retidos, usei `|j_m(tL)|≤1`, logo `|T̂_m(t)|≤2 sqrt(L(m+1/2))`.

O polo tem coeficientes `p_n=2 sqrt(L(n+1/2)) i_n(L/2)`. Da série positiva de `i_n`, segue o majorante `|p_n|≤2 sqrt(L(n+1/2)) (L/2)^n/(2n+1)!! · exp((L/2)^2/2)`. A mesma razão geométrica limita sua soma de cauda. Para os coeficientes retidos, Cauchy–Schwarz dá `Σ|p_m|≤sqrt(N) sqrt(2L) cosh(L/2)`, e cada coeficiente individual é no máximo `sqrt(2L) cosh(L/2)`.

Com `c=TM/π <3360`, esses majorantes fornecem o limite de Gershgorin para a cauda `ε_D` e os limites de soma por linha/coluna para o bloco de acoplamento `B`; o teste de Schur dá `||B||≤ε_B`. O replay Arb a 2048 bits encontrou:

- Par: `ε_D < 2,397 × 10^-205` e `ε_B < 2,215 × 10^-98`.
- Ímpar: `ε_D < 1,731 × 10^-206` e `ε_B < 5,959 × 10^-99`.

Os intervalos Arb completos, hashes e as somas por linha/coluna estão no arquivo `P239_ZHU_TAIL_COUPLING.json`.

## Fechamento da margem e reprodução

Como `β*−ε_D > 10^-18` nos dois setores, a cauda diagonal é mais coerciva que a margem finita comum `m₀`. Pelo limite de blocos `λ_min(R)≥min(λ_min(A),λ_min(D))−||B||`, obtemos `R > (10^-18−10^-90)I > 5×10^-19 I`. A desigualdade one-stroke então transfere o mesmo limite inferior para `Q`. A decomposição par/ímpar cobre funções reais gerais; a decomposição em partes real e imaginária estende o limite a funções complexas. [1]

Para reproduzir a auditoria:

```bash
python3 numeria/p237/p239_zhu_tail_coupling.py \
  --data numeria/p237/data \
  --p238 numeria/P238_SCHUR_MARGIN.json \
  --out numeria/P239_ZHU_TAIL_COUPLING.json \
  --bits 2048
python3 -m unittest numeria.test_p239_zhu_tail_coupling -v
```

O relatório é uma extensão computacional da janela P237, não uma prova formal em Lean. As desigualdades clássicas citadas acima e a redução one-stroke são os insumos analíticos; a integração GL60 finita permanece certificada pelos artefatos e verificadores P237/P238. O artigo consultado é um preprint arXiv, e o pacote-fonte arXiv obtido não continha o verificador numérico mencionado no texto. [1]

Elaborado por Manus AI para o projeto NUMERIA.

## References

[1]: https://arxiv.org/html/2608.24827v2 "Xuefeng Zhu, Weil positivity in compact windows: a finite reduction, certified two-sided bounds, and a Landau-Widom decay law, arXiv:2608.24827v2"
