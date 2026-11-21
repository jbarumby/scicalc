# NUMERIA–RIEMANN · ES8-P259 · Suzuki escalar e cauda por faixas

A P259 certifica positividade da forma de Weil numa janela racional **acima de `log(7)/2`** e verifica a função escalar de Suzuki em uma faixa contínua finita. Os dois resultados usam a fórmula explícita e aritmética dirigida. Esta branch parte da P258, no commit [9a46e7b](https://github.com/jbarumby/scicalc/commit/9a46e7bf6164edce5bc3b3cd91451a3420b0bbe6).

## Resultado para a forma de Weil

Para `L=0,972955074527656652552676371721590` e toda função complexa suave de suporte compacto em `[-L,L]`, Arb certifica `Q_W(f)≥2,42×10⁻²⁸ ‖f‖₂²`. O replay em C/MPFR certifica independentemente os 754 pivôs e uma margem superior a `2,41×10⁻²⁸`. A largura supera `log(7)/2` em aproximadamente `1,35×10⁻³⁴`. A identificação com a forma de Weil usa a fórmula de Guinand–Weil da P251 e incorpora o deslocamento de 7 pela estimativa de borda abaixo. [Certificado Arb](data/BANDED_WEIL_590.json) · [Replay MPFR](data/P259_MPFR.json) · [P251 Gate C](../NUMERIA_ES8P251/gate_c/README.md)

| Início do piso 1/2 | Setor par, bloco finito | Setor ímpar, bloco finito | Extensão ao operador |
| ---: | --- | --- | --- |
| 850 | Margem superior a `4,66×10⁻²⁸` | Rayleigh menor que `−4,77×10⁻²⁷` | O auxiliar tem direção negativa |
| 590 | Margem superior a `4,85×10⁻²⁸` | Margem superior a `1,93×10⁻²⁵` | Schur fecha nos dois setores |

A antiga testemunha ímpar da P258 fica positiva nos dois modelos. No corte 850 surge outra testemunha, com 57 coordenadas diádicas. O Gate B publicado no PR #66 foi reproduzido byte a byte e confirma essa obstrução específica. Antecipar o piso maior para 590 acrescenta `(G₈₅₀−G₅₉₀)/4` à forma, uma correção semidefinida positiva pela projeção sobre a faixa de frequências. [Gate B publicado](https://github.com/jbarumby/scicalc/pull/66) · [Certificado do corte 850](data/BANDED_WEIL_850.json)

### Símbolo e transferência modal

O auxiliar de quatro termos é `1/2 I`, mais os polos e o multiplicador retido até `T`. Seu multiplicador residual é `Ψ₄−1/2` até 504 e `−1/4` entre 504 e `T`. Isso equivale a piso absoluto `1/4` entre 504 e `T` e `1/2` a partir de `T`. A identidade coincide com a atualização por `(I−G_T)/4` do auxiliar P258. A integral ponderada até 504 vem do agregado P253 autenticado por SHA-256 e mantém a cota analítica GL60. As concentrações `G₅₀₄` e `G_T` usam as identidades de extremos de Bessel com Arb, com verificação de uma amostra GL60 em 850.

O piso `1/4` a partir de 504 foi regenerado em Arb. Para `t≥850`, Binet e a soma das amplitudes de primos dão `Ψ₄(t)≥β_antigo+log(850/504)+1/504−1/850>1/2`, pois `log(850/504)>346/677>1/2`. Os 286 subintervalos adicionais cobrem continuamente `[590,850]` e certificam o mesmo piso. [Auditoria da cauda no JSON](data/BANDED_WEIL_590.json)

O antigo teste que reservava uma parte da margem do complemento permitia usar apenas uma fração da atualização positiva no bloco retido. Isso deixava o ímpar inconclusivo. O argumento direto aplica a estimativa de cauda à representação do auxiliar com diagonal `1/2 I` e todo o multiplicador residual apoiado em `[-590,590]`. A cota P239 para a parte até 504 admite qualquer piso entre 0 e 1. Ela dá módulo inferior a 20, também válido na faixa residual constante.

Se `n` é o primeiro grau de Legendre descartado, a aproximação da exponencial `exp(itx)` por polinômios de grau menor que `n` na elipse de Bernstein de parâmetro 2 dá

    c = 6 exp(3LT/4) 2^(-n),    p = 12·2^(-n).

De fato, a cauda de Chebyshev em norma uniforme é limitada por `4M·2^(-n)` se a função tem módulo até `M` na elipse. O fator de conversão para L² é `√(2L)<3/2`. Para os polos, o majorante na elipse é `exp(5L/8)<2`. A projeção sobre os graus descartados anula o polinômio aproximante. As normas integrais dos vetores de frequência retidos e de polos retidos são menores que 2. Portanto

    ||B|| ≤ (40T/3)c+4p,     ||D−(1/2)I|| ≤ (20T/3)c²+2p².

O auditor encerra a exponencial em Arb e arredonda as cotas racionais para cima. Para uma cota inferior `a` do bloco retido, `d≤λ_min(D)` e `b²≥||B||²`, verifica `a,d>m` e `(a−m)(d−m)>b²`. Isso dá margem estrita `m` para o auxiliar inteiro. O replay C/MPFR usa centros diádicos, raios da montagem e a mesma cota analítica GL60. [Prova do Schur P254](../NUMERIA_ES8P254/README.md) · [Orçamento real P239](../numeria/README_ES8P239_ZHU_EXTENSION.md)

### Inclusão do primo 7

Com `Ω=3360`, Binet dá `Ψ₄(t)−1/2>log(20/3)−1/2>13/10+1/12` em `|t|≥Ω`. A prova de faixas separadas da P258, com consumo alto `13/10`, deixa reserva `1/12` e dá perda menor que `1839δ`, onde `δ=2L−log 7`. O auditor encerra `log 7` por séries racionais e verifica a margem no **próprio extremo racional**. Temos `2L<2<log 8`, de modo que apenas a potência prima 7 entra nessa janela. A extensão cofinal permanece por demonstrar. [Faixas de borda P258](../NUMERIA_ES8P258/README.md#limiar-do-primo-7)

## Resultado escalar de Suzuki

Supondo apenas a fórmula explícita, os **78.742 subintervalos** e a estimativa inicial cobrem todo o intervalo `[0,log(10⁶)]`. O auditor encontra **`Ψ(t)>0` para `t>0` nessa faixa**, com `Ψ(0)=0`. Também verificou a desigualdade aritmética abaixo para os **78.734 prefixos** de potências primas até um milhão. O resultado numérico da faixa escalar tem escopo distinto da margem de Weil. [Certificado escalar](data/SUZUKI_FINITE.json) · [Suzuki 2023](https://doi.org/10.1112/jlms.12785)

Com `λₙ=(4n+1)/2` e `c=(log π−ψ(1/4))/2`, Suzuki escreve sua função como `Ψ(t)=A(t)−P(t)`, onde

    A(t) = 4(e^(t/2)+e^(−t/2)−2)−ct + Σ (1−e^(−λₙt))/λₙ²,
    P(t) = Σ_{q≤e^t} (Λ(q)/√q)(t−log q) = ∫₀^(t/2) S_L(0) dL.

A integral decorre de cada potência prima entrar em `S_L(0)` com peso `2Λ(q)/√q` quando `L≥(log q)/2`. Em `0<t≤1/64`, os oito primeiros termos de Lerch e a desigualdade `1−e^(−u)≥u/(1+u)` dão `Ψ(t)≥m₀t` com `m₀>0`. Para `t≥1/64`, entre potências primas, a segunda derivada é `e^(t/2)+e^(−t/2)−e^(−t/2)/(1−e^(−2t))`. Ela é maior que `−33` e torna-se positiva a partir de `log 2`. Os tangentes com esse erro de curvatura cobrem continuamente cada intervalo, incluindo os contatos por continuidade. O resto de Lerch é majorado por uma série geométrica com razão `e^(−2t)`. [P256](../NUMERIA_ES8P256/README.md)

### Alvo aritmético para uma prova independente

Para uma potência prima `q`, sejam `A_q=Σ_{n≤q}Λ(n)/√n` e `B_q=Σ_{n≤q}Λ(n) log(n)/√n`. A função `P(t)` é o supremo das linhas `A_q t−B_q`, incluindo o prefixo vazio. O máximo de `A_qt−4e^(t/2)` é `2A_q(log(A_q/2)−1)`. Assim, uma forma equivalente da desigualdade eventual de Suzuki é

    B_q ≥ 2A_q(log(A_q/2)−1) para todas as potências primas q suficientemente grandes.

A condição implica `P(t)≤4e^(t/2)` eventualmente, pois as finitíssimas linhas de prefixos pequenos são dominadas para t grande. Na direção inversa, `A_q→∞` pela divergência de `Σ_p 1/p`, e o maximizador `2log(A_q/2)` acaba no intervalo eventual. O certificado confirma os prefixos até um milhão. Falta uma estimativa analítica que controle todos os prefixos grandes. [Suzuki 2025, Teorema 1](https://doi.org/10.1007/s11139-025-01238-9)

## O que foi conferido nas publicações

- Os três certificados da execução anterior foram reproduzidos byte a byte no commit `ab6f0f2`. A versão da banda 590 arquivada nesta branch acrescenta a transferência de Schur. Os 21 testes dos PRs #65–#66 e os 14 testes focados da P258 e desta branch passaram em Python normal e com `-O`. [Execução publicada](https://github.com/jbarumby/scicalc/actions/runs/44861182405)
- Os [PRs #67–#69](https://github.com/jbarumby/scicalc/pull/69) foram conferidos no commit `503e361`. Os 68 testes de P260 a P262 passaram tanto em Python normal quanto com `-O`. Foram reproduzidos os momentos theta e o desacordo entre a leitura literal de uma fórmula de KPS e a recorrência testada. Aquele resultado tem escopo finito.
- A atualização do MIT de 21 de novembro confirmou `con-leche` e `nanoda` na exportação principal de Alpöge. As novas desigualdades de Robin e Lagarias são consequências formais do mesmo teorema, por equivalências verificadas separadamente. O Comparator completo da prova principal aparecia em andamento. A P259 não usa esses teoremas em suas verificações. [Auditoria do MIT](https://math.mit.edu/~drew/grh/) · [Redução aritmética em Lean](https://math.mit.edu/~drew/grh/ArithSolution.lean)

## Reprodução

```bash
python -m pip install -r NUMERIA_ES8P249/requirements.txt
sha256sum -c NUMERIA_ES8P259/SHA256SUMS
python -m unittest NUMERIA_ES8P258.test_p258 NUMERIA_ES8P259.test_p259 -v
python -O -m unittest NUMERIA_ES8P258.test_p258 NUMERIA_ES8P259.test_p259 -v
python -m NUMERIA_ES8P259.suzuki --limit 1000000 --out SUZUKI_FINITE.json
python -m NUMERIA_ES8P259.banded_weil --high-cutoff 850 --out BANDED_WEIL_850.json
python -m NUMERIA_ES8P259.banded_weil --high-cutoff 590 --out BANDED_WEIL_590.json
python -m NUMERIA_ES8P259.mpfr_replay --out P259_MPFR.json
```

O replay C requer compilador C, GMP e `libmpfr.so.6` no Linux.
