# NUMERIA–RIEMANN ES8-P251 · Gate B

Os dois blocos de baixa frequência da janela `L=161/200`, com corte `T=309` e 259 modos de Legendre por paridade, são **positivos definidos** sob a cota analítica de quadratura da P237. A montagem foi refeita a partir do integrando em 618 painéis GL60, com 37.080 nós por setor. Foram certificados 518 pivôs em Arb e novamente em C/MPFR.

| Quantidade certificada | Par | Ímpar |
|---|---:|---:|
| Pivôs positivos do centro | 259/259 | 259/259 |
| Limite inferior estrito para o menor autovalor do bloco contínuo | `3,2373 × 10^-18` | `3,6300 × 10^-15` |
| Limite superior para a soma de raios por linha | `8,952 × 10^-75` | `8,814 × 10^-75` |
| Limite superior da contribuição analítica GL60 à norma matricial | `8,441 × 10^-48` | `8,439 × 10^-48` |

Os valores da tabela foram arredondados na direção conservadora. Os extremos exatos e os limites dos 518 pivôs Arb estão no certificado. [Certificado Arb](data/P251_GATE_B_CERT.json) · [Replay C/MPFR](data/P251_GATE_B_MPFR_AUDIT.json)

## Montagem

O cálculo preserva a janela e o valor fixo de `β = log(504/(2π))−1/504−A_L`. As potências primas são `2,3,4,5` com `A_L=Σ_q 2 log(p)/√q`. Na paridade `s∈{0,1}`, o índice de matriz `k` tem grau `n=2k+s` e amplitude real

```text
u_n(t) = (−1)^k √(2L(2n+1)) j_n(Lt)
```

O fator comum de fase no setor ímpar desaparece no produto bilinear. Cada painel `[p/2,(p+1)/2]` aplica os 60 nós `r_i` de Gauss–Legendre por `t_i=(2p+1+r_i)/4`. Os pesos do produto `u_n(t_i)u_m(t_i)` são `w_i(Ψ_L(t_i)−β)/(4π)`, com `Ψ_L(t)=Re ψ(1/4+it/2)−log π−Σ_q a_q cos(t log q)`. Por setor, a montagem acrescenta a diagonal `βI` e o polo `+2ppᵀ` no par ou `−2ppᵀ` no ímpar, com `p_n=√(2L(2n+1)) i_n(L/2)`. [Normalização e polos P237](../../numeria/p237/p237_ldl_full_arb.py)

Para cada nó, duas avaliações Arb nativas de Bessel esférico no extremo superior inicializam a recorrência descendente. Cada passo usa a identidade exata `j_(n−1)(x)=(2n+1)j_n(x)/x−j_(n+1)(x)`. As sementes e a recorrência usam 640 bits. A geração confere a interseção com `j_0(x)=sin(x)/x`. Os produtos em blocos usam 256 bits, e a soma dos oito lotes usa 384 bits. Os índices e os extremos dos intervalos são conferidos antes da soma. O arquivo agregado preserva 33.670 entradas triangulares por setor, antes de acrescentar `βI`, polos ou erro analítico. [Soma GL60 integral](data/P251_GATE_B_GL60.json.gz)

A montagem usa a fatoração `UᵀWU`. Cada lote mantém em memória até oito painéis. A opção de aplicação em vetor usa `Uᵀ(WUx)` em cada lote, com armazenamento proporcional ao número de nós do lote vezes o número de modos. A fatoração final usa as duas matrizes pequenas materializadas. A estrutura Toeplitz da P249 vale para a base de células uniformes e o termo de primos de espectro completo. O trecho retido em frequência na base de Legendre usa produtos de amplitudes dependentes de ambos os graus. A decomposição assinada P249/P250 entra na verificação do trecho de altas frequências com o mesmo β. [P249](../../NUMERIA_ES8P249/RELATORIO_ES8P249.md) · [P250](../../NUMERIA_ES8P250/README.md)

## Certificação

A cota analítica usa a elipse de Bernstein `E₃` de cada painel, `|Ψ_L−β|<50` no prolongamento complexo da P237 e `|Lt|<249` na faixa calculada. Em cada grau, toma-se a envoltória `a_n=ceil(√(2L(2n+1))) min(2, 2·249^n/(2n+1)!!)`. A exatidão até grau 119 e o resto de Chebyshev dão um erro por entrada inferior a `618·50 a_n a_m/3^120`. A soma máxima por linha limita a norma do erro simétrico. O erro analítico é computado em frações e aplicado uma vez, após somar a quadratura. [Derivação GL60 P237](../../numeria/p237/p237_gl60_analytic_certificate.py)

Seja `A₀` a matriz diádica exata de centros do bloco montado e `δ` a soma dos limites para a norma do erro de montagem e do erro GL60. A congruência LDL de `A₀` tem diagonal estritamente positiva. Uma inclusão dirigida da inversa limita `ν≥||A₀⁻¹||_F²`, logo o menor autovalor de qualquer matriz da caixa contínua é pelo menos `1/√ν−δ`. O replay independente recompila o kernel C/MPFR P242 autenticado por SHA-256, reconstrói os 518 pivôs e as inversas a 1.024 bits e verifica em frações a desigualdade de Weyl para uma margem racional. A etapa MPFR recebe os centros do agregado Arb. [Código C/MPFR P242](../../numeria/p242/p242_mpfr_interval_ldl.c) · [Replay](data/P251_GATE_B_MPFR_AUDIT.json)

O Gate A demonstra `H₀₀>2,8884×10^-10` para a contribuição retirada entre 309 e 504, suficiente para inviabilizar a transferência escalar da margem P238 de `10^-18`. O Gate B aplica Weyl ao erro de reconstrução da própria matriz em 309, cujo orçamento é inferior a `8,441×10^-48` em cada setor. [Gate A](../README.md)

## Auditoria independente

- Os novos replays dos painéis `[0,25)` e `[25,50)` contêm os **134.680 intervalos Arb arquivados** dos subtriângulos de 259 modos. Os quatro arquivos-fonte P236/P237 são vinculados por SHA-256. Em `[500,525)`, **20 entradas selecionadas** têm a bola histórica contida tanto no replay novo quanto no replay P245D em ponto fixo. Essas comparações são entre somas discretas GL60. [Auditoria P236/P237/P245](data/P251_GATE_B_HISTORICAL_AUDIT.json)
- Os testes dos painéis 0 e 617 usam os brackets racionais P245B autenticados por SHA-256, pesos pela identidade da derivada de Legendre e amplitudes escalares em forma fechada. Os testes adversariais cobrem lacunas, duplicatas, β incorreto, ordem do triângulo, pivôs e inversas inválidos e cotas de perturbação inconclusivas.
- Reexecutei a ação matriz-livre de dois vetores diádicos, um por paridade, em todos os 618 painéis. As **518 coordenadas** coincidem por interseção de intervalos com a ação do agregado. O módulo da diferença enclausurado ficou abaixo de `1,021 × 10^-74` por coordenada nos vetores executados. [Replay matriz-livre](data/P251_GATE_B_MATRIX_FREE.json)
- Foram aprovados **41 testes** em Python normal e com `-O`, sendo 13 do Gate B e 28 da P249, P250 e Gate A. O replay do Gate A também foi executado no commit-base deste Gate B.

## Complemento modal e alcance

A prova racional da P249 com os limites P243/P239 dá `||B||<10^-90` e `||D−βI||<10^-190` para os graus descartados a partir de 518 e 519. O limite inferior certificado para a cabeça é `m>3,2373×10^-18`. A desigualdade de Young permite subtrair `m/2` da cabeça e `2||B||²/m` da cauda. O teste exato dá uma margem **superior a `1,6186×10^-18` para o operador-modelo reduzido** na janela fixa.

A extensão usa a positividade do símbolo de P250 para `t≥309` e a cota analítica P239 para o complemento modal. Os resultados da quadratura usam as identidades clássicas de Bessel e digama, a estimativa de Bernstein e a exatidão de Gauss–Legendre. Uma afirmação para a forma de Weil depende da identificação Weil–Zhu na classe admissível. A extensão cofinal e a Hipótese de Riemann permanecem em aberto. [P239](../../numeria/README_ES8P239_ZHU_EXTENSION.md) · [Certificado Gate B](data/P251_GATE_B_CERT.json)

## Reprodução

Executar na raiz do repositório. O replay C/MPFR requer um compilador C, GMP e a biblioteca `libmpfr.so.6` (em Linux). As decisões de sinal usam Arb, MPFR e frações exatas.

```bash
python -m pip install -r NUMERIA_ES8P249/requirements.txt
python -m unittest NUMERIA_ES8P251.gate_b.test_finite_head -v
python -O -m unittest NUMERIA_ES8P251.gate_b.test_finite_head -v
sha256sum -c NUMERIA_ES8P251/gate_b/SHA256SUMS
python -m NUMERIA_ES8P251.gate_b.finite_head \
  --matrix NUMERIA_ES8P251/gate_b/data/P251_GATE_B_GL60.json.gz \
  --out /tmp/P251_GATE_B_CERT.json \
  --action-audit /tmp/P251_GATE_B_MATRIX_FREE.json
python -m NUMERIA_ES8P251.gate_b.mpfr_audit \
  --matrix NUMERIA_ES8P251/gate_b/data/P251_GATE_B_GL60.json.gz \
  --out /tmp/P251_GATE_B_MPFR_AUDIT.json
python -m NUMERIA_ES8P251.gate_b.historical_audit \
  --out /tmp/P251_GATE_B_HISTORICAL_AUDIT.json
```

A regeneração direta do agregado parte de todos os painéis e usa as mesmas divisões da execução arquivada.

```bash
work=$(mktemp -d)
for start in 0 80 160 240 320 400 480 560; do
  stop=$((start+80))
  if [ "$stop" -gt 618 ]; then stop=618; fi
  python -m NUMERIA_ES8P251.gate_b.finite_head_quadrature \
    --start "$start" --stop "$stop" --bits 256 --bessel-bits 640 --threads 1 \
    --out "$work/$start-$stop.json.gz"
done
python -m NUMERIA_ES8P251.gate_b.finite_head --chunks "$work/"*.json.gz \
  --matrix "$work/P251_GATE_B_GL60.json.gz" --out "$work/P251_GATE_B_CERT.json"
cmp NUMERIA_ES8P251/gate_b/data/P251_GATE_B_GL60.json.gz "$work/P251_GATE_B_GL60.json.gz"
cmp NUMERIA_ES8P251/gate_b/data/P251_GATE_B_CERT.json "$work/P251_GATE_B_CERT.json"
```

Os comandos `cmp` pressupõem a mesma implementação e as versões fixadas. As decisões matemáticas do replay se baseiam nos intervalos e nas cotas verificadas.
