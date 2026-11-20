# NUMERIA–RIEMANN · ES8-P252

A P252 amplia a janela certificada de meia-largura `0,805` para **`0,95`** e implementa um circuito reutilizável para a matriz Toeplitz do termo de primos em células. A certificação clássica da janela e a simulação quântica são resultados com escopos próprios.

## Janela (L=19/20)

Foram recalculados os **618 painéis GL60** até o corte 309, com 37.080 nós e 259 graus de Legendre por paridade. A janela continua abaixo de `log(7)/2`. Como `5<exp(2L)<7` e 6 não é potência prima, as quatro potências `2,3,4,5`, o símbolo `Ψ` e o β fixado em 504 permanecem os mesmos. Essa verificação permite reaproveitar o limite `|Ψ−β|<20` na banda retida e o certificado P250 de sinal em frequências `t≥309`. A quadratura, os argumentos de Bessel, a normalização de Legendre e os polos foram refeitos para o novo `L`. [P250](../NUMERIA_ES8P250/README.md) · [Método de montagem](../NUMERIA_ES8P251/gate_b/README.md)

| Setor | Pivôs LDL positivos | Margem estrita conservadora em Arb |
|---|---:|---:|
| Par | 259/259 | `6,82 × 10^-27` |
| Ímpar | 259/259 | `1,237 × 10^-23` |

Arb a 1.024 bits certificou a congruência LDL do centro de cada caixa, as cotas da inversa e a transferência de Weyl. O replay C/MPFR independente repetiu os 518 pivôs e o limite de Weyl sobre os mesmos dados de quadratura. Na elipse de Bernstein das meia-unidades, `|Lt|≤294` e `exp(|Im(Lt)|)<2`. A cota de erro GL60 por soma de linha ficou abaixo de `1,391×10^-47` em ambos os setores. Essa cota é aplicada uma vez, depois da soma dos painéis. [Certificado Arb](data/P252_WINDOW_CERT.json) · [Replay MPFR](data/P252_MPFR_AUDIT.json)

O teste racional do complemento dá `‖B‖<10^-54` e `‖D−βI‖<10^-117` nas duas paridades. O complemento de Schur fornece margem superior a `3,41×10^-27` para o operador-modelo. A identificação pela fórmula clássica de Guinand–Weil, feita no Gate C, dá

```text
Q_W(f) ≥ 3,41×10^-27 ‖f‖₂²
para f ∈ C_c^∞(ℝ;ℂ), supp(f) ⊆ [−19/20,19/20].
```

Por inclusão de suportes, a mesma cota vale em janelas menores. Os insumos analíticos são a fórmula explícita, as identidades de Bessel e Binet e o orçamento de Bernstein para GL60. A extensão cofinal e a Hipótese de Riemann permanecem em aberto. O próximo limiar aritmético nesta direção é `L=log(7)/2`, onde o suporte das correlações alcança o novo deslocamento `log(7)`. [Ponte de Weil no Gate C](../NUMERIA_ES8P251/gate_c/README.md)

## Circuito quântico para o termo de primos

Para `N` células uniformes, sejam `k_d≥0` os lags da matriz Toeplitz de translações simétricas da P249. Use um registrador de dados de tamanho `M`, potência de dois com `M≥2N−1`, e a translação cíclica `U_d|j⟩=|j+d mod M⟩`. A combinação linear de unitários é

```text
C = k₀ I + Σ_{d>0} k_d (U_d + U_−d)
α = k₀ + 2 Σ_{d>0} k_d
```

O circuito prepara o seletor com amplitudes `√(k_d/α)` para os deslocamentos assinados. Depois aplica o deslocamento condicionado e desfaz a preparação. O bloco com seletor zero é `C/α`. Projetar também o registrador de dados nas primeiras `N` células dá exatamente `K_N/α`, com `K_N[i,j]=k_|i−j|` no modelo ideal. O preenchimento garante que a diferença física dos índices seja representada sem colisão modular. A probabilidade da projeção conjunta para uma entrada normalizada `x` é `‖K_N x‖²/α²`. [Operador de translações da P249](../NUMERIA_ES8P249/RELATORIO_ES8P249.md) · [QASM](data/P252_PRIME_LCU.qasm)

Foram implementadas duas versões da seleção. A primeira usa soma reversível de constantes. A segunda aplica QFT, fases condicionadas por cada bit de frequência e a QFT inversa. Em frequência, `U_d` tem o fator `exp(2πi d k/M)`. As fases são multiplexadas pelos bits do seletor. As duas versões usam os mesmos registradores de dados e seleção.

| Simulação de 16 células, nove qubits | Soma reversível | QFT e fases |
|---|---:|---:|
| Portas CX compiladas | 6.764 | **241** |
| Profundidade compilada | 11.378 | **367** |

A comparação usou Qiskit 2.2.3, base lógica `u3/cx` e a mesma configuração de transpilação. A versão QFT reproduziu a ação Toeplitz por vetor de estado com diferença máxima inferior a `3,1×10^-15` na entrada complexa testada. A probabilidade de seletor zero dentro das células físicas foi de aproximadamente **1,90%** nessa entrada. O arquivo OpenQASM 2 preserva a fase global por uma sequência explícita de portas e foi reimportado na simulação. A contagem é de síntese lógica, anterior ao mapeamento para conectividade física. O circuito representa a parcela de primos na base de células. Os limites para Weil acima são obtidos pelos blocos intervalares de Legendre, incluindo as parcelas arquimediana e de polos. [Diagnóstico QFT](data/P252_PRIME_LCU.json) · [Diagnóstico soma reversível](data/P252_PRIME_LCU_REVERSIBLE.json)

## Reprodução

A partir da raiz do repositório:

```bash
python -m pip install -r NUMERIA_ES8P249/requirements.txt -r NUMERIA_ES8P252/requirements-quantum.txt
sha256sum -c NUMERIA_ES8P252/SHA256SUMS
python -m unittest NUMERIA_ES8P252.test_window_audit NUMERIA_ES8P252.test_prime_lcu -v
python -O -m unittest NUMERIA_ES8P252.test_window_audit NUMERIA_ES8P252.test_prime_lcu -v
python -m NUMERIA_ES8P252.window_audit \
  --matrix NUMERIA_ES8P252/data/P252_GL60.json.gz --out /tmp/P252_WINDOW_CERT.json
python -m NUMERIA_ES8P251.gate_b.mpfr_audit --half-width 19/20 \
  --matrix NUMERIA_ES8P252/data/P252_GL60.json.gz --out /tmp/P252_MPFR_AUDIT.json
python -m NUMERIA_ES8P252.prime_lcu --cells 16 \
  --out /tmp/P252_PRIME_LCU.json --qasm /tmp/P252_PRIME_LCU.qasm
python -m NUMERIA_ES8P252.prime_lcu --cells 16 --select reversible \
  --out /tmp/P252_PRIME_LCU_REVERSIBLE.json
```

Os nove testes da P252 e os 47 das etapas anteriores passaram em Python normal e com `-O`. O replay MPFR requer compilador C, GMP e `libmpfr.so.6` no Linux. Para regenerar os oito lotes originais de quadratura:

```bash
work=$(mktemp -d)
for start in 0 80 160 240 320 400 480 560; do
  stop=$((start+80))
  if [ "$stop" -gt 618 ]; then stop=618; fi
  python -m NUMERIA_ES8P251.gate_b.finite_head_quadrature \
    --start "$start" --stop "$stop" --half-width 19/20 \
    --bits 256 --bessel-bits 640 --threads 1 --out "$work/$start-$stop.json.gz"
done
python -m NUMERIA_ES8P252.window_audit --chunks "$work/"*.json.gz \
  --matrix "$work/P252_GL60.json.gz" --out "$work/P252_WINDOW_CERT.json"
cmp NUMERIA_ES8P252/data/P252_GL60.json.gz "$work/P252_GL60.json.gz"
cmp NUMERIA_ES8P252/data/P252_WINDOW_CERT.json "$work/P252_WINDOW_CERT.json"
```
