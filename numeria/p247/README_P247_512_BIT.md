# P247 — auditoria de precisão fixa a 512 bits

## Objetivo e alteração

Este marco habilita uma configuração opcional de 512 bits na auditoria de ponto fixo P245B/P245D. A precisão padrão permanece em 256 bits para preservar os resultados arquivados. A seleção é explícita por `NUMERIA_P245_FIXED_PRECISION=512`; valores diferentes de 256 ou 512 são rejeitados.

No modo 512-bit, o código aumenta o número de termos das séries de logaritmo e Bessel e da fórmula de Binet: `LOG_N=205`, `K=64` e `BERNOULLI_K=45`. Os restos dessas séries são calculados como racionais e devem ficar abaixo de uma unidade da grade ativa. O resto trigonométrico mantém seu limite analítico. Todos os intervalos são arredondados para fora com inteiros arbitrários.

## Resultados reproduzidos

Foram repetidos os blocos de painéis `[25,50]`, `[500,525]` e `[975,1000]`. Cada replay usa 25 painéis, 60 nós GL por painel, duas paridades e 13 pares selecionados por setor; registra 26 inclusões de somas de bloco contra as referências Arb. Os três JSONs contêm os intervalos, as frações exatas dos extremos das bolas Arb, os hashes das entradas e flags explícitas de escopo.

Nos pares selecionados que envolvem o índice 376, as caixas 512-bit foram mais estreitas que as de 256 bits em todos os blocos. A redução varia: no bloco 25–50 e no 500–525 foi aproximadamente `1,16 × 10^77` (cerca de `2^256`). No bloco terminal 975–1000, algumas entradas tiveram redução de aproximadamente `2^114` a `2^115`, pois a largura também depende dos intervalos propagados e dos majorantes analíticos, não apenas da grade de arredondamento. Assim, não se deve generalizar o fator `2^256` a toda a matriz.

O majorante analítico de alto grau continua sendo o majorante conservador zero-centrado do P245D. Esta alteração reduz a contribuição numérica de arredondamento, mas não troca o majorante por uma estimativa analítica de Bessel mais forte. Os brackets certificados das raízes GL continuam com 224 bits.

Os seis arquivos P237 compactados dos blocos foram incluídos em `references/`; seus SHA-256 são fixados nos replays e conferidos pelo teste. As raízes usam `numeria/p245/gate_c/P245B_ROOT_BRACKETS.json`, cujo SHA-256 também está em cada replay.

## Reprodução

Na raiz do repositório, execute:

```bash
export NUMERIA_P245_FIXED_PRECISION=512
export PYTHONPATH=numeria/p245/gate_c:numeria/p237
for start in 25 500 975; do
  stop=$((start + 25))
  python3 numeria/p245/gate_c/p245d_additional_chunk_replay.py \
    --references numeria/p247/references \
    --roots numeria/p245/gate_c/P245B_ROOT_BRACKETS.json \
    --start "$start" \
    --out "numeria/p247/data/P247D_REPLAY_${start}_${stop}_512.json"
done
python3 numeria/p247/test_p247_high_precision_replay.py
```

O teste verifica os três artefatos, recalcula com `fractions.Fraction` a contenção dos extremos Arb, valida os SHA-256 dos seis chunks e confirma a redução de largura frente aos certificados 256-bit arquivados. O replay não depende de Arb, MPFR, SciPy ou mpmath na etapa de certificação.

## Escopo lógico

Este resultado cobre apenas as inclusões selecionadas nesses três blocos. Não refaz todas as 142.506 entradas do triângulo de cada matriz, não verifica uma identidade geral entre a matriz e a forma de Weil/Zhu, não fecha a exaustão cofinal de janelas e não prova a Hipótese de Riemann. Os replays preservam explicitamente esses limites em campos booleanos falsos.
