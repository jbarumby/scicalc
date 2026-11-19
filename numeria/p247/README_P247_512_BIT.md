# P247 — auditoria de precisão fixa a 512 bits

## Objetivo e alteração

Este marco habilita uma configuração opcional de 512 bits na auditoria de ponto fixo P245B/P245D. A precisão padrão continua sendo 256 bits para preservar os resultados arquivados existentes. A seleção é explícita por `NUMERIA_P245_FIXED_PRECISION=512`; qualquer outro valor diferente de 256 ou 512 é rejeitado.

No modo de 512 bits, o módulo aumenta os comprimentos das séries de logaritmo e Bessel e os termos de Bernoulli usados na fórmula de Binet: `LOG_N=205`, `K=64` e `BERNOULLI_K=45`. Os restos são calculados como racionais e os testes internos exigem que fiquem abaixo de uma unidade da grade ativa. O resto de Taylor trigonométrico permanece certificado pelo limite analítico já existente. Os intervalos são arredondados para fora com inteiros arbitrários, não por ponto flutuante.

## Resultado reproduzido

Foi reexecutado o bloco de painéis `[25, 50]`, com 25 painéis, 60 nós GL por painel, duas paridades e 13 pares selecionados por setor. O replay arquivado registra 26 inclusões de somas de bloco contra as referências Arb, com aritmética inteira de 512 bits. O arquivo `data/P247D_REPLAY_25_50_512.json` contém todos os intervalos, frações exatas dos extremos das bolas Arb, hashes das entradas e declara explicitamente que o triângulo completo de 377 modos não foi refeito.

Para os pares selecionados envolvendo o modo de índice 376, o intervalo em ponto fixo ficou aproximadamente `1,16 × 10^77` vezes mais estreito que no artefato de 256 bits correspondente. Esse fator é essencialmente `2^256`. A interpretação é estritamente numérica: o majorante analítico de alto grau continua sendo o majorante conservador zero-centrado já usado no P245D; a mudança diminui a largura induzida pela grade de arredondamento, mas não substitui o majorante por uma estimativa analítica de Bessel mais forte. A precisão do certificado de raízes GL continua em 224 bits.

As duas entradas compactadas do bloco foram incluídas em `references/`. Seus SHA-256 são fixados no replay e verificados pelo teste local. Os brackets de raízes vêm de `numeria/p245/gate_c/P245B_ROOT_BRACKETS.json` e seu hash está registrado no JSON do resultado.

## Reprodução

Na raiz do repositório, execute:

```bash
NUMERIA_P245_FIXED_PRECISION=512 \
PYTHONPATH=numeria/p245/gate_c:numeria/p237 \
python3 numeria/p245/gate_c/p245d_additional_chunk_replay.py \
  --references numeria/p247/references \
  --roots numeria/p245/gate_c/P245B_ROOT_BRACKETS.json \
  --start 25 \
  --out /tmp/P247D_REPLAY_25_50_512.json
```

Valide o artefato arquivado e os hashes dos dados com:

```bash
python3 numeria/p247/test_p247_high_precision_replay.py
```

A execução é determinística para esses insumos e não requer Arb, MPFR, SciPy ou mpmath na etapa de certificação; o replay consome os brackets de raízes já fixados. O teste adicional recalcula em `fractions.Fraction` a contenção de cada intervalo Arb publicado e valida o SHA-256 dos dois arquivos de entrada.

## Escopo lógico

Este experimento reforça somente a auditabilidade de inclusões selecionadas no bloco 25–50. Ele não certifica todas as 142.506 entradas do triângulo de cada matriz, não prova uma identidade geral entre a matriz e a forma quadrática de Weil/Zhu, não fecha a exaustão cofinal de janelas e não prova a Hipótese de Riemann. O artefato preserva essas limitações em campos booleanos explícitos, todos falsos para as afirmações globais.
