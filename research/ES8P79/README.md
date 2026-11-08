# NUMERIA–RIEMANN ES8-P79 — Arb independente em três janelas

Branch: research/numeria-es8p79-three-window-replay-2026-11-07

## Escopo científico

Preserva as bases históricas de Legendre por partes:
- a = 0,30: 24 modos, 8 pares estratificados da matriz com 156 entradas independentes;
- a = 0,34: 26 modos, 8 pares estratificados da matriz com 182 entradas independentes;
- a = 0,40: 44 modos, blocos [11,22,11], 14 pares estratificados de 990 entradas; inclui primo 2.

O piloto tem 30 pares, divididos em 6 jobs do GitHub Actions.

## Reprodução sem acb.integral

replay.py usa python-flint==0.9.0, ctx.prec=256, Gauss–Legendre intervalar e erro baseado em uma elipse de Bernstein em 12 trechos. O termo da digama é majorado globalmente, sem 64 amostras do contorno. Não usa funções gaussianas. Os dados sample_reference.json são consultados somente ao comparar a nova inclusão com a caixa histórica. Todos os inteiros grandes são strings decimais exatas, não números IEEE-754.

verify_results.py exige os seis JSONs reais, parâmetros, precisão, completude e compatibilidade em aritmética racional. Ausência de evidência aciona status FAIL_CLOSED.

## Rodar

No GitHub: Actions → NUMERIA ES8-P79 three windows Arb sample.

O workflow é acionado por push na branch de pesquisa; um commit posterior à instalação do workflow foi previsto para disparo. Se a plataforma não registrar esse push, a opção workflow_dispatch pode depender da existência do workflow na branch padrão.

Com Python 3.11 instalado:

    cd research/ES8P79
    python -m unittest -v test_replay
    python -m pip install python-flint==0.9.0
    python replay.py --window a030 --shard 0 --out a030_sample_0.json

O job verify-artifacts reagrupa as seis saídas e valida as 30 entradas.

## Portões de aceitação

1. Seis jobs Arb terminam sem erros e emitem seis artefatos.
2. bits=256 em cada JSON e python-flint presente.
3. Os 30 pares e seus índices são completos, sem duplicações.
4. Nenhum par apresenta intervalos disjuntos.
5. O resultado consolidado é VERIFIED_SAMPLE_ARTIFACTS e foi conferido nos logs.
6. Um segundo auditor examina os limites analíticos de Bernstein/digama antes de afirmar que as caixas são rigorosas.

## Limitações

Não há neste commit reprodução das 1.328 entradas. Nenhuma conclusão universal sobre Weil ou a Hipótese de Riemann é pretendida. Um workflow publicado ou disparado não equivale a prova de conclusão dos jobs.

Os arquivos de referência são extratos com proveniência SHA-256 dos JSONs completos da ES8-P76; reproduzir os 1.328 pares requer publicar os dados e lançar jobs separados após o piloto.
