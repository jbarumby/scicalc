# NUMERIA–RIEMANN: reconciliação do CI — 2026-10-22

## Mudanças locais nesta branch

- Acrescentado o branch `numeria-riemann/weil-v2-adaptive-ci-fix` ao gatilho `push` do workflow, para que o push desta branch dispare a auditoria.
- Acrescentadas execuções adaptativas a 256 bits para `N=8` e `N=10`, cada uma seguida de `verify_certificate --compare` contra o relatório correspondente a 192 bits.
- Mantidos os passos existentes; não foi alterado o código matemático nem os parâmetros já registrados.

## Verificações executadas nesta sessão

- `git apply --check` do patch anterior: passou antes da alteração final do gatilho.
- `git diff --check`: passou.
- Parse YAML do workflow: passou.
- `python -m py_compile numeria_riemann/*.py`: passou.
- Verificador racional exato para os relatórios N=6, 8, 10 a 192 bits: todos `verified=true`.
- Verificador racional exato para N=6, 8, 10 a 256 bits, com comparação contra 192 bits: todos `verified=true`, `consistent_with_other_report=true`; os enclosures se intersectam.
- Checagem independente em ponto flutuante para N=6, 8, 10 a 192 bits: todos os valores ficaram dentro dos enclosures.
- Testes base sem Arb/FLINT: 8 testes, todos passaram.

## Bloqueio remanescente

`python-flint` não está instalado no ambiente desta sessão (`ModuleNotFoundError: No module named 'flint'`). Por isso, a suíte intervalar e o teste ponta a ponta com zeros não foram executados novamente nesta sessão. O arquivo `emulacao-ci.txt` fornecido pelo usuário relata 30 testes OK (1 skipped) e um teste ponta a ponta OK, mas seu commit indicado (`7e446af`) não corresponde ao HEAD original do bundle (`1782d9b...`). Esses resultados permanecem evidência fornecida, não uma reprodução independente nesta sessão.

O GitHub Actions real também continua pendente: a branch original não foi encontrada no remoto. Este commit prepara uma branch explícita para publicação, mas **não afirma que o código foi enviado nem que o CI remoto passou**.

## Escopo matemático

Mesmo com o CI remoto aprovado, os certificados dizem respeito às matrizes finitas da família testada. Não provam a hipótese de Riemann nem a positividade global do critério de Weil.
