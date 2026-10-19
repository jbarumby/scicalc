# NUMERIA-Q — Marco de reprodução independente de divsteps

**Branch:** `research/divsteps-reproduction-2026-10`  
**Estado atualizado:** 2026-10-19 — avaliador upstream reproduzido no CI; auditoria de confiança Lean documentada; comparação NUMERIA-Q ainda pendente.

## Objetivo

Estabelecer uma trilha rastreável para comparar o protótipo NUMERIA-Q com o circuito reversível de inversão modular publicado por Gu, Kemeny, Yelin e Lukin (2026). Não se reivindica descoberta, equivalência com o circuito publicado nem vantagem quântica nesta etapa.

## Referência primária e procedência

- Artigo: https://arxiv.org/abs/2610.22367
- Código, circuitos e avaliador: https://github.com/andigu/uncompute
- Tarefa: `tasks/modinv-point-add`
- Circuito: `tasks/modinv-point-add/circuits/divsteps/circuit.py`
- Prova: `lean/divsteps`
- Commit fixado pelo workflow: `64b5a13aa249f1360846657392e6c66817b26d0f`

O workflow registra o commit efetivo, versões de Python/uv e SHA-256 do circuito, contrato e lockfile. A fixação torna futuras execuções repetíveis; não prova retroativamente qual commit foi usado nas execuções históricas.

## O que está verificado

1. A execução CI de 18/10 concluiu os três testes clássicos então existentes.
2. O avaliador oficial aceitou o circuito em 100.000 disparos: `verdict OK`, 4 disparos falhos diante de orçamento 20, 935 qubits de pico, 4.006.673 Toffoli, 15.391.673 Clifford e custo Toffoli × qubits de 3.746.239.255.
3. O build Lean upstream concluiu 1.518 tarefas.
4. Uma referência independente Python para `hddivsteps` foi testada localmente: 5 testes passaram, incluindo casos de borda e 1.500 entradas aleatórias. Esses testes são validação algébrica amostral, não prova formal nem teste do circuito reversível.
5. Uma amostragem independente de 100.000 entradas uniformes (semente 20261018) observou mínimo 502, máximo 532 e média 516,83705 passos; nenhuma amostra excedeu 533. A amostragem não prova um limite universal.

Detalhes e limites: [relatório de execução](EXECUTION_RESULTS_2026-10-18.md) e [auditoria de 19/10](AUDIT_2026-10-19.md).

## Critérios de aceitação

- [x] Testes clássicos preexistentes passam no CI (3/3 na execução registrada).
- [x] Avaliador upstream emite `verdict OK` e registra as contagens.
- [x] Build Lean conclui; a ressalva de `native_decide` está documentada.
- [x] Reexecução local da referência Python `hddivsteps`: 5/5 testes nesta auditoria.
- [x] Workflow passa a fixar a revisão upstream e a guardar evidências de procedência nas próximas execuções.
- [ ] Reproduzir os quatro disparos falhos com entradas identificadas; os logs disponíveis registram contagens agregadas, não os valores das entradas que falharam.
- [ ] Auditar independentemente a lista de axiomas de todos os teoremas relevantes e avaliar uma alternativa sem `native_decide`.
- [ ] Comparar NUMERIA-Q e upstream com a mesma especificação e convenções de custo.

## Limites e transparência

O avaliador `OK` significa aceitação segundo a política do avaliador, não perfeição nem execução em hardware quântico. O build Lean não significa que os teoremas de limite sejam verificados exclusivamente pelo kernel: os resultados de limite de 590 passos dependem de `native_decide`, conforme a documentação upstream. Não há ainda comparação quantitativa de implementação própria NUMERIA-Q. Nenhum novo job foi enviado à IBM Quantum; `main` e as versões congeladas 05.19.1/05.20 não são alvo desta trilha.
