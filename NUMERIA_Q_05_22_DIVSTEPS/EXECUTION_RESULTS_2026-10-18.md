# Resultado da execução independente — divsteps

**Workflow:** https://github.com/jbarumby/scicalc/actions/runs/39131118399  
**Data:** 2026-10-18  
**Branch:** `research/divsteps-reproduction-2026-10`

## Resultado geral

O workflow concluiu com sucesso. Os testes clássicos do NUMERIA-Q passaram; o avaliador oficial upstream emitiu `verdict OK`; e a compilação Lean concluiu 1.518 tarefas. Isso valida a execução do circuito publicado sob o avaliador upstream, mas não prova novidade do NUMERIA-Q nem comparação de desempenho contra a implementação do NUMERIA-Q.

## 1. Testes clássicos NUMERIA-Q

- Resultado: sucesso.
- 3 testes aprovados.
- Estes testes cobrem a referência clássica de inversão modular, não a correção do circuito reversível.

## 2. Avaliador oficial do circuito divsteps

O log do avaliador upstream reporta:

| Métrica | Resultado |
|---|---:|
| Veredito | OK |
| Disparos de teste aleatórios | 100.000 |
| Disparos com saída incorreta | 4 |
| Disparos com sinal incorreto | 3 |
| Disparos falhos segundo o critério agregado | 4 |
| Orçamento permitido de falhas | 20 |
| Qubits de pico | 935 |
| Operações expandidas | 20.333.427 |
| Toffoli por disparo | 4.006.673 |
| Clifford por disparo | 15.391.673 |
| Custo Toffoli × qubits | 3.746.239.255 |
| Qubits sujos no fim | 0 |

A taxa agregada de falhas observada foi 4/100.000 = 0,004%, abaixo do orçamento de 20/100.000 = 0,02%. O avaliador aceitou o circuito. Não confundir esse resultado com fidelidade de estado quântico ou execução em hardware quântico: trata-se do avaliador de circuitos reversíveis do repositório upstream, usando seu modelo e teste aleatório.

Log completo e artefato: https://github.com/jbarumby/scicalc/actions/runs/39131118399/artifacts/12137999957

## 3. Prova Lean

- Resultado do build: `Build completed successfully (1518 jobs)`.
- Artefato com logs: https://github.com/jbarumby/scicalc/actions/runs/39131118399/artifacts/12138998509

### Ressalva de auditoria formal

O log de Lean informa que alguns teoremas dependem de axiomas, incluindo `example590_par._native.native_decide.ax_1_1`, além dos axiomas usuais `propext`, `Classical.choice` e `Quot.sound`. Portanto, o build compilou, mas antes de afirmar que todos os resultados são derivados sem axiomas adicionais, é necessário auditar essa dependência de `native_decide` e confirmar a intenção do projeto Lean upstream. Não se deve resumir o resultado como "toda a prova é axiom-free".

## Conclusão científica

- **Reprodução do avaliador upstream:** passou no ambiente GitHub Actions.
- **Compilação da prova Lean:** passou, com ressalva de dependências axiomáticas que precisam ser inspecionadas.
- **NUMERIA-Q ter descoberto uma melhoria:** não demonstrado.
- **Comparação quantitativa independente entre NUMERIA-Q e a arquitetura publicada:** ainda pendente.

Próximo passo: auditar o uso de `native_decide` nos arquivos Lean e comparar, em condições equivalentes, as arquiteturas do NUMERIA-Q com a referência upstream. Não enviar trabalhos à IBM Quantum para essa etapa.
