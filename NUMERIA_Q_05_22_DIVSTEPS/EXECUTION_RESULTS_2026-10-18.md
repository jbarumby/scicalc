# Resultado da execução independente — divsteps

**Workflow:** https://github.com/jbarumby/scicalc/actions/runs/39131118399  
**Data da execução:** 2026-10-18  
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

## 3. Auditoria da dependência `native_decide` em Lean

A documentação upstream de `lean/divsteps/README.md` esclarece expressamente que os teoremas de álgebra, teoria dos números, ponte entre recorrências e recuperação de ramos usam os axiomas usuais reportados pelo Mathlib (`propext`, `Classical.choice`, `Quot.sound`). Os teoremas de limite de 590 passos dependem adicionalmente do axioma `native_decide` gerado para `Hddivsteps.example590_par`.

A documentação descreve essa dependência como parte explícita da base de confiança: o certificado Booleano fechado é avaliado por compilação e runtime do Lean, e aceitar os teoremas que dependem dele implica confiar também nesse caminho de avaliação, além do kernel Lean e Mathlib. A documentação informa que `Main.lean` contém comandos `#print axioms` para os resultados públicos e que o projeto não contém `sorry`.

**Conclusão precisa:** o build Lean passou (1.518 tarefas), e a dependência adicional de `native_decide` é documentada intencionalmente pelo projeto upstream; não é, por si só, evidência de prova incorreta. Entretanto, a garantia formal dos teoremas de limite de 590 passos tem essa base de confiança computacional ampliada. Não é correto descrevê-los como verificados exclusivamente pelo kernel sem ressalva.

Referência upstream: https://github.com/andigu/uncompute/tree/main/lean/divsteps

## 4. Limites científicos

- **Reprodução do avaliador upstream:** passou no ambiente GitHub Actions.
- **Compilação Lean:** passou; a dependência `native_decide` e sua base de confiança foram esclarecidas pela documentação upstream.
- **NUMERIA-Q ter descoberto uma melhoria:** não demonstrado.
- **Comparação quantitativa independente entre NUMERIA-Q e a arquitetura publicada:** ainda pendente.
- **Execução em hardware quântico:** não realizada; nenhum job novo foi enviado à IBM Quantum.

Próximo passo: comparar implementações sob a mesma especificação, testes e convenções de custo, sem confundir reprodução upstream com resultado novo do NUMERIA-Q.
