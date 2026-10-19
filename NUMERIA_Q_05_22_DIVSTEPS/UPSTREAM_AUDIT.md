# Auditoria da referência upstream — divsteps

**Data da consulta:** 2026-10-18  
**Repositório consultado:** https://github.com/andigu/uncompute  
**Tarefa:** `tasks/modinv-point-add`  
**Artigo:** https://arxiv.org/abs/2610.22367

## Valores publicados no README upstream

| Construção | Toffoli | Qubits de pico | Toffoli × qubits |
|---|---:|---:|---:|
| Divsteps, configuração padrão (`j_store=5`) | 4.006.673 | 935 | 3.746.239.255 |
| Divsteps, variante `wide` | 3.293.278 | 1.191 | 3.922.294.098 |
| Divsteps, variante `fill` | 3.637.655 | 1.031 | 3.750.422.305 |
| Divsteps, variante `j5` | 4.317.926 | 927 | 4.002.717.402 |
| Divsteps, variante `j6` | 4.427.338 | 913 | 4.042.159.594 |
| Divsteps, variante `j7` | 4.649.514 | 904 | 4.203.160.656 |
| Divsteps, variante `j8` | 5.331.729 | 897 | 4.782.560.913 |
| Divsteps, variante `j9` | 7.866.779 | 893 | 7.025.033.647 |
| Euclidiano especializado | 20.521.949 | 843 | 17.300.003.007 |
| Baseline Luo et al. | 69.529.757 | 838 | 58.265.936.366 |

Fonte dos números: README da tarefa upstream. São resultados reportados pelos autores, não uma descoberta original do NUMERIA-Q.

Link direto: https://github.com/andigu/uncompute/blob/main/tasks/modinv-point-add/README.md

## Auditoria de confiança da prova Lean

A documentação pública de `lean/divsteps/README.md` declara que a prova usa os axiomas padrão reportados pelo Mathlib (`propext`, `Classical.choice`, `Quot.sound`) e que os teoremas de limite de 590 passos dependem adicionalmente de `native_decide`, gerado para `Hddivsteps.example590_par`. O projeto explica que essa dependência representa a avaliação compilada de um certificado Booleano fechado e requer confiança no compilador/runtime de Lean para esse cálculo, além do kernel e Mathlib.

O README upstream também afirma que `Main.lean` termina com comandos `#print axioms` para os teoremas públicos e que o projeto não contém `sorry`. Assim, o achado `native_decide` é uma dependência de confiança declarada pela própria referência, não evidência isolada de erro lógico. A ressalva formal permanece importante: não descrever os teoremas de 590 passos como verificados exclusivamente pelo kernel.

## O que foi verificado nesta etapa

- Os resultados de recurso foram lidos da documentação pública upstream.
- O avaliador upstream foi executado em GitHub Actions: 100.000 testes aleatórios, 4 falhas agregadas diante de orçamento 20, veredito `OK`.
- A compilação Lean foi concluída com sucesso (1.518 tarefas).
- Os testes clássicos independentes do NUMERIA-Q passaram (3/3).
- A documentação upstream esclarece explicitamente a dependência adicional de `native_decide`.

## O que ainda não foi verificado

- Comparação quantitativa independente entre implementações NUMERIA-Q e upstream sob especificação e custos equivalentes.
- Qualquer melhoria nova do NUMERIA-Q em relação ao circuito publicado.
- Execução em hardware quântico; não foi enviado nenhum novo job à IBM Quantum.

## Conclusão

A reprodução do avaliador upstream passou e a prova Lean compilou. A dependência `native_decide` está explicitada na base de confiança do upstream. A etapa de maior valor científico agora é a comparação controlada entre implementações, sem alegação de novidade antes dos resultados.
