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

Fonte dos números: README da tarefa upstream, consultado via GitHub. Eles são resultados reportados pelos autores, **não uma reprodução independente realizada pelo NUMERIA-Q**.

Link direto: https://github.com/andigu/uncompute/blob/main/tasks/modinv-point-add/README.md

## O que foi verificado nesta etapa

- O repositório upstream disponibiliza circuitos, avaliador, arquivos de resultados e uma prova Lean relacionada à recuperação de decisões de divsteps.
- Os números de recursos acima foram lidos da documentação pública upstream.
- A referência clássica independente do NUMERIA-Q passou em 3 testes locais.

## O que ainda não foi verificado

- Execução local do avaliador upstream e conferência do veredito certificado.
- Execução da prova Lean `lean/divsteps`.
- Comparação de implementações com a mesma especificação, mesmos testes e mesmas convenções de custo.
- Qualquer melhoria nova do NUMERIA-Q em relação ao circuito publicado.

## Bloqueio de execução

A tentativa de clonar o repositório upstream no ambiente de execução falhou porque o ambiente não conseguiu resolver o host `github.com` (falha de DNS). Portanto, não é correto afirmar que o circuito ou a prova Lean foram executados aqui. Os comandos reproduzíveis permanecem documentados em `README.md`.
