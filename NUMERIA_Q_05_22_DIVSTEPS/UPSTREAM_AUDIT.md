# Auditoria da referência upstream — divsteps

**Data da atualização:** 2026-10-19  
**Revisão fixada para auditoria:** `64b5a13aa249f1360846657392e6c66817b26d0f`  
**Repositório:** https://github.com/andigu/uncompute

## Resultados de recursos

| Construção | Toffoli | Qubits de pico | Toffoli × qubits |
|---|---:|---:|---:|
| Divsteps, padrão (`j_store=5`) | 4.006.673 | 935 | 3.746.239.255 |
| Divsteps, variante `wide` | 3.293.278 | 1.191 | 3.922.294.098 |
| Divsteps, variante `fill` | 3.637.655 | 1.031 | 3.750.422.305 |
| Divsteps, variante `j5` | 4.317.926 | 927 | 4.002.717.402 |
| Divsteps, variante `j6` | 4.427.338 | 913 | 4.042.159.594 |
| Divsteps, variante `j7` | 4.649.514 | 904 | 4.203.160.656 |
| Divsteps, variante `j8` | 5.331.729 | 897 | 4.782.560.913 |
| Divsteps, variante `j9` | 7.866.779 | 893 | 7.025.033.647 |
| Euclidiano especializado | 20.521.949 | 843 | 17.300.003.007 |
| Baseline Luo et al. | 69.529.757 | 838 | 58.265.936.366 |

São resultados reportados pelo upstream, não uma descoberta do NUMERIA-Q.

## Confiança formal e `native_decide`

O arquivo `Hddivsteps590/Certificate.lean` declara `example590_par` com `native_decide` sobre o certificado Booleano `certEmptyPar Mstar 590`. O teorema `example590` deriva desse resultado usando `certEmptyPar_eq`. A versão paralela usa `Task.spawn`/`Task.get`; o próprio código prova que ela é definicionalmente equivalente à especificação serial, descrevendo a paralelização como otimização de velocidade.

`Hddivsteps590/Main.lean` contém `#print axioms` para os teoremas públicos, incluindo `hddivsteps590_prime_inverse`, `zWalk590_bound` e `zWalk590`. Esses teoremas de limite dependem da prova do certificado e, portanto, herdam a dependência `native_decide`. A documentação upstream explicita que isso amplia a base de confiança para incluir o caminho de compilação/runtime usado na avaliação nativa, além do kernel e dos axiomas usuais de Lean/Mathlib.

**Verificado:** código fonte da revisão fixada e logs CI mostram build concluído.  
**Ainda em aberto:** uma auditoria independente de cada saída `#print axioms`; uma prova alternativa sem `native_decide`; qualquer alegação de que a dependência torna a prova incorreta (não há evidência para essa alegação).

## Avaliador e falhas

Duas execuções CI observadas retornaram `verdict OK`, com 100.000 disparos, quatro falhas agregadas, três sinais incorretos, zero qubits sujos, 4.006.673 Toffoli, 15.391.673 Clifford, 935 qubits e o mesmo hash de operações `f3436a16c61b3e91109b9cb65cc130d4c00cd07f1e4ab977248d36815c23fcf1`.

O README upstream informa que o circuito agenda 533 passos e que caminhadas raras mais longas são cobradas contra o orçamento de falhas. A amostra independente de 100.000 entradas uniformes com semente 20261018 teve máximo 532 e zero caminhadas acima de 533; isso é compatível com uma ocorrência rara não amostrada, mas não explica as quatro falhas do avaliador. Os logs acessíveis só mostram contagens agregadas, não as entradas falhas. A causa permanece **não confirmada**.

## Estado científico

- Reproduzir a avaliação do circuito publicado: sim, pelo avaliador upstream.
- Validar o build Lean: sim; ressalva de confiança `native_decide` explícita.
- Comparar custos NUMERIA-Q com upstream em condições equivalentes: pendente.
- Demonstrar novidade ou vantagem quântica: não demonstrado.
- Execução em hardware quântico: não realizada; nenhum novo job IBM Quantum enviado.
