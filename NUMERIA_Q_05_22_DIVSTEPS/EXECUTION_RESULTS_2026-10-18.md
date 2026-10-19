# Resultado da execução e auditoria independente — divsteps

**Branch:** `research/divsteps-reproduction-2026-10`  
**Atualização documental:** 2026-10-19  
**Execução CI principal:** https://github.com/jbarumby/scicalc/actions/runs/39131118399  
**Reexecução CI do avaliador:** https://github.com/jbarumby/scicalc/actions/runs/39180629421  
**Revisão upstream fixada para execuções futuras:** `andigu/uncompute@64b5a13aa249f1360846657392e6c66817b26d0f`

## 1. Resultados CI observados

Na execução 39131118399, os três jobs — testes clássicos, avaliador upstream e build Lean — terminaram com sucesso. O job do avaliador usou o executor de referência em Python puro e levou cerca de 16 min 35 s na execução de 18/10. Na execução 39180629421, os testes clássicos registrados continuaram em 3/3 e o avaliador levou cerca de 12 min 11 s, também com executor Python puro. Essas durações são tempos do job no CI, não tempos de uma reexecução local com kernel C.

O workflow antigo baixava `andigu/uncompute@main` e não registrava o SHA efetivamente usado. O workflow foi atualizado neste branch para fixar o commit e registrar procedência e artefatos estruturados nas execuções futuras. Isso não permite inferir retroativamente o commit exato das execuções históricas.

## 2. Avaliador oficial do circuito

| Métrica | Resultado observado |
|---|---:|
| Veredito | `OK` |
| Disparos aleatórios | 100.000 |
| Saídas incorretas | 4 |
| Sinais incorretos | 3 |
| Disparos falhos agregados | 4 |
| Orçamento de falhas | 20 |
| Qubits de pico | 935 |
| Operações expandidas | 20.333.427 |
| Toffoli por disparo | 4.006.673 |
| Clifford por disparo | 15.391.673 |
| Custo Toffoli × qubits | 3.746.239.255 |
| Qubits sujos no fim | 0 |

As categorias de sinal incorreto e saída incorreta se sobrepõem; não se devem somar 4 e 3 como sete falhas. O veredito confirma a aceitação pelo critério do avaliador, não perfeição do circuito nem execução em hardware quântico.

O log acessível registra apenas contagens agregadas. Não contém os valores das quatro entradas que falharam; portanto, a causa específica ainda não está demonstrada.

## 3. Investigação das quatro falhas e agenda de 533 passos

O README upstream diz que o circuito agenda 533 passos e atribui caminhadas raras mais longas ao orçamento de falhas. Isso torna o excesso de passos uma hipótese plausível, mas os logs agregados não permitem confirmá-la para as quatro falhas observadas.

Foi executada uma estimativa independente usando a recorrência `delta,f,g` em Python, semente fixa 20261018 e 100.000 entradas uniformes em `[1,p-1]`, com `p = 2^256 - 2^32 - 977`:

- mínimo observado: 502 passos;
- máximo observado: 532 passos;
- média: 516,83705 passos;
- entradas com mais de 533 passos: 0;
- entradas com mais de 590 passos: 0.

Essa amostra independente não é a mesma sequência de entradas usada pelo avaliador e não é prova do pior caso. Zero ocorrências acima de 533 não permite identificar a causa dos quatro disparos falhos. A próxima investigação deve preservar ou registrar as entradas exatas dos disparos falhos, se o avaliador suportar essa instrumentação, e comparar cada entrada com a referência clássica. Até lá, a causa fica marcada como **não confirmada**.

## 4. Referência Python de hddivsteps

A referência local usa a recorrência meio-delta e compara os inversos com `pow(x, -1, p)` ou verifica `x * inv ≡ 1 (mod p)`. Os cinco novos testes passaram localmente: primo correto, casos de borda, 1.500 entradas aleatórias, rejeição de entradas não canônicas e propriedade de inverso. A execução Monte Carlo de 100.000 amostras reproduziu média 516,83705 e máximo 532.

Limites: é uma referência clássica escrita separadamente; não simula nem prova o circuito reversível, não mede recursos quânticos e não substitui a prova Lean. O script de estimativa é informativo, não critério de aceitação.

## 5. Lean e `native_decide`

Na revisão upstream fixada, `Hddivsteps590/Certificate.lean` prova `example590_par` com `native_decide`; `example590` deriva do certificado por `certEmptyPar_eq`. `Main.lean` imprime as dependências de axiomas de `hddivsteps590_prime_inverse`, `zWalk590_bound`, `zWalk590` e de resultados de recuperação.

A documentação e o código confirmam que a dependência de `native_decide` é intencional e faz parte da base de confiança computacional declarada. A avaliação compilada do certificado implica confiar no caminho de compilação/runtime relevante, além do kernel e dos axiomas usuais de Lean/Mathlib. Não é evidência, por si só, de prova incorreta; também não é correto chamar os teoremas de limite de 590 passos de exclusivamente kernel-verificados.

O build CI concluiu 1.518 tarefas, mas esta auditoria ainda não fez uma execução independente, teorema por teorema, dos resultados `#print axioms` nem construiu uma alternativa kernel-only. Essas tarefas permanecem abertas.

Fontes upstream fixadas:
- [Certificate.lean](https://github.com/andigu/uncompute/blob/64b5a13aa249f1360846657392e6c66817b26d0f/lean/divsteps/Hddivsteps590/Certificate.lean)
- [Main.lean](https://github.com/andigu/uncompute/blob/64b5a13aa249f1360846657392e6c66817b26d0f/lean/divsteps/Hddivsteps590/Main.lean)
- [README Lean](https://github.com/andigu/uncompute/blob/64b5a13aa249f1360846657392e6c66817b26d0f/lean/divsteps/README.md)

## 6. Procedência, métricas e comparação

Baseline observado: 4.006.673 Toffoli, 15.391.673 Clifford e 935 qubits de pico; `h_ops = f3436a16c61b3e91109b9cb65cc130d4c00cd07f1e4ab977248d36815c23fcf1`. A conversão nominal `4.006.673 × 7 = 28.046.711 T` é apenas estimativa aritmética, não contagem T obtida por síntese. Qubits lógicos e profundidade T não foram reportados pelo avaliador.

A comparação de recursos com NUMERIA-Q permanece **N/D**: não há nesta trilha uma implementação NUMERIA-Q funcional e medida sob a mesma especificação, decomposição, tratamento de ancillas e convenção de custo. Nenhum delta percentual deve ser inferido.

## 7. Inconsistências e pendências

- O patch proposto fora do branch menciona 134 s no texto e 95,2637 s em `report.json`; esses números não são reconciliados aqui e não foram usados como duração verificada.
- A referência a reexecução com kernel C precisa de log/procedência verificáveis antes de ser declarada confirmada; os dois jobs CI observados usaram Python puro.
- Ainda faltam entradas de falha reproduzíveis e diagnóstico da causa.
- Ainda falta auditoria independente da base de axiomas e uma estratégia kernel-only, se exigida.
- Ainda falta comparação equivalente NUMERIA-Q vs. referência.

**Restrições preservadas:** nenhum merge; nenhuma alteração em `main` ou nas versões congeladas 05.19.1/05.20; nenhum novo job enviado à IBM Quantum.
