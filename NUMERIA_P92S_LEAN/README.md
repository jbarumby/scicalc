# NUMERIA P92-S — Formalização parcial em Lean 4 (núcleo R1–R7 da PR #107)

**Objetivo:** verificar no kernel Lean, sem \`sorry\` ou \`admit\`, os passos **algébricos e combinatórios** da reconstrução P92-R. A prova geral publicada de Preissmann continua matematicamente válida na literatura; aqui se formaliza a reconstrução NUMERIA. A presente etapa **NÃO** constitui formalização completa de Lemas 5/6, nem da Hipótese de Riemann.

## Ambiente

- Lean 4.35.0; Mathlib \`v4.35.0\`, iguais ao projeto \`NUMERIA_E8K_FORMAL\` deste repositório.
- Pacote Lake isolado em \`NUMERIA_P92S_LEAN\`.
- Build formal em GitHub Actions, exigindo compilação sem erros e rejeitando \`sorry\` e \`admit\`.
- Fonte principal: [NumeriaP92S.lean](NumeriaP92S.lean).

## Mapa dos resultados

| Enunciado Lean | Escopo formal exato | Situação |
|---|---|---|
| \`two_pole_partial_fractions\` | identidade racional de dois polos para \`u,h,u+h ≠ 0\` | prova algébrica por \`field_simp; ring\` |
| \`finite_two_pole_identity\` | aplicação somada sobre qualquer conjunto finito | usa somente a identidade anterior |
| \`pathCost_le_of_bellman\` | indução retroativa finita R1⇒R2 | **condicional:** assume \`base\` e \`step\`, que são a desigualdade de Bellman para o núcleo |
| \`boundary_accounting\` | cancelamento algébrico do truncamento \`f_X\` e erro \`(X-Ka)/X\` | **condicional:** assume o controle de \`f_X\` obtido via R2 e a soma das distâncias |
| \`sharp_bound_of_energy_and_compensator\` | combinação R6 + duas cotas de energia + identidade exata ⇒ constante \`-3\` de R7 | **condicional:** R6 e as duas cotas de energia são hipóteses, e a igualdade já sumarizada também |
| \`sharp_rhs_positive\` | positividade dimensionalmente correta do lado direito quando \`0<ds,dt≤h\` e \`K>3\` | prova algébrica independente |

**Interpretação:** a formalização confirma as **costuras lógicas explícitas** da reconstrução. Não prova ainda a desigualdade R1 para rampas convexas, a representação integral de misturas positivas, a comparação unilateral R2 em sua aplicação analítica, a etapa R5→R6 com a convergência pareada, nem a transferência infinita da cauda R4.

## Por que a formalização por partes é necessária

A auditoria tripla enviada confirmou com frações exatas e simbólico os componentes, mas não era uma prova de kernel para qualquer número real. O Lean obriga que os pressupostos não comprovados apareçam como parâmetros: neste arquivo não há uma forma de usar R6 ou Bellman sem fornecê-los no local de utilização.

O núcleo de frações parciais mantém os sinais corretos da P92-P, e \`sharp_rhs_positive\` corrige a escala do microlema C-iii da auditoria: a expressão natural é \`h³ · RHS\`, não \`h² · RHS\`.

## Reprodução

Da raiz do repositório:

\`\`\`bash
cd NUMERIA_P92S_LEAN
lake update
lake exe cache get
lake build NumeriaP92S
\`\`\`

O workflow \`.github/workflows/numeria-p92s-lean-formal-core.yml\` faz a instalação/compilação em GitHub Actions.

## Próximo gate: R1 não condicional

A próxima PR deverá trocar o parâmetro abstrato \`step\` de \`pathCost_le_of_bellman\` por uma prova em Lean da monotonicidade das somas para rampas \`(T-u)_+\` e da representação por medida positiva das funções convexas. Depois, formalizar R5→R6 como limite de somas pareadas absolutamente convergentes (sem assumir que \`Σ_s\` e \`Σ_t\` isoladamente convergem). Só então poderemos falar em formalização completa do Lema 5 e do Lema 6.

**Governança:** PR draft sobre \`numeria-p92r-lemma5-bellman-reconstruction\`. Não fazer merge automático. Não inferir RH.
