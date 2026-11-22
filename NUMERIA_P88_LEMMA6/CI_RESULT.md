# NUMERIA P88 — Auditoria do Lema 6: resultado e alcance

**CI aprovada:** [GitHub Actions #45047223939](https://github.com/jbarumby/scicalc/actions/runs/45047223939), `completed/success`, 22/11/2026.

**Artefato:** [numeria-p88-mixed-compensation-lattice-sharpness](https://github.com/jbarumby/scicalc/actions/runs/45047223939/artifacts/14764107518), id `14764107518`.

**SHA-256 certificado pelo runner:** `p88-mixed-compensation.json` = `000d1fe13956b2ae67a6c06d1d405fac75a6b20b4d9e88224ac7c8b48cb247b2`. O artefato também contém `p88-SHA256SUMS.txt` com o código fonte e as dependências P86/P87.

## Avanços matemáticos novos

1. **Reescrita EXATA de L6 como desigualdade de compensação:**

   `E_ab=Delta_a+Delta_b+2*C_ab-2*(delta_a+delta_b)/h² >= 0`,

   onde `Delta_a=pi²/(3*delta_a)-sum_(k!=a)delta_k/(x_k-x_a)²` (conforme P87) e `C_ab=sum_(k!=a,b) delta_k/[(x_k-x_a)(x_k-x_b)]`. A identidade `h²(RHS_L6-LHS_L6)=E_ab` está demonstrada termo a termo no README e verificada exatamente nos exemplos.
2. **Novo teorema geral INDEPENDENTE, mas mais fraco que L6:** Para todos os espaçamentos reais e toda dimensão finita,

   `S_ab<=2/h²*[pi²/3*(1/delta_a+1/delta_b)-(delta_a+delta_b)/h²]`.

   A prova usa `(u-v)²<=2(u²+v²)` e a cota diagonal P87.
3. **Igualdade ótima EXATA na rede uniforme infinita:** se `x_k=kd` e os dois nós estão separados por `m>=1` passos,
   
   `S_ab=d⁻³*(4*zeta(2)/m²-6/m⁴)=2*pi²/(3*m²*d³)-6/(m⁴*d³)`.

   A soma cruzada `sum_(k!=0,m)1/[k(k-m)] =2/m²` é deduzida das três somas harmônicas absolutamente convergentes; a igualdade na expressão do Lema 6 demonstra a **necessidade assintótica do termo negativo com coeficiente 3**. Para truncamentos simétricos, o resto omitido satisfaz `0<R<=2/(3*M³*d³)`, demonstrado por comparação integral. Os testes usam limites racionais rigorosos de π, obtidos por Machin/arctan.

## CI e verificações

- **8 testes unitários passaram** em Python normal e **8/8** sob `python -O`.
- **747 pares de índices** em configurações finitas explícitas, todas as desigualdades com margens **estritamente positivas** ao substituir π por limite racional inferior rigoroso.
- **588 pares com termo cruzado negativo** nos exemplos: usar `C_ab>=0` seria um erro e impediria uma prova válida.
- Sete testes de truncamento de redes uniformes, incluindo diferentes passos racionais e distâncias entre índices; limites de cauda por frações exatas.
- Os testes reproduzem o compensador e exigem `E_ab=h²(RHS_L6-S_ab)` em todos os casos, sem `float`.
- A CI proíbe declarar `independent_full_Preissmann_Lemma6_proved=true`, `RZ_S_epsilon_reproved=true` ou `RH_proved=true`.

## Pendente

**Ainda não foi provado** `E_ab>=0` para TODOS os espaçamentos reais. Equivalência algébrica + casos finitos + igualdade em rede uniforme NÃO são, isoladamente, uma prova independente completa do Lema 6. A estimativa afiada continua sendo **teorema externo publicado por Preissmann (1984)**.

A fonte primária (Lemme 6, pp. 108–110) mostra que seu argumento usa uniformização de extremos e **Lema 5** sobre somas recíprocas, cujas fórmulas na extração textual são parcialmente ilegíveis. É necessário reconstruir esse ponto com cautela: não preencher fórmulas por adivinhação. [Yangjit, Lema 3.2, (3.6)](https://arxiv.org/pdf/2203.14950) confirma o enunciado exato.

**Próxima prioridade P89:** reconstruir o Lema 5 de Preissmann e sua aplicação às três regiões de índices para provar a compensação `E_ab>=0` em dimensão arbitrária. Só então retirar a dependência desse lema da cadeia de `F_X` e `m0`. A estimativa de Ramaré–Zúñiga para `S_epsilon(X)` e a região `3e10<=X<3e12` também permanecem abertas.

**Sem repetição de S3**, sem alteração dos certificados P84/P85 ou da janela de Suzuki `t<=79`; não há prova global de RH. PR permanece **draft**, sem merge automático.
