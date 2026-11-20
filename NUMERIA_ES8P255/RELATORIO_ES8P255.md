# NUMERIA–RIEMANN · ES8-P255 · Três prioridades: RH formal, método e ponte NUMERIA

**Data de referência:** 20 de novembro de 2026. **Escopo:** auditoria documental/estática da prova Lean anunciada por Levent Alpöge e Claude, comparação matemática com o NUMERIA P249–P254 e priorização de verificações subsequentes. Não há uma prova independente produzida nesta etapa.

## 1. Cinco teoremas e definição padrão

Consultados diretamente no GitHub, na revisão do espelho `40e56bab` apontada pelo PR #58, os arquivos `src/Main.lean`, `src/Statements.lean` e `src/Supplements/Kernels/MainTheorems.lean`. Seus **SHA-256 conferem** com a auditoria registrada no próprio NUMERIA #58. A definição de `RiemannHypothesis` foi lida na Mathlib do commit `0df444a360eaa60ab8c11dca51a86af692955474`: para todo `s : ℂ`, se `riemannZeta s = 0`, `s` não é um inteiro par negativo e `s ≠ 1`, então `s.re=1/2`. A revisão Git blob do arquivo conferida pelo conector foi `3c92195fce20df7fda41cfd823b3266f83cf3070`.

| Constante Lean principal | Tipo identificado | Dependência imediata declarada |
|---|---|---|
| `riemannHypothesis` | `RiemannHypothesis` (definição da Mathlib) | `rhOfAll` |
| `riemannZeta_eq_zero_re` | `RiemannZetaZerosOnLine` | `rhOfAll` |
| `LFunction_ne_zero_of_half_lt_re` | `DirichletLNonvanishingRightOfHalf` | `allChi` |
| `LFunction_eq_zero_re` | `DirichletLZerosOnLine` | `claimAll.2.1` |
| `LFunction_eq_zero_re_of_isPrimitive` | `PrimitiveDirichletLZerosOnLine` | `critAll` |

A não anulação geral não é deduzida pela definição de RH: `allChi` é declarado como instância de `UpperTheory.FinalTheorems.everyCharacter`, e `rhOfAll` usa `Formalism.VW.rh_of_every allChi`. **O centro lógico a auditar é essa vasta cadeia de dependências**, incluindo a correção dos enunciados e a conformidade com os axiomas padrão.

**Metodologia e ressalvas:** o programa `audit_statements.py` autentica por SHA-256 três arquivos literais de origem e verifica o encadeamento sintático; opcionalmente, se o arquivo completo Mathlib for fornecido, confere seu hash de blob Git e a definição explícita. O relatório JSON registra `mathlib_file_locally_verified = null`, porque a verificação integral da Mathlib foi feita por consulta ao código autenticado no GitHub, não pela leitura local de um arquivo dessa dimensão. **Não invoca Lean, não exporta a solução completa e não reproduz o kernel.** O desafio MIT tem `sorryAx` de propósito; o verificador do desafio com esse axioma permitido não prova os teoremas.

O MIT (Andrew Sutherland, atualização 2026-11-20 19:45 UTC) relatou: 175.187 módulos compilados; verificação de dependências com os três axiomas usuais; e aceitação da exportação completa da solução pelo nanoda. O Comparator e outro replay estavam ainda pendentes naquele registro. A auditoria do NUMERIA #58 autenticou os 14 arquivos compactados e os tipos do desafio, mas também não reproduziu o kernel da solução completa localmente. **O Clay ainda classifica RH como não resolvida em sua página.**

Fontes: https://math.mit.edu/~drew/grh/ · https://alpo.ge/remark.pdf · https://github.com/jbarumby/scicalc/pull/58 · https://www.claymath.org/millennium/Riemann-Hypothesis/

## 2. Estratégia matemática do anúncio — o que foi efetivamente apresentado

A nota **não é uma demonstração matemática completa em texto corrido**; o autor a descreve como uma cronologia preliminar. A linha argumentativa anunciada é:

1. Método de **Stepanov** inspirado na hipótese de Riemann para corpos de funções: um determinante com multiplicidades força forte anulação.
2. **Guinand–Weil/Fourier** organiza uma forma com sinais e um teste de variância. O autor passa à soma ponderada
   
   `Σ_{n≤X} Λ(n)·(X/n−1)_+ − X log(X) + (1+γ)X`.
3. A multiplicidade, problemática sobre os inteiros, seria simulada por **levantamentos theta iterados** entre grupos unitários totalmente definidos. Para a soma de von Mangoldt de comprimento `X`, a cadeia proposta tem comprimento na ordem de `√X`.
4. O sumário do autor enumera estimativas para normas e níveis das elevações, divisibilidade, identidades seesaw, fórmulas de traço relativas, estruturas integrais sobre caracteres de Hecke e estimativas de covolume.
5. A passagem aos teoremas `allChi`/`rhOfAll` é codificada em Lean, mas a nota não especifica os lemas quantitativos completos suficientes para reconstruir os cálculos sem estudar os módulos formais.

**Não está demonstrado por esta análise** que a soma ou determinante denotado por `W(1_S)` na nota seja precisamente a forma quadrática e a normalização de Weil empregadas nas matrizes do NUMERIA. Não basta uma semelhança de notação.

Fonte primária: https://alpo.ge/remark.pdf (páginas 1, 3 e 4).

## 3. Comparação com NUMERIA P249–P254 e conservação de resultados

| Frente | NUMERIA existente | Relação com Alpöge | Trabalho seguinte necessário |
|---|---|---|---|
| Forma explícita | P251 Gate C relaciona Fourier–Mellin à forma de Weil em suporte limitado; P252 relata uma margem local até L=0,95 | A nota de Alpöge cita Guinand–Weil | Comparar **normalizações, transformadas e classes de funções** termo a termo, sem assumir identidade |
| Primos | P249 usa `S_L(t)=Σ_{p^k<e^{2L}}2log(p)/√(p^k)cos(tlog(p^k))`; P250 melhora o corte de cauda; P253 enfrenta o deslocamento de 7 | Alpöge usa soma **ponderada até X**, com `Λ(n)(X/n−1)_+` | Derivar, se possível, uma ponte explicita entre as duas janelas e a ponderação; não identificar expressões diferentes |
| Positividade | P251 B e P252 fornecem certificados finitos com hipóteses analíticas; P253 registra direções negativas em **modelos auxiliares**; P254 propõe piso 1/4 após 504 | Alpöge anuncia um argumento por determinantes e levantamentos theta | Investigar um enunciado **quantitativo de forma** ou norma que possa ser testado nas mesmas bases; RH formal, mesmo se validada, não fornece automaticamente margem coerciva explícita da matriz |
| Computação | Operadores Toeplitz/FFT e replays Arb/MPFR verificam integrandos | O argumento de Alpöge é automórfico, não uma simples matriz densa | Manter replays P249–P254 e não importar hipóteses automórficas sem provas/autenticação |
| Escopo | Provas/condicionais de janelas compactas e diagnósticos finitos | Teoremas Lean de RH e GRH para funções L de Dirichlet | Acompanhar validação integral em fontes primárias e produzir reprodução **autônoma** de um lema central, não apenas hashes |

**Preservação:** nenhum código, certificado, branch ou PR de P249–P254 foi modificado pela elaboração desta P255. A auditoria será publicada numa branch nova, derivada do último commit P254, e submetida em PR rascunho. Não mesclar sem revisão.

## Recomendação técnica imediata

A próxima pesquisa útil não é uma nova execução de matrizes gigantes. Escolha um único lema candidato no desenvolvimento Lean, relacionado a **somas de von Mangoldt ponderadas, não anulação de L-functions ou determinantes de positividade**, extraia seu tipo completo e **todas as hipóteses**, localize seu grafo de dependências e procure um teste ou demonstração legível em escala menor. Comparar esse lema ao controle de primos da P249–P253 sem substituir a justificativa pela RH já pressuposta.

Se o MIT concluir a auditoria integral, registrar o novo estado datado e revisar as conclusões, preservando o histórico anterior; não chamar um enunciado amplo “ERH” se as cinco declarações auditadas só tratam de Dirichlet.

## Reprodução local

```bash
python -m unittest NUMERIA_ES8P255_ALPOGE_PRIORIDADES.test_p255 -v
python -O -m unittest NUMERIA_ES8P255_ALPOGE_PRIORIDADES.test_p255 -v
python -m NUMERIA_ES8P255_ALPOGE_PRIORIDADES.audit_statements --out /tmp/P255_AUDIT.json
```

Foram executados 11 testes de identificação, proveniência e adulteração, ambos os modos normal e otimizado. A auditoria não equivale a verificação integral do Lean. Os 3 fontes menores em `fixtures/` são cópias verificadas por SHA-256 de arquivos públicos da revisão original. O relatório e auditor gerados são artefatos locais do NUMERIA, não da autoria de Alpöge.