# NUMERIA E8-K V18 — Auditoria independente de testemunhas

## Escopo e resultado

A V18 confere as testemunhas PUBLICADAS pela V17, sem reexecutar toda a Fase 4 ou os pivôs da V16. O script do usuário foi executado contra os artefatos V17 originais: **12 PASS / 1 FAIL** em `--strict`, devido ao caso ressonante circular k=5 (a mesma quadratura Simpson usada como valor esperado e observado). A revisão não disfarça esse problema: preserva a identificação da comparação histórica como prova inválida, e apresenta uma identidade analítica independente para substituí-la.

Após correção: **14 PASS / 0 FAIL / 0 SKIP**, com E7 e `--strict`. **13/13 mutações do relatório identificadas**; **10/10 casos de ressonância** conferidos numericamente por um caminho independente. Duas execuções do script corrigido geraram exatamente o mesmo JSON, byte a byte. Os SHA-256 de cinco artefatos foram verificados.

## Novo lema analítico para a ressonância

Defina `L=9/20`, `w_k=kπ/(2L)` e `f_k(t)=sin(w_k(t+L))/sqrt(L)` para `t∈[-L,L]`. Com a convenção `F_k(r)=∫ f_k(t)exp(i r t)dt`, substituindo `y=t+L` obtemos

`F_k(w_k) = exp(-i w_k L) / sqrt(L) ∫_0^{2L} sin(w_k y) exp(i w_k y)dy`.

A identidade `sin(z) exp(iz) = (exp(2iz)-1)/(2i)` é exata. Como `2w_k(2L)=2kπ`, a integral da parcela oscilatória é zero e a integral restante é `iL`. Portanto, para todo inteiro `k≥1`,

`F_k(w_k)=i sqrt(L) exp(-ikπ/2)` e, exatamente, `|F_k(w_k)|²=L=9/20`.

Este é um resultado analítico elementar, não uma estimativa numérica. A V18 também compara o valor complexo obtido com integração Simpson independente para `k=1..10`. Estas dez comparações são verificações float auxiliares, sem cota de erro rigorosa para a quadratura Simpson.

O k=5 publicado no V17 permanece rotulado como **testemunha circular histórica, não utilizada como prova**. A V18 registra separadamente a substituição por uma identidade fechada; não altera o relatório anterior.

## O que os PASS significam

- Consistência de fórmulas *declaradas* e de hashes da V17 com as fontes conhecidas;
- reprodução exata da envoltória de `ln 2` e `ln 3` e do cancelamento algébrico em parâmetros racionais de `X=π²`;
- quadraturas numéricas independentes de cancelamento dos momentos e da ressonância (sanity checks);
- sensibilidade de 13 controles por alterações reais do JSON;
- integridade do E7 original e da fonte B5; hash de 3 artefatos V17.

Os controles NÃO conferem a correção global da fórmula de Guinand–Weil nem a exatidão dos enclosures analíticos da Fase 4. A V18 não recalcula as matrizes LDU. O uso de SHA-256 sem assinatura não prova autoria ou autenticidade externa.

## Pendências científicas para a V19

1. Reexecutar a Fase 4 e ligar os seus dados brutos diretamente aos enclosures analíticos, incluindo erro de quadratura e cauda de integração.
2. Obter uma segunda demonstração, com majorantes rigorosos, dos restos de quadratura/cauda e dos intervalos usados na V16.
3. Submeter o módulo Lean à compilação, separando formalização do núcleo intervalar dos teoremas analíticos ainda externos.

**Status:** `V18_V17_WITNESSES_AUDITED_RESONANCE_INDEPENDENTLY_REDERIVED_GLOBAL_CERTIFICATION_PENDING`.

A Hipótese de Riemann **não** foi demonstrada.