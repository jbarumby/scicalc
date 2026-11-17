# NUMERIA — ES8-P227 — auditoria Arb cruzada

Este ramo contém uma **segunda avaliação** de 200 modos de Zhu para L=4/5, diferente da P226:

- python-flint/Arb **0.8.0** com **1024 bits**, contra 0.9.0 e 800 bits na P226;
- Bessel via `J_(n+1/2)` nativa, digama nativa, contra âncoras por série e Binet–Stirling;
- 14 lotes reais de 50 painéis GL60 (8 pares, 6 ímpares), 700 painéis, 42.000 nós e 281.400 intervalos triangulares;
- em auditoria externa aos runners, **281.400/281.400** pares de intervalos das P226/P227 se intersectaram com aritmética inteira exata;
- depois de recompor a matriz e adicionar uma vez o erro de Chebyshev–Bernstein, a conferência racional das 40.200 entradas de R ficou dentro de `10^-46` em relação às matrizes P193 e à testemunha P196;
- a ligação P207 permanece **condicionada** à demonstração dos lemas analíticos e da identidade da redução de Zhu.

## Reprodutibilidade

- Painéis iniciais e finais em Arb 0.8: `numeria/p227_replay_arb08.py`.
- Lotes completos: `numeria/p227_chunk_arb08.py`.
- CI: `.github/workflows/numeria-p227-arb08.yml`, `.github/workflows/numeria-p227-first-chunk.yml`, `.github/workflows/numeria-p227-remaining-13.yml`.
- Execuções: https://github.com/jbarumby/scicalc/actions/runs/43996361257, https://github.com/jbarumby/scicalc/actions/runs/43996705915, https://github.com/jbarumby/scicalc/actions/runs/43996977950.

Os comparadores inteiros, auditoria analítica e artefatos de execução estão documentados no pacote externo `NUMERIA_ES8P227.zip`. Eles não utilizam valores P193 para construir as novas quadraturas.

**Limite:** v0.8.0 e v0.9.0 compartilham a família Arb/FLINT; esta não é uma auditoria com biblioteca matemática inteiramente distinta. Nenhuma conclusão prova a Hipótese de Riemann ou positividade de Weil para janelas arbitrárias. **PR em rascunho; não fazer merge sem revisão matemática independente.**
