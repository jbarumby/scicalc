# NUMERIA ES8-P256 — Ponte entre somas de Mangoldt e o símbolo S_L

Esta etapa é uma continuação da P255. Dois objetos separados: (a) inspeção de **dependências declaradas** em módulos Lean de Alpöge; (b) **identidade matemática independente**, sem RH, entre a soma ponderada de Mangoldt e a integral de uma família de símbolos do NUMERIA.

A identidade nova é `W_L(t)=exp(L)∫_0^L cosh(L-u) S_u(t)du`; a prova está em [PROVA_P256.md](PROVA_P256.md). O efeito da nova potência prima `q=7` sobre a soma ponderada é `((X-7)/7)log7 cos(tlog7)` para `7<X<8`, enquanto no símbolo duro é um salto de `2log7/sqrt7` em `t=0`. Não se obteve margem de Weil ou RH.

## Reproduzir sem bibliotecas externas

```bash
python -m unittest NUMERIA_ES8P256.test_weighted_bridge -v
python -O -m unittest NUMERIA_ES8P256.test_weighted_bridge -v
python -m NUMERIA_ES8P256.weighted_bridge --out /tmp/P256_CERT.json
```

Para validar também as declarações dos 9 módulos Lean de origem, ponha os arquivos `src/...` no diretório informado por `--source-dir` (revisão `Deicyde/riemann-hypothesis@40e56bab`). O auditor valida o hash Git blob de cada um e a presença das hipóteses explícitas, mas **não invoca o Lean**.

```bash
python -m NUMERIA_ES8P256.weighted_bridge --source-dir /caminho/riemann-hypothesis/lean --out /tmp/P256_CHAIN.json
```

O commit e os blobs da origem estão fixados no script. O grafo textual **não é a árvore completa de dependências nem verifica os termos de prova**. Os resultados de kernel do MIT permanecem separados. Preserve P249–P255 e mantenha o PR em rascunho.