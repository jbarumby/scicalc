# NUMERIA–RIEMANN — ES8-P9: integração intervalar independente das 30 entradas

**Data da pesquisa:** 03/11/2026. **Situação:** soma nodal racional auditada integralmente e três testes matriciais finitos condicionais. **Não é prova da Hipótese de Riemann**.

## Objetivo

Continuação da ES8-P8, que certificou compatibilidade ponto a ponto para `768.000` intervalos. Esta fase de fato **integra** os nós com pesos de Newton–Cotes interpolatórios, alarga os intervalos pelo erro de Cauchy e a cauda de `r >= 160`, e reconstrói a matriz projetada em `N=6,8,10`.

## Resultado

- 1.600 células × 16 nós × 30 entradas = **768.000 parcelas** somadas de modo independente, em **32 lotes**.
- **16 pesos racionais** deduzidos por duas implementações algébricas independentes; momentos exatos de grau 0 a 15 verificados; pesos negativos tratados com reversão de limites.
- **30/30** intervalos de integrais nodais independentes intersectam os históricos.
- **30/30** intervalos de integrais incluindo o erro de interpolação intersectam os históricos.
- Erro de interpolação comum, exato: `4154891126177792/209922127267236328125`.
- A matriz projetada (arquimediano + primo 2, com polos eliminados) tem decomposição racional LDL positiva após subtração de `0,05 I` e de uma cota de perturbação, para **N=6,8,10**, **condicionada** às cotas analíticas e convenções explicitadas.

## Reprodução

Na pasta da fase, com Python 3.11+:

```bash
python -m pip install pytest
python check_manifest.py
python -m pytest -q test_es8p9.py
python es8p9_integrate.py --start 0 --end 1600 --batch 50 --terms 128 --digits 35 --out RESULTADOS_ES8P9_REPLAY.json
python validate_replay.py RESULTADOS_ES8P9_INTEGRAIS.json RESULTADOS_ES8P9_REPLAY.json
python es8p9_projected.py --input RESULTADOS_ES8P9_REPLAY.json --out RESULTADOS_ES8P9_PROJECAO_REPLAY.json
```

**Observação:** a execução `--start` / `--end` parcial não obtém status de conclusão da integral completa. A auditoria preserva por SHA-256 os arquivos históricos e os resumos integrais; o ZIP V16 original, enviado pelo usuário, foi aberto e validado localmente antes da exportação deste pacote.

## Arquivos

- `es8p9_integrate.py`: motor de soma nodal 16×1600×30, pesos assinados, intervalos, hashes de lote, cotejo histórico.
- `es8p9_projected.py`: propaga as entradas globais e primo 2 pelas projeções, testa LDL racional.
- `es8p8_full.py`, `es8p7_independent.py`: avaliador intervalar separado, preservado sem alterações da fase anterior.
- `es8p6_independent_audit.py`: cotas racionais de Cauchy e cauda, preservadas.
- `AUDITAR_FASE3_PROJECAO.py`: rotina de primo/projeção da fonte histórica; não é independente dessa parte aritmética.
- `HISTORICAL_PHASE4_REFERENCE.json`: 30 resultados e hashes extraídos dos ZIPs históricos autenticados.
- `RESULTADOS_ES8P9_INTEGRAIS.json`: 30 integrais, 32 hashes de lotes e comparações.
- `RESULTADOS_ES8P9_PROJECAO.json`: testes condicionais das três matrizes projetadas.
- `test_es8p9.py`, `validate_replay.py`: 17 testes e replay dos 32 lotes/30 entradas.
- `RELATORIO_ES8P9.md`: provas usadas, condições e limitações.
- `MANIFESTO_SHA256.json`, `check_manifest.py`: integridade dos arquivos de pesquisa.

## Limite lógico

A auditoria prova a consistência dos cálculos **dadas** as inclusões nodais, os lemas de erro e a fórmula explícita especializada. O algoritmo não prova da estaca zero a Guinand–Weil, não formaliza em Lean as identidades de digama, nem demonstra a positividade de Weil para todas as funções-teste ou a Hipótese de Riemann.
