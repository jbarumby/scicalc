# NUMERIA–RIEMANN ES8-P5: reconciliação histórica E7/V16 ↔ ES8-P3

**Data:** 2026-11-03. **Status:** `CONDITIONAL_HISTORICAL_INTERVAL_COMPATIBILITY_NOT_GLOBAL_WEIL_OR_RH`.

Esta etapa recupera, reabre e audita o **ZIP original V16** e seus dois ZIPs internos (E7 e Fase4). Reexecuta o script V16 sem alterações, checa **SHA-256 das fontes**, e compara por aritmética racional os limites históricos com a reconstrução ES8-P3 para os 30 pares de mesma paridade e as matrizes projetadas N=6,8,10.

## Reproduzir com os arquivos originais

Requisitos: Python 3.11+ e apenas biblioteca-padrão; `pytest` apenas para os testes.

```bash
python es8p5_reconcile.py \
  --v16-upload '/caminho/E8K_V16_GITHUB_COMPLETO(1).zip' \
  --extra-e7 '/caminho/E7_original (1).bin' \
  --extra-e7 '/caminho/E7(1).zip' \
  --extra-e7 '/caminho/NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3(1).zip' \
  --nested-e7 '/caminho/E7_DENTRO_DO_ZIP(1).zip' \
  --fixture-out FIXTURE_HISTORICO_ES8P5.json \
  --out RESULTADOS_ES8P5.json
python -m pytest -q test_es8p5_reconcile.py
```

A cópia local do ZIP E8K V16 completo (já fornecida pelo usuário) acompanha o **ZIP de entrega** deste estágio para facilitar a reprodução; não é adicionada ao GitHub em formato binário. Os outros quatro arquivos são duplicatas ou contêineres aninhados do mesmo E7, todos verificados.

## Reproduzir a comparação **sem o ZIP original**

```bash
python es8p5_reconcile.py --fixture FIXTURE_HISTORICO_ES8P5.json --M 64 --out RESULTADOS_ES8P5_REPLAY.json
python -m pytest -q test_es8p5_reconcile.py
```

**Distinção obrigatória:** este segundo modo compara matematicamente a implementação atual com campos racionais historicamente extraídos; **não reabre nem reautentica os ZIPs binários originais**. A autenticação integral foi executada no modo anterior e está documentada no relatório.

## Resultados

- E7 externo: 4 arquivos byte a byte iguais com SHA-256 `4acf99a2...02a7`; 1 contêiner ZIP aninhado com o mesmo E7.
- V16 íntegro: 7/7 entradas no SHA256SUMS, pacote E7: 32/32 hashes, Fase4: 30/30 hashes.
- V16 executado novamente: 12/12 gates e **JSON final idêntico byte a byte**.
- Arquimediano original com cauda R>160: **30/30** interseções de intervalos com ES8-P3.
- Primo 2: **30/30** interseções entre V16 e ES8-P2.
- Coeficientes de projeção que anulam os polos: **8/8** interseções.
- Matrizes projetadas combinadas: **67/67** interseções (N=6: 10; N=8: 21; N=10: 36).
- Sem cauda, **somente 2/30** intervalos históricos de integração em [0,160] interceptam os da integração infinita; não confundir essas quantidades.

## Limite científico

Interseção de intervalos é **teste de compatibilidade**, não prova autônoma da correção dos lemas que geraram os intervalos. Mesmo com matrizes finitas positivas sob as hipóteses da V16 e ES8-P3, a positividade universal de Weil e a Hipótese de Riemann **não são demonstradas**.

Detalhes e premissas em `RELATORIO_ES8P5.md`. O programa recusa uma alteração de hash dos arquivos originais e os testes negativos recusam troca de sinal, projeções ou intervalos.
