# NUMERIA–ES8-P4 — Cadeia de proveniência E7/V16 e gate de equivalência

**Escopo:** consolidar o controle de origem Git do ZIP E7, reexecutar as identidades algébricas de polos, a paridade e os certificados finitos da ES8-P3; classificar duas listas de referências (`14` e `15`).

## Reproduzir

```bash
python -m pip install pytest==8.4.2
python -m pytest -q test_es8p4_reconcile.py
python triage_es8p4.py
python es8p4_reconcile.py --M 64 --out RESULTADOS_ES8P4.json
# Opcional, quando ZIP E7 íntegro estiver acessível:
python es8p4_reconcile.py --e7 NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip
# A modalidade estrita REJEITA quando não foi provada a igualdade interna:
python es8p4_reconcile.py --strict   # saída 2: NOT ESTABLISHED
```

Os arquivos em `vendor/` são cópias não modificadas dos códigos ES8-P2/ES8-P3 para tornar o pacote autossuficiente. A cópia de ES8-P2 tem SHA-256 fixado e conferido. Os arquivos de origem não são editados.

**Resultado:** mesmo objeto Git SHA-1 do E7 em dois caminhos e dois commits; reexecução válida da matriz-modelo V17; `HISTORICAL_E7_V16_EQUIVALENCE=NOT_ESTABLISHED`. Identidade de objeto Git significa mesmos bytes armazenados pelo Git, **não** prova das equações internas nem comparação entrada a entrada com o E7.

Conferência crítica: a margem `1/20` do LDL é relativa às **coordenadas** da base projetada `B`. Como `B` não é ortonormal, ela não deve ser apresentada automaticamente como constante de coerção na norma física original; é uma cota na norma euclidiana dos coeficientes da base projetada.

Veja `RELATORIO_ES8P4.md` para detalhes, referências, limites e próximos passos. Não é prova de RH ou da positividade universal de Weil.
