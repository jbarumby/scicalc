# NUMERIA–RIEMANN / ES8-P10

Auditoria matemática independente da cauda da digama e melhoria racional de **seis milhões de vezes** no limite de erro da interpolação de Cauchy da Fase 4.

## Reproduzir

```bash
python -m pip install pytest
python -m pytest -q test_es8p10.py
python es8p10_lemmas.py --out RESULTADOS_ES8P10.json
python check_manifest.py
```

Os arquivos de código da ES8-P6/P9 são cópias da versão de entrada para execução offline e não devem ser atribuídos à ES8-P10. `RESULTADOS_ES8P9_INTEGRAIS.json` é um conjunto de **entradas autenticadas por SHA-256**, não um resultado recalculado nesta etapa. O cálculo completo dos 768 mil valores foi executado na ES8-P8/ES8-P9, por outro workflow.

**Escopo:** comprovação das desigualdades de resto da digama (via Laplace+coth) e da cota integrada de interpolação para as funções-teste, com transporte dos novos erros a 3 matrizes finitas. **Não** formaliza todas as identidades em Lean nem implica RH.

Ver `RELATORIO_ES8P10.md` para as demonstrações, entradas e limitações.
