# NUMERIA ES8-P15 — Teoremas gerais e auditoria de uniformidade

Esta etapa responde a três obrigações: decaimento controlado da cauda E7, passagem de somas finitas a séries infinitas, e critério de margem uniformemente positiva.

**Diferencie resultados:** o Lean pode demonstrar o limite do majorante e uma implicação para qualquer função limitada por ele. Não está provada aqui, no Lean, a estimativa da integral arquimediana verdadeira. Igualmente, o teorema `∀ n` exige hipóteses `∀ n` não fornecidas pelas três matrizes históricas.

Arquivos:
- `ES8P15.lean`: lemas e teoremas condicionais sem axiomas adicionais.
- `audit_es8p15.py`, `test_es8p15.py`: reprodução Python, contraexemplos, orçamento uniforme da cauda bruta, e testes.
- `RESULTADOS_ES8P15.json`: valores de diagnóstico e restrições explícitas.
- `RELATORIO_ES8P15.md`: demonstrações matemáticas e pendências.
- `numeria-es8p15-lean.yml`: compilação Lean e regressão Python.

Comandos:

```bash
cd research/2026-11-04/es8p15
python -m unittest -v test_es8p15.py
python audit_es8p15.py
```

No projeto Lean, `ES8P14.lean` precisa ser compilado para importar `ES8P15.lean`.