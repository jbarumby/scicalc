# NUMERIA ES8-P — auditoria de projeção e artigos (2026-11-03)

**Status:** projeção algébrica exata de polos; envelope racional projetado **condicional ao modelo V24**; diagnóstico numérico de primos; verificação de mismatch de paridade/semirreta V17–V24; **SEM PROVA RH e SEM POSITIVIDADE DO WEIL UNIVERSAL**.

## Executar

```bash
python -m pip install mpmath==1.3.0 sympy==1.14.0 pytest==8.4.2
python -m pytest -q test_es8p_projection.py
python es8p_projection.py --out RESULTADOS_ES8P.json
python triage_references.py --source REFERENCIAS_ORIGINAIS_9.csv --out REFERENCIAS_TRIADAS_ES8P.csv
```

A biblioteca `python-flint` não estava instalada no ambiente local; nenhuma integração de matriz com Arb foi declarada como certificada. O pipeline emprega inteiros e `fractions.Fraction` para as cotas do modelo e `sympy` para identidades exatas. A integração com `mpmath` é apenas diagnóstico.

## Arquivos

- `es8p_projection.py`: projeção exata e limites condicionais, fórmula fechada da convolução, teste de paridade
- `test_es8p_projection.py`: 18 testes de regressão
- `RESULTADOS_ES8P.json`: resultados reais da execução
- `RELATORIO_ES8P.md`: fórmulas, provas, limites, artigos e pendências
- `REFERENCIAS_ORIGINAIS_9.csv`: tabela original enviada, preservada
- `REFERENCIAS_TRIADAS_ES8P.csv`: 11 obras únicas, prioridade e ressalvas
- `triage_references.py`: algoritmo reproduzível de classificação do CSV
- `MANIFESTO_SHA256.json`: hashes e limitações

**Bloqueio principal:** incompatibilidade de paridade/parte real entre os integrandos arquimedianos expressos em V17 e V24. Antes de atualizar as matrizes originais é preciso verificar e corrigir o operador bilinear efetivamente usado em cada versão.
