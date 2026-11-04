# NUMERIA–RIEMANN ES8-P19

Esta fase é uma revisão **crítica** do arquivo Qwen recebido e uma formalização da **correção de dois polos**, não uma prova da Hipótese de Riemann.

## Entregáveis

- `QWEN_ORIGINAL_2026-11-04.tex` — fonte original preservada.
- `PARECER_QWEN_ES8P19.md` — erros lógicos e conclusão aproveitável.
- `ES8P19.lean` — identidade algébrica para todo `n` e transferência condicional a `c₀₀`, sem `fin_cases` na dimensão.
- `audit_es8p19.py` e `test_es8p19.py` — aritmética racional exata e testes adversariais.
- `RESULTADOS_ES8P19.json` — saída reproduzível.
- `RELATORIO_ES8P19.md` — demonstração e pendências científicas.
- `MANIFESTO_SHA256.json` — integridade das fontes e certificados.

## Executar localmente

```bash
python -m unittest -v test_es8p19.py
python audit_es8p19.py --out RESULTADOS_ES8P19_REPLAY.json
cmp RESULTADOS_ES8P19_REPLAY.json RESULTADOS_ES8P19.json
```

## Status

Python: 11/11 testes aprovados. Lean: fonte preparada para compilação em GitHub Actions; **não declarar aprovado antes da execução confirmar**. A forma infinita completa de Weil, sua semilimitação e o critério universal continuam abertos.
