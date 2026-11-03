# NUMERIA–ES8-P3 — Arquimediano projetado, LDL intervalar e estabilidade V24

**Status científico: `CONDITIONAL_FINITE_PROJECTED_WEIL_POSITIVITY_CERTIFICATE_NOT_RH`.**

A ES8-P3 fecha a parcela arquimediana **para o modelo matricial real V17** (`N=6,8,10`, `L=9/20`) por uma fórmula fechada em série absolutamente convergente, com cauda cercada por Euler–Maclaurin e frações racionais. Combina-a com a parcela prima certificada na ES8-P2; anula polos por projeção; demonstra LDL positivo para a matriz **deslocada por `0,05 I`** e aplica os erros condicionais de taper V24. A igualdade com os originais binários E7/V16 **ainda não foi revalidada**. Nenhuma prova de RH.

## Arquivos

- `es8p3_arch_projected.py`: série arquimediana, caudas Hurwitz–Euler–Maclaurin, projeção, LDL intervalar, propagação V24.
- `es8p2_prime_cert.py`: cópia imutável do módulo da etapa anterior para certificar os primos e a base de projeção.
- `test_es8p3_arch_projected.py`: 16 testes, inclusive controles negativos, quadratura diagnóstica e precisão independente de `mpmath`.
- `RESULTADOS_ES8P3.json`: `M=128` com intervalos racionais completos, incluindo extremos escalonados por `10^40`.
- `RESULTADOS_ES8P3_M64.json` e `RESULTADOS_ES8P3_M256.json`: validação por refinamento.
- `RELATORIO_ES8P3.md`: demonstrações e lista de premissas não verificadas.
- `TESTES_ES8P3.txt`: saída dos testes locais.
- `MANIFESTO_SHA256.json`: hashes dos arquivos reproduzíveis.

## Executar

```bash
python -m pip install mpmath==1.3.0 numpy pytest
python -m pytest -q test_es8p3_arch_projected.py
python es8p3_arch_projected.py --M 128 --delta 1/1000000 --out RESULTADOS_ES8P3_REPLAY.json
```

O certificado em si usa exclusivamente a biblioteca padrão (`fractions.Fraction`); `mpmath` e `numpy` são importados apenas pelos testes diagnósticos.

## Gates científicos

- **PASS**: normalizações **V17** e paridade; série arquimediana com cotas de cauda; inclusão racional de `A` e `Prime`; LDL de `M_N−0,05I`; norma de taper abaixo da margem.
- **CONDITIONAL**: aplicação da fórmula explícita clássica de Guinand–Weil e dos lemas V24 de suavização à família de funções da V17.
- **NOT_ESTABLISHED**: confronto byte a byte com os arquivos históricos E7/V16, prova da positividade universal de Weil ou da Hipótese de Riemann.

Leia o relatório antes de citar ou reutilizar os resultados.
