# NUMERIA ES8-M2 — Certificação de momentos μ₂ e μ₃

**Data:** 2026-11-03. **Status:** `CERTIFIED_FINITE_GAUSSIAN_MOMENTS_H1_AND_STIELTJES_H1`, condicionado às identidades analíticas clássicas explicitadas em `RELATORIO_ES8M2.md`. **Não é prova da Hipótese de Riemann nem da positividade universal de Weil.**

## O que foi feito

- Certificados para μ₀…μ₃ na fórmula global de Guinand–Weil, com a=1/50 e momentos normalizados por 400ᵏ;
- Limites arquimedianos por quadratura racional de ponto médio, cauda gaussiana e digama, incluindo pesos t⁴, t⁶;
- Soma de potências de primos com derivadas k=0…3, controle rigoroso de todo n>64;
- Certificação da positividade definida das matrizes 2×2 Hankel e Stieltjes deslocada, **somente** para essa família e parâmetro;
- 14 testes unitários, comparação com `mpmath` não certificada, mais refinamento M=4096/8192/16384;
- Triagem de 22 linhas bibliográficas em dois CSVs, reduzidas a 19 registros distintos.

## Reproduzir

Python 3.11+; `mpmath==1.3.0` (diagnóstico) e `pytest` (testes). Para a parte racional certificada, apenas a biblioteca padrão do Python é requerida.

```bash
python -m pip install mpmath==1.3.0 pytest==8.4.2
python -m pytest -q test_es8m2_certify.py
python es8m2_certify.py --M 16384 --T 48 --terms 192 --N 64 --out RESULTADOS_ES8M2_REPLAY.json --refs-out REFERENCIAS_TRIADAS_REPLAY.csv
```

O `GitHub Actions` compara os **extremos racionais exatos** reproduzidos com os registrados; a coluna de diagnóstico de ponto flutuante não participa dessa igualdade.

## Material

- `es8m2_certify.py` — cálculo racional, identidades, quadratura, derivados do primo e matrizes 2×2;
- `es7_archimedean.py`, `es6_weil_gaussian_baseline.py` — dependências versionadas da certificação;
- `test_es8m2_certify.py` — testes de polinômios, derivadas, sinais, determinantes, caudas, dados e parâmetros;
- `RESULTADOS_ES8M2.json` — certificado com extremos exatos em frações;
- `RELATORIO_ES8M2.md` — hipóteses, demonstração das cotas, conclusão e bibliografia;
- `REFERENCIAS_TRIADAS_ES8M2.csv` e `sources/` — arquivo bibliográfico deduplicado e originais;
- `MANIFESTO_SHA256.json` — hashes de integridade.

**Limitação expressa:** sem certificação de todos os zeros, sem identificação das matrizes projetadas E7/V24, sem generalização para todas as funções-teste. A positividade universal de Weil continua em aberto.
