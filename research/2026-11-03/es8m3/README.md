# NUMERIA–RIEMANN — ES8-M3: momentos até μ₅ e matrizes 3×3

**Data:** 2026-11-03. **Status:** duas matrizes gaussianas **finitas 3×3** positivas definidas por aritmética racional inclusiva, **sob a fórmula explícita clássica e as desigualdades analíticas documentadas**. **Não é prova de RH, positividade universal de Weil, nem certificação das matrizes E7/V24.**

## Arquivos

- `es8m3_certify.py`: intervalos racionais para os momentos `m_k=μ_k/400^k`, `k=0,...,5`; quadratura composta de Simpson com resto explícito de quarta derivada; soma de digama com correção de Euler–Maclaurin; derivadas do termo primo por recorrência polinomial; verificação dos três menores principais das duas matrizes.
- `es7_archimedean.py`, `es6_weil_gaussian_baseline.py`: versões-base **copiadas** da ES8-M2, sem alteração silenciosa das etapas publicadas.
- `RESULTADOS_ES8M3.json`: execução `M=8192`, `terms=96`, `T=48`, `N=64`.
- `RESULTADOS_ES8M3_REFINED.json`: execução independente de refinamento `M=16384`, mesmos demais parâmetros; os seis intervalos refinados estão contidos nos originais.
- `test_es8m3_certify.py`: testes algébricos, numéricos independentes, testes adversariais e reprodução real.
- `RELATORIO_ES8M3.md`: métodos e resultados, hipóteses, o que foi e não foi demonstrado.
- `MANIFESTO_SHA256.json`: integridade dos arquivos, exceto ele próprio.

## Reproduzir

```bash
python -m pip install mpmath==1.3.0 pytest==8.4.2
python es8m3_certify.py --M 8192 --terms 96 --T 48 --N 64 --out RESULTADOS_ES8M3.json
python es8m3_certify.py --M 16384 --terms 96 --T 48 --N 64 --out RESULTADOS_ES8M3_REFINED.json
python -m pytest -q test_es8m3_certify.py
```

**Limite:** o resultado pertence ao funcional gaussiano em `a=1/50`; nenhuma base ou projeção E7/V24 foi identificada com estas matrizes. O teste não envolve listas truncadas de zeros como premissa de certificação.
