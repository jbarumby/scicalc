# NUMERIA–RIEMANN — ES7: certificação da integral arquimediana gaussiana

**Data:** 03/11/2026. **Status:** intervalos racionais dirigidos para a parcela arquimediana da fórmula de Guinand–Weil e sua derivada em `a`. **Não** certifica a soma de todos os zeros, a equivalência das matrizes E7/V24 nem a Hipótese de Riemann.

## Conteúdo

- `es7_archimedean.py`: intervalos da integral e sua derivada; aritmética racional com escalonamento e arredondamento externo; majorantes explícitos para truncamento da série de digama, quadratura e cauda gaussiana.
- `test_es7_archimedean.py`: doze testes, incluindo avaliações independentes por `mpmath` usadas **somente como diagnósticos não certificados**.
- `RESULTADOS_ES7.json`: execução `M=8192`, `terms=192`, `T=48`, `a=1/50` e `a=1/25`.
- `es7_v24_gate.py`: verifica a presença de âncoras textuais de normalização na documentação V24 e ES6; **não** demonstra equivalência de matrizes. Precisa receber os caminhos dos arquivos originais do repositório.
- `RELATORIO_ES7.md`: deduções, intervalos, premissas e limitações.
- `MANIFESTO_SHA256.json`: hashes dos arquivos publicados.

## Executar

```bash
python -m pip install mpmath==1.3.0 pytest
python -m pytest -q test_es7_archimedean.py
python es7_archimedean.py --M 8192 --terms 192 --T 48 --out RESULTADOS_ES7.json
python es7_v24_gate.py --v24 ../../../v24/README_V24.md --es6 ../es6/es6_weil_gaussian.py --out COTEJO_ES7_V24.json
```

O comando de cotejo pressupõe estrutura de diretórios `research/2026-11-03/es6` e `v24/` da raiz do repositório; na raiz, prefira passar os caminhos absolutos ou os relativos `v24/README_V24.md`, `research/2026-11-03/es6/es6_weil_gaussian.py`.

**Critério de publicação:** registrar como `CERTIFIED_ARCH_GAUSSIAN_ONLY`; `NOT_ESTABLISHED` quanto à identificação das matrizes e da parte espectral global.
