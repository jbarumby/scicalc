# NUMERIA–ES3 — Auditoria Spin–Zeta (03/11/2026)

**Status:** objeções matemáticas verificadas + experimentos numéricos de falsificação de modelos; **nenhuma prova da Hipótese de Riemann**. Não é um cálculo adélico nem uma certificação da localização dos zeros da zeta.

## Conteúdo

- `es3_spin_zeta.py`: verificação simbólica exata dos 8 estados locais, modelo tensorial finito de `m=1..4` lugares, diagonalização do Hamiltoniano quadrático e comparação com ordenadas conhecidas dos zeros da zeta, somas parciais do critério de Li e contagem de Riemann–von Mangoldt.
- `test_es3_spin_zeta.py`: 8 testes automatizados independentes do ajuste aos zeros.
- `RESULTADOS_ES3.json`: saída efetivamente produzida com `mpmath` a 65 dígitos para os 20 primeiros zeros conhecidos da reta crítica.
- `RELATORIO_ES3.md`: hipóteses, deduções e limitações.
- `Qwen_latex_20261103_m86z47xnk.txt`: texto de origem enviado pelo usuário, preservado com indicação de que é especulativo e contém erros.
- `MANIFESTO_SHA256.json`: hashes dos arquivos de origem, código, resultados e testes.

## Reproduzir

Requer Python 3.11+ e as bibliotecas `numpy`, `mpmath`, `sympy`, `pytest`.

```bash
python -m pip install numpy mpmath sympy pytest
python es3_spin_zeta.py --out RESULTADOS_ES3.json --zeros 20 --train 8 --holdout 6
python -m pytest -q test_es3_spin_zeta.py
```

### O que está demonstrado

Oito verificações **algébricas exatas** (SymPy) das matrizes de spin; para a soma tensorial de `m` operadores locais, o espectro é `{-m,-m+2,...,m}` com multiplicidades binomiais; no regime `beta>2*abs(gamma)` o operador quadrático estável é unitariamente equivalente a um oscilador deslocado, com níveis `E_n = E_0+n*sqrt(beta**2-4*gamma**2)`. A sequência de zeros da zeta não obedece a essa contagem linear (fórmula de Riemann–von Mangoldt).

### O que NÃO está demonstrado

Não construímos nenhum operador em um espaço adélico, não houve identificação espectral com a zeta, os zeros da zeta foram usados apenas como **dados de referência numéricos**; somas parciais dos zeros críticos têm sinal não negativo automaticamente e **não certificam** os coeficientes completos de Li. A comparação numérica tampouco constitui verificação rigorosa da fórmula global de Weil.

Esta pesquisa deve permanecer em uma branch separada até revisão externa. Não substituir os resultados e as hipóteses das fases E7/V24 por esses modelos ilustrativos.
