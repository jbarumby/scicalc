# NUMERIA NL1/ES1 — pacote executável

**2026-11-03.** Newton 2D para três zeros conhecidos da zeta (numérico, NÃO certificado), Krawczyk exato para um modelo espectral 2×2 (certificado apenas neste modelo), análise de cotas V24 não projetadas, 44 referências classificadas e testes executáveis.

- `RELATORIO_NL1_ES1.md`: resultados, método, fronteiras matemáticas e próximos gates.
- `numeria_nl1.py`: solver, verificador racional, comparador de cotas V24, saída JSON.
- `test_numeria_nl1.py`: testes unitários.
- `RESULTADO_NL1_ES1.json`: execução já realizada.
- `EXTRATO_V24_COTAS.json`: extrato de quatro números da auditoria V24, não o original completo.
- `classificar_referencias.py` e `REFERENCIAS_CLASSIFICADAS.csv`: classificação e deduplicação bibliográfica.
- `MANIFESTO_SHA256.json`: integridade e procedência dos arquivos.

```bash
python -m pip install mpmath==1.3.0 pytest
python numeria_nl1.py --out RESULTADO_NL1_ES1.json --v24 EXTRATO_V24_COTAS.json
python -m pytest -q test_numeria_nl1.py
```

**Cuidado central:** não se apresentou prova da Hipótese de Riemann; a certificação Krawczyk não foi aplicada à zeta real ou à forma de Weil.
