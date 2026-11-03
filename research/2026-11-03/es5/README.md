# NUMERIA–ES5 — Derivadas logarítmicas e traço aritmético (03/11/2026)

**Escopo:** porta de entrada **incondicional** para a aritmética da função zeta, usando identidades conhecidas e um experimento numérico reproduzível. **Não é uma fórmula global de Weil, não certifica zeros nem prova a Hipótese de Riemann.**

## Componentes

- `es5_arithmetic_trace.py`: enclosures exatos (inteiros e frações) de duas representações independentes de `-ζ'(σ)/ζ(σ)` para σ inteiro ≥4; prova computacional da janela `log 2 < 0.9 < log 3`; traço gaussiano numérico; verificação da derivada logarítmica da zeta completada.
- `test_es5_arithmetic_trace.py`: 12 testes de regressão, inclusive sinais, fatoração, intervalos, refinamento, identidades e entradas inválidas.
- `RESULTADOS_ES5.json`: saída reproduzida com `N=4096`, `σ=4`, `dps=45`.
- `RELATORIO_ES5.md`: fórmulas, protocolos, comparação com ES4/V24, certificados, limitações.
- `TESTES_ES5.txt`: execução local dos testes.
- `MANIFESTO_SHA256.json`: hashes e contexto de reprodução.

## Executar

```bash
python -m pip install mpmath==1.3.0 pytest==8.4.2
python -m pytest -q test_es5_arithmetic_trace.py
python es5_arithmetic_trace.py --N 4096 --sigma 4 --dps 45 --out RESULTADOS_ES5.json
```

**Limitações:** a aritmética racional certifica os dois intervalos de log-derivada e a janela de suporte. A identidade gaussiana é analiticamente válida em `Re(s)>1`, mas suas integrais e somas finitas foram avaliadas com `mpmath`, não com Arb/FLINT. As cotas de cauda são analíticas; o erro numérico de quadratura não está coberto por elas. Nenhuma identificação do operador ES4 com a zeta foi construída.
