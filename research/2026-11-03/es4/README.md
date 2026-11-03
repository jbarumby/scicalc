# NUMERIA–ES4 — Contagem espectral não quadrática

Experimento de controle para a fase ES3/ES4: operador diagonal explícito de autoadjunção elementar e lei de contagem semelhante à de Riemann–von Mangoldt, **sem qualquer alegação de equivalência espectral com a função zeta**.

## Executar

```bash
python -m pip install mpmath==1.3.0 pytest==8.4.2
python -m pytest -q test_es4_spectral_count.py
python es4_spectral_count.py --zeros 25 --dps 60 --out RESULTADOS_ES4.json
```

## Arquivos

- `es4_spectral_count.py`: definição e inversão via Lambert W do termo suave de contagem, funções para traços finitos e somas tipo Li (adversariais).
- `test_es4_spectral_count.py`: 12 testes de aritmética numérica, identidade inversa, monotonicidade, contagem, modelo deslocado, erros e entradas inválidas.
- `RESULTADOS_ES4.json`: execução com 25 zeros conhecidos de referência via mpmath; **não certificado por intervalos**.
- `RELATORIO_ES4.md`: argumento matemático da construção e as lacunas da fórmula explícita de Weil.
- `MANIFESTO_SHA256.json`: hashes para os arquivos do pacote.

As `12` verificações testam o código, mas não substituem as provas analíticas expostas no relatório. A `main` do projeto é preservada; a nova frente deve permanecer em branch de pesquisa.
