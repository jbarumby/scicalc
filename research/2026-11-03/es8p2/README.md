# NUMERIA–ES8-P2 — projeção real de Weil e termo primo intervalar

Data 03/11/2026. **Resultado finito parcial**: intervalos racionais de `B^T Prime B` para modos seno compactos `N=6,8,10` e anulação exata de paridades opostas. Não é certificado de `B^T(G+P)B`, de RH, ou da positividade universal de Weil.

Arquivos:
- `es8p2_prime_cert.py`: cálculo com frações exatas, Taylor e limites dirigidos de `π`, `log 2`, `sqrt 2`; correção da convenção de semirreta da V24.
- `test_es8p2_prime_cert.py`: 18 testes, com checagens independentes de `mpmath` classificadas como diagnóstico, não certificados.
- `RESULTADOS_ES8P2.json`: matrizes com extremos inteiros escalonados exatos (`scale = 10**40`) e resumos decimais.
- `RELATORIO_ES8P2.md`: demonstrações, limites, referências bibliográficas e os gates ainda não satisfeitos.
- `References_LeapSpace (10|11|12).csv`: registros originais recebidos, imutados.
- `triage_es8p2.py` e `REFERENCIAS_TRIADAS_ES8P2.csv`: desduplicação dos 51 registros em 44 DOIs e avaliação de relevância.
- `RESUMO_REFERENCIAS_ES8P2.json`: contagens.
- `MANIFESTO_SHA256.json`: integridade dos arquivos acima.

Reprodução (Python 3.11+):

```sh
python -m pip install pytest==8.4.2 mpmath==1.3.0
python -m pytest -q test_es8p2_prime_cert.py
python es8p2_prime_cert.py --out RESULTADOS_ES8P2.json
python triage_es8p2.py --root .
```

**Importante**: os 18 testes não provam a validade da fórmula explícita de Weil desde o início, nem demonstram equivalência exata com os códigos-fonte integrais E7/V16. Para isso é necessário auditar a matriz completa, inclusive sua parcela arquimediana, no mesmo referencial.
