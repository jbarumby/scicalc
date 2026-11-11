# ES8-P151 — pacote reproduzível

## Conteúdo

- `TEOREMA_ES8P151.md`: definições, provas de erro `C¹`, correção exata de momentos em escala única, dimensão e transformação Toeplitz.
- `RELATORIO_ES8P151.md`: síntese científica, contagens e limites.
- `p151_dense_multiscale.py`: índices, contagens e cotas racionais; geração de JSON.
- `p151_toeplitz_exact.py`: congruência Toeplitz em frações racionais.
- `test_p151.py`: testes automatizados com biblioteca padrão.
- `P151_AUDIT.json`: contagens, cancelamentos e exemplos de orçamentos de erro.
- `*_EM_TXT.txt`: cópias idênticas dos programas `.py` em formato TXT para facilitar downloads.
- `MANIFEST_SHA256.json`: hashes dos arquivos do pacote, exceto ele próprio.

## Reprodução (Python 3.11+)

```bash
python -m unittest -v test_p151.py
python p151_dense_multiscale.py --N 8 --out P151_AUDIT.json
python p151_toeplitz_exact.py
```

As computações utilizam somente a biblioteca padrão. Os testes racionais não provam os lemas analíticos por si sós; suas demonstrações constam de `TEOREMA_ES8P151.md`.

**Não contém cálculo da Gram real de Weil nem prova da Hipótese de Riemann.**
