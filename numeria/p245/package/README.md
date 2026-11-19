# ES8-P245 Gate B — certificado GL60 por aritmética intervalar exata

Este pacote certifica **48 entradas selecionadas de somas GL60**, nos painéis 0, 503 e 1007, dois setores e oito pares de índices (0, 1, 16 e 376). **Não** certifica o erro analítico de quadratura ou as 142.506 entradas da matriz inteira, nem RH.

- `p245b_certified_gl60.py` — código independente de Arb, MPFR, SciPy e mpmath na execução com `--roots`. Usa inteiros Python, Fractions, intervalos dirigidos e restos analíticos.
- `P245B_ROOT_BRACKETS.json` — 60 propostas racionais de raízes GL60, verificadas novamente por sinais de Legendre exatos.
- `P245B_CERT.json` — 48 caixas GL60 e 18 limites adicionais de magnitudes de graus extremos, vinculados aos hashes das seis fontes P236.
- `pilot_sources/` — seis amostras Arb P236 originais, copiadas sem modificação do pacote P245.
- `test_p245b_certified_gl60.py` — 12 testes, inclusive adulterações adversariais.
- `RELATORIO_P245B.md` — prova das cotas de erro e delimitação de escopo.
- `SHA256SUMS_P245B.txt` — fingerprints dos arquivos, gerado depois do certificado.

**Comando:**

```bash
python3 p245b_certified_gl60.py --pilot-dir pilot_sources --roots P245B_ROOT_BRACKETS.json --out /tmp/P245B_REPLAY.json
python3 -m unittest -v test_p245b_certified_gl60
python3 -O -m unittest -q test_p245b_certified_gl60
```

A geração do certificado neste Gate B não utiliza valores numéricos aproximados para validar inclusões. Os inteiros das raízes foram sugeridos inicialmente por mpmath e posteriormente validados por avaliação exata de `P60`; o replay a partir de `--roots` é totalmente aritmético/racional.
