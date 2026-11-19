# P246 — auditoria dos novos certificados P237/P245D

**Reprodução (Python 3, sem dependências externas):**

```bash
python3 p246_verify_uploaded.py --data sources --out /tmp/P246_CHECK.json
python3 -m unittest -v test_p246_verify_uploaded
python3 -O -m unittest -q test_p246_verify_uploaded
```

A pasta `sources/` contém os nove arquivos originais enviados, preservados byte a byte. O relatório distingue a **verificação exata dos intervalos arquivados** da **reprodução da quadratura**, que exige os seis chunks-fonte não incluídos. Não houve publicação no GitHub.
