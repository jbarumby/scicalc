# NUMERIA NB1 — Estudo bibliográfico + certificado finito

1. `python nb1_gram_audit.py` — gera `RESULTADOS_NB1.json` com escolhas racionais em N=6,8,10 e limites de erro com frações exatas.
2. `python test_nb1_gram_audit.py` — cinco testes de consistência.
3. `RELATORIO_PESQUISA_NB1.md` — contexto bibliográfico, derivação e limites das conclusões.

Python + NumPy para a busca heurística inicial, mpmath para uma checagem de teste. A verificação de coeficientes fixos usa `fractions.Fraction` e limites racionais de logaritmo.

Não demonstra RH, não prova otimalidade e não integra diretamente o repositório GitHub.
