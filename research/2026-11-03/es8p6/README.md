# NUMERIA–RIEMANN — ES8-P6

Auditoria das **cotas analíticas** da Fase 4 e caudas E7, com `fractions.Fraction` e Python padrão. Os 30 intervalos históricos de 1.600 células cada foram autenticados; a nova demonstração cobre o majorante holomorfo, a interpolação de grau 15 e os dois envelopes de cauda (conservador e forte). **Não é prova de RH nem verificação integral de todos os valores nodais.**

## Execução completa, com os ZIPs recuperados do usuário

```bash
python es8p6_independent_audit.py \
  --phase4 /caminho/NUMERIA_FASE4_30_ENTRADAS_AUDITORIA_PROJECAO_2026-10-24.zip \
  --e7 /caminho/NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip \
  --out FIXTURE_FASE4_ES8P6.json
python es8p6_node_spotcheck.py
python -m pip install pytest
python -m pytest -q test_es8p6_independent.py
python check_manifest.py
```

O ZIP de Fase 4 estava dentro do arquivo original `E8K_V16_GITHUB_COMPLETO(1).zip`. O artefato de transferência para GitHub inclui apenas as **fixtures autenticadas**, os dois scripts nodais originais em `source_snapshot` com SHA-256 e o novo auditor; não publica os ZIPs binários históricos.

## Conferência em GitHub Actions

O workflow executa o replay **offline** dos 22 testes, verifica a fixture e o manifesto. A verificação dos 30 arquivos originais em ZIP depende da execução local do comando completo acima. A `main` não deve ser alterada.

## Conclusões condicionais

Majorante e erro de interpolação rederivados, cauda E7 (Proposições 5 e 6) rederivada, dois limites históricos de cauda conferidos, nove interseções de intervalos nodais verificadas. O código não reavalia todos os nós e não comprova o critério de Weil para todas as funções admissíveis.
