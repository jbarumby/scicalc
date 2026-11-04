# NUMERIA–RIEMANN ES8-P7 — nós e admissibilidade de funções-teste

Esta fase **não** demonstra RH nem a positividade universal de Weil. Ela realiza uma auditoria independente **amostral** de 450 entradas nodais da Fase 4 e fornece uma dedução analítica de admissibilidade das funções `H_ij` usadas na fórmula explícita clássica.

## Reproduzir

Requer Python 3.11 ou mais recente, `mpmath==1.3.0` e `pytest`:

```bash
python -m pip install mpmath==1.3.0 pytest
python es8p7_independent.py --out RESULTADOS_ES8P7.json
python -m pytest -q test_es8p7_independent.py
python check_manifest.py
```

O cálculo usa `fractions.Fraction` com somas positivas da digama e limites integrais, raízes por `isqrt`, reduções por `π` cercado pela fórmula de Machin, Taylor com restos rigorosos, intervalos externos sobre uma malha racional fixa. `mpmath` fornece apenas um *sanity check* independente; não produz certificados.

Os dois módulos `fase4_*_sem_arb.py` reproduzem **bytes idênticos** do ZIP histórico auditado ES8-P6; seus hashes são fixados nos testes. O módulo `es8p7_independent.py` executa a nova implementação e a confronta com a antiga. Os resultados no CSV são apenas metadados bibliográficos fornecidos pelo usuário; a quarta fonte suplementar (Nemes) vem de pesquisa externa.

## Gates

- 30 pares da mesma paridade por 15 nós de teste distribuídos em cinco células: 450 interseções intervalares.
- 15 checagens mpmath de alta precisão para `Ω(r)`; apenas diagnóstico.
- Admissibilidade de `H_ij`: inteiro, par, decaimento `O(|r|^-4)` em faixas horizontais, sem pressupor RH.
- **Não executado:** 30 × 1600 × 16 = 768.000 avaliações originais. Não presumir que as outras avaliações estejam corretas apenas porque a amostra passou.

Ver relatório técnico `RELATORIO_ES8P7.md` e resultados exatos `RESULTADOS_ES8P7.json`.
