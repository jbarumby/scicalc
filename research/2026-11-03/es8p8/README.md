# NUMERIA–RIEMANN — ES8-P8: auditoria exaustiva dos nós da Fase 4

**Data:** 03/11/2026. **Status:** 768.000 interseções intervalares verificadas; não é uma demonstração de RH ou da positividade universal de Weil.

## Objeto da auditoria

Malha racional: 1.600 células de largura `1/10`, 16 nós uniformes por célula, 30 pares `1 <= i <= j <= 10` de mesma paridade: `1600 × 16 × 30 = 768000` avaliações do integrando arquimediano real `q_ij(r)`. A fase ES8-P7 certificara amostralmente 450 interseções.

`es8p8_full.py` compara integralmente os intervalos produzidos por dois procedimentos com arredondamento **direcionado**:

- **Histórico:** `fase4_fixedpoint_sem_arb.py`, `10^-40`, 64 parcelas da série de digama e resto Euler–Maclaurin; arquivo original, fixado por SHA-256.
- **Novo:** `IndependentFixed`, `10^-35`, 128 parcelas da série de digama, cauda via integral e derivada; seno por Taylor com resto explícito, redução do argumento, sinc por Taylor / seno dividido e constante `Ω(0)` calculada por um terceiro módulo da ES8-P7, separado do avaliador histórico.

Cada uma das 32 partes de 50 células armazena um hash SHA-256 **da sequência ordenada exata de limites inteiros de ambos os métodos**, além de contagens e larguras máximas. O algoritmo não usa `float` nem `mpmath` para calcular intervalos. A escolha inteira da volta de `2π` é validada por uma condição de redução para `[-4,4]`.

## Reproduzir

Requer Python 3.11 ou superior e `pytest`, `mpmath` para os testes numéricos auxiliares. O cálculo integral não depende de `mpmath`.

```bash
python -m pip install pytest mpmath==1.3.0
python check_manifest.py
python -m pytest -q test_es8p8_full.py
python es8p8_full.py --start 0 --end 1600 --chunk 50 --terms 128 --out RESULTADOS_ES8P8_REPLAY.json
python validate_replay.py RESULTADOS_ES8P8.json RESULTADOS_ES8P8_REPLAY.json
```

A execução integral levou aproximadamente **150 segundos** no ambiente local de referência; o tempo pode variar. Para execução parcial: `python es8p8_full.py --start 0 --end 50 --chunk 50 --out parcial.json`. A execução parcial recebe status diferente da integral e não pode ser interpretada como conclusão de toda a malha.

## Conteúdo

- `es8p8_full.py`: motor intervalar independente + auditoria integral.
- `test_es8p8_full.py`: testes de integridade, sinais, seno, sinc, digama e controles negativos.
- `fase4_fixedpoint_sem_arb.py` e `fase4_nos_sem_arb.py`: fontes históricas autenticadas.
- `es8p7_independent.py`: biblioteca independente de constantes racionais (proveniência ES8-P7).
- `RESULTADOS_ES8P8.json`: 32 lotes, SHA-256 dos pares intervalares, todas as contagens, largura máxima e resultado final.
- `RESUMO_FASE4_1600CELULAS_D40.json`: resumo original da fase histórica, autenticado.
- `REFERENCIAS_ORIGINAIS_17.csv`, `REFERENCIAS_TRIADAS_ES8P8.csv`: trabalhos de Boyadzhiev e Kanemitsu et al.; registrados como fundamentação, não como terceiro avaliador executado.
- `PROVENIENCIA_ES8P8.json`, `RELATORIO_ES8P8.md`, `MANIFESTO_SHA256.json`: escopo dos hashes, demonstrações matemáticas, limites e integridade.

## Interpretação científica

`FULL_GRID_PAIRWISE_INTERVAL_INTERSECTIONS_CHECKED` é um teste de **compatibilidade de inclusões**, não igualdade de todos os limites e não demonstração em Lean. Os dois avaliadores compartilham a identidade clássica da digama e a mesma transformação Fourier/sinc; sua independência é aritmética e de controle dos restos, não de todo o fundamento analítico. A validade das inclusões permanece condicionada às hipóteses analíticas explicitadas, e a positividade da forma de Weil para uma classe infinita permanece em aberto.
