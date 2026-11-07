# NUMERIA–RIEMANN — ES8-P84
Auditoria de ordem fixa em 50 modos? **Não**: esta fase mantém as três matrizes históricas (24/26/44 modos), sem gaussianas.

## Objetivo
Recalcular as mesmas 30 entradas estratificadas da ES8-P79 com 12 subintervalos de r em [0,30], mas usando exatamente **110 nós de Gauss–Legendre** em todos os subintervalos em vez de ordem adaptativa. Isso verifica, de modo independente da seleção de n, a propagação do erro de quadratura analítica.

- Python 3.11, python-flint==0.9.0, ctx.prec=256
- parâmetro elíptico rho=34/25, semilargura por trecho L=5/4
- majorante global |f_ij(z)|<600 em cada elipse
- erro por trecho ≤8LM rho^(1-2·110)/(rho−1) < 9,5·10^-26
- erro somado dos 12 trechos <1,2·10^-24 por integral, antes do coeficiente de matriz.
- Não chama acb.integral; intervalos antigos não entram na integral, somente na comparação.

O arquivo fixed_order_replay.py contém a rotina, test_fixed.py os testes sem FLINT, verify_fixed.py exige os 6 arquivos Arb reais e confere 30 intervalos por frações racionais exatas.

A execução é iniciada por evento pull_request no PR #17, em rascunho. Não mesclar sem inspecionar os artefatos e auditar a identificação da forma de Weil.

## Limitações
O replay de ordem fixa não substitui provas de domínio/identificação nem de estabilidade espectral do complemento infinito. A positividade nas janelas compactas não demonstra a Hipótese de Riemann.
