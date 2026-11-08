# NUMERIA ES8-P115 — auditoria independente após 24/24 shards P111

Data: 2026-11-08

O run de origem da quadratura é https://github.com/jbarumby/scicalc/actions/runs/42509353671.

A verificação a seguir recupera todos os 24 artefatos P111 do run original,
reconstitui 14.365 entradas de intervalos Arb em 256 bits e confronta
cada uma com o resultado da P109 via frações exatas, sem confiar em um
simples número de pivôs fornecido pelo JSON.

O teste de positividade usa aritmética inteira FLINT para o polinômio
característico do centro deslocado. A desigualdade final só pode ser
aceita se a matriz deslocada for definida positiva e `shift>rho_Frobenius`.

O estado verde do preflight isolado NÃO é uma prova de positividade.
Caso o critério falhe, a auditoria deve ser marcada INCONCLUSIVA.

Este arquivo acompanha a execução da auditoria, não é testemunha matemática.
