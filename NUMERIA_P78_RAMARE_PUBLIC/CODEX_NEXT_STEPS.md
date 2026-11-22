# Retomada Codex — P78: substituir a dependência de comunicação privada

Branch inicial: `numeria-p78-ramare-public-source`, baseada em P77. Não alterar branches antigas.

1. Confirmar o SHA e o replay do certificado `NUMERIA_P78_RAMARE_PUBLIC/certificate.json` (frações exatas) e conferir no preprint Ramaré–Zúñiga Alterman **v1**, arXiv:2603.25961v1, os cinco limites citados. A fonte é um preprint público, não artigo necessariamente revisado.
2. Não atribuir à P78 prova de b4=0.54 em todo o Lema 6. As linhas 1–9 retêm b4=0.54 **apenas** na faixa X≥3e12, eta≤2delta/logX. Usar a tabela de constantes b4 alternativas nas linhas 10–16.
3. Criar um novo auditor Arb (não editar o P77 original) que aceite b4 por linha e reavalie `C3,C4,K3,L1,L2,A1,A2,B1` em 192/256 bits com a quadratura K1 certificada da P77.
4. Se os B1 do impresso não majorarem os recalculados, gerar uma tabela superior nova e adaptar a P75 em uma branch separada para refazer t∈[18,79]. Não chamar células com limites inferiores não positivos de contraexemplos.
5. Comparar com `audit/p75-journal-only-baseline` e verificar a integração da cauda, as simetrias e a identidade de Suzuki independentemente.
6. Apenas depois de reproduzir os gates, considerar t>79. Registrar JSON/SHA256/CI/PR em rascunho, sem mesclar PRs.