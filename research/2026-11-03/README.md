# NUMERIA–RIEMANN — pesquisa registrada em 03/11/2026

Esta pasta registra materiais das frentes NL1/ES1, NL2/ES2 e NB1, preservando a branch principal do projeto.

## Conteúdo
- `nl1_es1/`: código do método de Newton, Krawczyk racional para modelo 2×2, 10 testes, resultados numéricos e classificação bibliográfica.
- `nb1/`: código para certificado finito da aproximação de Nyman–Beurling, 5 testes e relatório.
- `inputs/NUMERIA_PESQUISA_NL2_ES2_2026-11-03.md`: levantamento e proposta metodológica, NÃO uma implementação executada.
- `inputs/NUMERIA_REVISAO_NYMAN_MELLIN_2026-11-03.tex`: revisão matemática; `inputs/deepseek_latex_20261103_c059d1.tex` é texto original recebido, conservado apenas para rastreabilidade, NÃO validado.
- `inputs/References_LeapSpace*.csv`: listas recebidas como insumo bibliográfico.

## Limites de validade
1. A certificação Krawczyk de NL1 refere-se somente a um sistema matricial ilustrativo 2×2; as raízes aproximadas da zeta NÃO foram certificadas por intervalos.
2. Os certificados NB1 são estimativas de erro para escolhas finitas de coeficientes, NÃO prova de convergência da aproximação ótima nem da Hipótese de Riemann.
3. As análises NL2/ES2 são propostas de pesquisa sem reproduções espectrais executadas.
4. O extrato da V24 é um dado auxiliar e NÃO uma prova de equivalência das matrizes projetadas.
5. Nenhum resultado contido nesta pasta prova a Hipótese de Riemann.

## Auditoria
Consulte `PROVENIENCIA_PUBLICACAO.json` para os SHA-256 de cada material de origem.
Testes locais reexecutados no registro: NL1/ES1: 10 aprovados; NB1: 5 aprovados.
O workflow da raiz `.github/workflows/numeria-nl1-es1-nb1.yml` verifica novamente os testes no GitHub Actions.
