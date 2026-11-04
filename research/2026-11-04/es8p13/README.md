# NUMERIA–RIEMANN · ES8-P13 · Certificados LDL concretos em Lean 4

Fonte: ES8-P9 (30 integrais, Fase 4) → ES8-P10 (cota melhorada de Cauchy) → ES8-P12 (teorema LDL geral). Data 4 de novembro de 2026.

**Status inicial:** 10 testes Python exatos aprovados; a compilação Lean será verificada em GitHub Actions.

## Conteúdo

- `ES8P13.lean`: prova de positividade das três matrizes racionais explicitamente fornecidas (`N=6,8,10`, dimensões `4,6,8`) por uma fatoração racional `UᵀDU`; contém a verificação por casos de todas as identidades das entradas e de todos os pivôs.
- `RESULTADOS_ES8P13.json`: entradas racionais simplificadas, pivôs exatos e certificados dos extremos intervalares originais de aritmética fixa em escala `10^90`.
- `scratch.py`: reconstrói os intervalos projetados a partir dos arquivos da ES8-P10.
- `generate_es8p13.py`: reexecuta a geração sem números de ponto flutuante para testes/certificados (float apenas em campos opcionais explicativos).
- `test_es8p13.py`: testes matemáticos positivos e controles de mutação negativa.
- `numeria-es8p13-lean.yml`: compila Lean 4.35.0 + Mathlib e reexecuta a reconstrução no GitHub.

## Margem condicional

Para cada `N`, seja `A_N` a matriz intervalar ES8-P10 e `C_N` a sua aproximação racional com denominador `10000`. Foi verificado em frações exatas:

\[
\max_i\sum_j \max\{|(A_N)_{ij,-}-(C_N)_{ij}|,\ |(A_N)_{ij,+}-(C_N)_{ij}|\}<\frac1{100}.
\]

No Lean, a fatoração exata prova \(C_N-\frac{3}{50}I\succ0\) para `N=6,8,10`. Junto à estimativa elementar para matrizes simétricas pelo máximo das somas absolutas das linhas, isso implica, **sob a validade das inclusões analíticas originais**, \(M_N-\frac1{20}I\succ0\).

**Limite formal:** o Lean comprova a positividade das matrizes racionais simplificadas. A checagem que relaciona seus coeficientes às inclusões de 90 dígitos e a cota de perturbação é feita em Python exato; os lemas analíticos para as integrais não foram totalmente mecanizados. Não é uma prova da positividade universal de Weil nem da Hipótese de Riemann.
