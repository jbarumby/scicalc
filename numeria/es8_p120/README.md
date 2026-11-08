# NUMERIA P120 — Ponte Lean condicional
O arquivo P120_PositiveBridge.lean expressa o passo puramente lógico: W=C+E hermitiana, forma central limitada por shift*||x||², erro >= -rho*||x||², com shift>rho, implica W.PosDef.

**Não** formaliza os 169 centros, a identidade de Weil, o algoritmo Bareiss, o polinômio FLINT nem as inclusões Arb/GL/Cauchy. Uma execução CI verde confirma somente a ponte condicional, e não a Hipótese de Riemann.

Referência de Mathlib: Matrix.PosDef.of_dotProduct_mulVec_pos.
