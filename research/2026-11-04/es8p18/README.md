# NUMERIA ES8-P18 — correção de polos em ℓ² e orçamento projetado melhorado

**2026-11-04.** Prossegue a ES8-P17 e reexamina o último LaTeX recebido (`deepseek_latex_20261104_dd6431.tex`); não houve novo anexo neste pedido.

Resultado novo: `||B_N||op² ≤ 6275281/3686400 <171/100` para todo `N≥3` a partir dos coeficientes históricos com `L=9/20` e `π²≥9`. Isso melhora a cota `21/10` anterior. O orçamento condicional da cauda arquimediana após projeção passa a `3857/506250 <1/100`. Provamos matematicamente que a **correção dos polos** se estende em ℓ² com erro de truncamento em norma `≤ (1609/800)/√N` quando a identidade é mantida em todas as coordenadas; as truncações de posto finito de `B_∞` convergem apenas fortemente.

`python -m unittest -v test_es8p18.py`; `python audit_es8p18.py --out RESULTADOS_ES8P18_REPLAY.json`; conferir o JSON.

**Lean:** `NUMERIA_E8K_FORMAL/ES8P18.lean`, compilado pelo workflow `numeria-es8p18-lean.yml` (o compilador é a condição para declarar os teoremas aprovados). O módulo prova implicações **sob hipóteses Cauchy–Schwarz e E7 explícitas**; não formaliza o operador ℓ² inteiro nem a identidade analítica E7.

Não demonstra forma infinita de Weil, margem positiva uniforme ou Hipótese de Riemann. A `main` permanece inalterada.
