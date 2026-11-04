# ES8-P17 — Norma da projeção, cauda projetada e crítica do LaTeX

A ES8-P17 incorpora `deepseek_latex_20261104_dd6431.tex` como **material de investigação**, nunca como prova. Reconstitui a fórmula histórica da matriz `B_N` em `es8p_projection.py` e `AUDITAR_FASE3_PROJECAO.py`, estabelece a cota dimension-free `||B_N||² ≤ 5148881/2560000 < 21/10`, e combina com o orçamento arquimediano bruto condicional da ES8-P15 para obter `||B_Nᵀ E_N B_N|| < 1/100` com `R(N)=200N²`.

**Escopo:** somente cauda arquimediana para `L=9/20`, com a estimativa integral E7 como premissa. A prova Lean somente verifica a parte racional e a implicação dadas as cotas como hipóteses. A derivação da norma de `B_N` é apresentada em papel e auditada com frações exatas em Python; não é ainda um teorema do Lean sobre `ℓ²`. Não há prova de positividade de Weil universal nem de RH.

Execução: `python -m unittest -v test_es8p17.py && python audit_es8p17.py --output RESULTADOS_ES8P17_REPLAY.json` e comparar o JSON ao certificado versionado. Compilação via `lake env lean ES8P17.lean` com Lean 4.35.0/Mathlib no GitHub Actions.
