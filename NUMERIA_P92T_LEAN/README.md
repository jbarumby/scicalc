# NUMERIA–RIEMANN P92-T — Bellman global de rampas formalizado em Lean 4

**Status verificado:** [GitHub Actions #45308852409](https://github.com/jbarumby/scicalc/actions/runs/45308852409) — `success`; Lean 4.35.0, Mathlib v4.35.0; os três módulos obrigatórios compilam sem `sorry` nem `admit`. Branch do [PR #109 (draft)](https://github.com/jbarumby/scicalc/pull/109) baseada na [PR #108, P92-S](https://github.com/jbarumby/scicalc/pull/108).

A P92-T **formaliza universalmente a desigualdade de Bellman para cada rampa convexa** e **para combinações finitas com coeficientes não negativos**. Ela ainda NÃO formaliza R1 para toda função convexa, nem os Lemas 5/6 completos de Preissmann ou a Hipótese de Riemann.

## 1. O novo teorema global

Para `L≥0` e `h>0`, seja
\[
n_h=\lceil L/h\rceil,\qquad
U_0(L,h)=h\sum_{j=0}^{n_h-1}(L-jh)_+,\qquad
U_1(L,h)=U_0(L,h)-hL.
\]
O Lean provou, **sem restrição de suporte comum**, que, se `0<a≤b`,
\[
\boxed{U_0(L,a)\le U_0(L,b),\qquad U_1(L,b)\le U_1(L,a).}
\]
Em consequência, para todos `L≥0` e `h,b>0`,
\[
\boxed{\min(h,b)L+U_1(L,b)\le U_0(L,h).}
\]

Para a rampa original `f_T(w)=(T-w)_+=\max(T-w,0)`, com `T,w\in\mathbb R`, defina o potencial real
\[
U_T(w,h)=U_0((T-w)_+,h).
\]
O teorema Lean `bellman_for_all_ramps` prova exatamente:
\[
\boxed{\min(h,b)f_T(w)+U_T(w+b,b)\le U_T(w,h),\qquad h,b>0.}
\]
A identidade necessária `U_1(L,b)=U_0((L-b)_+,b)` é **também formalizada** em `fullRamp_shift` usando `Nat.ceil_add_one`, inclusive os casos de fronteira.

## 2. Como a prova cruza todos os limiares

A prova não enumera casos racionais nem amostra valores de `h`.

1. `NumeriaP92T.lean` prova a fórmula polinomial exata da soma real finita, as duas monotonias em cada célula de suporte, a nulidade dos termos posteriores ao suporte e a continuidade em `L=nh`.
2. `rampPoly_mono_across_cells` e `finiteRamp_mono_across_cells` colam uma quantidade **arbitrária** de células por indução em `m≥n`: cortam em `q=L/m`, passam pelas igualdades de fronteira e usam a hipótese indutiva.
3. `GlobalRamp.lean` escolhe os suportes exatos `n_h=Nat.ceil(L/h)`, prova as duas desigualdades globais com `Nat.ceil_mono` e obtém `fullRamp_bellman`.
4. `ShiftRamp.lean` prova `finiteRamp0_shift`, `fullRamp_shift` e `bellman_for_all_ramps` com os argumentos explícitos `(T,w,h,b)`.
5. `bellman_finite_ramp_mixture` multiplica cada desigualdade por um peso não negativo e soma em qualquer `Finset`. Portanto uma família **finita** de rampas com pesos não negativos também satisfaz Bellman.

Todos os enunciados foram compilados em Lean 4.35.0/Mathlib v4.35.0. O workflow `.github/workflows/numeria-p92t-lean-bellman.yml` rejeita `sorry`/`admit` em **todos os três módulos** e exige `lake build` dos três alvos.

## 3. Como reproduzir

```bash
cd NUMERIA_P92T_LEAN
lake update
lake exe cache get
lake build
```

Fontes:

- [NumeriaP92T.lean](NumeriaP92T.lean) — células, continuidade e colagem finita arbitrária;
- [GlobalRamp.lean](GlobalRamp.lean) — suporte `Nat.ceil` e Bellman global;
- [ShiftRamp.lean](ShiftRamp.lean) — potencial deslocado e misturas finitas;
- [GLUING_RAMPAS.md](GLUING_RAMPAS.md) — reconstrução matemática em papel, originalmente escrita antes da colagem global verificada.

## 4. Pendências científicas e próximo gate P92-U

O resultado da P92-T **não** basta para remover o parâmetro de Bellman abstrato de todos os teoremas da P92-S, porque as funções `u^{-2}`, `f_X(u)=(1/u-1/X)_+` e `φ_H(u)=1/[u²(u+H)²]` são **misturas integrais**, não combinações finitas de rampas.

Próximas obrigações, em ordem:

1. Formalizar a extensão de `bellman_for_all_ramps` a integrais de rampas com medidas positivas (Tonelli e domínio de integrabilidade). Incluir a **massa atômica `X^{-2}` em `T=X`** para `f_X`.
2. Especializar o Bellman integral a `f_X` e `u^{-2}` e conectar diretamente os teoremas `pathCost_le_of_bellman` e `boundary_accounting` da P92-S, eliminando as premissas agora demonstráveis.
3. Formalizar a comparação R4 das caudas reais e o limite pareado absolutamente convergente R5→R6; **não** tratar as duas somas assinadas individuais como convergentes.
4. Somente então compor o Lema 6 afiado R7 no Lean, confrontando ainda as convenções com o artigo original de Preissmann.

**Governança:** PR #109 permanece em rascunho, sem merge automático. Os Lemas 5/6 gerais são resultados publicados de Preissmann, mas a reconstrução NUMERIA em Lean ainda é **parcial**. Não há prova da Hipótese de Riemann.
