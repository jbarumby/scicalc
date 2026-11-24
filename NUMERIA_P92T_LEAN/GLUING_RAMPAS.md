# P92-T — Roteiro rigoroso para colar as células de suporte da rampa

**Status atualizado:** o argumento de colagem finita, a escolha automática dos suportes por `Nat.ceil`, a identidade da rampa deslocada e o Lema R1 global **para rampas e combinações positivas finitas** agora estão formalizados e compilados em Lean: [NumeriaP92T.lean](NumeriaP92T.lean), [GlobalRamp.lean](GlobalRamp.lean) e [ShiftRamp.lean](ShiftRamp.lean); [CI #45308852409](https://github.com/jbarumby/scicalc/actions/runs/45308852409). Este documento preserva a demonstração em papel; **a extensão a misturas integrais e os Lemas 5–6 completos ainda NÃO estão formalizados**.

## 1. Potenciais para a rampa

Sejam `L>=0` e `h>0`. Defina
\[
V_0(L,h)=h\sum_{j=0}^{\infty}(L-jh)_+,
\qquad
V_1(L,h)=h\sum_{j=1}^{\infty}(L-jh)_+
=V_0(L,h)-hL.
\]
Por `h>0` e `L<∞`, essas somas têm número finito de termos não nulos. Para `L=0`, são nulas. Para `L>0`, ponha
\[
n=n(h)=\left\lceil\frac Lh\right\rceil\ge1.
\]
Então `(n-1)h<L<=nh`; os termos de índices `0,...,n-1` são positivos, e os de índices `j>=n` são nulos, inclusive no caso `L=nh`. Segue
\[
V_0(L,h)=hnL-\frac{h^2n(n-1)}2=:P_0(L,h,n)
\]
e
\[
V_1(L,h)=h(n-1)L-\frac{h^2n(n-1)}2=:P_1(L,h,n).
\]

## 2. Monotonia dentro da célula

Na célula `L/n<=h<L/(n-1)` (o lado direito é omitido para `n=1`), `n(h)=n` permanece constante. Portanto
\[
\frac{\partial}{\partial h}P_0(L,h,n)
=n[L-(n-1)h]\ge0,
\]
\[
\frac{\partial}{\partial h}P_1(L,h,n)
=(n-1)[L-nh]\le0.
\]
O arquivo Lean demonstra uma versão **inteiramente algébrica** de ambas as monotonias entre dois passos `h_1<=h_2` que pertençam à mesma célula. Não depende de diferenciação automática.

## 3. Continuidade em cada nó

Na fronteira `h=L/n` com `n>=1` ocorre uma mudança de número de parcelas. Porém:
\[
P_0(L,h,n+1)-P_0(L,h,n)=h(L-nh)=0
\quad\text{quando }L=nh.
\]
Também `P_1(L,h,n+1)-P_1(L,h,n)=0`, pois ambos se obtêm de `P_0` retirando a mesma parcela `hL`.

Essa igualdade vale até com `h,L` simbólicos, sob a única hipótese `L=nh`; inclui a condição de igualdade do `Nat.ceil`. As versões correspondentes para a soma com `max` também constam do módulo Lean em desenvolvimento.

## 4. Colagem finita de células

Fixe quaisquer `0<h_1<h_2`. Os nós que podem ocorrer no intervalo compacto `[h_1,h_2]` são
\[
\{L/k:k\in\mathbb N_{\ge1},\ h_1\le L/k\le h_2\}.
\]
Esse conjunto é **finito**, porque `k<=L/h_1`. Ordene os nós e subdivida `[h_1,h_2]` pelos pontos. Em cada subintervalo, `n(h)` é constante e `V_0` é não decrescente, `V_1` não crescente (passo 2). Nas junções, os valores de cada lado coincidem (passo 3). Encadeando as desigualdades ao longo da partição finita obtemos
\[
\boxed{V_0(L,h_1)\le V_0(L,h_2),\qquad
V_1(L,h_1)\ge V_1(L,h_2).}
\]
Essa é a **justificativa matemática global**. A indução que atravessa qualquer número de células foi implementada em `rampPoly_mono_across_cells` e `finiteRamp_mono_across_cells`; a passagem aos suportes exatos com `Nat.ceil` foi implementada em `GlobalRamp.lean` e verificada no kernel Lean.

## 5. Bellman para uma rampa, deduzido da colagem

Seja `f_T(u)=(T-u)_+` e `w>0`, com `L=(T-w)_+`. Então `f_T(w)=L`,
\[
U_0(w,h)=V_0(L,h),\qquad U_0(w+b,b)=V_1(L,b).
\]
Se `b<=h`, usando `V_0(L,b)\le V_0(L,h)`:
\[
\min(h,b)L+V_1(L,b)=V_0(L,b)\le V_0(L,h).
\]
Se `h<=b`, usando `V_1(L,b)\le V_1(L,h)`:
\[
\min(h,b)L+V_1(L,b)
=hL+V_1(L,b)\le hL+V_1(L,h)=V_0(L,h).
\]
Logo a desigualdade R1 vale para cada rampa, **em matemática ordinária**. Os módulos Lean agora verificam também a colagem de todas as células, a identidade `U_1(L,b)=U_0((L-b)_+,b)` e o enunciado `bellman_for_all_ramps` na forma deslocada original.

## 6. Extensão positiva — gate posterior

Se `f(u)=\int_{(0,\infty)}(T-u)_+\,d\mu(T)` com `μ` positiva e ambos os potenciais finitos, então Tonelli permite integrar a desigualdade da rampa sobre T:
\[
\min(h,b)f(w)+U_0(w+b,b)\le U_0(w,h).
\]

- Para `f(u)=u^{-2}`, pode-se escolher `d\mu(T)=6T^{-4}dT`.
- Para `f_X(u)=(1/u-1/X)_+`, a medida é `2T^{-3}1_{T<X}dT+X^{-2}\delta_X`. **A massa atômica na fronteira `T=X` é indispensável.**
- Para `phi_H(u)=1/[u^2(u+H)^2]`, a densidade positiva é `phi_H''(T)` e os termos de fronteira se anulam no infinito.

A passagem por `tsum`/integrais não foi formalizada ainda em Lean. Esse é o passo seguinte à colagem.

## 7. Lista de obrigações verificáveis para o próximo commit

1. **Concluído no Lean:** definir `V_0(L,h)` via `Nat.ceil (L/h)` e provar a identidade com as somas de rampas e com `P_0(L,h,n(h))`.
2. **Concluído no Lean por abordagem alternativa:** controlar toda a cadeia finita de nós por indução em seus índices, sem enumerar um conjunto explícito de pontos de transição.
3. **Concluído no Lean:** indução de colagem e Lema R1 global de cada rampa, sem hipóteses de célula.
4. **Em aberto:** levantar R1 a misturas integrais positivas por Tonelli, com domínios de integrabilidade explícitos. Misturas positivas FINITAS já foram formalizadas.
5. Remover a hipótese `step` da cadeia P92-S para os núcleos realmente cobertos.

**Governança:** prova em papel/Lean parcial explicitamente separados; PR #109 em rascunho. Não alegar formalização completa de Preissmann 5/6 ou RH.
