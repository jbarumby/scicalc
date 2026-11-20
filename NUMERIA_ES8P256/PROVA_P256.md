# NUMERIA–RIEMANN · ES8-P256 — Ponte exata entre von Mangoldt ponderada e o símbolo de primos

**Escopo:** identidade elementar independente da RH; rastreio estático de lemas Lean de Alpöge; nenhum uso de cálculo intervalar como substituto de prova funcional, nenhum novo limite coercivo para a forma de Weil.

## 1. O elo `allChi` visto nos fontes originais

Na revisão fixa `Deicyde/riemann-hypothesis@40e56bab`, os módulos `Main.lean`, `MainConclusion.lean`, `EveryCharacter.lean`, `RiemannHypothesisOfEveryCharacter.lean`, `RowSrc.lean` e `EveryCharacterOfLiteratureOfVWin.lean` exibem a cadeia:

```
riemannHypothesis : RiemannHypothesis  := rhOfAll
  rhOfAll                         := Formalism.VW.rh_of_every allChi
  allChi : UpperTheory.EveryCharacter
    := UpperTheory.FinalTheorems.everyCharacter
  FinalTheorems.everyCharacter
    := everyCharacter_of_row Gal.Int.galoisLocalGlobalStrong_recParameter
  everyCharacter_of_row (h : GalStar.Row0b ...) : EveryCharacter
    := every_character_of_literature_of_vWin O(h) (vwin_of_row h)
  every_character_of_literature_of_vWin
    (O : OpenLiteratureStatements) (hV : Top.VWin) : EveryCharacter
    := Lit.every_character_of_theta_integrals_local_matching ...
```

`RowSrc.lean` obtém o parâmetro Galois via `galoisRepresentations_of_thetaLift ... tLiftK`. **As hipóteses `O` e `hV` aparecem explicitamente no lema intermediário**; são preenchidas em módulos posteriores, não podem ser descartadas por omissão. Os imports transitivos são extensos; o rastreio de declarações NÃO verifica integralmente o termo Lean. Em `ClaimStrip.lean` há também uma ponte elementar: `factor_ne_zero` garante que o fator de Euler associado a um primo não se anula quando `Re(s)>0`, sendo invocado por `prim_zero` ao reduzir zeros de caracteres imprimitivos a primitivos.

Um fortalecimento quantitativo matematicamente imediato — **derivado aqui**, não enunciado como tal no módulo Lean — é: para `p` primo e `σ=Re(s)>0`, como `|χ(p)|≤1`, 

\[
 |1-\chi(p)p^{-s}|\ge 1-|\chi(p)|p^{-\sigma}\ge1-p^{-\sigma}>0.
\]

Esse limite é útil localmente e **não prova a não anulação de uma função L inteira em `1/2<σ≤1`**, onde o produto de Euler não converge absolutamente da mesma maneira. O passo global profundo continua localizado em `everyCharacter`, via resultados automórficos.

## 2. Identidade ponderada exata (nova ponte independente do NUMERIA)

Defina a função de von Mangoldt `Λ(q)=log p` para `q=p^k` (k≥1) e `Λ(q)=0` para outros inteiros. Para `L>0`, `X=e^{2L}` e `t∈R`, o símbolo usado no NUMERIA é

\[
 S_u(t)= 2\sum_{2\le q<e^{2u}} \frac{\Lambda(q)}{\sqrt q}\cos(t\log q).
\]

Defina a versão oscilatória da soma ponderada do anúncio,

\[
 W_L(t)=\sum_{2\le q<e^{2L}}\Lambda(q)\Big(\frac{e^{2L}}q-1\Big)\cos(t\log q).
\]

Então, **exatamente, sem RH e sem erro numérico**,

\[
 \boxed{W_L(t)=e^L\int_0^L\cosh(L-u)\,S_u(t)\,du.}\tag{P256-1}
\]

**Prova.** A soma sobre `q<e^{2L}` é finita; podemos trocar soma e integral. Um termo com `q` entra em `S_u(t)` precisamente para `u > a_q=(log q)/2`. Sua integral é

\[
 e^L\frac{2\Lambda(q)}{\sqrt q}\cos(t\log q)
 \int_{a_q}^L\cosh(L-u)du
 =e^L\frac{2\Lambda(q)}{\sqrt q}\sinh(L-a_q)\cos(t\log q)
 =\Lambda(q)\Big(\frac{e^{2L}}q-1\Big)\cos(t\log q).
\]

A última igualdade decorre de `2 sinh v = e^v-e^{-v}` e `e^{a_q}=√q`. Se `q=e^{2L}`, o termo tem peso zero; as convenções `<` ou `≤` não alteram a identidade.

Para `t=0`, a soma `W_L(0)` é precisamente a parcela ponderada positiva mencionada na nota de Alpöge, **antes** da subtração dos termos principais `−X log X+(1+γ)X`. Além disso, com `ψ(v)=Σ_{q≤v}Λ(q)` (Chebyshev), temos a identidade de Abel

\[
 \boxed{W_L(0)=X\int_1^X \frac{\psi(v)}{v^2}\,dv.}\tag{P256-2}
\]

Pois `X∫_q^X v^{-2}dv=X/q−1` para cada `q≤X`.

Se `P_L(z)=Σ_{q<X} Λ(q)q^{-1/2}e^{iz\log q}`, então outra identidade exata é

\[
 \boxed{W_L(0)=XP_L(i/2)-P_L(-i/2).}\tag{P256-3}
\]

Ela exprime a ponderação por dois valores complexos de uma **soma finita inteira**, não por uma continuação analítica de ζ ou do operador de Weil.

## 3. O limiar do primo 7

Para `7<X<8`, nenhum novo termo de Mangoldt além de `q=7` entra relativamente à faixa `5<X<7`, pois `Λ(6)=0`. Exatamente,

\[
 W_X(t)-W_X^{\{2,3,4,5\}}(t)=\frac{X-7}{7}\log7\cos(t\log7).
\]

Assim, ao tomar `X=7+ε`, a amplitude da contribuição ponderada cresce como `ε log7/7` e **tende a zero** quando `ε→0+`. Em contrapartida, o símbolo duro 

\[
 S_L(t)=2\sum_{q<e^{2L}}\frac{\Lambda(q)}{\sqrt q}\cos(t\log q)
\]

adquire imediatamente o termo de amplitude `2 log7/√7` ao ultrapassar o limiar. A soma ponderada `W_X(t)` é contínua em `X`; sua derivada (com `t` fixo) sofre salto `log7/7 cos(t log7)`.

Isso sugere **testar pesos de suporte suaves** em modelos auxiliares para atenuar a ativação do primo 7. **Não elimina a obstrução da P253:** substituir o símbolo duro por um peso ponderado muda o operador. Não há nesta etapa uma prova de que a positividade do novo modelo implique a positividade original de Weil, nem um ganho da margem da P252.

## 4. Auditoria exata em Python

O arquivo `weighted_bridge.py` trata `log(p)` e `cos(t log q)` como símbolos formais indexados por primos/potências primas e verifica em frações racionais as identidades de coeficientes de Abel, Mellin finita e ativação. Para `X=6`, verifica-se

\[
 W(6)=\tfrac52\log2 + \log3 + \tfrac15\log5.
\]

Testes abrangem os limiares 4, 6, 7 e 8, valores racionais estritos de `X`, denominadores exatos e alterações adversariais. O programa não avalia numericamente a integral do item (P256-1), cuja prova é analítica, nem reivindica compilar Lean ou certificar globalmente RH.

## 5. Próxima etapa sugerida — P257

Definir uma família de símbolos de teste com taper `w_L(q)` que se anula em `q=e^{2L}`; construir sua forma no espaço de suporte e **demonstrar uma comparação operacional com a forma de Weil original** (ou registrar rigorosamente por que falha). Só depois comparar margens P252 e atravessar `log7/2`. Em paralelo, escolher um lema quantitativo realmente global do bloco automórfico de Alpöge para tradução legível; o limite de Euler acima é útil mas claramente insuficiente para RH.