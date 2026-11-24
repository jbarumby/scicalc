# NUMERIA P92-U — Bellman e misturas integrais positivas em Lean 4

**Base:** [PR #109 — P92-T](https://github.com/jbarumby/scicalc/pull/109). Objetivo: estender a desigualdade de Bellman global para rampas `f_T(u)=(T-u)_+`, já verificada pelo kernel Lean na P92-T, para **integrais de rampas com respeito a medidas positivas**, e em seguida especializar a representação de `f_X(u)=(1/u-1/X)_+`.

## 1. Teorema de integração

Se `μ` é qualquer medida positiva boreliana em `ℝ`, `h,b>0`, e as três funções de `T` abaixo são integráveis, a P92-U demonstra em Lean:

\[
\boxed{
\min(h,b)\int (T-w)_+\,d\mu(T)
+\int U_T(w+b,b)\,d\mu(T)
\le \int U_T(w,h)\,d\mu(T).
}
\]

`U_T(w,h)=h\sum_{j\ge0}(T-w-jh)_+` é o potencial da rampa com suporte finito exato, construído na P92-T.

O argumento **não assume Bellman como hipótese**: aplica o teorema já formalizado `bellman_for_all_ramps` ponto a ponto, integra a desigualdade por `integral_mono`, e distribui a integral pelos dois termos usando `integral_add` e `integral_const_mul`.

**Hipóteses que não desaparecem:** integrabilidade do integrando da rampa, do potencial deslocado e do potencial original. O teorema é universal sobre medidas e passos, mas **não afirma** que qualquer medida produz uma mistura integrável.

## 2. Por que a etapa é relevante para o Lema 5

Para `0<u≤X`, a função truncada `f_X(u)=(1/u-1/X)_+` admite a representação positiva

\[
f_X(u)=\int_u^X \frac{2(T-u)}{T^3}\,dT+\frac{X-u}{X^2}.
\]

O último termo é a contribuição do **átomo de massa `X^{-2}` em `T=X`**. Sem esse átomo a integral fica aquém de `f_X` pelo termo `(X-u)/X^2`.

A prova de Bellman para as misturas integrais é um componente; para aplicá-la de verdade a `f_X` ainda é necessário formalizar a identidade integral acima **com seus domínios**, além das condições de integrabilidade da medida `2 T^{-3} \mathbf 1_{(0,X)}dT+X^{-2}\delta_X` quando integrada contra as rampas relevantes.

## 3. Certificado parcial do átomo de fronteira

O módulo adicional TruncatedKernel.lean demonstra no Lean **a identidade algébrica dos extremos da primitiva** da densidade 2(T-u)/T³:

\[
\left[-\frac2T+\frac u{T^2}\right]_{u}^{X}
+\frac{X-u}{X^2}=\frac1u-\frac1X.
\]

Ele também formaliza, pelo lema integral_dirac da Mathlib, que
\[
\frac1{X^2}\int (T-u)_+\,d\delta_X(T)=\frac{X-u}{X^2}
\quad(0<u\le X),
\]
e comprova a não negatividade da contribuição atômica. A igualdade para a integral da **densidade contínua** ainda não foi formalizada: a identidade da primitiva e a avaliação pelo teorema fundamental do cálculo são obrigações diferentes.

**Porta seguinte:** provar a avaliação efetiva ∫_u^X 2(T-u)/T³dT e demonstrar que as medidas resultantes tornam integráveis os três potenciais de bellman_integral_ramps. Não confundir a igualdade de extremos da primitiva com uma prova formal do Teorema Fundamental do Cálculo.

## 4. Troca finita integral–soma e suporte uniforme (Lean)

O módulo `FiniteInterchange.lean` demonstra a igualdade de Bochner
\[
\int h\sum_{j\in S}(T-w-jh)_+\,d\mu(T)
=
h\sum_{j\in S}\int(T-w-jh)_+\,d\mu(T),
\]
para todo conjunto finito `S` e sempre que as rampas individuais são integráveis. Não reivindica intercâmbio de série infinita com integral.

Também demonstra que, **para todos os parâmetros de rampa `T≤X`**, um passo `h>0` e
\[
j\ge\left\lceil\frac{X-w}{h}\right\rceil
\]
implicam \((T-w-jh)_+=0\). O limite é uniforme em `T`; consequentemente, uma mistura concentrada em `T≤X` tem apenas um número finito de parcelas potenciais não nulas para cada `w,h`.

**Importante:** ainda falta demonstrar em Lean que a medida concreta para `f_X` (densidade `2/T³` abaixo de `X` mais o átomo em `X`) satisfaz as hipóteses de integrabilidade e que sua integral contínua tem a primitiva indicada. Sem esses dois resultados, o Lema 5 quantitativo ainda não está formalizado por completo.

### Estado de CI

- [Execução #45309959995](https://github.com/jbarumby/scicalc/actions/runs/45309959995): `IntegralRamp.lean` aceito pelo Lean.
- [Execução #45310189271](https://github.com/jbarumby/scicalc/actions/runs/45310189271): `TruncatedKernel.lean` aceito pelo Lean.
- [Execução #45310589835](https://github.com/jbarumby/scicalc/actions/runs/45310589835): todos os seis alvos Lean compilados, incluindo `FiniteInterchange.lean`, e verificação sem `sorry`/`admit` aprovada.

## 5. Escopo honesto

| Resultado | Estado |
|---|---|
| Bellman global para rampas individuais | Formalizado P92-T |
| Bellman para misturas finitas positivas | Formalizado P92-T |
| Integração positiva de Bellman para medidas arbitrárias **sob integrabilidade explícita** | **Compilado no Lean — CI #45310589835** |
| Extremos algébricos da primitiva e átomo Dirac | **Compilado no Lean**; o teorema fundamental do cálculo continua pendente |
| Troca de somas finitas com integrais e limite uniforme de suporte | **Compilado no Lean — CI #45310589835** |
| Lema 5 quantitativo para configurações reais arbitrárias no Lean | Pendente |
| Compensação infinita de dois índices R5→R6 | Pendente |
| Lema 6 geral no Lean ou Hipótese de Riemann | **Não demonstrados** |

## 6. Build e governança

Código: `NUMERIA_P92T_LEAN/IntegralRamp.lean`, `NUMERIA_P92T_LEAN/TruncatedKernel.lean` e `NUMERIA_P92T_LEAN/FiniteInterchange.lean`, incluídos nos alvos obrigatórios da biblioteca Lake. Lean 4.35.0 + Mathlib v4.35.0.

```bash
cd NUMERIA_P92T_LEAN
lake update
lake exe cache get
lake build
```

Workflow dedicado: `.github/workflows/numeria-p92u-lean-integral.yml`; exige compilação e proíbe `sorry`/`admit`.

PR P92-U em rascunho, sem merge. A demonstração geral publicada dos Lemas 5/6 de Preissmann não constitui prova da Hipótese de Riemann.
