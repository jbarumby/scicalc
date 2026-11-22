# P89 — Parecer técnico dos seis arquivos P87 enviados

**Base revisada:** `p87_teorema_preissmann.pdf` (10 pp), `.tex`, `p87_verificacao.py`, `p87_ci.yml`, `p87_capa.html`, `p87_SHA256SUMS.txt`, fornecidos em 22/11/2026. Conferência **local** de SHA-256 confirma os cinco arquivos abrangidos pelo manifesto:

| Arquivo | SHA-256 recebido e verificado localmente |
|---|---|
| `p87_teorema_preissmann.pdf` | `616de178ba794135deafe3dccd041c50f60f6553a370c6b0bd93e99362096b05` |
| `p87_teorema_preissmann.tex` | `d16771dab750ef4f00e00bf841ddcf5f179e412ef6392914ea819a0e7b9d706d` |
| `p87_verificacao.py` | `5089f7016827dffb5a6202cd783ca70ab099ba5a87db37e954e2192a1b8ff14a` |
| `p87_ci.yml` | `d6fccbdfaee9e8bf1bd5116f69d38ebba9e31049a87a00aac52b51360ed2be2a` |
| `p87_capa.html` | `a1d420bdf90d56009d9fdb2574b1344b50262ddd1a6041ed8d5dee3153558135` |

O próprio `p87_SHA256SUMS.txt` referencia também `p87_ledger.json` e `p87_execucao.log`, **não enviados nesta rodada**. Seus conteúdos/resultados não podem ser auditados apenas a partir do manifesto.

## 1. Erro matemático na nota, Lema 5.3, pp. 5–6

O PDF e o TeX afirmam, a partir da mudança de variável `u=t s`, que

```
t ψ_1(t) = integral_0^infinity u exp(-u)/(1-exp(-u/t)) du
```

**FALTA o fator `1/t`.** A identidade correta é

\[
\boxed{t\psi_1(t)=\int_0^\infty \frac{u e^{-u}}{t(1-e^{-u/t})}\,du.}
\]

Além disso, para `t2>t1` e `u>0`,
`e^{-u/t2}>e^{-u/t1}` implica **`1-e^{-u/t2}<1-e^{-u/t1}`**, oposto à desigualdade no texto. As duas falhas invalidam o raciocínio de monotonicidade **como redigido**, embora a conclusão seja verdadeira.

### Reparo rigoroso, sem cálculos numéricos

Defina \(\mathcal D_u(t)=t(1-e^{-u/t})\). Então
\[
\mathcal D'_u(t)=1-(1+u/t)e^{-u/t}>0,
\]
pois \(e^v>1+v\) para \(v>0\). Logo, \(\mathcal D_u(t)\) cresce estritamente com `t>0`; consequentemente o integrando positivo
\(\frac{u e^{-u}}{\mathcal D_u(t)}\)
decresce estritamente. A comparação por integração prova que
\(t\mapsto t\psi_1(t)\) é estritamente decrescente. Como
\(\Phi^*(r)=t\psi_1(t)\) para \(t=1/r\), \(\Phi^*\) cresce estritamente.

Isso **corrige o Lema 5.3** mantendo intactos o teorema de Bellman e sua conclusão, desde que os demais passos do documento sejam preservados.

**Correção secundária de notação:** na Proposição 5.2, ao reindexar a soma que representa \(\frac1{1+r}\Phi^*(r/(1+r))\), o termo \(\sum r(1+r)/[1+(k+1)r]^2\) corresponde à \(\Phi^*\) *antes* da multiplicação por \(1/(1+r)\), e não depois. A identidade final \(\Phi^*(r)=r+V(r)\) está correta.

## 2. Erro de direção do certificado B5

No arquivo `p87_verificacao.py`, a função `block_b5(pi_u)` calcula
`bound23 = pi_u**2/3` e `bound26 = pi_u**2/6`, onde `pi_u` é uma **cota superior** rigorosa de \(\pi\). Da desigualdade `Q<pi_u²/3` não se segue `Q<pi²/3`: o primeiro alvo é *maior* que o segundo.

**Patch conceitual obrigatório:**

```python
# Em main(), substituir:
r["B5"] = block_b5(pi_u)
# por:
r["B5"] = block_b5(pi_l)
# e renomear o parâmetro de block_b5 para pi_lower;
# os valores bound23, bound26 passam a ser cotas INFERIORES rigorosas.
```

Também **preservar os resultados exatos no JSON**:
```python
"one_sided_min_margin_exact": str(min_margin26),
"Qi_min_margin_exact": str(min_margin23),
# e em families:
"margin_exact": str(marg)
```
Em vez de serializar a única margem certificadora como `float`, a CI deve validar `Fraction(...)>0` para todas as margens. `float` pode ser usado só para apresentação.

## 3. Separação rigorosa de gates

- A demonstração matemática universal da energia `Q_i<=pi²/3` é conhecida via P87 (Shan) e continua válida; o Lema 5.3 de Bellman pode ser reparado como acima, tornando uma via alternativa útil.
- B0: cotas de π por Machin têm origem exata, mas a comparação com `mpmath` em 40 dígitos é diagnóstico.
- B1/B2: varreduras de `mpmath` de monotonicidade **não substituem** a demonstração analítica no TeX.
- B3/B4/B6/B7: iterações float/`numpy`, otimização coordenada, autovalores `float64` e testes finitos são **diagnósticos**; não constituem certificado de norma para todos os N.
- B5: desigualdades de conjuntos racionais podem ser realmente certificadas com `Fraction` e **cota inferior** correta de π, mas a implementação enviada precisa do reparo.
- `p87_ci.yml`: a aprovação eventual da CI significa que as condições programadas passaram; o código presente chama `float(...)` nas verificações do B5. Também aponta para branch `numeria-p87-preissmann-energy-lemma`, diferente da cadeia atual P87/P88.
- A hipótese espectral `(H)` dos Lemas 2–6 de Preissmann **permanece importada**, como a própria nota reconhece. O título "Teorema geral de Preissmann" descreve a montagem **condicional a (H)**, não uma prova independente completa do limite `pi*m0`.

## 4. Próxima etapa: P89 (Lema 5, caudas exteriores)

A prova original de Preissmann usa o Lema 5 para controlar somas de recíprocos após uniformizar os pontos exteriores da dupla `a,b` (Lema 6, pp. 108–110). A nova P89 demonstra uma **subdesigualdade universal explícita para essas duas caudas** via o Lema de convexidade da P87, sem importar o Lema 5 ou o Lema 6. Veja [README.md](README.md).

O **trecho central** `a<k<b`, cujo termo cruzado é negativo, segue necessitando uma prova variacional afiada. Não chamar a P89 de reconstrução completa do Lema 5, de Lema 6 ou da constante `m0`.

**Fontes:** [Preissmann (1984)](https://doi.org/10.5169/seals-53823); [Yangjit (2022/2023)](https://arxiv.org/abs/2203.14950); [P87](https://github.com/jbarumby/scicalc/pull/87); [P88](https://github.com/jbarumby/scicalc/pull/88).
