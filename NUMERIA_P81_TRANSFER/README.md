# NUMERIA P81 — transferência pública de S_epsilon até C3/C4 (condicional)

**Status:** auditoria algébrica independente do fator ζ, assumindo os teoremas publicados de Ramaré–Zúñiga Alterman e a desigualdade prévia para F_X citada em Chourasiya–Simonič v2. **Não** é prova independente de todas essas hipóteses, nem prova global de RH.

## Motivação: novos documentos recebidos

- `p78_s2_passagem_completa (1).pdf` — SHA-256 `69baebd7d0d79a00ed76109f3a5aacdc1be2e79dd66489900c82d813566a2b43` (7 páginas).
- `p78_s2_cadeia_verificacao.json` — SHA-256 `0e821f2f285b3d15fe127d2728836cbc3af5b044ca916bb9209741becd73ce83`.

O JSON registra concordância em ponto flutuante a cerca de 8×10⁻¹⁴ entre as fórmulas reconstruídas e o protótipo. Isso é uma **checagem diagnóstica**, não um certificado Arb. O PDF é útil por separar ε₂=2δ/log X e ε₁=δ/log X, mas a sua hipótese H1 **não enuncia a folga adicional necessária** e sua "prova" da monotonicidade de G usa um O(ε²) sem cota uniforme.

O PDF ainda descarta `exp(γ ε₁)` por um "a fortiori" que, isoladamente, tem **sentido incorreto**. Não é legítimo eliminar um fator maior que 1 de um limite superior.

## Proposição P81 (transferência dos dois ramos)

**Hipóteses importadas:**

1. Para cada uma das 16 linhas e cada uma das 3 zonas de altura `X>=3e12`, existe uma cota publicada `0<=S_ε(X)<=U` válida uniformemente para `0<=ε<=ε_cap`, conforme as fontes v1/v3 auditadas parametricamente na P80.
2. Escreva `l` para o limite inferior racional certificado de `log X` na zona (`28,36,75`) e `η₂=2δ/l`, `η₁=δ/l`. A escolha da fonte da P80 requer `η₂<=ε_cap`; assim, ambas as escalas estão cobertas.
3. O número racional `b4` selecionado na P80 satisfaz **estritamente** `b4 > (1+η₂)*U` em todas as zonas.

**Resultado, para os dois valores ε_j = jδ/log X, j=1,2**:

```
ζ(1+ε) = 1 + ∑_{n≥2} n^(-1-ε)
       ≤ 1 + ∫_1^∞ x^(-1-ε)dx = 1+1/ε.

ε_j ζ(1+ε_j) S_{ε_j}(X)
  ≤ (1+ε_j) U ≤ (1+η₂) U < b4.

Λ_X(1+ε_j) := ζ(1+ε_j) S_{ε_j}(X) < b4/ε_j.
```

Todas as comparações numéricas/de faixa para essa proposição são realizadas **com `fractions.Fraction`, sem float**, em `transfer.py`: 16 linhas × 3 zonas × 2 ramos = **96 gates**. A etapa seguinte (transferência para C3/C4) é uma implicação algébrica das expressões do artigo original, e não uma reconstrução pela regressão de L ou m0.

## Derivação de C3 e C4

Chourasiya–Simonič v2, Lema 6, usa a estimativa preliminar (derivada de KLN) para `σ=1+δ/log X`:

```
F_X(σ,0,T,1/2)
 ≤ π*m0*Λ_X(2σ−1)
   + ((T+π*m0)/X) * [ b3*X^(−2(σ−1))
                      + Λ_X(σ)/(5*(5X)^(σ−1)) ].
```

Aqui `2σ−1=1+ε₂`, `σ=1+ε₁`. Aplicando a proposição:

- **C3:** `πm0 Λ_X(1+ε₂) < πm0*b4*log X/(2δ) ≤ C3(X0,δ)*log X` porque `γ>0` e `X≥X0`.
- **C4, primeiro termo:** `Λ_X(1+ε₁)/(5(5X)^ε₁) < [b4*log X/(5δ)]*e^(−δ)*5^(−ε₁) ≤ [b4/(5δ e^δ)]*log X`. Aqui a redução é justificada por **`5^(−ε₁)≤1`**, sem suprimir indevidamente `exp(γε₁)`.
- **C4, residual:** `b3*X^(−2δ/log X) = b3*e^(−2δ) ≤ [b3*e^(−2δ)/log X0]*log X`.

Obtém-se a mesma forma `C3=(πm0*b4/(2δ))*exp(2δγ/log X0)`, `C4=b4/(5δ e^δ)+b3 e^(−2δ)/log X0` **para X≥X0=3×10¹²**, sob as duas entradas externas acima. É suficiente para a faixa de aplicação com `X=hT` e `hT0>=3×10¹²`, conforme a parametrização verificada nos PRs anteriores; a hipótese exata `X=hT` deve continuar checada nas futuras utilizações.

**Em particular, não é preciso provar a desigualdade mais forte `ζ(1+ε)≤exp(γε)/ε` no intervalo `0<ε≤1/32`** para estabelecer a transferência pública da P80: a margem `(1+η₂)` já foi contabilizada na construção do b4.

## O que NÃO está demonstrado por P81

- Prova independente dos teoremas de Ramaré–Zúñiga v1/v3 para `S_ε(X)`. São estimativas externas citadas, não rederivadas.
- Prova independente da estimativa prévia de `F_X` (KLN/Chourasiya–Simonič), nem seus fatos de somas de Dirichlet.
- Extensão ao domínio original do Lema 6 `3×10¹⁰≤X<3×10¹²`. Os coeficientes públicos P80 não foram certificados ali.
- Modificação da tabela `B1`, nova janela da função de Suzuki ou prova da Hipótese de Riemann.

## Fontes

- [Chourasiya–Simonič, v2 (2025), Lema 6](https://arxiv.org/html/2507.15184v2).
- [Ramaré–Zúñiga Alterman, v1](https://arxiv.org/html/2603.25961v1) e [v3](https://arxiv.org/html/2603.25961v3).
- [P80 / fonte auditada](https://github.com/jbarumby/scicalc/pull/80) e [P79 / B1 e Suzuki](https://github.com/jbarumby/scicalc/pull/79).

**Próximo gate:** reproduzir criticamente a estimativa KLN para F_X e a identidade dos coeficientes de Dirichlet, sem assumir que a cadeia completa está demonstrada só porque as constantes fecham numericamente.
