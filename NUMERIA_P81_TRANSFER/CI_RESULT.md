# P81 — Resultado da auditoria da passagem Sε → ζSε → C₃/C₄

**Execução:** [GitHub Actions #45034330569](https://github.com/jbarumby/scicalc/actions/runs/45034330569) — job `conditional-transfer`: **completed/success** em 22/11/2026.

**Artefato:** [numeria-p81-conditional-lemma6-transfer](https://github.com/jbarumby/scicalc/actions/runs/45034330569/artifacts/14758604876) (ID `14758604876`). Inclui `p81-transfer.json`, com as 96 comparações racionais, parâmetros de fonte e margens exatas, e `p81-sha256sums.txt`.

## Gates e margem

- 16 linhas completas da P78/P80;
- 3 regiões de altura: `3×10¹²≤X<10¹⁶`, `10¹⁶≤X<10³³` e `X≥10³³`;
- 2 ramos por região: `ε₂=2δ/log X` e `ε₁=δ/log X`;
- **96/96** limites estritos `b4>(1+ε_j)U` confirmados por `fractions.Fraction`, sem float;
- menor folga positiva exata: `17/2800000` (aprox. `6.0714285714×10⁻⁶`);
- seis testes unitários passaram sob Python normal e `-O`, incluindo perturbações adversariais de b4.

## Correção da justificativa enviada

O PDF `p78_s2_passagem_completa (1).pdf` tem progresso real ao separar os dois ramos e identificar a região válida, mas as alegações de demonstração integral ultrapassam a evidência:

1. O controle de G(ε) baseado em `O(ε²)` sem uma cota uniforme não demonstra o sinal de G' em todo o intervalo.
2. A hipótese H1 não enuncia que a P80 já absorveu o fator `1+ε` na cota b4. Esse fator é **necessário** se usarmos somente o limite trivial `ζ(1+ε)≤1+1/ε`.
3. Não se pode simplesmente apagar um `exp(γε)>1` de um majorante. A fórmula correta de CS contém `(5X)^{-ε₁}=5^{-ε₁}e^{-δ}`, e `5^{-ε₁}≤1` justifica C4 **quando** a cota de `ζ(1+ε₁) S_{ε₁}` já foi estabelecida.
4. A identidade `L=log T0` está fixada pelo parâmetro do auditor; um ajuste mpmath concordante em aproximadamente 10⁻¹³ não é uma demonstração de identidade exata.

**Alternativa:** em cada zona, o teorema externo fornece `0≤S_ε≤U` para `ε≤ε_cap`, a P80 verifica `η₂=2δ/log_floor≤ε_cap` e `b4>(1+η₂)U`. Portanto, pela comparação integral para ζ, **ambos os ramos satisfazem** `ε_j ζ(1+ε_j)S_{ε_j}<b4`. As formas de C3 e C4 seguem por desigualdades algébricas do Lema 6 de CS v2, respeitando `X≥3×10¹²`.

## Distinção essencial

Esta auditoria é **rigorosa na aritmética racional**, mas **condicional nas hipóteses analíticas importadas**. Não prova independentemente os teoremas de Ramaré–Zúñiga sobre Sε, nem a desigualdade preliminar de F_X de KLN/Chourasiya–Simonič, nem o Lema 6 para `3×10¹⁰≤X<3×10¹²`. Não produz nova extensão t para Suzuki, e não prova a Hipótese de Riemann.

**Próxima prioridade:** verificar a desigualdade preliminar de `F_X` e a identidade `Λ_X(τ)=ζ(τ)S_{τ−1}(X)` diretamente com as definições de KLN, e auditar se todos os consumidores do Teorema 1 usam `X=hT≥3×10¹²`.
