# NUMERIA P82 — Resultado da auditoria F_X

**GitHub Actions:** [run #45035493798](https://github.com/jbarumby/scicalc/actions/runs/45035493798) — job `fx-audit` **completed/success**, sem falhas.

**Artefato:** [numeria-p82-fx-arithmetic-audit](https://github.com/jbarumby/scicalc/actions/runs/45035493798/artifacts/14759970358), ID `14759970358`: `p82-fx-audit.json` e `p82-sha256sums.txt`.

**SHA-256 do certificado JSON:** `34a3f74c7cca548d3d3c7e535384694673cb196393b369ab823221fcc18e1677`.

## Testes executados

- **9 testes unitários**, cada um aprovado em Python normal e `-O`.
- **10 identidades racionais de convolução truncada**, incluindo a contribuição `n=1` responsável pelo termo `−1` em `ζ(τ)S_{τ−1}(X)-1`.
- **6 reagrupamentos racionais da forma KLN**, preservando `πm0(n+1)` e produzindo `πm0∑λ²/n^{2σ−1}+(T+πm0)∑λ²/n^{2σ}`.
- **6 verificações de decomposição `X<n≤5X` e `n>5X`** com expoentes inteiros, além de derivação geral por monotonicidade de potências.
- **3 testes diagnósticos de integrais finitas** de séries de Dirichlet, calculados pela fórmula trigonométrica explícita com `mpmath` (não são prova do teorema analítico).

## Conclusão auditável

Para `X≥3×10¹²`, `σ=1+δ/log X`, `0<δ≤0.303` e `T≥0`, a cadeia usada por CS v2 leva à desigualdade

```
F_X(σ,0,T,1/2)
  ≤ log X * [ C3(X0,δ) + C4(X0,δ)*(T+πm0)/X ],

C3 = πm0*b4*exp(2γδ/log X0)/(2δ),
C4 = b4/(5δ*expδ) + 0.605*exp(-2δ)/log X0,
X0=3×10¹².
```

A conclusão é **condicional** aos seguintes teoremas externos:

1. KLN (2018), Lema 3.6 (3.14): desigualdade de valor médio das séries de Dirichlet (m0 exato).
2. KLN (2018), Lema 3.4 (3.8): cota de soma curta `0.605/X` para `X<n≤5X`, atribuída a Ramaré.
3. Ramaré–Zúñiga Alterman (2026), cotas para `S_ε(X)` auditadas na P80/P81.

**Diferença exata:** `Σλ²/n^τ = ζ(τ)S_{τ−1}(X) − 1` porque `ζM_X` tem coeficiente `1` em `n=1`. A fórmula de CS em que `Σλ²/n^τ≤ζS` continua correta; a P81 também é conservadora.

**Localização de domínio:** a prova do Corolário 1 de Chourasiya–Simonič usa `h=1`, `X=T` e `T≥3×10¹²`, dentro da faixa pública. O Lema 6 na íntegra afirma também `X≥3×10¹⁰`; a nossa versão pública **não** cobre a parte adicional.

Os PDFs fornecidos são bons documentos auxiliares. Sua contribuição matemática principal é a separação dos dois ramos e a monotonicidade de B1. Não adotar as alegações de "demonstração ponta a ponta" baseadas em ajuste numérico ou resto `O(ε²)` não majorado.

**Status final:** identidades e álgebra de F_X auditadas, cadeia condicional fechada, três entradas analíticas externas ainda abertas; sem nova extensão de Suzuki além de t=79, sem prova global de RH.

**Próxima prioridade:** auditar de forma independente o Lema 3.6/3.14 de KLN (ou procurar prova do limite Montgomery–Vaughan–Ramaré com constante m0) e a origem da cota local b3=0.605, antes de fortalecer a reivindicação de independência.
