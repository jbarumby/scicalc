# P87 — Resultado matemático e reprodução no GitHub

**CI de cálculo confirmada:** [GitHub Actions #45046458986](https://github.com/jbarumby/scicalc/actions/runs/45046458986), **completed/success**, em 22/11/2026. Reexecução do HEAD com este relatório será monitorada separadamente.

**Artefato:** [numeria-p87-shan-lemma1-lemma6-ledger](https://github.com/jbarumby/scicalc/actions/runs/45046458986/artifacts/14764046381), contendo `p87-shan.json`, `p87-lemma6.json`, `p87-SHA256SUMS.txt`.

**SHA-256:**
- `p87-shan.json`: `55e1f60ffaf1bfdb1a46100aa3dce9a3c5ad85baff29505a220fc520a7f`.
- `p87-lemma6.json`: `c990336911f9b1de0496a23bbe6eca240426ba52deadb1a04127e366fede0d5e`.

## 1. Lema 1 — demonstração universal concluída por argumento de Shan/Yangjit

**Teorema (todos os N finitos e todos os espaçamentos reais):**
[
Q_i=sum_{j
e i}rac{delta_idelta_j}{(x_i-x_j)^2}lerac{pi^2}{3}.
]

A demonstração **completa** está no [README, §§1–2](README.md). Usa os Lemas 2.3–2.5 de Yangjit (§2.2, [arXiv:2203.14950](https://arxiv.org/pdf/2203.14950)), por método de Shan; reconstruímos as estimativas de convexidade e telescopagem **sem presumir o Lema 1 de Preissmann**. A constante é ótima no limite de grandes grades uniformes centradas.

A prova da P86 `Q_i<=4` foi aprimorada para `pi²/3≈3.289868` para **toda dimensão**, não apenas as 9 configurações finitas do certificado.

**Testes de reprodução (apenas apoio à prova escrita):**
- 9/9 testes aprovados em Python normal e 9/9 sob `-O`.
- 9 conjuntos de pontos, incluindo grade uniforme de 41 pontos, com normalização dos dois lados, ascensão dos descensos e comparação convexa em `Fraction`.
- Entradas adversariais com descensos e pontos repetidos; cada execução registra as etapas de correção.
- Para grandes espaçamentos normalizados, os certificados de exemplo usam a cota racional da cauda de Basel depois do termo 256, evitando um somatório de milhões de frações; **isso não é uma hipótese adicionada à demonstração**, apenas um majorante rigoroso para diagnosticar amostras.

## 2. Lema 6 — fórmula auditada, demonstração independente ainda aberta

O enunciado de três denominadores foi verificado contra [Yangjit, Lema 3.2, equação (3.6), p. 9](https://arxiv.org/pdf/2203.14950), que cita explicitamente Preissmann Lema 6 (1984).

**223 pares de índices** de conjuntos finitos passaram pelos testes em frações exatas, com margem estritamente positiva mesmo após substituir `pi` por um limite **inferior racional certificado** via arctan/Machin. **Isso não demonstra o Lema 6 para qualquer N**; a condição `Preissmann_Lemma6_all_N_proved_here=false` é obrigatória na CI.

## 3. Ponte espectral `m0` — álgebra fechada condicionalmente ao Lema 6

Usando o novo Lema 1, a estimativa `sigma=4` de Shan e a hipótese do Lema 6:

[
U^2le Vleft(rac{pi^4}{45}V+rac{2pi^2}{3}Uight).
]

A solução da desigualdade quadrática produz `U/V<=(π²/3)(1+sqrt(6/5))` e, por conseguinte,

[
|B|_2^2le pi^2left(1+rac{2}{3}sqrt{rac65}ight).
]

A identidade foi confirmada **exatamente no anel \(\mathbf Q[\sqrt{6/5}]\)**, sem pontos flutuantes. Esta parte é uma dedução algébrica válida **se** o Lema 6 for aceito, não uma prova independente do Lema 6.

## 4. Deduplicação S3

Os cinco uploads da nota S3 e dos DOCs são **iguais por SHA-256** aos que já foram examinados na P84. Não houve nova adoção da TABLE nem nova certificação de `B1`. A tabela S2 robusta e o replay Suzuki `18<=t<=79` preservam seu status de **certificação computacional finita, condicional aos insumos analíticos**.

## Conclusão e próxima prioridade P88

**P87 concluiu a cota diagonal \(\pi^2/3\) em dimensão arbitrária e a álgebra da combinação espectral condicional.** O **único elo faltante para uma rederivação independente completa do teorema de Preissmann nesse percurso é o Lema 6 uniforme** (a estimativa de três denominadores e seu termo negativo). Ele será o objetivo exclusivo da P88; os testes finitos são evidência de consistência, não substitutos da demonstração.

Ramaré–Zúñiga `S_epsilon` e a faixa `3e10<=X<3e12` também continuam dependentes de insumo analítico externo. Não há prova global da Hipótese de Riemann.
