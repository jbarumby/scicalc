# NUMERIA P89 — Resultado auditado: caudas exteriores e errata da nota P87

**Workflow aprovado:** [GitHub Actions #45050636080](https://github.com/jbarumby/scicalc/actions/runs/45050636080), **completed/success**, 22/11/2026.

**Artefato:** [numeria-p89-exterior-lemma5-and-p87-errata](https://github.com/jbarumby/scicalc/actions/runs/45050636080/artifacts/14765093913), ID `14765093913`, com `p89-exterior.json` e `p89-SHA256SUMS.txt`.

SHA-256 do JSON: `6448c00e395e4891d443f6e2730caf1933dac52c5885edc6222a3e58f5cc99c0`.

## 1. Novo resultado geral: caudas exteriores controladas sem Preissmann Lema 5/6

Para quaisquer pontos reais distintos `x_1<...<x_N`, índices `a<b`, distância `h=x_b-x_a`, `δ_i` = intervalo mais próximo:

```
sum_(k<a) δ_k / [(x_a-x_k)² (x_b-x_k)²]
   <= δ_a sum_(j>=1) 1 / [(jδ_a)²(h+jδ_a)²]

sum_(k>b) δ_k / [(x_k-x_a)² (x_k-x_b)²]
   <= δ_b sum_(j>=1) 1 / [(jδ_b)²(h+jδ_b)²].
```

A prova universal completa está no [README.md](README.md), usando a desigualdade unilateral de Shan da P87 com o núcleo `φ_h(t)=1/[t²(t+h)²]`, positivo, decrescente e convexo. O argumento não usa nem a prova variacional original do Lema 5 nem o Lema 6 de Preissmann para *essas duas parcelas exteriores*.

Também obtivemos o resto superior racional
`δ Σ_{j>K}φ_h(jδ)≤1/[3(Kδ)³]`, para `K>=1`, por comparação integral.

**Divisão de sinais provada analiticamente:** na soma cruzada `C_ab=Σδ_k/[(x_k-x_a)(x_k-x_b)]`, as parcelas de índices externos são **não negativas**, enquanto as internas `a<k<b` são **não positivas**. O termo interno afiado é o gate não resolvido para `E_ab>=0` e, consequentemente, para uma prova independente do Lema 6 completo.

## 2. Verificações executadas

- **7 testes unitários**, todos aprovados em Python normal e **7** em `python -O`.
- **30 pares de índices** de configurações racionais explícitas (grades de até 31 pontos, aglomerados extremos, padrões geométricos), auditados por `fractions.Fraction`.
- Em cada par: igualdade exata de separação esquerda+centro+direita, controle racional rigoroso das duas caudas, conservação da soma cruzada e sinais exterior/interior.
- Os casos numéricos são regressão; **a prova de todo N e todos os reais** é o argumento de Shan no README, não uma extrapolação.

## 3. Seis uploads da nota P87: conteúdo e erros identificados

O usuário forneceu PDF (10 páginas), `.tex`, `.py`, `.yml`, `.html`, `SHA256SUMS.txt`. Os cinco hashes referenciados pelos arquivos presentes foram conferidos **localmente** e registrados no [P87_UPLOAD_ERRATA.md](P87_UPLOAD_ERRATA.md); `p87_ledger.json` e `p87_execucao.log` referidos no manifesto não estavam entre os seis arquivos desta rodada.

### Falha matemática no Lema 5.3, corrigida pela P89

A nota inclui uma representação de `t ψ₁(t)` sem o fator obrigatório `1/t` após a mudança `u=ts` e afirma um sinal de comparação exponencial invertido. A identidade correta é:

```
t ψ₁(t) = ∫_0^∞ [u exp(-u)]/[t (1-exp(-u/t))] du.
```

A monotonicidade desejada segue da derivada positiva
`d/dt [ t(1-exp(-u/t)) ] = 1-(1+u/t) exp(-u/t) > 0`, em virtude de `e^v>1+v`. A conclusão do Lema 5.3 é válida depois desta correção, mas a prova da nota original não era válida como redigida.

### Falha de direção no certificado B5, corrigida no parecer

O `p87_verificacao.py` original testa `Q < π_upper²/3` e `S_uni < π_upper²/6` com uma cota **superior** de π. Isso **não** certifica as desigualdades desejadas com π. A CI adequada deve usar `π_lower` rigorosa de Machin, e registrar margens `Fraction` exatas em JSON, sem reduzi-las apenas a `float`.

Os 24 gates do script original **não são todos provas intervalares ou racionais**: B3/B4/B6/B7 usam malhas e autovalores float, classificados como diagnóstico.

## 4. Posição científica e próximo gate

**P89 demonstrou um controle uniforme de ambas as caudas exteriores.** Não reproduziu inteiramente o Lema 5 na sua formulação variacional, nem a cota afiada do Lema 6. Permanecem externas:
- compensação `E_ab>=0` para quaisquer espaçamentos, sobretudo o trecho **interior** `a<k<b`;
- desigualdade ponderada uniforme de Preissmann como insumo da cota `F_X`;
- limites Ramaré–Zúñiga para `S_epsilon(X)` e faixa `3×10^10 <= X < 3×10^12`.

**Próxima fase P90:** reconstruir especificamente o controle misto **interior** de Lema 5/6, com termo corretivo `-3(δ_a+δ_b)/h⁴`, mantendo ambas as faixas exteriores já fechadas. Os resultados P84/P85 e Suzuki até `t=79` são preservados como **finitos e condicionais**.

PR em rascunho, sem merge automático. **Não há prova global da Hipótese de Riemann.**
