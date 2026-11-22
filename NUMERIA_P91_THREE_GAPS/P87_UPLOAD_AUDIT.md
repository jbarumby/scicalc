# P91 — Auditoria complementar dos novos uploads P87 (22/11/2026)

**Natureza:** auditoria de proveniência e correção matemática. A prova P91 para três lacunas já estava fechada e não depende desses uploads; não é uma reinterpretação do Lema 5/6 geral.

## 1. SHA-256 verificado localmente

| Material recebido | SHA-256 | Conclusão |
|---|---|---|
| p87_teorema_preissmann (2).tex | ac31e3443e223128d585952951714551be9b378aa3f3224b3d639520155303a1 | TeX |
| p87_teorema_preissmann (3).tex | ac31e3443e223128d585952951714551be9b378aa3f3224b3d639520155303a1 | **idêntico bit a bit** ao (2) |
| p87_teorema_preissmann (1).pdf | a5f82234583b389d17a32268973b0b23c9f1e1d94c9f799843b69617c6b41dd2 | 9 páginas; o Lema 5.3 continua com dois erros na demonstração |
| p87_ledger (1).json | f593dc80743a5510c064aaa6ad3d4b2149f9950389c2294e16e425121b728da7 | mesmo hash do ledger já auditado; 24 gates, zero falhas reportadas |

O PDF possui hash diferente do PDF anteriormente revisado em P89; não supor identidade entre as versões PDF. O ledger marca expressamente a hipótese dos termos cruzados como importada de Preissmann, não provada aqui.

## 2. Lema 5.3: correção da representação de Laplace

No TeX e no PDF novos, para \(u=ts\), falta o fator \(1/t\) em:

\[
\boxed{t\psi_1(t)=\int_0^\infty \frac{u e^{-u}}{t(1-e^{-u/t})}\,du.}
\]

Há também uma inversão: se \(t_2>t_1\) então \(e^{-u/t_2}>e^{-u/t_1}\), portanto \(1-e^{-u/t_2}<1-e^{-u/t_1}\), e não o contrário. A redação entregue não demonstra corretamente a monotonicidade.

O reparo já documentado na [errata P89](../NUMERIA_P89_LEMMA5/P87_UPLOAD_ERRATA.md) é definir \(F_u(t)=t(1-e^{-u/t})\). Para \(u,t>0\),

\[
F_u'(t)=1-(1+u/t)e^{-u/t}>0
\]

pela desigualdade \(e^v>1+v\). Logo o integrando positivo \(u e^{-u}/F_u(t)\) decresce estritamente com \(t\), como desejado.

**Segundo reparo independente:** usando Tonelli,

\[
\Phi^*(r)=\sum_{k=0}^\infty\frac{r}{(1+kr)^2}
=\int_0^\infty u e^{-u}\frac r{1-e^{-ru}}\,du.
\]
A função \(r/(1-e^{-ru})\) tem derivada
\[
\frac{1-(1+ru)e^{-ru}}{(1-e^{-ru})^2}>0.
\]
Portanto \(\Phi^*(r)\) é estritamente crescente. O resultado de Bellman pode usar a monotonicidade reparada. Isso não prova a hipótese espectral cruzada, que continua importada.

## 3. Ledger B5: sentido da cota de pi

Comparação **exata** usando Fraction nos racionais do ledger enviado:

- o valor B5.bound_pi2_3 é **igual** a B0.pi_upper²/3;
- o valor B5.bound_pi2_3 **não é igual** a B0.pi_lower²/3.

A desigualdade com cota **superior** de pi é insuficiente para concluir a desigualdade com pi exato. Reexecutar o bloco B5 com **pi_lower** e exportar as margens exatas como Fraction, não somente float. Os 24/24 gates reportados não corrigem essa direção de prova. B3/B4/B6/B7 também são diagnósticos numéricos, não provas universais.

## 4. Relação com P91/P92

A auditoria **não altera** o resultado P91: três lacunas (dois pontos interiores), quaisquer pontos exteriores, 447 coeficientes positivos em nove triângulos, certificado JSON reproduzido com SHA-256:

74a71549ef55ebadfe7cebf8dac5c85f906b6486c9b97d4880836ce5a8d18516

A estratégia de supersolução/Bellman na P87 é conceitualmente interessante, mas atua sobre o funcional unilateral de energia, **diferente** do funcional de dois índices e dos termos interiores da P91. Não presumir que o Bellman de P87 resolve o caso \(m\ge4\).

**Prioridade matemática seguinte:** quatro lacunas interiores (\(b-a=4\)), sob análise dos cones de ordenação e da compensação variacional do Lema 5. Não alegar Lema 5/6 completo, constante espectral geral independente, Ramaré–Zúñiga pendente, nem RH global. Manter PR draft.
