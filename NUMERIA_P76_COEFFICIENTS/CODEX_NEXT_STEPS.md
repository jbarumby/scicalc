# Continuação Codex: P76-B1

Comece na branch `numeria-p76-cs-v2-coefficients`, baseada no PR #75.

1. Confirme o CI do P76-A: 16 desigualdades `B₂,B₃` por Arb/FLINT, 192/256 bits; confirme os SHA dos artefatos.
2. Estude Chourasiya–Simonič v2, Teorema 1 e Eqs. (35)–(40). Derive as contas de `K₁,K₂,K₃,L₀,L₁,L₂,A₁,A₂` em Arb, com **rigor na quadratura da integral de K₁**. Não use `mpmath.quad` como certificador.
3. Refaça cada `B₁` do Teorema 1 até a Tabela 1 e prove que o decimal publicado é uma **cota superior**; se falhar, obtenha nova cota com arredondamento rigoroso para fora e recertifique toda a faixa `18≤t≤79` da P75.
4. Mantenha distintas a cadeia do preprint v2 e a cadeia apenas com artigos publicados (`audit/p75-journal-only-baseline`).
5. Preserve todos os JSONs, os hashes e as branches anteriores. Só abra PR subsequente se a auditoria incluir uma nova evidência verificável. Não declare prova global de RH.