# NUMERIA–E8K — registro imutável da versão auditada

**Data de fixação:** 2026-11-01  
**Repositório:** `jbarumby/scicalc`  
**PR:** [#8](https://github.com/jbarumby/scicalc/pull/8) — aberta, não mesclada  
**Branch:** `recovery/e8k-v15-restore-v13-v14-2026-11-01`  
**Commit auditado (HEAD):** `fe6bdb7463172ea7107cd050410fe0f94196de29`  
**Tree SHA:** `4a50bd64248c01a49b09fc2c9bf540726e2402a2`

## Toolchain Lean

- Arquivo `NUMERIA_E8K_FORMAL/lean-toolchain`: `leanprover/lean4:v4.35.0`
- Lean: **4.35.0** (versão solicitada pelo toolchain).
- Mathlib: revisão Git `56a005855ea17e0b0958f73ec801ae140b0a1679`, entrada `v4.35.0) no manifesto.
- SHA-256 do manifesto `NUMERIA_E8K_FORMAL/lake-manifest.json`: o manifesto fixa as dependências transitivas por revisões Git; o blob SHA do GitHub é `ca9698950cd2408ab702f623d9d3081eb9f7267d`.
- SHA-256 do conteúdo Lean, conforme `SHA256SUMS.txt`: `b6397d5a3bd5b73d3516cd90d39cd478f93af8ae2234fa3e4cbad0ff7d255fb2`.
- SHA-256 do conteúdo `NUMERIA_E8K_FORMAL/lean-toolchain), conforme manifesto: `ca5f1f601d34c49079e38d7e7041eb0f2a9c227fbc1df0f404820069a052330c`.
- SHA-256 do conteúdo `NUMERIA_E8K_FORMAL/lakefile.toml), conforme manifesto: `f11b514659b8d3b8c933b38f306fe967f96818cd244cc21c517527d8d1f494b0`.

> O hash de um blob Git (SHA-1) não é o mesmo que o SHA-256 do conteúdo. Os valores acima são rotulados por tipo para evitar confusão.

## Execuções finais no mesmo commit

Todas as três execuções abaixo reportaram `completed / success` e apontam para o mesmo HEAD SHA.

| Auditoria | Run | Artefato(s) | Digest SHA-256 do ZIP |
|---|---:|---|---|
| V13 — cotas racionais + build Lean | [41322092455](https://github.com/jbarumby/scicalc/actions/runs/41322092455) | [e8k-v13-auditoria-racional](https://github.com/jbarumby/scicalc/actions/runs/41322092455), ID `13089038592`; [e8k-v13-lean-diagnostics](https://github.com/jbarumby/scicalc/actions/runs/41322092455), ID `13089479372` | `d5bb68255aaead38437122c68eb6a2a16ff9870bd66fe3144407f5cae871a01f`; `fe528e04b2c15fba41537dacb12abd286993b9254b3128c4088afad096dbc476` |
| V14 — replay/testemunhas | [41322092472](https://github.com/jbarumby/scicalc/actions/runs/41322092472) | `num-e8k-v14-testemunhas`, ID `13089279676` | `969d7cdbaefc0ae8217c38aba94db29d1f9e096949a5eef898d0a11726627e80` |
| V15 — bloco arquimediano 2×2 | [41322092440](https://github.com/jbarumby/scicalc/actions/runs/41322092440) | `num-e8k-v15-certificado-condicional`, ID `13089329496` | `5aeb01657adce6ec1857cf3cdca70c8d7cc4b30f507cdbd76595be11c8244e0d` |

### Resultado técnico reproduzível

- V13: jobs de cotas racionais e formalização Lean concluídos com sucesso; o log do job Lean registra `lake build` bem-sucedido (3138 jobs).
- V14: replay integral de 15.754 pontos, 11 somas de banda e 128 oráculos racionais, conforme os relatórios versionados.
- V15: 11/11 gates aprovados; manifesto com 27 entradas verificadas; 7/7 controles negativos detectados; o JSON foi reexecutado e comparado byte a byte no runner.

## Limites do registro

Esta é a fixação de uma **versão de software e de resultados de auditoria**, não uma certificação matemática global. A compilação Lean cobre o arquivo-alvo e as dependências resolvidas pelo build; os scripts Python conferem dados e consequências racionais sob premissas analíticas declaradas. Nenhum desses resultados demonstra a fórmula explícita global de Weil, a positividade da matriz de Weil em dimensão arbitrária ou a Hipótese de Riemann.

Para reproduzir: fazer checkout do commit acima, conferir `sha256sum -c SHA256SUMS.txt`, entrar em `NUMERIA_E8K_FORMAL` e executar `lake build`; depois executar os três workflows via GitHub Actions. Não editar este registro retroativamente: uma nova versão deve receber novo commit e novo registro.
