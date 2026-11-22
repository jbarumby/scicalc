# P84 — comparação com a Seção 6 da nota S3 (sem replicação de trabalho)

**Fonte examinada:** `p78_s3_nota_adocao_repo_certified.pdf`, 6 páginas, especialmente o **Runbook de adoção**, Seção 6 (pp. 4–5 da numeração impressa). Esta comparação é uma revisão dos procedimentos; não é um certificado matemático adicional.

## Identidade de arquivos recentemente reenviados

Os hashes SHA-256 dos arquivos de entrada repetidos são os mesmos dos materiais S3 já recebidos e arquivados nas fases anteriores:

| Nome e versões locais | SHA-256 | Resultado |
|---|---|---|
| `p78_s3_doc1_certificado.json`, `... (1).json`, `...(1).json` | `093dd70b4fc6ddf0385094cd112699034f2fad2c36a3c51758c7f7c044652137` | idênticos |
| `p78_s3_doc2_certificado.json`, `...(1).json` | `ac356d2997e93945e8ce0751dd9e2f00054b45b9284a800d3ed690b57d8c6846` | idênticos |
| `p78_s3_adocao_table.json` e duas variantes com `(1)` | `35a8bd405b7a368587a842318448f8af3513b406e78d77882c8a39eb0d065fb2` | idênticos |
| `p78_s3_downstream_scan.json` e `...(1).json` | `f5c399a8a773af7e8817cc757ee43114d780d9b3df022d44e07e004052e91547` | idênticos |
| `p78_s3_adocao_table.py` e `...(1).py` | `1610f923b74da5e7b4253006e882e35f52c955084650331cde249d045c1b9479` | idênticos |
| `p78_s3_docs_execucao.log` e `...(1).log` | `a52d7c82537295adfd0be103321d1e54e7eeb39c1241d023827be1ef90fe9b34` | idênticos |

O `p78_s3_certificado_oficial (1).json` **não é byte a byte idêntico** ao primeiro JSON: SHA-256 `6c0be960a3824310aa73b04bbb2f5bf1eb54a688b845ab8391bb542b48fdef1a` contra `0bcc21e80cc1eb99f3085fce265861f83cc7bdd61ad075186cd6eb3002ad9d56`. Comparação estrutural: **a única diferença é `generated_utc`** (15:14:03Z vs. 15:05:08Z).

## Correspondência detalhada: Seção 6 versus P84

| Etapa da nota | Situação P84 (GitHub Actions #45038163289 e #45038444472) |
|---|---|
| 1. Copiar kit S2/S3 e workflow original | **Substituído por uma implementação própria auditável**, `NUMERIA_P84_S3_REPO/` e workflow `numeria-p84-s3-repo.yml`, preservando as entradas originais sem executá-las como autoridades |
| 2. Rodar `p78_s3_verifica_tudo.py --variant robust --dyadic 192,256,384` em modo REPO | **Não executado literalmente**: substituído por auditor independente `audit_s3.py`, três precisões, dados da tabela real e quadratura efetiva P79 |
| 3. Rodar `p78_s3_recert_docs.py` | **Substituído** por `build_docs.py` que gera `DOC1_repo_interval.json` e `DOC2_repo_interval.json` com endpoints Arb diádicos, não `float` envelopado |
| 4. Editar sete `b1` da TABLE publicada | **Equivalência funcional, implementação versionada:** `TABLE_S2_PUBLIC_ROBUST` existe no módulo `table.py`, e `TABLE_V2_S2_PUBLIC_ROBUST` no `cs_v2.py`, com a tabela original intacta |
| 5. Rodar `p78_s2_recertify.py --dyadic 192,256,384` pós-edição | **Executado por equivalente mais forte**: 16 linhas × 3 precisões no auditor Arb original adaptado da P79, `K1` por `acb.integral`, `K2` pela fórmula analítica, margens e hashes |
| 6. Rodar `p78_s2_downstream_scan.py` | **Substituído por exame real do checkout:** 717 arquivos, 36 referências detectadas, 0 potenciais consumidores sem classificação dentro das regras do scanner (antes eram 17 arquivos de `/tmp/kit_p78`) |
| 7. Regenerar `certificate.json` + SHA256 | **Artefatos novos, não sobrescritos:** relatórios de 192/256/384 bits, DOC1/DOC2, `p84-SHA256SUMS.txt`. O arquivo histórico `certificate.json` é preservado |
| 8. Acionar `NUMERIA-P78-S3-certificacao`, robust + G + minimal | **Não executado literalmente**: a CI P84 completa passou para a versão **robusta**, incluindo Suzuki até `t=79`. Não foram recertificadas as alegações extras de `G`, 28 valores ζ e limiares b4* provenientes do kit STANDALONE |

**Conclusão operacional:** repetir os comandos do kit original é desnecessário para as partes fechadas e pode dar um certificado enganoso se o marcador `REPO` for tomado como prova de quadratura. A P84 já fornece a grade intervalar e o replay robusto. **Manter os dois executáveis históricos apenas como evidência/diagnóstico.**

## Falha de escopo no script original

O `p78_s3_verifica_tudo.py` recebido define:
```python
repo = find_repo()
mode = 'REPO' if repo else 'STANDALONE'
```
Isso **somente detecta a existência do diretório** `NUMERIA_P76_COEFFICIENTS`. A implementação de `b1_upper_interval()` usa `alpha_b4_linear`, `beta_b4_linear` e um coeficiente `D0` previamente ajustados por `float`, inflados com `H_ENV=1e-12`. O próprio programa serializa endpoints com `float(lo)` e `float(hi)`. Portanto, suas saídas não são uma prova intervalar independente do modelo original, mesmo se o rótulo for `mode=REPO`.

O script `p78_s3_adocao_table.py` fornecido consolida valores e declara `gate_pos_adocao=PASS` com dados anteriores; **não modifica a TABLE nem executa a quadratura**. A proposta de inserir suas linhas literalmente na TABLE publicada sobrescreveria um insumo histórico. A P84 adota corretamente uma **variante versionada**.

## Escopo matemático que continua aberto

- Os relatórios P84 DOC1/DOC2 são computacionalmente rigorosos **nas quantidades efetivamente reconstruídas**; não encerram as reivindicações sobre `G'(ε)<0` em todo o domínio, os limiares `b4*` nem os alegados 28 gates `arb.zeta`.
- Os teoremas de Montgomery–Vaughan/Preissmann (valor médio) e Ramaré–Zúñiga (limite de `S_ε(X)`) continuam entradas analíticas **importadas**.
- A faixa `3×10¹⁰ ≤ X < 3×10¹²` permanece **fora da cobertura pública** dessa cadeia. A afirmação da nota de que testar os extremos basta não é um certificado de uniformidade no interior; em particular, `ε = cδ/\log X` **cresce** quando X diminui, e a margem `b4 − (1+ε)U`, mantendo U fixo, **diminui**, não aumenta.
- A P78 possui certificado com 16 linhas e a P84 refez os gates de B1 das linhas 1–9. Isso não prova de maneira independente todos os teoremas analíticos externos nessas linhas.
- A positividade até `t=79` é **finita e condicional**; nenhuma prova global da Hipótese de Riemann.

## Evidência

- [PR #84](https://github.com/jbarumby/scicalc/pull/84)
- [GitHub Actions #45038444472 (success)](https://github.com/jbarumby/scicalc/actions/runs/45038444472)
- [Relatório de certificação e downstream](CI_RESULT.md)
- [Artefato #14759814197 com os certificados anteriores da grade](https://github.com/jbarumby/scicalc/actions/runs/45038163289/artifacts/14759814197)

**Orientação:** não reexecutar nem integrar o mesmo kit como nova etapa. Deduplicar fontes, manter a P84 em PR rascunho e alocar esforço apenas aos gates ainda não comprovados.
