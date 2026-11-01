# NUMERIA E8-K V18 — Instalação no GitHub

**Arquivos contidos**

- `v17/RESULTADO_E8K_V17.json`, `v17/PARECER_TECNICO_V17.md`, `v17/AUDITAR_E8K_V17_NORMALIZACAO.py` — originais preservados.
- `v18/AUDITAR_E8K_V18_ORIGINAL_USUARIO.py` — cópia literal do novo código enviado pelo usuário.
- `v18/AUDITAR_E8K_V18_TESTEMUNHAS.py` — revisão executável; revisa o caso ressonante `k=5`, valida o SHA E7, endurece `--strict` e estabiliza o JSON.
- `v18/RESULTADO_E8K_V18.json` — relatório reproduzido em modo estrito; 14 PASS.
- `v18/SHA256SUMS_V18.txt` — cinco hashes de arquivos (caminhos relativos à pasta `v18`).
- `v18/PARECER_TECNICO_V18.md` — escopo, resultados e lacunas.
- `.github/workflows/auditoria-e8k-v18-testemunhas.yml` — CI independente (não altera os workflows V14/V17).
- `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip` — somente no pacote completo.

Extraia o ZIP na raiz da sua **branch de auditoria** do `scicalc` sem alterar caminhos. Para o ZIP incremental, mantenha na raiz o E7 original com SHA-256 `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7`.

Execute:

```bash
python v18/AUDITAR_E8K_V18_TESTEMUNHAS.py \
  --v17-report v17/RESULTADO_E8K_V17.json \
  --v17-parecer v17/PARECER_TECNICO_V17.md \
  --v17-source v17/AUDITAR_E8K_V17_NORMALIZACAO.py \
  --e7 NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip \
  --out v18/RESULTADO_E8K_V18.json --strict --emit-sums
cd v18 && sha256sum -c SHA256SUMS_V18.txt
```

Na aba Actions, execute `auditoria-e8k-v18-testemunhas`. Os resultados no ZIP foram obtidos **localmente**, não no GitHub remoto.

Se o E7 falhar ao baixar, use o pacote completo. Ele já inclui o E7 dentro da raiz, sem necessidade de baixar o ZIP E7 separadamente.

**Cuidado:** o arquivo Lean das fases anteriores não é substituído nem recompilado nesta fase. `PASS` aqui não equivale a prova da fórmula de Weil ou da Hipótese de Riemann.