# Instruções E8-K V15

**A V15 é incremental**: incorpora a V14 e acrescenta o certificado de sinal e a robustez da matriz arquimediana 2×2. Não substitui os arquivos históricos da V14.

## Arquivos do pacote

- `v15/VERIFICAR_E8K_V15_BLOCO.py` — verificador stdlib / `Fraction` (novo).
- `v15/RESULTADO_E8K_V15.json` — resultado determinístico da execução local (novo).
- `v15/PARECER_TECNICO_V15.md` — escopo, desigualdades e limites matemáticos (novo).
- `.github/workflows/auditoria-e8k-v15.yml` — GitHub Actions (novo).
- `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip` — original incluído no pacote completo.
- `SHA256SUMS.txt` — manifesto original de 27 itens V14 sem alteração.
- `v15/SHA256SUMS_V15.txt` — SHA de documentos V15 e âncoras dos dados anteriores.

## Execução

Extraia o ZIP completo para uma pasta e, dentro dela:

```bash
sha256sum -c SHA256SUMS.txt
python v15/VERIFICAR_E8K_V15_BLOCO.py --root . --out v15/RESULTADO_E8K_V15.json
sha256sum -c v15/SHA256SUMS_V15.txt
```

O manifesto V15 usa o JSON pré-gerado determinístico; deve corresponder ao resultado reexecutado. Não edite os JSONs manualmente.

No GitHub, copie o conteúdo extraído para a raiz da branch de auditoria, faça commit e execute `auditoria-e8k-v15-bloco-arquimediano` em Actions. A execução local foi realizada; a execução no GitHub não foi confirmada.

## O que a V15 não prova

Ela não prova a fórmula global de Weil, não prova RH, não formaliza em Lean as cotas analíticas e não verifica todas as matrizes do projeto. O resultado de negativo definido é exclusivamente do bloco arquimediano 2×2, condicionado às inclusões originais.