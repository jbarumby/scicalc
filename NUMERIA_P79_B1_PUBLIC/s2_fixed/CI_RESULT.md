# P79-S2 — correções de execução e recertificação do patch

**Run GitHub Actions:** [45029602798](https://github.com/jbarumby/scicalc/actions/runs/45029602798) — job `recertify` concluído com **success**, em 22/11/2026.

**Artefato:** [numeria-p79-s2-fixed-recert](https://github.com/jbarumby/scicalc/actions/runs/45029602798/artifacts/14757670084) (ID `14757670084`). Contém os certificados-diag­nóstico `s2-robust-certificate.json`, `s2-minimal-certificate.json`, `s2-original-table-diagnostic.json`, varredura downstream, logs e hashes.

## Três erros corrigidos

1. **ZeroDivisionError em `--proposed` ausente:** as primeiras nove linhas preservam `b4=0.54`, então `log(B1_pub/B1_0.54)/log(K3_pub/K3_0.54)=log(1)/log(1)=0/0`. O programa agora identifica o caso de parâmetros idênticos e desativa esse diagnóstico não aplicável, mantendo o gate Arb. Testado com 16 linhas.
2. **TypeError em `min(dict)`:** o resumo tentava obter `min()` dos dicionários de resultados por precisão. Agora compara **margens inferiores diádicas exatas** convertidas para `Fraction`; registra o mínimo como racional exato, além do `float` informativo.
3. **Precisão de entrada:** o `Fraction` para `arb` dos parâmetros b4 ocorria fora de `ctx.workprec(prec)`. Agora os racionais são transformados em Arb **dentro** do contexto de 192/256 bits.

O patch passou a conferir `b4_public` e `delta` contra o certificado de fonte pública da P78 e a testar o alvo `1/1000` da variante robusta com frações racionais exatas, não com floats.

## Gates executados

- **Pré-adoção robusta, 192 e 256 bits:** sete linhas (10–16), todas com `(B1_novo-B1_Arb_upper).lower()>0`; todos os limites robustos têm margem inferior racional `>=0.001`.
- **Pré-adoção minimal, 192 bits e `python -O`:** sete linhas, todas com margens inferiores estritamente positivas.
- **Pós-adoção sobre a tabela histórica não modificada, 192 bits:** execução **sem erro de programação**, com resultado matemático esperado: linhas 1–9 passam, linhas **10–16 falham** com `S2 GATE FAIL: linhas [10,11,12,13,14,15,16]`.
- **Varredura downstream:** 323 arquivos `.py/.json` inspecionados. Os resultados são **heurísticos**, não uma instrução de editar automaticamente os consumidores.

### Variante robusta (2 decimais)

| linha | b4 público | B1 robusto | margem Arb inferior aproximada |
|---|---:|---:|---:|
| 10 | 0.557 | 17.49 | 0.005441334 |
| 11 | 0.557 | 19.67 | 0.010748839 |
| 12 | 0.558 | 22.25 | 0.010384122 |
| 13 | 0.559 | 25.48 | 0.007350429 |
| 14 | 0.560 | 29.89 | 0.008783929 |
| 15 | 0.561 | 36.59 | 0.003831118 |
| 16 | 0.606 | 51.65 | 0.001904444 |

Estas cotas são superiores às da P79 em seis casas decimais e portanto **podem ser usadas como majorantes alternativos**, condicionados à mesma cadeia de hipóteses analíticas. Porém a cobertura finita P79 até `t=79` **foi executada com a tabela de seis casas, não com a variante S2 robusta**. Não transferir esse resultado automaticamente.

Os arquivos originais fornecidos foram preservados em `prototype_input/`, sem alterar a `TABLE` original P76 nem os certificados anteriores. O projeto continua condicionado à justificação externa do Lema 6 e suas fontes; nenhuma conclusão prova RH global.
