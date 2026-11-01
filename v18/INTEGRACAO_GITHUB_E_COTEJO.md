# V18: incorporacao e cotejo no GitHub

Esta branch recupera as fontes originais V17, o script V18 enviado pelo usuario,
a revisao V18 executavel, as testemunhas, os hashes fixados e um cotejador
independente com testes de adulteracao.

## Execucao

```bash
cd v18 && sha256sum -c SHA256SUMS_V18.txt && cd ..
python v18/COTEJAR_V17_V18.py
```

O cotejador exige 14/14 gates, 13/13 mutacoes, 10 ressonancias e
replay byte a byte do JSON V18. Executa controles negativos de
ausencia do E7, alteracao da normalizacao do primo, testemunha de
suporte e hash da fonte.

O workflow `auditoria-e8k-v18-testemunhas.yml` roda no mesmo commit e
publica o artefato de diagnostico, mesmo em caso de falha.

## Limites

V18 NAO recalcula o LDL da Fase 4/V16; as cotas analiticas de
quadratura/cauda e a aplicabilidade global de Guinand-Weil nao sao
certificadas por estes testes. A ressonancia da V17 (k=5) obtida
Simpson-vs-Simpson e desconsiderada como evidencia independente;
a V18 utiliza identidade analitica e compara por quadratura float
apenas como sanity check. RH nao demonstrada.
