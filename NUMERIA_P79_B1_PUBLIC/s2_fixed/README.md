# P79-S2 — correção do kit de recertificação

**Originais intactos:** `NUMERIA_P79_B1_PUBLIC/prototype_input/`.

**Corrigidos:** `p78_s2_recertify_fixed.py`, `p78_s2_downstream_scan_fixed.py`.

## Bug confirmado

Em `--proposed`, as linhas 10–16 usam `b4>0.54` e passam pela razão numérica dos logaritmos. No modo sem patch, linhas 1–9 mantêm `b4=0.54`. A expressão original `log(B1new/B1old)/log(K3new/K3old)` torna-se `0/0`, gerando `ZeroDivisionError`. A versão corrigida detecta exatamente `bpub==B54`, marca o teste de expoente como não aplicável nessa linha e mantém o **gate Arb** independente.

Correção de precisão: os racionais `b4` são criados **dentro** de `ctx.workprec(prec)`; grava-se também o endpoint superior de B1 e a margem inferior como **diádicos exatos** (`man_exp`). A variante `robust` exige margem `>=1/1000` por frações exatas, não por `float`.

A `TABLE` original da P76 **não é modificada**; a P79 calcula os coeficientes novos num overlay independente. O modo sem patch é intencionalmente **inconclusivo/falha de gate** para as sete linhas antigas, não uma falha de execução.

## Reproduzir

A partir da raiz do repo:

```bash
python -m pip install python-flint==0.9.0 mpmath==1.3.0
python -m NUMERIA_P79_B1_PUBLIC.s2_fixed.p78_s2_recertify_fixed \
  --proposed NUMERIA_P79_B1_PUBLIC/prototype_input/p78_s2_table_patch.json \
  --variant robust --dyadic 192,256 \
  --out /tmp/p79-s2-robust.json
python -m NUMERIA_P79_B1_PUBLIC.s2_fixed.p78_s2_downstream_scan_fixed \
  --root . --saida /tmp/p79-s2-scan.json
```

A variante `minimal` também pode ser reavaliada; nenhuma variante de duas casas substitui automaticamente a tabela de seis casas da P79 cuja cobertura `18<=t<=79` passou com **512 bandas**, sob hipóteses analíticas externas.

O arquivo `p78_s2_table_patch.json` foi criado originalmente a partir de protótipos em ponto flutuante. Esta execução com Arb certifica os limites escolhidos *condicionalmente à cadeia analítica*; não autentica sozinha a fonte do Lema 6, e não prova RH.
