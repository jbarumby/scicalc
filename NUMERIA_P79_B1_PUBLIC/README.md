# NUMERIA–RIEMANN P79 — Recalcular B1 com b4 de fonte pública e testar Suzuki

**Branch**: `numeria-p79-b1-public-replay`, baseada na P78, sem sobrescrever P75–P78.

## Diagnóstico recebido

Os arquivos `prototype_input/p78_b1_sensitivity.py` e `prototype_input/p78_b1_sensitivity_results.json` foram preservados **como protótipos não certificadores**, sem edição. Indicam que as linhas 10–16 não podem continuar usando os coeficientes B1 impressos no preprint CS v2 após substituir `b4=0.54` pelos limites da P78: os aumentos variam aproximadamente de 1,64% a 12,15%. Não é uma testemunha negativa de Suzuki.

## Novo gate formal computacional

1. O script `audit_public.py` copia as fórmulas do auditor P77 (`K1,K2,K3,L0,L1,L2,A1,A2,B1`), **substituindo apenas as duas ocorrências de b4 em C3 e C4** por um parâmetro por linha extraído do certificado racional P78. `K1` permanece a quadratura com `acb.integral` e majorante explícito da cauda.
2. Lê e autentica os 16 `b4` de `NUMERIA_P78_RAMARE_PUBLIC/certificate.json`, conferindo linha, delta, limite racional e as nove linhas com 0,54. Recusa alterações incompatíveis.
3. Enclausura cada B1 com Arb, guarda o endpoint superior **diádico exato**, e produz a menor grade decimal de 6 casas acima do endpoint, acrescida de uma unidade da última casa para folga estrita. Nos nove primeiros intervalos, preserva o coeficiente publicado **apenas se** continuar sendo majorante estrito.
4. Compara a nova tabela entre execuções separadas de 192 e 256 bits.
5. Substitui a tabela P75 **apenas em memória**, reexecuta 500 zeros, 256 bandas e a cobertura de `18≤t≤79` em FLINT/Arb. O resultado `inconclusive` é válido como diagnóstico, nunca apresentado como positividade.

CI arquiva os JSONs integrais, intervalo de primeira falha se houver, hashes SHA-256 de código e de artefatos, e compara as tabelas entre duas precisões.

## Limitações matemáticas

A validade deste replay é **condicional** às desigualdades da fonte Ramaré–Zúñiga Alterman **v1** usadas no P78, às hipóteses que as transportam ao Lema 6, às fórmulas publicadas de Chourasiya–Simonič **v2** e a todos os outros ingredientes de Suzuki/P74/P75. A fonte Ramaré–Zúñiga já possui versões posteriores (v3 em julho de 2026), que **não** foram automaticamente incorporadas à cadeia P78 v1. As consequências para os coeficientes e as condições de fronteira exigem revisão matemática.

A P79 não prova RH global nem transfere a conclusão até 79 sem uma execução completa e sem revisar as hipóteses externas.

### Reproduzir

```bash
python -m pip install python-flint==0.9.0 mpmath==1.3.0
python -m unittest NUMERIA_P79_B1_PUBLIC.test_public -v
python -O -m unittest NUMERIA_P79_B1_PUBLIC.test_public -v
python -m NUMERIA_P79_B1_PUBLIC.audit_public --precision 192 --out p79-b1-192.json
python -m NUMERIA_P79_B1_PUBLIC.audit_public --precision 256 --out p79-b1-256.json --cover-out p79-suzuki-cover.json
```

**Referências:** [Ramaré–Zúñiga v1](https://arxiv.org/abs/2603.25961v1) · [Chourasiya–Simonič v2](https://arxiv.org/abs/2507.15184v2) · [P77](https://github.com/jbarumby/scicalc/pull/77) · [P78](https://github.com/jbarumby/scicalc/pull/78).
