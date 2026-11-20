# NUMERIA–RIEMANN ES8-P255 — auditoria dos enunciados de Alpöge

Esta etapa executa três prioridades: (1) conferir os cinco enunciados Lean de RH/GRH de Dirichlet contra os arquivos originais e a Mathlib fixada; (2) analisar criticamente a estratégia Stepanov–theta–Guinand–Weil divulgada na nota de quatro páginas; (3) comparar com P249–P254 e preservar suas referências.

## Reprodução sem rede

No diretório raiz do repositório:

```bash
python -m unittest NUMERIA_ES8P255.test_p255 -v
python -O -m unittest NUMERIA_ES8P255.test_p255 -v
python -m NUMERIA_ES8P255.audit_statements --out /tmp/P255_AUDIT.json
```

O diretório `fixtures/` contém apenas **três arquivos pequenos** do código-fonte original do espelho fixado por hash, já autenticados no NUMERIA [PR #58](https://github.com/jbarumby/scicalc/pull/58). O auditor compara hashes SHA-256 e cinco pares de nomes/tipos e seus pontos de entrada. Opcionalmente, `--mathlib /caminho/Mathlib/NumberTheory/LSeries/RiemannZeta.lean` confere também o hash de objeto Git e a definição de RH da revisão Mathlib `0df444a360eaa60ab8c11dca51a86af692955474`.

**Atenção:** verificação estática de declarações, não compilação Lean nem verificação independente do corpo da solução. O `sorryAx` é parte do desafio com provas omitidas, não da solução completa. A auditoria do MIT em 20/11/2026 19:45 UTC relata a compilação de 175.187 módulos, as verificações de axiomas e aceitação da exportação integral pelo nanoda; outros replays ainda estavam em andamento nesse registro. O Clay ainda mantinha o status `Unsolved`.

## Resultados e próximos passos

Leia [relatório matemático detalhado](RELATORIO_ES8P255.md) e [snapshot SHA de P249–P254](P255_GITHUB_SNAPSHOT.json).

O código identifica `allChi : UpperTheory.EveryCharacter` como ponto de entrada da não anulação para todos os caracteres e `rhOfAll : RiemannHypothesis` como a passagem final via equação funcional. A nova hipótese a examinar **não é uma identidade já conhecida** entre o determinante de Stepanov/`W(1_S)` da nota e as matrizes (Q_W) do NUMERIA; ainda é necessário demonstrar essa ponte. O próximo trabalho útil é decompor um único lema quantitativo do desenvolvimento Lean para uma demonstração legível e verificar suas hipóteses, sem confundir RH formal com um certificado independente do NUMERIA.

Fontes primárias: [anúncio Alpöge](https://alpo.ge/remark.pdf) · [reprodução MIT](https://math.mit.edu/~drew/grh/) · [PR #58](https://github.com/jbarumby/scicalc/pull/58) · [Clay](https://www.claymath.org/millennium/Riemann-Hypothesis/).
