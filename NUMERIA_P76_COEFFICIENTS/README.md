# NUMERIA–RIEMANN P76 — Auditoria independente dos coeficientes de Chourasiya–Simonič

**Branch proposta:** `numeria-p76-cs-v2-coefficients`. **Base:** PR #75 (`numeria-p75-cs-v2-local-curvature`).

## Escopo e resultado

Esta etapa **não amplia** a janela finita de Suzuki; ela audita as dependências da P75. Verifica **as 16 linhas de B₂ e B₃ da Tabela 1**, de [Chourasiya–Simonič v2, arXiv:2507.15184](https://arxiv.org/html/2507.15184v2), diretamente a partir das fórmulas do **Teorema 1**, equação (10), e das constantes dos lemas 5, 7 e 8, equações (33) e (34).

Em particular, para `H=T₀=H₀=3·10¹²`, `h=1`, e cada linha `(σ₁,σ₂,d)` da Tabela 1:

- `A₃ = z/(π log 2)`, `A₄ = 3z/(2 log 2)`, `z=(3−2σ₁)/(2d)+1/log H`;
- `A₅ = C₆ z/(2π)+C₅/(2πd)`;
- `B₂_calc=A₃+A₄` e `B₃_calc=A₃ log H+A₅`.

O código reavalia `C₅` e `C₆` a partir de suas fórmulas **com Arb**; utiliza `arb.const_euler()`, `arb.pi()`, `arb.const_log2()`, a expressão exata de `m₀`, e converte cada decimal publicado em fração racional exata antes de convertê-lo em intervalo. Em cada linha exige `lower(B₂_pub − B₂_calc)>0` e `lower(B₃_pub − B₃_calc)>0`. Em caso contrário encerra com erro.

**Não confundir com a reconstituição de toda a Tabela 1:** a verificação de `B₁` continua em aberto, pois exige `K₁` (uma integral imprópria com erros controlados), `K₂`, `K₃`, `L₀`, `L₁`, `L₂`, `A₁` e `A₂` das equações (35)–(40) e do Teorema 1. Portanto, a auditoria **não autentica ainda a P75 até t=79** e não prova RH.

## Evidências

- `audit_b23.py`: cálculo intervalar independentemente das linhas de B₂/B₃ em 192 e 256 bits;
- `table.py`: transcrição literal da Tabela 1 de v2, incluindo a mudança **d=0,28 → d=0,29** no 13º intervalo (`σ₁=0,875`);
- `test_p76.py`: testes incluindo manipulação adversarial para reduzir artificialmente B₂ ou B₃;
- workflow `numeria-p76-cs-v2-coefficients.yml`: reproduz as 16 linhas em Arb e arquiva os JSONs dos dois níveis de precisão.

O ensaio exploratório local com `mpmath` (sem Arb disponível) indicou que os 16 valores de B₂ e B₃ ficam abaixo dos respectivos decimais publicados. **Esta indicação não é certificação rigorosa**; o gate é o arquivo de saída de Arb no CI.

## Execução

```bash
python -m pip install python-flint==0.9.0 mpmath==1.3.0
python -m unittest NUMERIA_P76_COEFFICIENTS.test_p76 -v
python -O -m unittest NUMERIA_P76_COEFFICIENTS.test_p76 -v
python -m NUMERIA_P76_COEFFICIENTS.audit_b23 --precision 192 --out /tmp/p76_b23_192.json
python -m NUMERIA_P76_COEFFICIENTS.audit_b23 --precision 256 --out /tmp/p76_b23_256.json
```

## Próximo gate P76-B1

1. Implementar enclosures intervalares para as integrais impróprias de `K₁`, equação (37), decompondo intervalo compacto e cauda analítica majorada; comparar método distinto e malhas refinadas.
2. Implementar `a₁(A₀)`, `a₂(A₀)` e `C₁…C₄` a partir do texto original, sem valores estimados ou arredondamento implícito.
3. Recalcular `A₁` e `A₂` e verificar **todas as 16 desigualdades B₁_published ≥ B₁_calculated** com margem Arb explícita, inclusive o caso `σ₁=1/2` que o próprio corolário trata separadamente.
4. Se algum coeficiente impresso não majorar a fórmula, construir uma tabela substituta com **arredondamento para cima** e refazer a P75, comparando até `t=79`. Não usar mera inflação de 0,1% como demonstração desse fato.
5. Reproduzir a cadeia *journal-only* da branch `audit/p75-journal-only-baseline` independentemente e não alterar/mesclar PRs #71–#75 automaticamente.

**Fontes:** [arXiv v2](https://arxiv.org/html/2507.15184v2) · [PR #75](https://github.com/jbarumby/scicalc/pull/75) · [auditoria paralela](https://github.com/jbarumby/scicalc/tree/audit/p75-journal-only-baseline).