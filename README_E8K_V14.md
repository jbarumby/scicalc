# NUMERIA E8-K V14 — replay, testemunhas e controle independente

## O que muda com a especificação V13 recebida

A V14 implementa o schema **`NUMERIA.E8K.V13.FULL_POINT_REPLAY.1`** com `phase=E8K_V14`, preservando o significado original da especificação e acrescentando rastreabilidade:

1. Publica os hashes SHA-256 e tamanhos **dos dois B3 originais** associados à fonte B5. Não substitui seus bytes por B3 corrigidos do E7.
2. Publica os **dois extremos racionais exatos da soma de cada uma das 11 bandas**; a soma dos extremos é novamente conferida a partir de **todos os 15.754 registros JSONL**.
3. Quantifica o desvio B3/B5 **somente para os intervalos de SOMA FINITA**. A `nested_fraction` representa zero ou um comparativo de intervalo agregado (amostra de tamanho um), **não** a fração de nós individuais, que não é recuperável dos relatórios B3 disponíveis. Há `scope=aggregate_finite_sum_only` e `per_point_metrics=NOT_AVAILABLE_NO_B3_NODE_INTERVALS`.
4. Calcula independentemente 128 integrandos em 11 bandas, usando Pi/Machin, gamma/Euler–Maclaurin N=128, digama por série, log/atanh, seno/cosseno/sinc por Taylor com aritmética racional dirigida. Publica todos os pares de intervalos, índices, argumentos, sucesso e SHA do JSON de detalhes. **Não usa as rotinas B5 como oráculo.** O resto analítico da digama ainda depende de lemas clássicos enunciados.
5. Publica aproximações **arredondadas para fora** para intervalos e cotas; conserva os valores exatos em strings racionais. Toda a verificação usa `Fraction`, nunca decimais aproximados.
6. SHA256SUMS inclui arquivos JSONL, `MANIFESTO_SHA256.json`, `b5/RESULTADO_B5.json`, B3 originais, E7 original, códigos e resultados V12/V13/V14. Não há hash autorreferente do próprio `SHA256SUMS` nem do relatório posterior de conferência.
7. O conferidor V14 usa gates bloqueantes, falha fechada, recalcula os 128 oráculos, e verifica **5 mutações isoladas** com o gate-alvo aprovado na linha de base.

## Resultados e escopo

- V12: todos os nós gerados pela fonte B5 original de SHA fixado e confrontados com registros E7.
- V13: recomputação racional dos 11 limites de quadratura e das duas caudas, **sob premissas analíticas declaradas**, sem formalizá-las todas.
- V14: reconciliação independente de JSONL e somas + 128 oráculos racionais de integrandos.
- Lean: código formal V13 incluído, **compilação ainda não verificada**. O workflow `auditoria-e8k-v13` é separado e pode revelar falhas.
- **Não prova** a fórmula global de Weil, validade analítica de todos os majorantes, positividade para dimensões arbitrárias ou a Hipótese de Riemann.

## Executar no GitHub Actions

Copie **o conteúdo descompactado** para a raiz de `jbarumby/scicalc` em uma branch de auditoria. Os caminhos devem incluir `.github/workflows/auditoria-e8k-v14.yml`, E7 ZIP e `proveniencia_b3/*.json`. Faça commit/push. Abra `Actions > auditoria-e8k-v14-testemunhas > Run workflow`.

## Executar localmente

```bash
python VERIFICAR_E8K_V12.py --zip NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip --provenance-dir proveniencia_b3 --strict-provenance --out REPLAY_E8K_V12_INTEGRAL.json
python VERIFICAR_E8K_V13.py --zip NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip --out RESULTADO_E8K_V13.json
python GERAR_E8K_V14.py --root .
python GERAR_SHA256_V14.py --root .
python CONFERIR_E8K_V14.py --root . --selftest --out CONFERENCIA_E8K_V14.json
```

O último comando deve sair com exit 0 somente se todos os gates passarem. Não é necessário instalar bibliotecas Python de terceiros.

## Limitações documentadas da especificação

A especificação solicitou uma fração de pontos B3↔B5 aninhados; como não há intervalos B3 individuais exportados, o relatório não pode inventar esses percentuais. A comparação de intervalos de soma agregada foi disponibilizada com escopo explícito. A indicação do PDF de usar o conferidor V12 com schema V13 requer adaptação: a implementação entregue é `CONFERIR_E8K_V14.py`, própria para o schema atualizado.

Hash isolado representa **integridade relativa, não autenticação por terceiros**. Para autenticação, assinar o manifesto fora do próprio pacote.

## Entradas preservadas

A especificação PDF e o código do conferidor proposto foram conservados em `referencias/`, sem substituir a revisão nova. O CSV/JSONL de pontos e o E7 original estão incluídos no ZIP completo.
