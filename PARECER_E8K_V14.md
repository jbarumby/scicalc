# Parecer técnico — NUMERIA E8-K V14

**Data da execução:** 28/10/2026 (sessão atual).
**Objeto:** incorporação da especificação V13 de replay integral com testemunhas, melhoria do conferidor V12 e início da auditoria V14.

## Resultado reproduzido

- Reexecução B5 com fonte de SHA-256 fixado: **15.754/15.754** nós (G11 9.524; G22 6.230); 11 faixas, proveniência B3 preservada, sem lacunas.
- Conferência independente dos JSONL: pesos contêm exatamente `integrand * passo`; malhas de ponto médio corretas; todas as somas racionais das bandas coincidem com o B5 e seus totais.
- Oráculo racional em **128 integrandos**, estratificados por faixa: **128/128 inclusões** (B5 contém o intervalo independente). Testemunhas completas em `oraculo/AMOSTRA_128.json`, SHA registrado no relatório V14.
- Recomputação V13 das cotas publicadas de quadratura/cauda: **10/10 comparações racionais**, ainda subordinadas às premissas analíticas que o próprio relatório enumera.
- Conferidor V14: **9/9 gates PASS**, inclusive **5/5 mutações isoladas detectadas**. Reexecução do conferidor sobre o pacote ZIP extraído produziu relatório idêntico, byte a byte.
- `SHA256SUMS.txt` conferiu **26 arquivos**; o ZIP E7 original permaneceu inalterado.

## Mapeamento das exigências da especificação

| Exigência V13 | Implementação V14 | Escopo real |
|---|---|---|
| Hash original de cada B3 | `original_b3_input_sha256` e bytes publicados | Hash confrontado diretamente com B5/E7 |
| Soma por faixa com testemunha | `bands[].sum = [lower, upper]` racional exato | Todos os 11 intervalos e 15.754 weighted somados |
| Desvio B3/B5 | `b5_b3_metrics` | **Intervalos agregados** das somas, não pontos individuais |
| Oráculo independente amostral | `independent_oracle_sample` | 128 integrandos, com séries racionais e arredondamento dirigido |
| Decimais outward | `rounding=outward` e campos `*_approx` | Exibidores validados contra racionais exatos |
| SHA256SUMS e warning | Manifesto com E7/B3/JSONL/resultados | Integridade relativa, não assinatura/autenticidade |
| Testes negativos | Verificador V14 `--selftest` | 5 alterações isoladas exigem falha do gate responsável |

## Compatibilidade e nuance sobre B3/B5

O relatório B5 publica intervalos das **somas finitas**, mas não limites B3 individuais para cada um dos 15.754 nós. Assim, o campo `nested_fraction` em V14 mede **um único par de intervalos agregados** (igual a 1 para G11 e G22); **não** representa 100% dos nós aninhados. Para ambos, a soma B3 está contida na soma B5, embora a afirmação inversa `b5_enclosed_in_b3` seja falsa. `disjoint_count = 0` e `max_gap = 0` nos dois agregados. Conclusões por nó requerem publicação de novos dados B3.

## Dependências matemáticas em aberto

O oráculo racional depende das identidades clássicas para seno, cosseno, log, pi e digama com seus restos. A conferência V13 usa majorantes de segundas derivadas e da cauda baseados em premissas analíticas declaradas. O workflow Lean V13 foi preservado, mas **não compilado neste ambiente**. Nenhuma dessas verificações demonstra a validade global da fórmula de Weil, positividade em todas as dimensões, ou a Hipótese de Riemann.

## Próximo gate científico

Formalizar ou demonstrar independentemente os majorantes analíticos da quadratura em cada faixa e a hipótese da Proposição 6 na cauda além de R=160. Tratar este trabalho como etapa separada da integridade/replay de dados.
