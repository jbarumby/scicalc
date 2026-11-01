# NUMERIA–E8K — comparação V17/V18 e plano de padronização

**Data:** 2026-11-01  
**Repositório:** `jbarumby/scicalc`  
**PR:** [#8](https://github.com/jbarumby/scicalc/pull/8) — manter aberta; não mesclar automaticamente  
**Branch-alvo:** `recovery/e8k-v15-restore-v13-v14-2026-11-01`

## 1. Resumo executivo

A V17 concentra a normalização e validação do relatório anterior. A V18 acrescenta verificações independentes de testemunhas, incluindo a substituição explícita de uma comparação circular Simpson-vs-Simpson por uma identidade analítica exata acompanhada de quadraturas numéricas independentes. O teste local V18 executado em 2026-11-01 reportou 14 PASS, 0 FAIL, 13/13 mutações detectadas e manifesto V18 com 5/5 arquivos conferidos. Esses números são resultados locais; não devem ser descritos como execução do GitHub Actions até que o workflow V18 seja executado no repositório.

## 2. Evolução V17 → V18

| Área | V17 | V18 / evolução observada | Próxima melhoria |
|---|---|---|---|
| Responsabilidade principal | Normalização do relatório V17 | Auditoria de testemunhas e consistência das evidências V17 | Separar lógica comum de normalização, hash, schema e controles negativos |
| Organização dos scripts | `AUDITAR_E8K_V17_NORMALIZACAO.py` | `AUDITAR_E8K_V18_ORIGINAL_USUARIO.py` e `AUDITAR_E8K_V18_TESTEMUNHAS.py` | Criar biblioteca comum e manter dois perfis explícitos de CLI, evitando um script monolítico |
| Testemunhas aritméticas | Testes de sanidade e relatório normalizado | Reprodução bit a bit de ln(2)/ln(3); comparação racional exata de log₂ < 9/10 < log₃ | Documentar origem, domínio e critério de aceitação de cada testemunha |
| Cancelamento de polos | Evidência anterior registrada | 8 amostras de X por k para identidade Cₖ + βCₚ = 0; 8 colunas | Manter teste de identidade separado dos testes numéricos |
| Ressonância | Evidência histórica V17 baseada em Simpson versus Simpson excluído | V18 identifica a limitação histórica e usa identidade analítica exata + 10 quadraturas float independentes | Preservar a evidência antiga como histórica, sem promovê-la a prova independente |
| Robustez | 9 gates estruturais, 6 negativos e 25 verificações de sanidade reportados pela V18 ao auditar a estrutura V17 | 13/13 mutações detectadas pelo self-test V18 | Expandir mutações para hash incorreto, campo ausente, fórmula alterada e artefato trocado |
| Integridade | Manifesto de três artefatos V17 ancorados por hash | Manifesto V18 e verificação de seis constantes SHA no relatório/código | Uma política de hashes comum para V17 e V18, sem sobrescrever o manifesto raiz |
| CI | Workflow V17 existe no pacote, mas precisa ser integrado e executado na branch-alvo | Workflow V18 existe no pacote, mas ainda não foi publicado/executado no branch-alvo no momento deste registro | Publicar artefatos e workflows, executar ambos no mesmo HEAD e comparar hashes |

## 3. Arquitetura recomendada para reduzir duplicação

Não recomendo fundir imediatamente todas as rotinas num único arquivo enorme. Recomendo uma biblioteca comum com dois perfis explícitos, preservando os comandos e relatórios históricos.

Estrutura-alvo:

```text
e8k_audit/
  __init__.py
  common.py          # leitura estrita, schema, canonicalização e erros
  integrity.py       # SHA-256, manifesto e ancoragem de proveniência
  mutations.py       # mutações negativas comuns e contagem de detecção
  witnesses.py       # testemunhas racionais/analíticas e quadraturas
  v17.py             # regras específicas de normalização V17
  v18.py             # regras específicas de testemunhas V18
v17/
  AUDITAR_E8K_V17_NORMALIZACAO.py  # wrapper de compatibilidade
v18/
  AUDITAR_E8K_V18_ORIGINAL_USUARIO.py
  AUDITAR_E8K_V18_TESTEMUNHAS.py
```

### Regras de migração

1. Primeiro extrair apenas utilitários puros compartilhados; não alterar critérios matemáticos nesta etapa.
2. Manter os pontos de entrada antigos como wrappers de compatibilidade e preservar os argumentos CLI.
3. Criar testes de regressão que comparem os relatórios JSON anteriores e posteriores, excluindo apenas campos comprovadamente voláteis.
4. Não remover a distinção entre modo original do usuário e modo testemunhas.
5. Fazer a migração em commit separado e exigir igualdade byte a byte quando o JSON for declarado determinístico.
6. Não declarar redução de duplicação concluída antes de medir linhas duplicadas e executar os testes.

## 4. Padronização de integridade e CI

- Cada versão deve ter um manifesto próprio, por exemplo `v17/SHA256SUMS_V17.txt` e `v18/SHA256SUMS_V18.txt`.
- O manifesto do pacote completo deve ter nome distinto, como `SHA256SUMS_PACOTE_V18.txt`; não substituir o `SHA256SUMS.txt` raiz.
- Verificar os manifestos com caminhos relativos à raiz correta, sem ambiguidade sobre o diretório de execução.
- Em CI, executar `sha256sum -c` antes dos scripts, executar duas vezes os geradores determinísticos e comparar os resultados com `cmp`.
- Usar `permissions: contents: read`, versões fixadas/seguras de ações e `if: always()` para publicar logs e relatórios de diagnóstico mesmo em falhas.
- Publicar como artefatos os JSONs, pareceres, manifestos e logs; registrar run ID, commit SHA, nome do artefato e digest do ZIP.
- Rodar V17 e V18 no mesmo commit e conferir a cadeia de hashes; um relatório local não substitui a execução remota.
- Fazer o workflow falhar se arquivos necessários estiverem ausentes, se o manifesto não conferir, se uma mutação esperada não for detectada ou se a saída determinística divergir.

## 5. Plano de execução ordenado

1. **Integrar os artefatos V17/V18 do pacote** na branch da PR #8, preservando os arquivos originais e registrando hashes.
2. **Publicar os workflows V17 e V18** com permissões somente leitura, dependências e caminhos explícitos.
3. **Executar V17 e V18 no GitHub Actions** e verificar jobs, artefatos e digests no mesmo HEAD.
4. **Adicionar testes de regressão da padronização** antes de mover código compartilhado.
5. **Extrair a biblioteca comum** em mudança isolada; provar que wrappers mantêm os resultados.
6. **Consolidar os pareceres** com este changelog como índice, sem substituir os pareceres versionados originais.
7. **Retomar a lacuna matemática V16**: a positividade algébrica de (1-(w/R)^2) sob hipóteses declaradas não demonstra, por si só, convergência absoluta nem a desigualdade integral da cauda. Recuperar o enunciado original da Proposição 6, formalizar hipóteses, majorante pontual e comparação integral.
8. **Manter separados os níveis de evidência**: prova Lean do que está efetivamente formalizado, verificações racionais, testes float e hipóteses analíticas externas.

## 6. Critérios de aceite

- [ ] Arquivos V17/V18 presentes na branch da PR #8 e ancorados por SHA-256.
- [ ] Workflows V17 e V18 passam em runner limpo no mesmo commit.
- [ ] Manifestos por versão e manifesto de pacote não conflitam.
- [ ] Duas execuções geram os mesmos relatórios quando a saída é declarada determinística.
- [ ] Testes de mutação cobrem corrupção de hash, alteração de fórmula, ausência de campo e troca de arquivo.
- [ ] Wrappers antigos preservam CLI e resultado; regressões documentadas.
- [ ] Parecer comparativo diferencia evidência exata, numérica e hipótese.
- [ ] Nenhuma conclusão afirma prova global de Weil ou da Hipótese de Riemann sem o argumento completo.

## 7. Limites científicos

Os testes aprovados demonstram que os verificadores implementados passaram nos casos e dados testados. Eles não demonstram automaticamente a validade de todas as premissas analíticas, a identidade global de Weil, a positividade em dimensão arbitrária ou a Hipótese de Riemann. A auditoria da cauda deve permanecer marcada como condicional até que a convergência e a desigualdade integral estejam provadas ou explicitamente assumidas.
