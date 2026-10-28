# NUMERIA E8-K — resultados locais reproduzíveis

Data: 2026-10-28. Ambiente disponível nesta sessão: Python 3.13.

## Integridade
- ZIP E7 SHA-256: `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7` — confere.
- Fonte `fontes/AUDITOR_PONTOS_B5.py` SHA-256: `a1753d3010e34f29371a9706cd6bb05b3b0c8491ca3cafb44bf98b7cd53d3df9` — confere.
- Os 14 hashes de arquivos listados em `MANIFESTO_SHA256.json` conferem com o conteúdo do ZIP E8-K.
- `py_compile` dos scripts Python do pacote: passou.

## Execuções
1. `ci/verify_e7.py`: passou, ZIP e fonte B5 verificados.
2. `ci/lean_draft_status.py`: executou; detectou `sorry` nas linhas 171, 172, 197, 223, 227 e 232. O rascunho não foi compilado.
3. `VERIFICAR_E8K_V10_UNIFICADO_CORRIGIDO.py --zip ... --out ... --trials 100`: passou com status `CHECKS_EXECUTED_AND_PASSED_WITH_SCOPE_LIMITS`; 8/8 checagens AST, 16/16 casos racionais de borda, 900 operações racionais exatas, identidade de Machin aprovada, concordância formula gamma/B5 aprovada, 3/3 mutações negativas detectadas. Esta execução usou 100 trials, não os 1800 definidos no workflow.
4. `ci/v11_gate.py independent`: padrão conhecido confirmado, 1 PASS e 5 SKIP; não é aprovação completa.
5. `ci/v11_gate.py sourced`: padrão conhecido confirmado, 3 PASS e 3 SKIP; três gates continuam pendentes.
6. `pytest -q -x testes_property_e8k.py`: bloqueado na coleta porque o ambiente não possui o pacote `hypothesis` (`ModuleNotFoundError`). Não foi aprovado.

## Limitações e próximos passos
- Os testes property-based precisam ser repetidos após instalar `pytest` e `hypothesis`.
- Executar o V10 com os 1800 trials definidos no workflow.
- O ZIP E7 foi usado localmente nesta execução, mas não foi adicionado ao repositório.
- O workflow no pacote ainda não foi ativado na raiz de `.github/workflows/`; a integração no GitHub continua parcial e não deve ser mesclada.
- Os gates não certificam a quadratura contínua, a cauda global de Weil ou a Hipótese de Riemann.
