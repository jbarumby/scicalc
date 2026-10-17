# NUMERIA-Q 05.21 — status de validação

O protótipo local foi relatado como tendo 8 testes aprovados e 3.200 casos exaustivos verificados nas instâncias pequenas. Esses números são registros da execução local informada; não foram reexecutados pelo GitHub nesta publicação.

## Distinção de evidências

- **Fato do projeto:** duas famílias implementadas no protótipo local: mintermos e ANF.
- **Resultado local relatado:** equivalência exaustiva nos primos pequenos incluídos no relatório local.
- **Não reproduzido aqui:** resultados quantitativos atribuídos à Anthropic/Claude.
- **Não demonstrado:** vantagem quântica, descoberta algorítmica inédita ou desempenho em hardware.

## Próxima etapa sugerida

Implementar uma terceira família baseada em divsteps reversível com histórico explícito e, depois, testar se as decisões podem ser reconstruídas a partir de registradores intermediários. Manter 05.19.1 e 05.20 congeladas e executar tudo em branch isolada.
