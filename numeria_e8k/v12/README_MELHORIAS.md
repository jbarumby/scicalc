# NUMERIA E8-K V12 — conferidor de evidências, V2

## Conclusão da execução local (28/10/2026)

O relatório fornecido afirma **8/8 gates PASS**, incluindo **6/6 mutações isoladas detectadas**, e resultado JSON idêntico no modo `--root`. Esses resultados são alegações registradas no pacote e devem ser reproduzidos pelo CI após todos os artefatos necessários serem incorporados.

## O que a V2 confere

- Integridade dos hashes declarados em SHA256SUMS.
- Integridade do ZIP E7, código B5 e manifesto interno.
- Proveniência dos dois B3 originais.
- Estrutura de 15.754 pontos, 11 faixas, intervalos ordenados, limites e contiguidades.
- JSONL, malhas de pontos médios, intervalos e produtos ponderados por inclusão dirigida.
- Somas exatas por banda e reconstrução das somas finitas e compostas.
- Testes negativos com mutações isoladas, exigindo que o gate-alvo previamente PASS seja o que detecta a mutação.

## Limites científicos

A V2 verifica consistência de relatórios contra evidências exportadas do mesmo B5; não prova que os integrandos estejam matematicamente corretos, nem que as cotas globais da quadratura ou a cauda integral para r>160 sejam válidas. Hashes sem assinatura digital não autenticam a origem perante terceiros.

**Próxima etapa: V13**, fundamentar com desigualdades analíticas independentes as cotas de quadratura e cauda.