# NUMERIA E8-K V14 — replay, testemunhas e controle independente

## Pacote de origem

Pacote recebido nesta conversa: `NUMERIA_E8K_V14_TODOS_ARQUIVOS.zip`.

- SHA-256 do ZIP recebido: `794ac407aa2278639fd2182a65f1ea36f02fa37044676ba90b8f71cdc5c29d12`
- Integridade estrutural ZIP: passou (`testzip` sem erro).
- Entradas: 29.
- Tamanho descompactado total: aproximadamente 10,0 MB.

O manifesto original do pacote foi preservado em [SHA256SUMS_ORIGINAL_PACKAGE.txt](SHA256SUMS_ORIGINAL_PACKAGE.txt). Não confundir a publicação desse manifesto com o upload dos próprios arquivos listados.

## O que implementa

1. Replay de 15.754 registros B5, com hashes fixados e proveniência.
2. Recomposição racional das somas exatas de 11 bandas e reconciliação dos JSONL.
3. Métricas B3/B5 explicitamente limitadas às somas finitas agregadas; não inventa métricas por ponto quando não há intervalos B3 individuais.
4. Oráculo racional independente para 128 integrandos em 11 bandas.
5. Arredondamento para fora e verificação por aritmética racional.
6. Conferidor fail-closed com gates e testes de mutação.
7. Código Lean V13 incluído no pacote, mas a compilação Lean permanece não verificada.

## Validação local executada sobre o ZIP completo

- `py_compile` dos sete scripts Python: passou.
- `VERIFICAR_E8K_V12.py`: código de saída 0.
- `VERIFICAR_E8K_V13.py`: código de saída 0.
- `GERAR_E8K_V14.py`: código de saída 0.
- `CONFERIR_E8K_V14.py --selftest`: **9/9 gates PASS**, zero FAIL.
- `sha256sum -c SHA256SUMS.txt`: todos os 27 itens listados passaram.

Esses testes foram executados na cópia local descompactada, não no GitHub Actions. Os arquivos completos ainda precisam ser enviados à raiz do repositório preservando exatamente os caminhos e bytes; consulte [STATUS_INCORPORACAO_V14.md](STATUS_INCORPORACAO_V14.md).

## Limites científicos

Esses resultados verificam apenas as propriedades computacionais e as premissas declaradas no pacote. Não comprovam a fórmula global de Weil, todos os majorantes analíticos, positividade em dimensões arbitrárias nem a Hipótese de Riemann.
