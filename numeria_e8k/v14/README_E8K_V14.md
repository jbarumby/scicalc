# NUMERIA E8-K V14 — replay, testemunhas e controle independente

Este diretório prepara a incorporação do pacote **E8K V14** recebido em `E8K_V14_GITHUB_COMPLETO.zip`.

## Escopo declarado pelo pacote

- Replay da fonte B5 com 15.754 registros (G11: 9.524; G22: 6.230).
- Recomputação racional das somas em 11 bandas e reconciliação dos registros JSONL.
- Oráculo racional independente para uma amostra de 128 integrandos.
- Conferidor fail-closed e autotestes de mutação.
- Formalização Lean associada à revisão V13.
- Relatórios de resultados, proveniência e hashes SHA-256.

## Arquivos de referência no pacote

O ZIP contém os scripts `VERIFICAR_E8K_V12.py`, `VERIFICAR_E8K_V13.py`, `GERAR_E8K_V14.py`, `GERAR_SHA256_V14.py`, `ORACULO_RACIONAL_V14.py` e `CONFERIR_E8K_V14.py`; relatórios JSON; `b5/*.jsonl`; evidências B3; a fonte Lean; o ZIP E7 original; e dois workflows GitHub Actions.

## Comandos locais

```bash
python VERIFICAR_E8K_V12.py --zip NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip --provenance-dir proveniencia_b3 --strict-provenance --out REPLAY_E8K_V12_INTEGRAL.json
python VERIFICAR_E8K_V13.py --zip NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip --out RESULTADO_E8K_V13.json
python GERAR_E8K_V14.py --root .
python GERAR_SHA256_V14.py --root .
python CONFERIR_E8K_V14.py --root . --selftest --out CONFERENCIA_E8K_V14.json
```

## Limites científicos

Os resultados de replay, hashes e testes de mutação verificam aspectos computacionais definidos pelo pacote. Não demonstram, por si sós, a fórmula global de Weil, todas as premissas analíticas, positividade para dimensões arbitrárias ou a Hipótese de Riemann. O relatório do pacote também declara que a compilação Lean não foi verificada.

## Integridade do pacote

O arquivo ZIP recebido foi verificado localmente quanto à integridade estrutural (CRC do ZIP). SHA-256 do ZIP recebido:

`043183f0313f7e62ac106402181327f9b7b82e3769d84d42ebc9e75ed6a12c3a`

Esse hash identifica o ZIP enviado nesta conversa; não é uma assinatura externa nem substitui a conferência de `SHA256SUMS.txt` depois que os arquivos forem extraídos e colocados nos caminhos finais do repositório.
