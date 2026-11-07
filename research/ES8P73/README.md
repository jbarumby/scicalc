# NUMERIA–RIEMANN — ES8-P73 (replay Arb em 50 modos)

**Escopo:** 50 modos Legendre por partes, largura a=9/20, R=30, partição [14,22,14], incluindo primo 2. Sem gaussianas.

O objetivo é **reproduzir de novo**, com python-flint==0.9.0, 19 entradas selecionadas de uma matriz histórica de 1275 pares e verificar sua compatibilidade por aritmética racional exata. Os centros históricos não são usados na integração: apenas na comparação de saída.

## Arquivos
- replay_arb.py: integrador com três segmentos [0,5],[5,15],[15,30], teste de polos em ±i/2 e comparação com as caixas históricas.
- historical_sample.json: 19 intervalos selecionados, copiados do replay histórico ES8-P72; identifica SHA-256 do arquivo-fonte completo.
- test_replay.py: teste prévio de completude, raios, ordenação e cobertura.
- .github/workflows/numeria-es8p73.yml: workflow de push e disparo manual, com quatro jobs de integração independentes.

## Executar
Na pasta research/ES8P73:

    python -m pip install python-flint==0.9.0
    python -m unittest -v test_replay.py
    python replay_arb.py --historical historical_sample.json --shard 0 --out replay_0.json --bits 256

Repita shard 0,1,2,3 para todos os 19 pares. Cada execução recusa sobrescrever seu JSON.

O workflow anexa os resultados JSON de cada shard, quando estes são criados.

## Portas para aceitar resultados
1. Instalação e importação de python-flint confirmadas no runner.
2. Testes prévios passam; os 19 pares estão distribuídos nos quatro shards.
3. Testes de polo da digama em ±i/2 geram bolas não finitas.
4. Todas as integrais terminam com bolas finitas; a aritmética é intervalar de 256 bits.
5. Nenhuma bola nova é DISJOINT das históricas.
6. Os artefatos são baixados, autenticados e inspecionados antes de declarar a reprodução concluída.

**Limites:** concordância ou inclusão entre duas execuções Arb não prova, sozinha, que a fórmula de Weil foi implementada corretamente. Não é prova da Hipótese de Riemann. Os arquivos gerados só devem receber status de execução comprovada após conferir os logs e artefatos reais.
