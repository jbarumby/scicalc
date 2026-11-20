# NUMERIA–RIEMANN ES8-P251 — transferência de corte: Gate A

Continua a P249/P250 no corte espectral de 309 com o mesmo beta de 504. **Resultado rigoroso condicionado às identidades clássicas de Binet e ao certificado P250:** a contribuição positiva H da faixa [309,504] na entrada diagonal do modo zero satisfaz **H00 > 2.8884327859588365e-10**. Como a margem conservadora finita P238 é 1e-18, a simples transferência por Weyl `lambda_min(A309) >= 1e-18-||H||` **não fecha**. Isto **não prova** que o bloco reduzido seja indefinido.

O script `p251_transfer_gate.py` reexecuta integralmente P250, exige seu SHA-256 determinístico, calcula quatro intervalos de seno com arredondamento dirigido P245B e verifica em frações exatas a cota integral para H00. As funções do pacote P250 estão na branch-base deste PR.

## Execução
```bash
python -m NUMERIA_ES8P251.p251_transfer_gate
python -m unittest NUMERIA_ES8P251.test_p251_transfer_gate -v
python -O -m unittest NUMERIA_ES8P251.test_p251_transfer_gate -v
```

## Próximo gate

Montar ou pré-condicionar rigorosamente o bloco **259×259 por setor**, 618 painéis GL60, evitando assumir que a antiga margem de 377 modos passa ao corte 309. Uma exploração exploratória em float64 é inconclusiva a 1e-15, e não fornece prova de positividade nem negatividade. O certificado detalhado e o diagnóstico estão disponíveis no pacote de reprodução P251 Gate A entregue no ChatGPT. A identificação Weil–Zhu e RH permanecem em aberto.
