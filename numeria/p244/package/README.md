# ES8-P244 — auditoria P242 e lemas P237–P241

Resultados efetivamente reproduzidos neste ambiente, sem alegar RH.

- `P244_REPLAY_CERT.json`: C/MPFR reexecutado a 1024 bits nos dois setores, saídas idênticas às originais, 754 pivôs e margens de Weyl.
- `P244_ANALYTIC_CHECK.json`: desigualdades racionais para o majorante E3, GL60, Binet e simbolização da cauda.
- `P244_NUMERICAL_SPOTCHECK.json`: 12 testes da montagem em float64; apenas diagnósticos.
- `p244_replay_fail_closed.py`: fluxo estrito que recusa adulterações de ν e reexecuta o cálculo C, a partir de fontes P237 autenticadas.
- `p244_analytic_gates.py`: verificação racional independente dos orçamentos matemáticos elementares.
- `test_p244.py`: 7 testes normais e otimizados; inclui adulteração adversarial.
- Scripts P242, rastreamentos MPFR originais e certificado P242, incluídos para reprodução.

Pré-requisitos para replay integral: Python 3, GCC, libMPFR.so.6, libGMP, os cinco arquivos P237 com hashes do relatório e os dois rastreamentos. Exemplo, estando na pasta de scripts:

```sh
python p244_analytic_gates.py --out P244_ANALYTIC_CHECK.json
python -m unittest -v test_p244
python -O -m unittest -q test_p244
python p244_replay_fail_closed.py \
  --data /caminho/arquivos-P237 \
  --traces /caminho/rastreamentos-P242 \
  --cert /caminho/P242_AUDIT_CERT.json \
  --work /tmp/p244-replay \
  --out P244_REPLAY_CERT.json
```

O verificador exato original P242 **não deve ser usado sozinho** para aceitar dados não confiáveis. O novo wrapper autentica e reexecuta as saídas primeiro. Mesmo assim, a correção do algoritmo C, dos lemas clássicos importados e da origem das integrais P237 não está formalizada em Lean.
