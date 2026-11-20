# Log de execução — marco divsteps

**Data:** 2026-10-18  
**Branch:** `research/divsteps-reproduction-2026-10`

## Teste executado localmente

Comando:
```bash
python -m unittest discover -s /tmp/numeria_divsteps/tests -v
```

Resultado observado:
```
test_all_nonzero_residues_small_primes ... ok
test_out_of_range_non_coprime_inputs_rejected ... ok
test_zero_has_no_inverse ... ok

Ran 3 tests in 0.000s
OK
```

Interpretação: a referência clássica independente de inverso modular passou nos três testes. Isso não testa o circuito reversível de divsteps.

## Tentativa de executar o repositório upstream

Comando:
```bash
git clone --depth 1 https://github.com/andigu/uncompute.git /tmp/uncompute
```

Resultado:
```
fatal: unable to access 'https://github.com/andigu/uncompute.git/':
Could not resolve host: github.com
```

## Prova Lean

Verificação de ferramentas locais:
- `lean`: não encontrado
- `lake`: não encontrado

A prova upstream não foi executada. Para seguir, é necessário um ambiente com acesso de rede e Lean 4/Lake conforme a versão indicada pelo repositório upstream. O README da prova recomenda `cd lean/divsteps && lake exe cache get && lake build`.

## Conclusão

- [x] Referência clássica: 3 testes passaram.
- [ ] Avaliador oficial do circuito divsteps: não executado por bloqueio de DNS.
- [ ] Prova Lean: não executada; ferramentas ausentes.
- [ ] Comparação de recursos sob protocolo comum: pendente.
- [ ] Descoberta nova: não alegada.

Não houve submissão de trabalhos quânticos, alteração de `main` ou modificação das versões congeladas 05.19.1 e 05.20.
