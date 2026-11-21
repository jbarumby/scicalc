# NUMERIA P75 — execução Arb com 256 faixas (replay completo)

**GitHub Actions:** https://github.com/jbarumby/scicalc/actions/runs/44947512003  
**Artefato completo JSON:** https://github.com/jbarumby/scicalc/actions/runs/44947512003/artifacts/14723274201  
**SHA-256 de `p75-cs2-79.json`:** `292cd84348575ef0f1328d58d82a0334f82e526ab9f7c38e52d05ff6616ea1cb`  
**Data da execução:** 2026-11-21

O workflow reproduziu os dez testes em Python normal e `-O`, verificou os hashes do núcleo e executou a cobertura **determinística** com `python-flint==0.9.0` a 256 bits:

| Propriedade | Resultado |
| --- | --- |
| `begin_t`, `end_t` | 18, 79 |
| Zeros isolados e auditados | 500 |
| Faixas de densidade | 256 |
| Malha de entrada | 1/64 |
| Status matemático do programa | `conditional_finite_positive` |
| Células aceitas | 4.284 |
| Células tentadas | 4.664 |
| Avaliações de extremos | 4.285 |
| Máxima profundidade de subdivisão | 3 |
| Mínima margem inferior | `46075230549864325426580486543243240430512775527604831320105092563618417550869 * 2^(-273)` |
| Mínima margem em decimal | > 0,00000303583875460407 |
| Célula limitante | [2643/64, 661/16] |
| Célula inconclusiva | nenhuma |
| SHA-256 das células aceitas | `73d2a88ff77a0ee345e3d7dac9189e7ea725ec6c312ff9eb87189cedb1e45722` |
| SHA-256 dos intervalos de zeros auditados | `482002a69ea975efa119d9a8346a89303d9825697031598fa384c5e8e01d2896` |
| Sobreposição com P71 | `true`: t=18 abaixo de log(10^8) |

**Atenção às distinções:**

- Uma tentativa anterior com **128 faixas** concluiu `inconclusive` perto de `t=78.62913`. A primeira célula sem margem positiva é uma limitação do majorante, não uma observação de violação de RH. A tentativa de 256 faixas fechou as células até 79.
- O verificador fornece uma implicação rigorosa **se** as hipóteses de Suzuki, Platt–Trudgian, Bellotti–Wong e **Chourasiya–Simonič arXiv:2507.15184v2** forem aplicadas exatamente como declarado. O último ainda é **preprint**: o audit técnico externo de suas condições e constantes permanece pendente.
- O programa conserva a procedência dos 500 zeros por isolamento via Arb, sinais de Hardy Z e contagem de zeros; a auditoria científica independente da cadeia completa continua necessária.
- Uma faixa finita de positividade não demonstra a Hipótese de Riemann global.

## Reexecução

Da raiz do repositório:

```bash
python -m pip install python-flint==0.9.0 mpmath==1.3.0
sha256sum -c NUMERIA_P75_CS2/SHA256SUMS
python -m unittest NUMERIA_P75_CS2.test_math NUMERIA_P75_CS2.test_flint -v
python -O -m unittest NUMERIA_P75_CS2.test_math NUMERIA_P75_CS2.test_flint -v
python -m NUMERIA_P75_CS2.certify --start 18 --end 79 \
    --zeros 500 --precision 256 --segments 256 --mesh 1/64 \
    --out /tmp/p75-cs2-79.json
sha256sum /tmp/p75-cs2-79.json
```

**Próxima prioridade:** reprodução em um segundo ambiente, prova formal dos fatores e da integral por faixas, e comparação da cadeia do preprint com uma cadeia exclusivamente de periódicos; só então buscar extensão além de t=79.
