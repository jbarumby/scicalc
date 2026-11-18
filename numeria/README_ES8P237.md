# ES8-P237 — certificado finito completo de 377 modos por paridade

**Status:** certificado rigoroso de positividade das duas matrizes finitas, condicionado à construção reduzida documentada no runner P237. **Não** é uma prova de positividade global da forma de Weil nem da Hipótese de Riemann.

## Resultado

Foram agregados os 1008 painéis por setor para `L = 161/200 = 0,805`, com `T_sharp = 504`, e montadas as duas matrizes simétricas completas, cada uma de dimensão 377. Cada arquivo de quadratura contém 71.253 entradas triangulares; a auditoria-sombra revalidou as 142.506 somas inteiras das duas paridades e vinculou os bytes das matrizes aos hashes das 41 faixas de painéis de cada setor.

O LDL intervalar direto das matrizes finais, mesmo a 2048 bits, permaneceu inconclusivo nos índices 60 (par) e 66 (ímpar). Isso é sobrelargura intervalar e **não** é evidência de pivô negativo. Um diagnóstico separado em ponto flutuante mostrou pivôs centrais positivos, mas não foi usado como prova. O certificado substituto abaixo usa somente aritmética Arb com arredondamento dirigido:

1. Para cada matriz final de bolas Arb `R`, toma-se a matriz central dyádica `C` formada pelos pontos médios das entradas.
2. Um LDL intervalar a 2048 bits certifica os 377 pivôs estritamente positivos de `C`, separadamente em cada paridade; logo, `C` é SPD.
3. A inversa de `C` é encerrada em Arb e obtém-se um limite superior rigoroso `ν` para `||C⁻¹||_F²`. A soma das larguras das entradas em cada linha fornece `δ ≥ ||R-C||∞`. Como a perturbação e `C` são simétricas, `||R-C||₂ ≤ δ`. Verifica-se `δ² ν < 1`, o que implica `δ < 1/||C⁻¹||_F ≤ λ_min(C)` e, pela desigualdade de Weyl, `λ_min(R) > 0` para toda matriz coberta pela caixa de intervalos.

| Setor | Pivôs LDL positivos do centro | `δ` (limite superior) | `δ² ν` (limite superior) | Conclusão para toda a caixa |
|---|---:|---:|---:|---|
| Par | 377/377 | aproximadamente `3,6612 × 10⁻⁴⁷` | aproximadamente `5,8392 × 10⁻⁵⁹` | SPD |
| Ímpar | 377/377 | aproximadamente `3,6576 × 10⁻⁴⁷` | aproximadamente `4,9288 × 10⁻⁶⁵` | SPD |

Os valores decimais são apenas apresentação; os enclosures Arb completos estão serializados no certificado. A matriz final inclui a quadratura GL60, o majorante analítico de erro adicionado uma única vez por entrada, o termo diagonal `β` e os termos de polo, segundo `p237_ldl_full_arb.py`. O replay testa de forma fail-closed a cobertura e a auditoria dos dois setores antes de montar as matrizes.

## Artefatos e reprodução

- `p237/data/P237_FULL_C_even.json.gz` e `P237_FULL_C_odd.json.gz`: quadraturas completas agregadas.
- `p237/data/P237_SHADOW_FULL_AUDIT.json`: auditoria independente das somas inteiras e hashes das fontes.
- `p237/data/P237_FINITE_LDL_2048.json.gz`: duas matrizes finais completas em caixas Arb; registra a tentativa de LDL direto inconclusiva, mas preserva todas as entradas para o replay alternativo.
- `p237/data/P237_CENTER_WEYL_2048.json.gz`: certificado completo dos 754 pivôs do centro e das duas transferências Weyl.
- Os scripts auxiliares em `p237/` reconstroem a matriz final a partir dos `FULL_C` e repetem o certificado.

Com Python e `python-flint` (Arb) instalados:

```bash
cd numeria/p237
python3 p237_ldl_full_arb.py \
  --folder data \
  --shadow-audit data/P237_SHADOW_FULL_AUDIT.json \
  --out /tmp/P237_FINITE_LDL_2048.replay.json.gz \
  --bits 2048
python3 p237_center_perturbation.py \
  --input /tmp/P237_FINITE_LDL_2048.replay.json.gz \
  --out /tmp/P237_CENTER_WEYL_2048.replay.json.gz \
  --bits 2048
```

O primeiro replay pode novamente marcar o LDL direto como inconclusivo em 60/66; isso é esperado. O segundo é o critério decisivo centro–Weyl.

## Proveniência

| Arquivo | SHA-256 |
|---|---|
| `P237_FULL_C_even.json.gz` | `b5419ac9af938915e646fa34e6990efb78c5f64aca62519720569320c8cd8c80` |
| `P237_FULL_C_odd.json.gz` | `0232b9ca923654ff89330efa0d1a8945a37085e9bbd4ca7ce218fdc7879f195c` |
| `P237_SHADOW_FULL_AUDIT.json` | `e51a28fe0180c1a6215be23b0f268d570bed3c748210dcca0a0d9ff3f59bea26` |
| `P237_FINITE_LDL_2048.json.gz` | `2695768a3cca11237dd63cb40c0b6c34427b390329e5ff2ec089e31111c7898a` |
| `P237_CENTER_WEYL_2048.json.gz` | `272caafee875c2d5301782ed56bcb4d895df98ee0551a2e878d99577f44156d3` |

## Limite matemático e próximo passo P207

O resultado é uma afirmação sobre duas compressões finitas de 377 modos no parâmetro fixo `L=0,805`. Ele não mostra que a matriz reduzida coincide, para todas as funções admissíveis, com a forma global de Weil; não limita a cauda projetada e todos os acoplamentos com o complemento infinito; não prova a uniformidade em `L`, dimensão ou suporte; e não resolve a identificação/normalização completa da fórmula explícita.

A nota `referencias/TAIL_PROOF_ENVIADO_2026-11-01.md` fornece uma estimativa de cauda arquimediana por entrada para uma base senoidal sob suas precondições, mas declara expressamente que não deriva a forma de Weil completa nem trata, por si só, o operador reduzido de Zhu, o termo de primos e a transferência global. `NUMERIA_E8K_FORMAL/TailBoundV16.lean` prova apenas uma condição algébrica racional do denominador, não essa estimativa analítica. A próxima ponte deve construir cotas verificáveis para o bloco fora dos 377 modos e seu acoplamento — por exemplo, via um complemento de Schur ou uma desigualdade uniforme de forma quadrática — e derivá-las da mesma normalização e dos mesmos dados analíticos usados em P237. Até essa ponte ser provada, a positividade global de Weil e a Hipótese de Riemann permanecem em aberto.
