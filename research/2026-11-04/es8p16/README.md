# NUMERIA–RIEMANN · ES8-P16 · Auditoria de extensão infinita

**Data:** 2026-11-04. **Estado:** ferramentas de análise finita auditadas por frações exatas; Lean enviado para compilação no CI. Não é prova de RH.

## Origem

Analisamos `deepseek_latex_20261104_525471.tex` (SHA-256 `3e3cec45b5a8d2fce3a1f38c400edd72675d60216bfe88ee5d650142189f8792`) confrontando suas afirmações com a fonte Lean exata de ES8-P14 (`ES8P14Concrete.lean`). A fonte foi recuperada do pacote ES8-P14, e seu hash é mantido no JSON.

## Resultados novos

- Centros e envelopes de erro são aninhados: 4→6 e 6→8. Os extremos `lower` e `upper` também coincidem exatamente nos cantos correspondentes (208 comparações racionais no total).
- O aninhamento dessas **três** seções não determina uma extensão infinita única. O módulo Lean fornece duas extensões distintas com o mesmo prefixo 8×8.
- Pivôs LDL todos iguais a 1 não implicam `λ_min≥1/20`: `A_t = [[1,t],[t,t²+1]]`, fatorado por `U_tᵀ U_t`, no vetor `(-t,1)` tem quociente de Rayleigh `1/(t²+1)`. Para `t=10`, `1/101<1/20`.
- Somar o shift `3/50` aos *pivôs individuais* da fatoração LDL é inválido: o segundo pivô muda de maneira não aditiva (ver JSON exato).
- A projeção `B` amplifica o orçamento da forma pelo fator de norma ao quadrado: se o erro bruto é ≤δ `||y||²` e `y=Bx`, para controlar o erro projetado por `δ κ ||x||²` é preciso verificar `||Bx||²≤κ||x||²`. Sem hipótese sobre `B`, o orçamento bruto é insuficiente.
- Para as seções finitas de tamanho 4, 6, 8, as somas quadráticas dos envelopes de erro (norma de Frobenius ao quadrado) são inferiores a `1/10000`. **Não** é um limite de Hilbert–Schmidt em dimensões infinitas; essa seria uma nova hipótese a provar.

## Arquivos

- `ES8P16.lean`: provas propostas de nestedness finita, não unicidade da extensão, contraexemplo espectral e propagação de erros por projeção.
- `audit_es8p16.py`: auditoria independente das 12 matrizes de `ES8P14Concrete.lean` por `fractions.Fraction`, usando a fonte original preservada em `sources/`.
- `test_es8p16.py`: testes e mutações adversariais.
- `RESULTADOS_ES8P16.json`: resultados exatos; nenhum cálculo `float` determina um resultado.
- `RELATORIO_ES8P16.md`: relatório técnico e correções necessárias ao LaTeX.
- `PROVENIENCIA_ORIGINAL.sha256`: SHA-256 do documento recebido.
- `numeria-es8p16-lean.yml`: workflow Lean 4.35.0 + Mathlib + auditoria Python.

## Executar

```sh
python audit_es8p16.py --source sources/ES8P14Concrete.lean --out RESULTADOS_ES8P16.json
python -m unittest -v test_es8p16.py
```

A prova Lean deve ser compilada no ambiente do repositório (Lean 4.35.0+Mathlib). **Até existir resultado remoto de sucesso, não a classificar como mecanicamente verificada.**

## Próxima hipótese verificável

Construir `B_N` exata e estimar `sup_N ||B_N||²`, somada a uma cota de erro total que inclua projeção, interpolação, avaliações nodais, primos e cauda. Para uma extensão `M_∞` de fato bem definida, provar `sup_N ||C_N||` ou closability da forma. Para a margem `1/20`, requer-se um argumento estrutural que forneça `C_N ≥ 3/50 I` em **todas as dimensões**, não apenas pivôs positivos.
