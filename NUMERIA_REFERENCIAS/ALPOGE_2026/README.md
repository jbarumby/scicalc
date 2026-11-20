# Reprodução da formalização de Riemann de Alpöge e Claude

Levent Alpöge publicou em 20 de novembro de 2026 um anúncio e a formalização Lean de uma demonstração da Hipótese de Riemann e da sua generalização para funções L de Dirichlet. A verificação independente de Andrew Sutherland no MIT registrou compilação, auditoria de axiomas e aceitação da exportação integral pelo verificador Rust `nanoda`. Este relatório reproduz a autenticação do código-fonte e a conferência dos tipos dos cinco teoremas no desafio. [Anúncio](https://alpo.ge/remark.pdf) · [Verificação do MIT](https://math.mit.edu/~drew/grh/)

## Reprodução local

Baixamos as 14 partes da página de Alpöge. A concatenação tem **239.275.554 bytes** e SHA-256 `9631b4f020199e580ba86178cb56c6e9b09f690c2369a80dc4334aa93f49c7fc`, igual ao publicado pelo MIT e pelo espelho Deicyde. Percorremos o arquivo sem executar o código e contamos **175.187 arquivos Lean**, com **17.374.215 linhas**. A busca textual encontrou zero ocorrências de `sorry`, `admit`, comandos `axiom` e `native_decide`. Encontrou uma ocorrência de `addDecl` em metaprogramação. A dependência transitiva de axiomas foi verificada no ambiente compilado pelo MIT. [Resultado local](RESULTADO_LOCAL.json) · [Proveniência do MIT](https://math.mit.edu/~drew/grh/) · [Notas do espelho](https://github.com/Deicyde/riemann-hypothesis/blob/main/verification/NOTES.md)

Também baixamos a exportação do desafio do MIT, verificamos o SHA-256 descompactado e lemos o tipo de cada um dos cinco teoremas diretamente nos registros NDJSON. O desafio contém `sorryAx` para definir as obrigações a conferir pelo Comparator. A exportação da solução é um arquivo diferente. Use Bash, Python 3.11 ou superior e `zstd`. [Hash e exportações do MIT](https://math.mit.edu/~drew/grh/)

```bash
set -euo pipefail
mkdir -p original
for parte in {01..14}; do
  curl -fLSs --retry 2 \
    "https://alpo.ge/lean.tar.gz.part${parte}" \
    -o "original/lean.tar.gz.part${parte}"
done
cat original/lean.tar.gz.part{01..14} > original/lean.tar.gz
printf '%s  %s\n' \
  9631b4f020199e580ba86178cb56c6e9b09f690c2369a80dc4334aa93f49c7fc \
  original/lean.tar.gz | sha256sum -c -
gzip -t original/lean.tar.gz
curl -fLSs https://math.mit.edu/~drew/grh/challenge.ndjson.zst -o original/challenge.ndjson.zst
zstd -d --long=27 original/challenge.ndjson.zst
python3 source_audit.py original/lean.tar.gz --challenge original/challenge.ndjson \
  --out resultado-local.json
cmp RESULTADO_LOCAL.json resultado-local.json
```

Executamos o `nanoda_bin` 0.4.17 da distribuição oficial de Lean 4.36.0 no desafio. Com os três axiomas permitidos para a solução, o teste encerrou com código 1 e identificou `sorryAx`. Permitindo também esse marcador apenas para examinar o desafio, o núcleo verificou **53.722 declarações sem erro de tipo** e encerrou com código 0. A rotina de exibição registrou `Unable to print axioms` nessa execução. Use `jq` e ponha no PATH o executável `nanoda_bin` do Lean 4.36.0 para reproduzir os dois testes. [Lean 4.36.0](https://github.com/leanprover/lean4/releases/tag/v4.36.0) · [Enunciados do desafio](https://github.com/Deicyde/riemann-hypothesis/blob/main/src/Challenge.lean)

```bash
cat > original/nanoda-strict.json <<'JSON'
{
  "use_stdin": true,
  "num_threads": 3,
  "permitted_axioms": ["propext", "Classical.choice", "Quot.sound"],
  "unpermitted_axiom_hard_error": true,
  "nat_extension": true,
  "string_extension": true,
  "print_success_message": true
}
JSON

code=0
nanoda_bin original/nanoda-strict.json < original/challenge.ndjson \
  > original/nanoda-strict.log 2>&1 || code=$?
test "$code" -eq 1
grep -F 'unpermitted axiom "sorryAx"' original/nanoda-strict.log

jq '.permitted_axioms += ["sorryAx"]' original/nanoda-strict.json > original/nanoda-challenge.json
nanoda_bin original/nanoda-challenge.json < original/challenge.ndjson
```

## Enunciados verificados

`src/Main.lean` declara `theorem riemannHypothesis : RiemannHypothesis := rhOfAll`. O tipo vem da Mathlib oficial fixada no commit `0df444a360eaa60ab8c11dca51a86af692955474`, com Lean `v4.33.1`. A prova invoca o resultado de não anulação para todos os caracteres de Dirichlet e aplica a redução pela equação funcional. [Código-fonte](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab/src/Main.lean) · [Mathlib fixada](https://github.com/leanprover-community/mathlib4/blob/0df444a360eaa60ab8c11dca51a86af692955474/Mathlib/NumberTheory/LSeries/RiemannZeta.lean) · [Redução](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab/src/Supplements/Kernels/MainTheorems.lean)

| Declaração | Enunciado |
| --- | --- |
| `riemannHypothesis` | Todo zero complexo de ζ, excetuados os inteiros pares negativos e `s = 1`, tem parte real `1/2`. Usa a definição da Mathlib. |
| `riemannZeta_eq_zero_re` | O mesmo enunciado escrito explicitamente. |
| `LFunction_ne_zero_of_half_lt_re` | Para todo módulo positivo, todo caráter de Dirichlet complexo χ e todo `s` com `Re(s) > 1/2`, vale `L(χ,s) ≠ 0`, exceto no par caráter principal e `s = 1`. |
| `LFunction_eq_zero_re` | Para todo módulo positivo e todo caráter de Dirichlet complexo χ, qualquer zero de `L(χ,s)` com `0 < Re(s) < 1` tem `Re(s) = 1/2`. |
| `LFunction_eq_zero_re_of_isPrimitive` | O enunciado anterior para caracteres primitivos. |

[Definições Lean dos enunciados](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab/src/Statements.lean)

No passo final, `rh_of_every` especializa o teorema de não anulação ao caráter de módulo 1, cuja função L é ζ. Para um suposto zero com `Re(s) < 1/2`, `rh_of_right` usa a equação funcional para transferi-lo a `1−s`, onde a parte real é maior que `1/2`. A prova trata a não anulação dos fatores gama e cosseno, excluindo `s=0` e os zeros triviais. Essa redução depende da não anulação geral provada no desenvolvimento, que é o passo substancial. Na nota, Alpöge descreve a inspiração em Stepanov, o uso de cadeias de levantamentos theta e a fórmula de Guinand–Weil. [Redução pela equação funcional](https://github.com/Deicyde/riemann-hypothesis/blob/40e56bab/src/Formalism/GammaFactorNonvanishing/RhOfRight.lean) · [Nota do autor](https://alpo.ge/remark.pdf)

## Verificações independentes e alcance

Sutherland relata que os **175.187 módulos** foram compilados e novamente examinados por `leanchecker`. O log de `#print axioms` registra exatamente `propext`, `Classical.choice` e `Quot.sound` para cada um dos cinco teoremas. `nanoda` aceitou a exportação completa de **61.570.595.788 bytes**, cujo SHA-256 é `5ea74b912147c5c93cc2b6d70cc727cf908df2f62ac6a20480a3fa02f6be244e`. A saída registra código 0. O MIT publicou também a exportação do desafio, construída a partir dos enunciados e da Mathlib oficial. [Auditoria de axiomas](https://math.mit.edu/~drew/grh/audit.log) · [Log do nanoda](https://math.mit.edu/~drew/grh/nanoda.log) · [Hashes](https://math.mit.edu/~drew/grh/SHA256SUMS)

Na atualização consultada do MIT, de 20 de novembro às 16h45 em Brasília, o Comparator de ponta a ponta continuava em execução. A comparação dos enunciados e o teste de axiomas tinham passado em uma execução. Posteriormente, Szymon Teżewski registrou a aceitação dos cinco teoremas com uma versão mais recente do `con-leche` em modo `--verified`. O MIT descreve a compilação original da sua auditoria em máquinas dedicadas e sem sandbox, e uma comparação de 77 módulos com a saída de `lake build`. [Estado e metodologia do MIT](https://math.mit.edu/~drew/grh/) · [Resultado do con-leche](https://leanprover.zulipchat.com/#narrow/channel/113488-general/topic/Discussion.3A.20Alp.C3.B6ge.20Proves.20the.20Riemann.20Hypothesis/near/642338051)

O teste local de instalação executou Lean 4.33.1. A tentativa de compilar o `Statements.lean` original parou no primeiro `import Mathlib`, pois o cache da Mathlib retornou HTTP 403 neste ambiente. O replay integral da solução requer aproximadamente 100 GB de memória e 70 GB de disco segundo o MIT, acima dos recursos locais. Na consulta, o Clay Mathematics Institute ainda classificava o problema como não resolvido. [Recursos e comandos para replay](https://math.mit.edu/~drew/grh/) · [Clay Mathematics Institute](https://www.claymath.org/millennium/riemann-hypothesis/)
