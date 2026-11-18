# NUMERIA ES8-P231 — ponte finita P230 → prefixo de 4 modos

## Resultado

P231 testa uma interface concreta entre os artefatos finitos P227/P228/P230 e a ponte quadrática condicional do P120. O escopo é deliberadamente pequeno: apenas os quatro primeiros modos pares, índices `[0,1,2,3]` (graus `[0,2,4,6]`), no setor `even`, `T*=200`.

Do centro racional P193, cada entrada do bloco 4×4 foi arredondada exatamente para uma fração decimal com denominador `10^40`. O extrator recompõe as caixas Arb P227, acrescenta os mesmos termos beta/polos e o erro de quadratura independente P229, e verifica com `Fraction` que **os 16 intervalos do prefixo estão a menos de `10^-40` do centro arredondado**. A maior distância reportada é aproximadamente `4.9977031301347021 × 10^-41`; o maior erro de arredondamento do centro é aproximadamente `4.9977031298379128 × 10^-41`.

A proveniência é ligada por SHA-256: P231 confere a hash comprimida da matriz P193 com a hash registrada em P230 e recalcula, a partir dos artefatos P227, a cota de soma de linha dos 200 modos. A cota resultante é aproximadamente `5.9357849352196727 × 10^-49` e deve ser exatamente igual ao campo P230. A margem inferior de autovalor finita P230 também é positiva (racional armazenada no JSON; aproximadamente `10^-20`).

Para a matriz central arredondada `C₄`, o arquivo Lean registra uma fatoração LDL racional exata de `C₄ − 10^-20 I`. Os quatro pivôs são estritamente positivos, aproximadamente `0.07239918`, `0.00212209`, `0.00394246` e `0.00173125`; tanto a identidade matricial quanto os sinais são verificados exatamente. O teorema `center_lower4` deriva daí

\[
 x^T C_4 x \ge 10^{-20}\sum_i x_i^2.
\]

A instância `P231_box_positive` aplica a ponte finita no estilo P120 e mostra que qualquer matriz racional hermitiana 4×4 cujas entradas satisfaçam `|Wᵢⱼ − C₄ᵢⱼ| ≤ 10^-40` é positiva definida. O orçamento quadrático usado é `rho = 4×10^-40`, estritamente menor que o deslocamento `shift = 10^-20`. O extrator demonstra que as caixas P227 do prefixo estão contidas nessa vizinhança.

## O que isso não prova

O módulo Lean exige como hipótese que `W` seja hermitiana e que suas entradas estejam na caixa uniforme. A checagem racional externa conecta as caixas finitas arquivadas P227 a `C₄`; **não prova** que as integrais/coeficientes matemáticos de Weil estejam dentro delas. Essa inclusão continua condicional às fórmulas e majorantes P227/P229. A hermiticidade do objeto analítico não foi derivada no Lean deste estágio.

P231 certifica somente uma afirmação finita sobre um prefixo de quatro modos. Não fornece controle dos acoplamentos para modos fora do prefixo, não resolve a cauda P207, não estabelece margem uniforme para todos os modos, não identifica a forma global de Weil, não prova positividade global e não prova a Hipótese de Riemann.

## Arquivos

- `P231_PREFIX4_DATA.json`: centros, caixas racionais, hashes, fatoração LDL e resultados exatos do prefixo.
- `p231_extract_prefix.py`: recompõe as caixas P227/P229 para o prefixo, verifica a hash e a cota P230 e gera o JSON. O replay integral exige os chunks P227 e entradas P193 usados anteriormente; esses diretórios temporários não estão incluídos neste repositório.
- `p231_emit_lean.py`: gera o módulo Lean com os racionais exatos do JSON e incorpora as hashes de proveniência.
- `test_p231_prefix.py`: quatro testes exatos de proveniência, inclusão das 16 caixas, simetria, fatoração LDL e margem.
- `../NUMERIA_E8K_FORMAL/P231_P230Prefix4.lean`: prova Lean da fatoração, margem central, erro por caixa e conclusão condicional `Matrix.PosDef`.
- `README_ES8P231_BRANCH_SYNTHESIS.md`: auditoria anterior dos 12 ramos; registra por que este é o próximo gate local.

## Reprodução

Com as entradas upstream disponíveis localmente:

```bash
python3 numeria/p231_extract_prefix.py \
  --chunks /caminho/para/P227all \
  --p193 /caminho/para/P193/inputs \
  --p230 numeria/P230_FINITE_SPD.json \
  --out numeria/P231_PREFIX4_DATA.json
python3 numeria/p231_emit_lean.py \
  --data numeria/P231_PREFIX4_DATA.json \
  --out NUMERIA_E8K_FORMAL/P231_P230Prefix4.lean
python3 -m unittest -v numeria/test_p231_prefix.py
```

Para compilar a prova, com Lean 4.35.0 e Mathlib conforme `NUMERIA_E8K_FORMAL/lean-toolchain` e `lakefile.toml`:

```bash
cd NUMERIA_E8K_FORMAL
lake env lean P231_P230Prefix4.lean
```

O arquivo Lean inclui localmente o lema condicional de matriz empregado no P120; a branch-base P227 não rastreia o arquivo Lean P120. O módulo é autônomo quanto a essa ponte, mas depende de Mathlib.

## Estado de validação

O extrator foi executado contra os artefatos locais; o SHA do centro P193 e a cota de linha P230 coincidem. Os quatro testes Python exatos passaram (`Ran 4 tests ... OK`), e os scripts passaram por `py_compile`. A prova Lean foi compilada com Lean 4.35.0 e o cache Mathlib direcionado: `lake env lean P231_P230Prefix4.lean` terminou com código zero, sem erros ou avisos. O módulo gerado não contém `sorry`, `admit` ou declarações de axiomas adicionais.
