# Parecer técnico — Qwen, "Forma de Weil como Perturbação de Posto Finito"

Fonte avaliada: `QWEN_ORIGINAL_2026-11-04.tex` (cópia byte a byte do documento fornecido), seções 1–6. Data da auditoria: 2026-11-04.

## O que o texto demonstra corretamente

1. A relação `center_n = matrix_n + (3/50) I` foi demonstrada para `n=4,6,8` no arquivo publicado `NUMERIA_E8K_FORMAL/ES8P14Concrete.lean`, e os certificados racionais positivos permanecem válidos nas caixas intervalares sob hipóteses analíticas externas.
2. Para cada dimensão FINITA fixa, toda diferença de matrizes tem posto no máximo `n`. Isso é correto, mas não distingue a forma de Weil de qualquer outra matriz finita.
3. Os limites de soma por linha (`0,001937`, `0,004379`, `0,008333`) dizem respeito SOMENTE aos casos verificados, não à família para todo `N`.

## Erros lógicos que devem ser corrigidos

**E1. Posto finito no limite não segue de posto finito em cada dimensão.**
As matrizes `(1/1000) I_n` têm posto `n` e norma `1/1000`; o operador limite `(1/1000) I` em `ℓ²` tem posto infinito e é não compacto. Por isso o item 3 da suposta prova principal é inválido como argumento para `K_∞`.

**E2. Erro pequeno em norma não implica posto finito.**
`||R||<1/100` controla uma perturbação limitada, NÃO sua dimensão de imagem. A conclusão em linhas 98–104 do LaTeX é válida para a matriz `K_N` apenas porque **todo o espaço é finito**; não pela propriedade do erro.

**E3. O deslocamento `(3/50) I` é posto `n` na seção `n`, mas tem posto infinito em `ℓ²`.**
A decomposição do centro em `S_n+(3/50)I_n` NÃO define uma decomposição global da forma de Weil como positiva + posto finito. O texto não constrói `Q_0` nem uma família coerente de `K_∞`.

**E4. As matrizes `N=6,8,10` não determinam todas as entradas de uma forma infinita.**
Pode-se acrescentar a uma extensão um coeficiente diagonal negativo depois desses prefixos sem mudar as três matrizes existentes. Sem positividade de TODAS as seções, não se pode aplicar o critério de positividade da forma em `c₀₀`.

**E5. Bibliografia e terminologia.**
O texto atribui a `Baez, J. C.; Duarte, A.` um título sem dados bibliográficos verificáveis. Não o utilizamos como fonte. A equivalência com RH é tradicionalmente descrita pelo *critério de Weil*, não demonstrada no documento por uma relação com um autor específico.

## Correção matemática aproveitável (ES8-P19)

Se `H` é um espaço de Hilbert, `A` é operador **limitado** em `ℂ²⊕H` com decomposição em blocos, `C:H→ℂ²` é limitado e `B(x)=(-Cx,x)`, então

```
B* A B − A_HH = C* (A_00 C − A_0H) − A_H0 C.
```

A imagem da correção está contida em `ran(C*) + ran(A_H0)`, subespaço de dimensão **no máximo 4**. Este resultado é genuinamente uniforme em relação à dimensão, mas diz respeito **somente à correção introduzida pelo cancelamento de dois polos**. Não implica que `A_HH` seja positivo, semibounded ou a forma verdadeira de Weil. Para um `A` não limitado é necessário construir o domínio e justificar cada bloco.

O arquivo `ES8P19.lean` formaliza, sem `fin_cases` no tamanho `n`, a identidade exata que expressa a correção da compressão finita como quatro produtos exteriores. O arquivo `audit_es8p19.py` reproduz a identidade com frações exatas, e encontra posto efetivo 4 para dimensões `n≥4` em um exemplo racional simétrico não trivial.

## Classificação

- **Aceito:** identidades de deslocamento da ES8-P14 no escopo finito e as cotas por linha já verificadas.
- **Corrigido:** a única afirmação genuína de baixo posto é a correção DOIS-POLOS `B*A B−A_HH`, sob hipótese de blocos bem definidos.
- **Rejeitado como prova:** inferência de baixo posto da forma completa a partir de matrizes finitas e erros pequenos.
- **Não estabelecido:** forma infinita efetiva de Weil, continuidade/fechamento, margem uniforme, positividade universal ou RH.
