# Relatório técnico — NUMERIA ES8-P151

**Data:** 10/11/2026  
**Objeto:** ampliar a família da P150 para centros positivos/negativos e escalas variáveis, com densidade construtiva e controle explícito de momentos.

## Conclusão principal

A família `B_d=Λ_d*φ_d`, com `d=2^(-m)`, e translações `j d` para `|jd|≤m` fornece, dentro de um **único nível** `m`, uma base admissível explícita de dimensão `2m·2^m-1`. Cada elemento da base é uma combinação de **três** translações contíguas com coeficiente central `-2cosh(d/2)`. Os dois momentos exponenciais cancelam simbolicamente.

Para cada função admissível `f ∈ C_c²([-L,L])`, há aproximantes dentro desse núcleo, com suporte fixo e erro `C¹` limitado por

`(1+48A·2^A){ d ||f'||∞ + (3d/2+d²/8)||f''||∞ }`, `A=max(ceil(L)+1,3)` e `d=2^(-m)`, para `m≥max(3,ceil(L))`.

O detalhe não trivial é que os dois corretores podem ser gerados **dentro do mesmo nível**: suas matrizes de momentos têm inversas limitadas independentemente de `m` porque possuem uma representação exata por dois transladados de uma mesma função par, positiva e não nula.

## Custo de gerar a família

| N | P147: geradores acumulados | P151: geradores acumulados | P151: geradores de uma escala | P151: base admissível de uma escala |
|---:|---:|---:|---:|---:|
| 4 | 559.238 | 244 | 129 | 127 |
| 6 | 214.748.360 | 1.518 | 769 | 767 |
| 8 | 73.300.775.178 | 8.168 | 4.097 | 4.095 |

A grande redução decorre da interpolação por chapéus e posterior suavização, em vez de somas de Riemann de passo `η^4`. Os espaços acumulados P151 são crescentes; a prova de densidade mais forte já funciona usando apenas os espaços individuais de cada nível.

## Redução de Gram por translação

Se a forma bilinear de Weil na normalização P111 for invariante por translação simultânea, a matriz bruta de cada escala é Toeplitz. Após aplicar o filtro tridiagonal de momentos, a matriz restrita continua Toeplitz; sua função de lag é uma combinação de cinco lags originais com coeficiente `c=2cosh(d/2)` (fórmula T4 no teorema).

Essa estrutura reduz os valores independentes necessários para uma matriz restrita de dimensão 4.095 de mais de 8 milhões (armazenamento genérico) para aproximadamente **4 mil lags**, antes de verificações analíticas e intervalares. A redução não transforma os lags em certificados positivos.

## Validação executada

- 19 testes de programação aprovados: índices e contagens inteiros, orçamento de erro racional, fator conservador de projeção, anulamentos exatos em polinômios de Laurent e identidades Toeplitz por congruência com racionais.
- Anulamento simbólico conferido em todos os `8.152` triplos dos níveis 1 a 8 na configuração acumulada de janela `N=8`.
- Cálculos da forma de Weil real: **não realizados** nesta etapa.
- Formalização em Lean/Arb: **não realizada** nesta etapa.

## Próxima fase recomendada — ES8-P152

1. Especificar uma fórmula explicitamente validada para os **lags** da forma de Weil sobre `B_d` e confirmar sinais/polos frente à referência clássica.
2. Calcular uma primeira janela pequena com controle de arredondamento e integral arquimediana.
3. Utilizar a transformação Toeplitz exata para calcular a Gram admissível; certificar PSD em cada janela verificada, sem inferir o resultado universal.
4. Investigar um princípio analítico que trate **todos** os níveis; essa é a condição que ainda falta para a aplicação do critério de Weil.

**Estado científico:** progresso construtivo em densidade e compressão de matrizes; positividade universal e RH não demonstradas. Nenhuma alteração no GitHub.
