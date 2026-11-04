# NUMERIA–RIEMANN — ES8-P14 (04/11/2026)

**Objetivo:** construir e formalizar a ponte de perturbação entre as caixas intervalares de 90 dígitos reconstruídas em ES8-P13 e as matrizes racionais cujas fatorações LDL foram compiladas em Lean.

## Arquivos

- `ES8P14.lean`: teorema matemático geral de limite inferior de uma forma quadrática sob erros de entradas simétricos e orçamento por linha, usando a desigualdade para pares de coordenadas; critério da margem `3/50 - 1/100 = 1/20`.
- `ES8P14Concrete.lean`: matrizes racionais de dimensão 4, 6 e 8 com os extremos exatos de escala `10^90`, envoltórias arredondadas para cima em escala `10^6`, provas de simetria, erros e margem para qualquer matriz **racional** em cada caixa.
- `RESULTADOS_ES8P13.json`: fonte imutável dos extremos e dos centros racionais, copiada do certificado ES8-P13.
- `RESULTADOS_ES8P14.json`: metadados da reconstrução e erros por linha.
- `generate_es8p14.py`: gera novamente as definições Lean a partir do JSON exato; não usa ponto flutuante para decisões.
- `test_es8p14.py`: controles por frações exatas e testes de mutação adversarial.
- `numeria-es8p14-lean.yml`: compilação em GitHub Actions com Lean 4.35.0 e Mathlib; reexecuta os testes e regenera os arquivos.

## Cotas por linha

| N | dimensão | cota racional de erro por linha | orçamento |
|:-:|:-:|--:|--:|
| 6 | 4 | 0,001937 | < 0,01 |
| 8 | 6 | 0,004379 | < 0,01 |
| 10 | 8 | 0,008333 | < 0,01 |

**Dependências e limites:** as inclusões da ES8-P10 identificam os valores das entradas **apenas se forem válidos os lemas analíticos externos** da fórmula explícita, cauda arquimediana, digama e interpolação. As provas Lean verificam os extremos racionais e a implicação **dentro de \(\mathbb Q\)**. A extensão para matrizes reais e a verificação mecanizada das identidades analíticas exigem mais trabalho. Isso não demonstra positividade universal nem RH.

## Execução local

```sh
python -m unittest -v test_es8p14.py
python generate_es8p14.py
```

## Compilação Lean

No GitHub Actions, `lake env lean` compila `ES8P14.lean`, `ES8P13.lean` e `ES8P14Concrete.lean` (importando os módulos anteriores). O status deve ser considerado **pendente** até que o workflow remoto termine.
