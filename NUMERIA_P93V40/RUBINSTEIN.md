# Auditoria do guia sobre a prova de Alpöge

A publicação da página Mathematics, intitulada *The Alpöge–Claude proof of the Riemann hypothesis explained for the Riemann zeta function (Part 2)*, apresenta imagens da primeira versão do guia publicado por Michael O. Rubinstein. A primeira imagem corresponde à página impressa 80, abertura da seção III.5. A expressão para $V(j,t,p)$ recebe o número (III.76), tal como na primeira versão. Na versão 2 ela passa a (III.77). A segunda imagem também confere com a numeração da primeira versão. O PDF integral que corresponde às imagens tem 164 páginas, com data de 21 de novembro de 2026. A versão 2, datada de 22 de novembro, tem 265. Foi possível conferir o início do recorte; a navegação pelas demais imagens passou a exigir login. A extensão de 77 páginas informada para o recorte ficou sem verificação integral. [Post do Facebook](https://www.facebook.com/100064026345553/posts/1645378550939674/) · [Primeira versão](https://math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper.pdf) · [Versão 2](https://math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf)

| Fonte | SHA-256 do PDF obtido na Waterloo |
| --- | --- |
| Primeira versão | `667d88b3c375adf0466678094f452ca966aef0f4d368f6340d2bd5328b7fb706` |
| Versão 2 | `b49de5037af55c6df66843db710332a358c4eb1d2799581c673cd2ea17e8ea40` |

## Resultados conferidos

O replay racional cobre 6.138 instâncias do passo de otimização dos raios no Lema 6 de §III.5 e 3.267 instâncias da fórmula para órbitas de dois lugares, mediante enumeração independente de vértices do programa linear. Os exemplos 5-ádicos e de raios também conferem. Esses cálculos verificam instâncias finitas do problema abstrato. A identificação dos núcleos teta da construção com as séries que satisfazem a hipótese de integralidade (B) e a aplicação de `LinkG` pertencem à cadeia de argumentos referenciada pelo guia. [Replay racional](rubinstein_replay.json)

Na seção I.1, o objeto analítico é

$$
G(x)=\sum_{n\le x}\Lambda(n)(x/n-1)-x\log x+(1+\gamma)x.
$$

Um crivo inteiro e Arb a 128 bits reproduziram os quatro valores tabulados, com bolas inteiramente contidas nas células de arredondamento a três decimais.

| $x$ | Valor a três decimais |
| ---: | ---: |
| 10.000 | 1,749 |
| 100.000 | 3,671 |
| 1.000.000 | 11,040 |
| 2.000.000 | −12,268 |

Os testes de adulteração rejeitaram uma alteração em um contador do JSON e rejeitaram o PDF da versão 2 quando solicitado o SHA-256 da primeira versão.

Com $\Re s>1$, a transformada de Mellin de $(x/n-1)_+$ é $n^{-s}/[s(s-1)]$. A soma absolutamente convergente dá $-(\zeta'/\zeta)(s)/[s(s-1)]$, e os dois termos principais de $G$ cancelam a parte polar em $s=1$. A Parte II explica por que, para cada $\varepsilon>0$, um limite superior de ordem $x^{1/2+\varepsilon}$ uniforme em janelas $[Y,Y^{3/2}]$ arbitrariamente altas exclui zeros à direita da reta crítica. Na cadeia do guia, essa cota vem da desigualdade aritmética (M) e de duas estimativas para caracteres de Dirichlet. A versão 2 apresenta treze hipóteses clássicas que sustentam a redução de (M). [Versão 2, §§I.1, II.1 e III.11](https://math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf)

Em 21 de novembro de 2026, a Lean FRO e a Mathlib Initiative informaram a reprodução da compilação da biblioteca de Alpöge, a coincidência do enunciado com a Hipótese de Riemann de Mathlib e a validação do arquivo de prova por verificadores independentes. Essa verificação externa e o escopo condicional da exposição escrita são fatos sobre artefatos distintos. [Comunicado conjunto](https://leanprover-community.github.io/blog/posts/riemann-hypothesis/)

## Relação com a V.40

O parâmetro $x$ da janela aritmética controla o comprimento da soma de von Mangoldt. Na V.40, $a$ controla o suporte físico, $t$ é a variável de Fourier e $2a-\log16$ determina a interação do novo termo primo. A verificação formal externa não fornece os recintos de $A_0-B$, um limite coercivo para esse minorante ou a transferência através de 16. O cruzamento local é avaliado pelos certificados dirigidos da própria V.40. [Relatório da V.40](README.md)

Para reproduzir a auditoria combinatória a partir da raiz do repositório

```sh
python -O -m NUMERIA_P93V40.rubinstein_replay --verify NUMERIA_P93V40/rubinstein_replay.json
```
