# ES8-P238 — auditoria da ponte P207 e alvo quantitativo de Schur

**Status:** auditoria de lacunas e implicação suficiente condicional. Não certifica a cauda de Weil nem positividade global.

## Resultado finito reutilizado de P237

Recomputando em Arb a partir do artefato P237 versionado, os limites inferiores transferidos para os dois blocos intervalares são aproximadamente `4,7912 × 10⁻¹⁸` (par) e `5,2098 × 10⁻¹⁵` (ímpar). Portanto, `m₀ = 10⁻¹⁸` é um limite comum conservador para ambos. A verificação reproduzível está em `p238_schur_bridge_audit.py`, e o resultado estruturado em `P238_SCHUR_MARGIN.json`.

Essa margem vale para as matrizes finitas P237 e suas caixas Arb. Usá-la numa forma global exige, antes, provar que tais blocos são as compressões principais da forma de Weil normalizada, com diferença controlada por `ε_bridge`.

## Lema de transferência necessário

Considere uma decomposição ortogonal real ou complexa `H = H₀ ⊕ H₁`, onde `dim H₀=377` por setor, e um operador/forma simétrica

```text
W = [[A, B], [B*, D]].
```

Uma hipótese suficiente para positividade estrita é:

```text
A >= (m₀ - ε_bridge) I,   D >= η I,   η > 0,
||B|| <= β,   ε_bridge + β²/η < m₀.
```

De fato, para `x ∈ H₀`, `y ∈ H₁`,

```text
<W(x,y),(x,y)> >= (m₀-ε_bridge)||x||² - 2β||x||||y|| + η||y||²,
```

e a forma escalar do lado direito é estritamente positiva quando o determinante `(m₀-ε_bridge)η - β²` é positivo. Equivalentemente, sob as hipóteses de domínio/invertibilidade pertinentes, basta o complemento de Schur `A - B D⁻¹ B*` ser positivo. Com `ε_bridge=0` e o limite comum P237, a meta normalizada é

```text
β / sqrt(η) < 10⁻⁹.
```

Se não houver uma cota uniforme `D >= ηI` com `η>0`, essa versão não se aplica. Seria necessário um controle relativo, por exemplo de `||D^{-1/2}B*||`, com domínio, faixa e eventual núcleo do operador do complemento tratados explicitamente.

## O que a auditoria dos materiais encontrou

1. **ES8-P15**, branch `origin/research/numeria-es8p15-tail-tsum-uniformity-2026-11-04`, HEAD `939d5bb`: `ES8P15.lean` e `ES8P15Budget.lean` provam implicações condicionais e um orçamento racional para uma cauda arquimediana **bruta, não projetada**, no setup `L=9/20`, `R(n)=200n²`. O orçamento `16/(3645n)+2/(30375n³)<1/100` é exato em `ℚ`, mas as premissas analíticas da cota E7 não são derivadas, e a margem uniforme depende de hipóteses para todo `n` não instanciadas pelas matrizes P237. Esse parâmetro `L` também difere de P237 (`L=161/200`).
2. **E8K V19**, branch `origin/audit/e8k-v19-analytic-tail-2026-11-01`, HEAD `b13120c`: replay racional de cotas P5/P6 e testes cobrem caudas arquimedianas por entrada para uma base senoidal, sob precondições como `R >= max(1,2w_max)`. P5/P6 não fornecem a coerção do bloco infinito `D` nem a norma do acoplamento `B`. As identidades analíticas e as cotas da digama/integrais permanecem parcialmente assumidas nos materiais; o arquivo `TailBoundV16.lean` formaliza somente `1-(w/R)^2 >= 3/4` sob suas hipóteses racionais.
3. **Ligação com P207/P237:** buscas nas duas branches por `P207`, `A_R`, `Zhu`, `P237`, `377` e `Schur` não acharam uma instanciação explícita. A árvore principal também não contém um certificado de bloco complementar ou de acoplamento que possa ser aplicado às matrizes P237. Isso descreve o conteúdo auditado, não exclui uma relação matemática ainda não registrada.

## Lacunas concretas a fechar

- Definir o espaço funcional global, o domínio da forma de Weil e a projeção cujo bloco principal é a matriz P237; provar simetria/hermiticidade e a identidade com a normalização explícita.
- Controlar `ε_bridge`: quadratura e arredondamentos já estão na caixa P237, mas ainda é preciso derivar que o modelo finito coincide com a compressão pretendida, incluindo primos, polos e restrições da redução de Zhu.
- Provar coerção para o complemento (`D >= ηI`) ou uma estimativa relativa equivalente; uma cota de cauda entrada a entrada não basta para isso em dimensão infinita.
- Derivar um limite de operador para `B` incluindo termos arquimedianos projetados, primos e demais contribuições; então verificar `ε_bridge + β²/η < 10⁻¹⁸` (ou a condição relativa correspondente) com aritmética rigorosa.
- Só depois justificar passagem de vetores de suporte finito ao domínio completo por densidade/continuidade, sem perda da margem.

Assim, P237 acrescenta uma margem finita explícita que torna o alvo de acoplamento mensurável; P15/V19 oferecem lemas úteis para componentes arquimedianas sob hipóteses, mas não fecham o elo P207. A positividade global de Weil e a Hipótese de Riemann continuam em aberto.
