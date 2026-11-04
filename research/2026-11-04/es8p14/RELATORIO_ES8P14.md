# ES8-P14 — Ponte entre caixas intervalares e certificados LDL

**Data:** 04/11/2026. **Escopo:** matrizes projetadas de dimensões 4, 6, 8. Não é uma prova da Hipótese de Riemann.

## 1. Teorema estrutural

Sejam \(C,M\in\mathbb{Q}^{n\times n}\), e \(e_{ij}=e_{ji}\ge 0\) tais que \(|M_{ij}-C_{ij}|\le e_{ij}\), \(\sum_j e_{ij}\le\delta\). Para qualquer vetor racional \(x\),

\[
x^\top Mx\ge x^\top Cx-\delta\sum_i x_i^2.
\]

A prova é inteiramente algébrica: para quaisquer \(a,e,u,v\in\mathbb Q\), com \(-a\le e\le a\), as quantidades \((a+e)(u+v)^2\) e \((a-e)(u-v)^2\) são não negativas. Sua soma implica
\[
  uev\ge -\frac{a}{2}(u^2+v^2).
\]
Somando sobre \(i,j\), a simetria \(e_{ij}=e_{ji}\) transforma o limite da soma em \(\sum_i(\sum_j e_{ij})x_i^2\le\delta\|x\|_2^2\).

Se, adicionalmente, \(C-\frac{3}{50}I\succ0\) e \(\delta\le\frac1{100}\), conclui-se
\[
  x^\top Mx > \left(\frac{3}{50}-\frac1{100}\right)\|x\|_2^2
 =\frac1{20}\|x\|_2^2,
\]
para todo \(x\ne0\). Este é o objetivo das definições e teoremas `positive_from_LDL_and_box` e `box_positive{4,6,8}`.

## 2. Proveniência e verificação dos extremos

O arquivo da ES8-P13 fornece 2 extremos inteiros por entrada em escala \(10^{90}\) e a matriz racional central com denominador \(10^4\). A ES8-P14 calcula, com `fractions.Fraction`, um limite \(e_{ij}\) arredondado **para cima** em incrementos de \(10^{-6}\) para ambos os extremos. Não utiliza `float` nas decisões matemáticas. As versões geradas em Lean incorporam *todos os extremos*, provam que cada um fica na envoltória, a simetria da envoltória e a soma das linhas, com os seguintes resultados:

| Matriz | Máximo de soma dos erros por linha | Margem restante acima de \(0.05\) |
|---|---:|---:|
| 4×4 | 0.001937 | 0.008063 |
| 6×6 | 0.004379 | 0.005621 |
| 8×8 | 0.008333 | 0.001667 |

O teorema LDL da ES8-P13 fornece \(C_N-0.06I\succ0\); os intervalos extremos representam os dados auditados da ES8-P10, não valores especiais comprovados em Lean.

## 3. Condições não resolvidas

1. **Escopo algébrico:** o enunciado mecanizado aplica-se a matrizes e vetores racionais. É possível ampliar o teorema para \(\mathbb R\) via conversão das identidades e o mesmo argumento algébrico, mas isso não foi automatizado nesta versão.
2. **Inclusão analítica:** a afirmação de que as integrais e funções zeta realmente pertencem às caixas ainda utiliza a justificação externa de ES8-P10. O Lean não verificou essas identidades.
3. **Positividade universal:** a propriedade das matrizes N=6,8,10, mesmo com caixas demonstradas, não implica positividade da forma de Weil sobre todos os suportes e todas as funções.

## 4. Hipótese de pesquisa, não teorema

Em vez de simplesmente elevar a dimensão mantendo \(L=0.45\), propõe-se testar uma **estratégia dual de testemunhas negativas**: estudar, assumindo a existência hipotética de um zero fora da linha crítica, se e com que crescimento de suporte uma função-teste admissível criaria uma forma de Weil negativa. Esse método deve buscar falsificações da hipótese estrutural proposta e observar o que acontece quando \(L\) cruza limiares \(\frac12\log p\), incorporando mais primos. Não há afirmação de descoberta ou solução de RH.

## 5. Estado da reprodução

Python exato: 6 testes aprovados localmente. O arquivo Lean foi criado, e sua compilação remota deve ser consultada separadamente antes de classificá-lo como formalmente verificado.
