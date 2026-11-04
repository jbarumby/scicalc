# NUMERIA–RIEMANN — ES8-P19: ponte algébrica para a forma global

**Data:** 2026-11-04. **Base:** ES8-P18, commit `2aec1bb46d3a568317fa37a615d44062e32660ed`.

**Classificação:** `ALGEBRA_DE_BLOCOS_E_FORMA_C00; RH_NAO_DEMONSTRADA`.

## 1. Origem e decisão

Foi analisado o documento `QWEN_ORIGINAL_2026-11-04.tex`, entregue pelo usuário. Seu teorema de "perturbação de posto finito da forma de Weil" NÃO é válido como teorema de operadores infinitos: o argumento utiliza apenas fatos triviais de dimensão finita e confunde limitação de norma com posto. O arquivo `PARECER_QWEN_ES8P19.md` registra a análise ponto a ponto.

## 2. Um teorema algébrico correto — correção de dois polos

Considere uma decomposição `V=P⊕H`, com `dim P=2`, e a matriz de blocos de um operador limitado `A`:

\[
A=\begin{pmatrix}A_{00}&A_{0h}\\ A_{h0}&A_{hh}\end{pmatrix}.
\]

A projeção do projeto é a inclusão modificada `B:H→P⊕H`, `B(x)=(-Cx,x)`, com `C:H→P`. Expansão exata:

\[
B^* A B=A_{hh}-C^*A_{0h}-A_{h0}C+C^*A_{00}C.
\]

Rearranjo:

\[
\boxed{B^*AB-A_{hh}= C^*(A_{00}C-A_{0h})-A_{h0}C.}
\]

Logo

\[
\operatorname{ran}(B^*AB-A_{hh}) \subseteq \operatorname{ran}(C^*)+\operatorname{ran}(A_{h0}),\qquad
\boxed{\operatorname{rank}(B^*AB-A_{hh})\le4}.
\]

A hipótese de **A limitado** (ou uma formulação de formas com domínios apropriados) é essencial para interpretar esses produtos como operadores limitados em H. Esta propriedade vale para blocos finitos independentemente de qualquer limite. O Lean da fase prova a identidade dos quatro produtos exteriores para **todo tamanho `n`**, sem enumerar dimensões específicas; o argumento linear-algébrico de posto do operador infinito permanece apresentado no relatório, não formalizado em Mathlib.

O resultado NÃO diz que a **forma completa de Weil** seja uma perturbação de posto ≤4 de uma **forma positiva**. A compressão de referência `A_hh` pode ser indefinida e o operador arquimediano não limitado precisa de domínio próprio.

## 3. Forma quadrática definida inicialmente em `c₀₀`

Para uma matriz infinita hermitiana/coerente `a(i,j)`, define-se sem ambiguidade uma forma algébrica em sequências com suporte finito:

\[
q_a(x)=\sum_{i\in\mathrm{supp}(x)}\sum_{j\in\mathrm{supp}(x)}\bar x_i a(i,j)x_j.
\]

É uma soma FINITA; não requer `tsum` nem converge automaticamente a uma forma fechada em `ℓ²`. Em Lean usamos coeficientes reais racionais (sem conjugação) e `Finsupp`:

- Se `q_a` é não negativo em **toda seção finita** (`∀ s : Finset ℕ`), então é não negativo em `c₀₀`.
- Se cada seção satisfaz `q_a≥γ||x||²` com o **mesmo γ**, o mesmo vale em `c₀₀`.
- Esses teoremas deixam as hipóteses universais explícitas. A verificação apenas das dimensões 4, 6 e 8 NÃO cumpre essas hipóteses.

Para transferir ao completamento é necessário um resultado de densidade **no domínio da forma**, continuidade ou fechamento e uma semilimitação apropriada; não é lícito apenas usar a densidade de `c₀₀` em `ℓ²` para uma forma não limitada.

## 4. Auditoria racional independente

O código `audit_es8p19.py` avalia uma família de blocos simétricos não triviais com frações exatas em dimensões `1,2,3,4,8,10,16,32`:

- Identidade das compressões verificada entrada por entrada.
- Correção fatorada em quatro produtos exteriores. Posto é `min(4,n)` nas amostras, demonstrando que a cota 4 é ótima em geral.
- A invariância à simetria e mutações adversariais são verificadas.
- Contraexemplo: `ε I` em `ℓ²`, com `ε=1/1000`, tem norma pequena e posto infinito.
- Contraexemplo: `diag(1/(k+1))` tem todas as seções positivas, mas nenhum limite espectral uniforme estritamente positivo.
- Contraexemplo: prefixos positivos até oito índices são compatíveis com uma entrada diagonal negativa no índice 9.

**Todos esses testes usam modelos de álgebra racional, não integrais certificadas de Weil.**

## 5. Não confundir as três afirmações

| Afirmação | Estado |
|---|---|
| Correção de polos `B* A B−A_hh` possui posto ≤4, para A limitado com blocos | Demonstrada matematicamente e identidade finita codificada em Lean |
| Os centros ES8-P13/P14 são positivos para três tamanhos | Provado previamente, condicional à identificação analítica |
| Toda a forma de Weil é positiva em `c₀₀` | NÃO demonstrado; faltam todos os coeficientes e todas as seções |
| A forma em `c₀₀` é fechável/semilimitada e sua extensão é a forma de Weil | NÃO demonstrado |
| A Hipótese de Riemann | NÃO demonstrada |

## 6. Próxima obrigação concreta — ES8-P20

1. Definir e calcular os coeficientes arquimedianos e de primos para **todo par de modos** sob uma normalização verificada, em vez de postular a matriz infinita.
2. Fixar o domínio de formas (inicialmente a imagem por Fourier de `c₀₀`); demonstrar que as restrições a `Finset.range N` recuperam literalmente os certificados históricos.
3. Estabelecer uma desigualdade de semilimitação/fechabilidade no domínio apropriado; investigar se a parte primo é limitada relativa à parte arquimediana para L fixo e suporte crescente.
4. Não utilizar *norma pequena*, *compacidade* ou *posto finito de uma correção* como substituto da positividade da referência.

## 7. Observação sobre compilação

`ES8P19.lean` é um arquivo de provas proposto com hipóteses explícitas e sem lacunas textuais. Só receberá classificação `LEAN_COMPILADO` após retorno de sucesso do workflow de GitHub Actions. Os resultados Python já foram executados localmente.
