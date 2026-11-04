# NUMERIA–RIEMANN — ES8-P18: extensão infinita da correção dos polos e constante melhorada

Data: 2026-11-04. Base: ES8-P17 (`eca8245ddf980e80a4e67944d08692786d422eb2`). Fontes: fórmula original `es8p_projection.py`, `ES8P17.lean`, `RELATORIO_ES8P17.md` e último LaTeX encaminhado `deepseek_latex_20261104_dd6431.tex` (já revisado na ES8-P17). **Nenhum novo arquivo foi anexado na solicitação atual**.

## 1. Resultado principal: melhora de \(\|B_N\|^2\)

O cancelamento dos polos (restrito a `L=9/20`) é descrito por

\[
 c_k(X)=\frac{k}{p}\frac{p^2X+L^2}{k^2X+L^2},\qquad X=\pi^2\ge 9,
 \quad p=\begin{cases}1 & k\text{ ímpar},\\2&k\text{ par}.\end{cases}
\]

Para `k≥3`, `p<k` e `X≥9`, temos o *majorante racional rigoroso*

\[
0<c_k(X)\le \frac{a_p}{k},\qquad
 a_1=\frac{409}{400},\quad a_2=\frac{1609}{800}.
\]

**Paridade ímpar.** Para `k=2m+1`, `m≥1`:

\[
\frac{1}{(2m+1)^2}\le\frac{1}{4m(m+1)},\quad
\sum_{m\ge1}\frac{1}{4m(m+1)}=\frac14.
\]

Portanto \(\sum_{k\ge3,\;k\,\mathrm{ímpar}}c_k^2\le 167281/640000\).

**Paridade par.** Para `k=2m`, `m≥2`, separamos `m=2` e `m=3` e aplicamos
\(m^{-2}\le[m(m-1)]^{-1}\) a `m≥4`:

\[
\sum_{m\ge2}\frac{1}{(2m)^2}\le
\frac1{16}+\frac1{36}+\sum_{m\ge4}\frac1{4m(m-1)}
=\frac1{16}+\frac1{36}+\frac1{12}=\frac{25}{144}.
\]

Assim \(\sum_{k\ge4,\;k\,\mathrm{par}}c_k^2\le 2588881/3686400\).

As duas linhas de cancelamento são ortogonais entre si por terem suportes disjuntos. Para `B_N=[-C_N;I]`, vale a identidade de Gram

\[
B_N^*B_N=I+C_N^*C_N,\qquad
\|B_N\|_{\mathrm{op}}^2=1+\max(\|c_{\mathrm{ímpar},N}\|_2^2,\|c_{\mathrm{par},N}\|_2^2).
\]

Consequentemente, para **todo** `N≥3`, usando as mesmas fórmulas de coeficientes:

\[
\boxed{\|B_N\|_{\mathrm{op}}^2\le\frac{6275281}{3686400}<\frac{171}{100}=1.71.}
\]

A cota anterior da ES8-P17 era `5148881/2560000≈2.01128` (substituída nesta fase por `≈1.70228`). A nova cota **não** pressupõe positividade de uma matriz de Weil.

## 2. Consequência para a cauda arquimediana

A partir da cota **condicional** da ES8-P15 para `R(N)=200 N²`,

\[
\|E_N^{\mathrm{cauda}}\|_{\mathrm{op}}\le\delta_N\le
\frac{406}{91125},
\]

segue a estimativa após projeção:

\[
\|B_N^*E_N^{\mathrm{cauda}}B_N\|_{\mathrm{op}}
\le\frac{171}{100}\,\delta_N
\le \boxed{\frac{3857}{506250}}\approx0.007619
 <\frac1{100}.
\]

Folga racional exata: \(1/100-3857/506250=2411/1012500\). A cota **depende** da validade analítica da estimativa E7 da integral, que *não* foi formalizada aqui. Os demais erros (interpolação, digama, primos, polos, normalização) não estão cobertos.

## 3. Construção efetiva de um operador de correção em \(\ell^2\)

Seja \(H=\ell^2(\{3,4,\ldots\};\mathbb C)\). As sequências `c_odd`, `c_even`, definidas acima, pertencem a `H`, como provado pelas somas majorantes. Logo os funcionais

\[
\lambda_o(x)=\sum_{k\ \mathrm{ímpar}}c_kx_k,\qquad
\lambda_e(x)=\sum_{k\ \mathrm{par}}c_kx_k
\]

convergem absolutamente por Cauchy–Schwarz e são limitados. A aplicação

\[
B_\infty:H\longrightarrow\mathbb C^2\oplus H,\qquad
x\longmapsto(-\lambda_o(x),-\lambda_e(x),x)
\]

é linear limitada, com \(\|B_\infty\|^2\le 6275281/3686400\). Ela preserva a norma inferior: \(\|B_\infty x\|\ge\|x\|\). A diferença \(B_\infty^*B_\infty-I\) é de posto no máximo dois, portanto é compacta e Hilbert–Schmidt. **O operador inteiro \(B_\infty\) não é Hilbert–Schmidt**, porque contém a identidade em espaço infinito-dimensional.

Defina `C_≤N` pela soma somente dos modos `k≤N`, mas mantenha a identidade **em todas as coordenadas** no operador estendido `B_≤N^ext=(-C_≤N,I)`. A desigualdade

\[
\sum_{k>N}c_k^2\le \left(\frac{1609}{800}\right)^2\sum_{k>N}\frac1{k^2}
\le\frac{(1609/800)^2}{N}
\]

usa telescopagem \(1/k^2\le1/[k(k-1)]\). Consequentemente

\[
\boxed{\|B_\infty-B_{\le N}^{\mathrm{ext}}\|_{\mathrm{op}}
\le\frac{1609}{800\sqrt N}.}
\]

Essa é uma **convergência em norma de operador da correção de polos**, não das truncações de posto finito do operador completo. Se `P_N` projeta as coordenadas para as `N` primeiras, então \(\|B_\infty-B_\infty P_N\|\ge1\) em todo `N` (teste com vetor-base na cauda). Já `B_\infty P_N→ B_\infty` **fortemente** por densidade dos suportes finitos e continuidade de `B_\infty`.

## 4. O que foi formalizado e o que permanece em aberto

- **Python com frações exatas:** constantes, limites em diferentes `N`, monotonicidade em `X=π²≥9`, telescopagem finita e testes adversariais.
- **Lean:** constantes exatas, parcela da telescopagem, conversão da desigualdade de Cauchy–Schwarz (assumida explicitamente) em orçamento uniforme, e propagação do orçamento condicional de cauda. O teorema completo de `B_∞` como operador em `ℓ²`, incluindo somabilidade analítica, **ainda não foi formalizado**.
- **Demonstração matemática no relatório:** coeficientes `ℓ²` e o operador infinito de **correção de polos**. Não estabelece uma forma quadrática infinita de Weil nem uma margem espectral uniforme.

### Próxima obrigação decisiva

Para avançar a um teorema sobre a forma de Weil, é preciso **definir uma forma sesquilinear global com domínio comum**, provar que suas restrições reproduzem exatamente todas as entradas `G+P` e controlar os termos de primos para suportes crescentes. Sem essa identidade, uma matriz infinita qualquer com prefixos positivos não fornece o critério universal. A positividade para todas as funções-teste e a Hipótese de Riemann seguem **em aberto**.
