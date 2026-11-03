# NUMERIA–RIEMANN — ES8-P3: bloco arquimediano projetado com caudas racionais e LDL intervalar

**Data:** 03/11/2026  
**Alvo:** modelo matricial real da V17, `L=9/20`, modos `k=1,...,N`, `N=6,8,10`, projeção de anulação dos dois polos (`dim=N−2`).  
**Status:** **CERTIFICAÇÃO FINITA CONDICIONAL** das matrizes *reconstruídas nesta etapa*, sob a fórmula explícita de Guinand–Weil com as convenções V17 e os majorantes analíticos de suavização V24. **Não equivale a uma auditoria binária integral dos arquivos históricos E7/V16, nem a uma prova de RH.**

## 1. Problema resolvido nesta fase

ES8-P2 certificou em frações a parcela dos primos, mas não as entradas arquimedianas. Nesta fase calculamos, para a **mesma** base e projeção,

\[
A_{ij}=\frac{1}{\pi}\int_{0}^{\infty}\Omega(r)\operatorname{Re}(F_i(r)\overline{F_j(r)})\,dr,
\quad \Omega(r)=\operatorname{Re}\psi(1/4+ir/2)-\log\pi,
\]

em que

\[
f_k(x)=\frac{\sin(k\pi(x+L)/(2L))}{\sqrt L}\,\mathbf1_{[-L,L]}(x),
\quad F_k(r)=\int_{-L}^{L}f_k(x)e^{irx}\,dx,
\quad\omega_k=\frac{k\pi}{2L}.
\]

A integração de um produto oscilatório por `mpmath` foi substituída por **uma identidade de séries absolutamente convergentes, com limites racionais externos**. Isso também resolve matematicamente a simetrização encontrada na ES8-P: a entrada real entre modos de paridades opostas é exatamente zero. A matriz-modelo associada usa `A + Prime + Polo` com a convenção V17 e a expressão de `Prime` da ES8-P2.

## 2. Lema novo: núcleo exponencial e fórmula fechada

Na expansão clássica da digama, escreva `q_n=n+1/4`. Então

\[
\Omega(r)=\Omega(0)+\sum_{n\ge0}
\frac{r^2}{q_n(4q_n^2+r^2)},\quad
\Omega(0)=-\gamma-\frac\pi2-3\log2-\log\pi.
\]

Use a identidade, com nossa transformada de Fourier,

\[
\frac{1}{2\pi}\int_{\mathbb R}
\frac{4q}{4q^2+r^2}\,e^{ir(x-y)}\,dr=e^{-2q|x-y|}.
\]

Como

\[
\frac{r^2}{q(4q^2+r^2)}=\frac1q-\frac{4q}{4q^2+r^2},
\]

cada termo da expansão pode ser convertido no núcleo integral `exp(−2q|x−y|)`. A avaliação direta da convolução contra os senos dá, para `q>0`,

\[
\begin{aligned}
K_{ij}(q)&=\iint_{[-L,L]^2}f_i(x)f_j(y)e^{-2q|x-y|}\,dx\,dy\\
&=\delta_{ij}\frac{4q}{4q^2+\omega_i^2}
+\mathbf1_{i\equiv j\ (2)}
\frac{2\omega_i\omega_j}{L}
\frac{1-(-1)^i e^{-4Lq}}{(4q^2+\omega_i^2)(4q^2+\omega_j^2)}.
\end{aligned}
\]

Daí a identidade **exata**:

\[
\boxed{\begin{aligned}
A_{ij}={}&\delta_{ij}\Omega(0)
+\delta_{ij}\sum_{n\ge0}\frac{\omega_i^2}{q_n(4q_n^2+\omega_i^2)}\\
&-\mathbf1_{i\equiv j\ (2)}
\sum_{n\ge0}\frac{2\omega_i\omega_j}{L}
\frac{1-(-1)^i e^{-4Lq_n}}
{(4q_n^2+\omega_i^2)(4q_n^2+\omega_j^2)}.
\end{aligned}}
\]

As duas séries convergem absolutamente (`O(n⁻³)` na primeira e `O(n⁻⁴)` na segunda), inclusive para todas as entradas de `N≤10`. A passagem da integral à soma usa convergência dominada: após o pareamento original a parcela espectral é `O(r²/q³)`, majorável por `r² |F_i F_j|`, que é integrável para estes modos (pertencem a `H^1_0([-L,L])`).

**Nota:** este resultado avalia a **parcela arquimediana do modelo V17**; a fórmula explícita de Guinand–Weil que o interpreta como contribuição da zeta é um teorema externo e a especialização histórica E7/V16 continua sendo um gate.

## 3. Certificação das constantes e da cauda sem integração oscilatória

As rotinas operam com `fractions.Fraction`, escala fixa `10^40`, arredondamento externo a cada operação e intervalos para `π`, `ln 2`, `ln π`, `γ` e `exp(−x)`:

- `π`: identidade de Machin e resto de arctangente alternada, com a orientação de intervalo **corrigida na ES8-P2**.
- `γ`: `H_256−ln 256−1/(512)` mais oito termos Bernoulli da fórmula de Euler–Maclaurin. O resto fica amplificado conservadoramente pela **magnitude do último termo incluído**, `|B16|/(16·256^16)`, sem depender do sinal de um termo seguinte.
- `ζ(s,q)=Σ_{n≥0}(q+n)^{-s}` nas caudas: Euler–Maclaurin até `B6`, com resto absoluto limitado pelo **módulo da última correção** `|B6| (s)_5/[6!·q^{s+5}]`. A desigualdade vem do resto integral com o polinômio de Bernoulli periódico e do seu máximo absoluto.
- Frações racionais: usamos a expansão geométrica em `x_i/q²`, `x_i=ω_i²/4`, com alternância e erro limitado pelo primeiro termo omitido. Para todos os `i,j≤10` e `M≥32`, `(x_i+x_j)/q²<1`.
- Para a parcela `exp(−4Lq)`: Taylor alternado dirigido para argumentos pequenos e, na cauda, a desigualdade conservadora `0≤exp(−x)≤2^(-floor x)`.

O programa exporta todos os extremos da **matriz completa**, não apenas apresentações com 18 casas. As avaliações `mpmath` aparecem **somente nos testes diagnósticos**, nunca na construção dos certificados.

## 4. Projeção, paridade e LDL intervalar

Para `k=3,...,N`, seja `p=1` para `k` ímpar e `p=2` para `k` par. A coluna de projeção é

\[
b_k=f_k-\frac{k}{p}\frac{p^2\pi^2+L^2}{k^2\pi^2+L^2}f_p.
\]

As duas avaliações nos polos `s=0,1` cancelam-se algebricamente em `b_k`; os coeficientes correspondem à fórmula auditada na V17. A matriz real **reconstruída aqui** é

\[
M_N=B^\top(A+\mathrm{Prime})B,
\qquad \mathrm{Prime}_{ij}=-\sqrt2\log2\,C_{ij}^{\rm sym}(\log2).
\]

Paridades opostas dão zeros exatos em ambas as parcelas: blocos `2+2`, `3+3` e `4+4` para `N=6,8,10`. A decomposição LDL sem pivoteamento é efetuada **sobre matrizes de intervalos racionais**, e não sobre autovalores arredondados. Mais forte: todas as decomposições de

\[
\boxed{M_N-\frac1{20}I}\tag{ES8-P3}
\]

possuem limites inferiores de pivôs estritamente positivos. Consequentemente, sob as identidades das entradas adotadas, **cada** `M_N` satisfaz `λ_min(M_N)>0,05`.

| N | dimensão | menor limite inferior dos pivôs de `M_N−0,05I` | limite superior do erro projetado com `δ=10⁻⁶` |
|---:|---:|---:|---:|
| 6 | 4×4 | `0,18162803187279` | `0,000022602484143709` |
| 8 | 6×6 | `0,18162803187279` | `0,000051980451443969` |
| 10 | 8×8 | `0,18162803187279` | `0,000099852324942144` |

A última coluna é a aplicação **condicional**, à mesma projeção `B`, do envelope V24 para as diferenças de suavização: contribuição arquimediana em `[0,∞)`, termos primos e dois termos de polos. O limite da norma espectral decorre da soma absoluta máxima por linha. Como `erro < 0,0001 < 0,05`, as matrizes **suavizadas correspondentes** também permanecem positivas definidas **se** as hipóteses analíticas do envelope V24 e a especialização de Guinand–Weil usada na V17 estiverem satisfeitas.

**Cuidado:** o menor pivô é **não** o menor autovalor. A cota `λ_min>0,05` vem de ter provado LDL positivo para a matriz *deslocada*, `M_N−0,05I`, não de tomar o menor pivô como autovalor.

## 5. Reprodutibilidade e refinamento

Três configurações independentes foram executadas no mesmo algoritmo, mudando o corte exato das séries:

| Corte M | Largura máxima de uma entrada arquimediana projetada (N=10) | LDL da matriz deslocada |
|---:|---:|---|
| 64 | `< 3,240×10⁻8` | aprovado |
| 128 | `< 7,40×10⁻12` | aprovado |
| 256 | `< 1,0×10⁻14` | aprovado |

Os três JSONs estão incluídos. A largura mais fina não é consequência de aumentar a precisão flutuante: resulta da redução analítica de caudas e dos arredondamentos racionais. Os testes `pytest` contêm comparação numérica independente com a digama, núcleo exponencial, Hurwitz zeta, controles adversariais, anulação de polos, sensibilidade do taper e LDL intervalar.

## 6. Cotejo com os documentos E7/V16/V17/V24

| Fonte auditada | Resultado | Pendência |
|---|---|---|
| `v17/RESULTADO_E8K_V17.json`, `v17/AUDITAR_E8K_V17_NORMALIZACAO.py` | Mesma convenção `F_k=∫f_k e^{irx}`, mesma metade positiva `1/π`, mesmo `H_ij=Re(F_i conj F_j)`, `−2 ln2/√2` e projeção de polos | Continua condicionada à fórmula externa de Guinand–Weil |
| `v24/README_V24.md`, `v24/AUDITAR_V24_CONTINUIDADE.py` | Envelope de erro reconduzido à matriz **real simetrizada**, usando `|Re z|≤|z|`; mesmas cotas de `e_k`, `D'_k`, polos e primos | Lemas de derivadas/cortes não foram formalizados em Lean nesta etapa |
| `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip` (SHA256 histórico `4acf99a2...02a7`) e `V16` | V17 reporta SHA e testes anteriores de especialização; registro histórico localizado no GitHub | **Não houve replay byte a byte dos originais E7/V16 nesta fase**; preservar `NOT_ESTABLISHED` para identidade literal dos resultados históricos |

Em particular, o módulo E7 original dentro do repositório é um ZIP binário e seu conteúdo não foi integralmente reaberto no ambiente desta fase. Há forte correspondência documental das convenções **V17**, mas não uma revalidação das matrizes históricas originais. Não reutilizar pivôs históricos como se tivessem sido recalculados.

## 7. Não extrapolação

O resultado é **finito**, para `L=9/20`, `N=6,8,10` e `δ=10^-6` com a base e projeção V17. Não demonstra que a positividade vale para todo `N`, `L`, toda função-teste ou o espaço de Weil. O critério universal de Weil equivale a RH, e essa universalidade **não** é consequência de três matrizes positivas. A fórmula explícita e os fatos de Euler–Maclaurin utilizados são entradas analíticas conhecidas, não provas novas da zeta.

## 8. Referências

- [NIST DLMF, §5.7](https://dlmf.nist.gov/5.7) — expansão da digama e valores especiais.
- [NIST DLMF, §25.11(iii)](https://dlmf.nist.gov/25.11.iii) — representação da zeta de Hurwitz por Euler–Maclaurin e restos de Bernoulli.
- [Bombieri–Lagarias, *Complements to Li's Criterion* (1999)](https://websites.umich.edu/~lagarias/doc/bombieri.pdf) — conexão entre fórmula explícita e positividade de Weil.
- [V17, normalização histórica](https://github.com/jbarumby/scicalc/blob/main/v17/RESULTADO_E8K_V17.json).
- [V24, continuidade e suavização](https://github.com/jbarumby/scicalc/blob/main/v24/README_V24.md).
- [ES8-P2](https://github.com/jbarumby/scicalc/tree/43722ffe48ee6b9bb78c32f0890d4c6518330d15/research/2026-11-03/es8p2) — parcela de primos em intervalos racionais.

**Conclusão:** certificação finita de autovalor inferior e estabilidade sob suavização **dentro do modelo V17**, com condições externas declaradas; não uma nova prova da Hipótese de Riemann.
