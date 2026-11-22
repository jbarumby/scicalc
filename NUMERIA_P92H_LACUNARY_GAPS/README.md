# NUMERIA–RIEMANN P92-H — Limite de Δ para lacunas com crescimento pelo menos geométrico

**Base:** [P92-G, PR #96](https://github.com/jbarumby/scicalc/pull/96). Nova prova analítica para **razões variáveis**, sem a identidade Toeplitz exata de uma progressão geométrica. Esta nota pertence à trilha espectral P88/P92-D–G e não resolve o orçamento interior de cinco lacunas da P93 geométrica.

## 1. Teorema demonstrado (todos os parâmetros reais, dimensão arbitrária)

Sejam `N>=2` e pontos reais `x_0<...<x_(N-1)`. Defina as lacunas sucessivas `g_j=x_(j+1)-x_j>0`, `j=0,...,N-2`, e os parâmetros de vizinho mais próximo

\[
\delta_i=\min_{j\ne i}|x_i-x_j|.
\]

Se as lacunas satisfazem **para todo índice consecutivo**

\[
\boxed{g_{j+1}\ge\frac{64}{25}\,g_j,}\qquad j=0,\dots,N-3,\tag{H0}
\]

a matriz real simétrica, com diagonal zero,

\[
\Delta_{ij}=\frac{\sqrt{\delta_i\delta_j}
 (\sqrt{\delta_i}-\sqrt{\delta_j})^2}{2(x_i-x_j)^2}\quad(i\ne j)
\]

satisfaz o **teorema P92-H.1**:

\[
\boxed{\|\Delta\|_2\le\frac{12793}{13824}\approx0.9254195601851852
<\frac{15}{16}=0.9375.}\tag{H1}
\]

É válido para TODOS os `N`, todas as lacunas **reais** positivas que cumprem (H0), e razões consecutivas que podem variar sem limite superior. A cota não depende de `N`, nem do contraste entre escalas, nem da uniformidade das razões. É uma condição **suficiente**, não se afirma que a razão 64/25 ou a constante 12793/13824 sejam ótimas.

Isso demonstra uma classe mais ampla que a P92-G (progressões geométricas exatas), mas **não cobre todas as razões q>1** com o limiar 15/16.

## 2. Estrutura do núcleo e estimativa escalar

Como as lacunas crescem, cada ponto tem os seguintes parâmetros:

\[
\delta_0=\delta_1=g_0,\qquad \delta_i=g_{i-1}\ (i\ge2).\tag{H2}
\]

Para `1<=i<j<=N-1`, ponha `k=j-i` e

\[
s=\sqrt{\frac{\delta_i}{\delta_j}}\in(0,1].
\]

De (H0) e (H2), segue `s<=(5/8)^k`. Além disso,

\[
x_j-x_i=\sum_{\ell=i}^{j-1}g_\ell\ge g_{j-1}=\delta_j.
\]

Substituindo o denominador maior que `delta_j`, obtemos a estimativa, **sem qualquer hipótese de razão constante**:

\[
0\le\Delta_{ij}\le\frac12s(1-s)^2
\le\frac12\min\left\{\frac4{27},\left(\frac58\right)^k\right\}.\tag{H3}
\]

A desigualdade para o máximo `4/27` vem da função `phi(s)=s(1-s)^2`: sua derivada é `(1-s)(1-3s)` e seu máximo em `[0,1]` ocorre em `s=1/3`; também `phi(s)<=s`.

## 3. A linha âncora e a norma espectral

A igualdade `delta_0=delta_1=g_0` implica

\[
\Delta_{01}=0,\qquad 0\le\Delta_{0j}\le\Delta_{1j}\quad(j\ge2),\tag{H4}
\]

porque os numeradores são idênticos e `x_j-x_0>x_j-x_1`.

Defina `m_k=min(4/27,(5/8)^k)`. Por (H3), as somas de linha do bloco de índices `1,...,N-1` são no máximo

\[
S:=\sum_{k\ge1}m_k
\le\sum_{k=1}^{3}\frac4{27}
+\sum_{k=4}^{\infty}\left(\frac58\right)^k
=\frac49+\frac{625}{1536}.\tag{H5}
\]

Pela comparação com a linha 1, a soma da linha 0 é no máximo `S/2`. Cada outra linha recebe no máximo uma entrada adicional da linha 0, limitada por `m_(j-1)/2<=2/27`. Logo

\[
\max_i\sum_{j\ne i}\Delta_{ij}\le S+\frac2{27}
\le\frac49+\frac{625}{1536}+\frac2{27}
=\boxed{\frac{12793}{13824}}.\tag{H6}
\]

A matriz `Delta` é real simétrica e tem entradas não negativas. Pelo limite de Schur/Gershgorin, `||Delta||_2` não excede a maior soma de linha. Como

\[
\frac{15}{16}-\frac{12793}{13824}=\frac{167}{13824}>0,
\]

segue (H1), **uniformemente em N**, como queríamos.

## 4. Corolário para todo crescimento mínimo q>1 (cota dependente de q)

Mesmo sem a restrição `q>=64/25`, o argumento fornece uma **estimativa dimension-free para toda q>1**: suponha `g_(j+1)>=q*g_j`, escreva `t=q^(-1/2)` e escolha qualquer inteiro `J>=0`. Então

\[
\boxed{\|\Delta\|_2\le \frac{4J}{27}+\frac{t^{J+1}}{1-t}+\frac2{27}.}\tag{H7}
\]

Isso decorre de usar `4/27` para os primeiros `J` valores de `k` e `t^k` na cauda, com o controle da âncora de (H4). O limite é válido para **todo real q>1 e todo N**, mas pode ser grande para `q` próximo de 1, e **não** certifica o limiar 15/16 nesse regime. Tomar `J=3,t<=5/8` dá (H1).

**Não confundir:** a P92-G obtém 15/16 até para q próximo de 1 **somente quando a razão é constante**. A P92-H generaliza para razões variáveis sob a condição quantitativa (H0).

## 5. Certificados racionais independentes

O programa [certificate.py](certificate.py) emprega apenas `fractions.Fraction` e testa famílias de lacunas `g_j=r_j^2`, com razões de raízes consecutivas `r_(j+1)/r_j>=8/5`. Assim, as entradas da matriz são frações **exatas**, sem aproximação de radicais.

- **71 configurações** com razões alternadas, crescimentos próximos do limite, saltos grandes e dimensões até 28.
- **35 configurações** tiveram fatorações LDL exatas tanto de `(12793/13824)I-Delta` quanto de `(12793/13824)I+Delta`: **70 verificações de definitude positiva** com frações racionais, independentes das somas de linha.
- Todos os exemplos verificam identidades de vizinho mais próximo, limite par-a-par (H3), dominação (H4) e cota de Schur (H6).
- O JSON preserva os valores exatos das razões e margens, os pivôs racionais LDL e limitações do escopo. Essas verificações são **regressões**, e não substituem a demonstração analítica acima, que vale para parâmetros reais arbitrários.

Reprodução da raiz do repositório:

    python -m unittest -v NUMERIA_P92H_LACUNARY_GAPS.test_p92h
    python -O -m unittest -v NUMERIA_P92H_LACUNARY_GAPS.test_p92h
    python -m NUMERIA_P92H_LACUNARY_GAPS.certificate --out p92h-lacunary-certificate.json

A CI executa dois conjuntos de 8 testes, verifica o JSON com `Fraction`, gera `SHA256SUMS` e disponibiliza tudo como artefato.

## 6. Alcance, fronteiras e próxima prioridade

- **Novo teorema efetivamente demonstrado:** cota uniforme <15/16 para **razões não constantes** com `g_(j+1)/g_j>=64/25`, sem limite do número de pontos nem das escalas, e cota dependente de `q` para todo `q>1`.
- **Não demonstrado:** limite 15/16 para todo `q>1` em razões variáveis, nenhum limite universal para configurações completamente arbitrárias sem condição de crescimento, Lemas 5–6 gerais de Preissmann, P93 geométrica de cinco lacunas, ou Hipótese de Riemann.
- **Próxima direção:** tratar blocos de crescimento lento (`q` próximo de 1) e acoplamentos entre trechos lacunários; a soma direta das constantes por bloco não é uma prova global. Procurar estimativas espectrais que levem em conta sobreposição e ortogonalidade, sem extrapolar amostras.

A PR permanece **draft**; não fazer merge automático.