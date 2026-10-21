# NUMERIA–RIEMANN: prova da cota da cauda arquimediana

*Tradução para português de [TAIL_PROOF.md](TAIL_PROOF.md). Em caso de divergência, vale a versão em inglês.*

**Data:** 2026-10-21. Este arquivo prova o enclosure analítico da cauda usado em `weil_interval_audit.py::_tail_bound`. Ele prova **somente** esse lema. A dedução da própria forma quadrática de Weil (normalização da fórmula explícita, o fator $1/\pi$, o termo de primos $-2\Lambda(n)/\sqrt n$, as restrições de polo) vem do protótipo e é conferida numericamente, não provada aqui (ver §6).

## 0. Contexto

Sejam $L>0$, $w_k = k\pi/(2L)$, $g_k(x)=\sin(w_kx)/\sqrt L$ em $[0,2L]$ e $h_k(t)=g_k(t+L)$ em $[-L,L]$. Escrevemos
$$F_k(r)=\int_{-L}^{L}h_k(t)\,e^{irt}\,dt,\qquad \Omega(r)=\operatorname{Re}\psi\!\left(\tfrac14+\tfrac{ir}{2}\right)-\log\pi,$$
e $G_{ij}=\frac1\pi\int_0^\infty \operatorname{Re}\!\big(F_i\overline{F_j}\big)\,\Omega\,dr$. O código calcula $\int_0^R$ por quadratura certificada e envolve $\int_R^\infty$ pela cota abaixo. Hipóteses verificadas em aritmética de bolas: $R\ge 2w_{\max}=N\pi/L$ e $R\ge 1$.

## 1. Lema 1 (forma fechada de $F_k$)

$$F_k(r)=i\,e^{-ik\pi/2}A_k(r),\qquad A_k(r)=\frac{2\sqrt L\,w_k}{w_k+r}\,\operatorname{sinc}\!\big(L(r-w_k)\big),\quad \operatorname{sinc}x=\frac{\sin x}{x}.$$

*Prova.* Com $x=t+L$, $F_k(r)=\frac{e^{-irL}}{\sqrt L}\int_0^{2L}\sin(w_kx)e^{irx}dx$. Escreva $\sin=(e^{iwx}-e^{-iwx})/2i$ e use $e^{\pm 2iw_kL}=e^{\pm ik\pi}=(-1)^k$. Com $c=(-1)^ke^{2irL}$,
$$\int_0^{2L}\sin(w_kx)e^{irx}dx=\frac{c-1}{2i\cdot i}\Big(\frac1{r+w_k}-\frac1{r-w_k}\Big)=\frac{(c-1)\,w_k}{r^2-w_k^2}.$$
Além disso, $e^{-irL}(c-1)=(-1)^ke^{irL}-e^{-irL}=2i\,e^{ik\pi/2}\sin\!\big(L(r+w_k)\big)$, e $\sin(L(r+w_k))=(-1)^k\sin(L(r-w_k))$, porque $2Lw_k=k\pi$. Logo
$F_k(r)=\frac{2i\,w_k}{\sqrt L}\,e^{ik\pi/2}e^{ik\pi}\frac{\sin(L(r-w_k))}{(r-w_k)(r+w_k)}=i e^{3ik\pi/2}A_k(r)=ie^{-ik\pi/2}A_k(r)$. $\square$

*Consequência.* $\operatorname{Re}(F_i\overline{F_j})=\cos\!\big(\tfrac{(i-j)\pi}2\big)A_iA_j$: é **identicamente nula** quando $i-j$ é ímpar, e vale $(-1)^{(i-j)/2}A_iA_j$ quando $i-j$ é par. Isso justifica os zeros exatos (inclusive na cauda) e a função `_phase`.

## 2. Lema 2 (decaimento de $A_k$)

Se $r\ge 2w_k$, então $|A_k(r)|\le \dfrac{4w_k}{\sqrt L\,r^2}$.

*Prova.* $|\operatorname{sinc}x|\le 1/|x|$. Para $r\ge 2w_k$: $r-w_k\ge r/2$, logo $|\operatorname{sinc}(L(r-w_k))|\le \frac{2}{Lr}$, e $w_k+r\ge r$. Assim $|A_k|\le \frac{2\sqrt L\,w_k}{r}\cdot\frac{2}{Lr}$. $\square$

## 3. Lema 3 (série para $\Omega$)

Para $r\ge 0$, com $y_n=n+\tfrac12$ e $S(r)=\sum_{n\ge0}\dfrac{r^2}{y_n(y_n^2+r^2)}$,
$$\operatorname{Re}\psi\!\left(\tfrac14+\tfrac{ir}2\right)=-\gamma-3\log2+S(r)-\frac{\pi}{2\cosh\pi r}.$$

*Prova.* De $\psi(z)=-\gamma+\sum_{n\ge0}\big(\tfrac1{n+1}-\tfrac1{n+z}\big)$ e $\psi(\tfrac12)=-\gamma-2\log2$: $\operatorname{Re}\psi(\tfrac12+ir)=\psi(\tfrac12)+\sum_n\big(\tfrac1{y_n}-\tfrac{y_n}{y_n^2+r^2}\big)=-\gamma-2\log2+S(r)$.
Seja $z=\tfrac14+\tfrac{ir}2$. A fórmula de duplicação dá $\psi(\tfrac12+ir)=\tfrac12\psi(z)+\tfrac12\psi(z+\tfrac12)+\log2$. Como $z+\tfrac12=\overline{1-z}$, a reflexão dá $\operatorname{Re}\psi(z+\tfrac12)=\operatorname{Re}\psi(z)+\pi\operatorname{Re}\cot\pi z$, e $\operatorname{Re}\cot(a+ib)=\frac{\sin2a}{\cosh2b-\cos2a}$, com $a=\tfrac\pi4$ e $b=\tfrac{\pi r}2$, vale $1/\cosh\pi r$. Então $\operatorname{Re}\psi(\tfrac12+ir)=\operatorname{Re}\psi(z)+\frac{\pi}{2\cosh\pi r}+\log2$. Combine as duas identidades. $\square$

## 4. Lema 4 (cota para $\Omega$)

Para $r\ge1$: $\;-3{,}14\le\Omega(r)\le\log r-0{,}99$, e portanto $|\Omega(r)|\le\log r+3{,}14\le 8+\log(3r)$.

*Prova.* A função $f(y)=r^2/(y(y^2+r^2))$ é positiva e decrescente em $y>0$, com $\int f\,dy=-\tfrac12\log(1+r^2/y^2)$. A comparação série–integral dá
$$\tfrac12\log(1+4r^2)=\int_{1/2}^\infty f\;\le\; S(r)\;\le\; f(\tfrac12)+\int_{1/2}^\infty f\;<\;2+\tfrac12\log(1+4r^2),$$
pois $f(\tfrac12)=2r^2/(r^2+\tfrac14)<2$. Para $r\ge1$: $\tfrac12\log5\le\tfrac12\log(1+4r^2)\le\log r+\tfrac12\log5$ (já que $1+4r^2\le5r^2$), com $0{,}8047<\tfrac12\log5<0{,}8048$. Além disso, $0<\frac{\pi}{2\cosh\pi r}\le\frac{\pi}{2\cosh\pi}<0{,}1356$ e $\gamma+3\log2+\log\pi\in(3{,}8013;\,3{,}8014)$. Pelo Lema 3, $\Omega=-(\gamma+3\log2+\log\pi)+S-\frac{\pi}{2\cosh\pi r}$, então
$$\Omega\ge-3{,}8014+0{,}8047-0{,}1356>-3{,}14,\qquad \Omega<-3{,}8013+2+\log r+0{,}8048<\log r-0{,}99.$$
Como $\log r\ge0$, $|\Omega|\le\max(3{,}14,\log r)\le\log r+3{,}14$. $\square$

(A constante $8$ no código é, portanto, muito conservadora; ver §7.)

## 5. Proposição (cauda de cada entrada)

Seja $R\ge\max(1,2w_{\max})$ com $w_{\max}=N\pi/(2L)$. Para $1\le i,j\le N$ com $i-j$ par,
$$\Big|\tfrac1\pi\!\int_R^\infty\!\operatorname{Re}(F_i\overline{F_j})\,\Omega\,dr\Big|\le\frac{16\,w_iw_j}{\pi L}\cdot\frac{8+\log(3R)+\tfrac13}{3R^3},$$
e a cauda é exatamente $0$ quando $i-j$ é ímpar.

*Prova.* Para $r\ge R\ge 2w_{\max}$, o Lema 2 dá $|A_iA_j|\le 16w_iw_j/(Lr^4)$ e o Lema 4 dá $|\Omega|\le 8+\log 3r$. Então
$$\int_R^\infty\frac{8+\log3r}{r^4}dr=\frac{8}{3R^3}+\Big[-\frac{\log 3r}{3r^3}\Big]_R^\infty+\int_R^\infty\frac{dr}{3r^4}=\frac{8+\log3R}{3R^3}+\frac1{9R^3}=\frac{8+\log 3R+\tfrac13}{3R^3}.$$
O caso ímpar é consequência do Lema 1. $\square$

A cauda entra como raio simétrico $\pm$(cota) em cada entrada de $G$ antes da projeção no subespaço restrito, e portanto se propaga pelo produto intervalar $B^{\mathsf T}WB$.

## 6. O que isto não prova

- A forma quadrática da fórmula explícita (normalização, sinais, termo de primos, restrições de polo) não é deduzida aqui. Ela é conferida numericamente: `test_quadratic_form_equals_sum_over_real_zeta_zeros` (`NUMERIA_SLOW=1`) compara $c^{\mathsf T}Ac$ mais a cauda arquimediana com $\sum_\gamma|F_c(\gamma)|^2$ sobre zeros reais da zeta, com concordância de cerca de $10^{-8}$. Isso é evidência, não prova.
- A correção da avaliação do Arb é assumida (aritmética de bolas com arredondamento para fora); as pré-condições $R\ge1$ e $R\ge2w_{\max}$ são verificadas com cotas de bola no código.
- A nota anterior em `WEIL_INTERVAL_AUDIT_CORRECTED_RESULTS_2026-10-20.md` cita outra cota, $|\Omega|\le\log r+\log\pi+\tfrac4{3r}$. Ela não é usada nem verificada aqui; o Lema 4 é a cota em que o código se apoia.
- Esta prova é nova e precisa de revisão independente.

## 7. Refinamento opcional

O Lema 4 permite o mesmo argumento com $8$ trocado por $3{,}14$: o fator da cauda passa a ser $(3{,}14+\log3R+\tfrac13)/(3R^3)$, cerca de $1{,}7$ vez menor em $R=160$. Como a cauda agora domina a largura restante dos enclosures, aumentar `--R` (a cauda cai como $R^{-3}\log R$) é a alavanca mais barata.
