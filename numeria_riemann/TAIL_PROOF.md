# NUMERIA–RIEMANN: proof of the archimedean tail bound

*Versão em português: [TAIL_PROOF_PT.md](TAIL_PROOF_PT.md).*

**Date:** 2026-10-21. This file proves the analytic tail enclosure used in `weil_interval_audit.py::_tail_bound`. It proves **only** that lemma. The derivation of the Weil quadratic form itself (explicit-formula normalization, the factor $1/\pi$, the prime term $-2\Lambda(n)/\sqrt n$, the pole constraints) is taken from the prototype and is cross-checked numerically, not proved here (see §6).

## 0. Setting

Let $L>0$, $w_k = k\pi/(2L)$, $g_k(x)=\sin(w_kx)/\sqrt L$ on $[0,2L]$, and $h_k(t)=g_k(t+L)$ on $[-L,L]$. Write
$$F_k(r)=\int_{-L}^{L}h_k(t)\,e^{irt}\,dt,\qquad \Omega(r)=\operatorname{Re}\psi\!\left(\tfrac14+\tfrac{ir}{2}\right)-\log\pi,$$
and $G_{ij}=\frac1\pi\int_0^\infty \operatorname{Re}\!\big(F_i\overline{F_j}\big)\,\Omega\,dr$. The code computes $\int_0^R$ with certified quadrature and encloses $\int_R^\infty$ by the bound below. Assumptions checked in ball arithmetic: $R\ge 2w_{\max}=N\pi/L$ and $R\ge 1$.

## 1. Lemma 1 (closed form of $F_k$)

$$F_k(r)=i\,e^{-ik\pi/2}A_k(r),\qquad A_k(r)=\frac{2\sqrt L\,w_k}{w_k+r}\,\operatorname{sinc}\!\big(L(r-w_k)\big),\quad \operatorname{sinc}x=\frac{\sin x}{x}.$$

*Proof.* Substituting $x=t+L$, $F_k(r)=\frac{e^{-irL}}{\sqrt L}\int_0^{2L}\sin(w_kx)e^{irx}dx$. Write $\sin=(e^{iwx}-e^{-iwx})/2i$ and use $e^{\pm 2iw_kL}=e^{\pm ik\pi}=(-1)^k$. With $c=(-1)^ke^{2irL}$,
$$\int_0^{2L}\sin(w_kx)e^{irx}dx=\frac{c-1}{2i\cdot i}\Big(\frac1{r+w_k}-\frac1{r-w_k}\Big)=\frac{(c-1)\,w_k}{r^2-w_k^2}.$$
Moreover $e^{-irL}(c-1)=(-1)^ke^{irL}-e^{-irL}=2i\,e^{ik\pi/2}\sin\!\big(L(r+w_k)\big)$, and $\sin(L(r+w_k))=(-1)^k\sin(L(r-w_k))$ because $2Lw_k=k\pi$. Hence
$F_k(r)=\frac{2i\,w_k}{\sqrt L}\,e^{ik\pi/2}e^{ik\pi}\frac{\sin(L(r-w_k))}{(r-w_k)(r+w_k)}=i e^{3ik\pi/2}A_k(r)=ie^{-ik\pi/2}A_k(r)$. $\square$

*Consequence.* $\operatorname{Re}(F_i\overline{F_j})=\cos\!\big(\tfrac{(i-j)\pi}2\big)A_iA_j$: it is **identically zero** for $i-j$ odd, and $(-1)^{(i-j)/2}A_iA_j$ for $i-j$ even. This justifies the exact zeros (also in the tail) and `_phase`.

## 2. Lemma 2 (decay of $A_k$)

If $r\ge 2w_k$, then $|A_k(r)|\le \dfrac{4w_k}{\sqrt L\,r^2}$.

*Proof.* $|\operatorname{sinc}x|\le 1/|x|$. For $r\ge 2w_k$: $r-w_k\ge r/2$, so $|\operatorname{sinc}(L(r-w_k))|\le \frac{2}{Lr}$, and $w_k+r\ge r$. Thus $|A_k|\le \frac{2\sqrt L\,w_k}{r}\cdot\frac{2}{Lr}$. $\square$

## 3. Lemma 3 (series for $\Omega$)

For $r\ge 0$, with $y_n=n+\tfrac12$ and $S(r)=\sum_{n\ge0}\dfrac{r^2}{y_n(y_n^2+r^2)}$,
$$\operatorname{Re}\psi\!\left(\tfrac14+\tfrac{ir}2\right)=-\gamma-3\log2+S(r)-\frac{\pi}{2\cosh\pi r}.$$

*Proof.* From $\psi(z)=-\gamma+\sum_{n\ge0}\big(\tfrac1{n+1}-\tfrac1{n+z}\big)$ and $\psi(\tfrac12)=-\gamma-2\log2$: $\operatorname{Re}\psi(\tfrac12+ir)=\psi(\tfrac12)+\sum_n\big(\tfrac1{y_n}-\tfrac{y_n}{y_n^2+r^2}\big)=-\gamma-2\log2+S(r)$.
Let $z=\tfrac14+\tfrac{ir}2$. Duplication gives $\psi(\tfrac12+ir)=\tfrac12\psi(z)+\tfrac12\psi(z+\tfrac12)+\log2$. Since $z+\tfrac12=\overline{1-z}$, reflection gives $\operatorname{Re}\psi(z+\tfrac12)=\operatorname{Re}\psi(z)+\pi\operatorname{Re}\cot\pi z$, and $\operatorname{Re}\cot(a+ib)=\frac{\sin2a}{\cosh2b-\cos2a}$ with $a=\tfrac\pi4$, $b=\tfrac{\pi r}2$ equals $1/\cosh\pi r$. So $\operatorname{Re}\psi(\tfrac12+ir)=\operatorname{Re}\psi(z)+\frac{\pi}{2\cosh\pi r}+\log2$. Combine. $\square$

## 4. Lemma 4 (bound for $\Omega$)

For $r\ge1$: $\;-3.14\le\Omega(r)\le\log r-0.99$, hence $|\Omega(r)|\le\log r+3.14\le 8+\log(3r)$.

*Proof.* $f(y)=r^2/(y(y^2+r^2))$ is positive and decreasing in $y>0$, with $\int f\,dy=-\tfrac12\log(1+r^2/y^2)$. Integral comparison gives
$$\tfrac12\log(1+4r^2)=\int_{1/2}^\infty f\;\le\; S(r)\;\le\; f(\tfrac12)+\int_{1/2}^\infty f\;<\;2+\tfrac12\log(1+4r^2),$$
because $f(\tfrac12)=2r^2/(r^2+\tfrac14)<2$. For $r\ge1$: $\tfrac12\log5\le\tfrac12\log(1+4r^2)\le\log r+\tfrac12\log5$ (as $1+4r^2\le5r^2$), with $\tfrac12\log5<0.8048$ and $>0.8047$. Also $0<\frac{\pi}{2\cosh\pi r}\le\frac{\pi}{2\cosh\pi}<0.1356$, and $\gamma+3\log2+\log\pi\in(3.8013,3.8014)$. By Lemma 3, $\Omega=-(\gamma+3\log2+\log\pi)+S-\frac{\pi}{2\cosh\pi r}$, so
$$\Omega\ge-3.8014+0.8047-0.1356>-3.14,\qquad \Omega<-3.8013+2+\log r+0.8048<\log r-0.99.$$
Since $\log r\ge0$, $|\Omega|\le\max(3.14,\log r)\le\log r+3.14$. $\square$

(The constant $8$ in the code is therefore very conservative; see §7.)

## 5. Proposition (tail of each entry)

Let $R\ge\max(1,2w_{\max})$ with $w_{\max}=N\pi/(2L)$. For $1\le i,j\le N$ with $i-j$ even,
$$\Big|\tfrac1\pi\!\int_R^\infty\!\operatorname{Re}(F_i\overline{F_j})\,\Omega\,dr\Big|\le\frac{16\,w_iw_j}{\pi L}\cdot\frac{8+\log(3R)+\tfrac13}{3R^3},$$
and the tail is exactly $0$ when $i-j$ is odd.

*Proof.* For $r\ge R\ge 2w_{\max}$, Lemma 2 gives $|A_iA_j|\le 16w_iw_j/(Lr^4)$; Lemma 4 gives $|\Omega|\le 8+\log 3r$. Then
$$\int_R^\infty\frac{8+\log3r}{r^4}dr=\frac{8}{3R^3}+\Big[-\frac{\log 3r}{3r^3}\Big]_R^\infty+\int_R^\infty\frac{dr}{3r^4}=\frac{8+\log3R}{3R^3}+\frac1{9R^3}=\frac{8+\log 3R+\tfrac13}{3R^3}.$$
The odd case is Lemma 1's consequence. $\square$

The tail is added as a symmetric radius $\pm$(bound) to each entry of $G$ before projection onto the constrained subspace, so it propagates through the interval product $B^{\mathsf T}WB$.

## 6. What this does not prove

- The explicit-formula quadratic form (normalization, signs, the prime term, the pole constraints) is not derived here. It is cross-checked numerically: `test_quadratic_form_equals_sum_over_real_zeta_zeros` (`NUMERIA_SLOW=1`) matches $c^{\mathsf T}Ac$ plus the archimedean tail to $\sum_\gamma|F_c(\gamma)|^2$ over real zeta zeros to about $10^{-8}$. That is evidence, not proof.
- Correctness of Arb evaluation is assumed (outward-rounded ball arithmetic); the preconditions $R\ge1$ and $R\ge2w_{\max}$ are checked with ball bounds in the code.
- The earlier note in `WEIL_INTERVAL_AUDIT_CORRECTED_RESULTS_2026-10-20.md` quotes a different bound, $|\Omega|\le\log r+\log\pi+\tfrac4{3r}$. It is not used and not verified here; Lemma 4 is the bound the code relies on.

## 7. Optional tightening

Lemma 4 implies the same argument with $8$ replaced by $3.14$: the tail factor becomes $(3.14+\log3R+\tfrac13)/(3R^3)$, about $1.7\times$ smaller at $R=160$. Since the tail now dominates the remaining enclosure width, raising `--R` (the tail scales as $R^{-3}\log R$) is the cheaper lever.
