# NUMERIA–RIEMANN · P259 — prova racional para uma testemunha ímpar

## 1. Premissas autenticadas

Entrada: certificado original P258_GL60.json do [GitHub Actions P258](https://github.com/jbarumby/scicalc/actions/runs/44858303395/artifacts/14688468811). SHA-256: **1aa765357cfe1fa46f864d5b57fc407157375e75b3fa8b23b77f0db5c2e4985c**.

A P258 forneceu um vetor ímpar de 29 coordenadas diádicas e um quociente de Rayleigh **intervalar exato** para o auxiliar truncado de quatro primos (corte 504, piso 1/4), na meia-largura L = 0,972955074527656652552676371721590. Os extremos do Rayleigh são próximos de −1,73007525684909 × 10⁻²⁴. O auditor da P259 lê o **extremo inferior** original por Fraction, sem confiar nos decimais resumidos.

## 2. Piso válido por faixas

O limite de Binet usado nas P249/P250 fornece, para t ≥ 504,

\[
\Psi_4(t)\ge \beta+\log(t/504)+1/504-1/t,\qquad \beta>0.
\]

No ponto 850, uma desigualdade logarítmica elementar garante

\[
\log(850/504)>\frac{2(850-504)}{850+504}=\frac{346}{677}>\frac12.
\]

O piso 1/4, anteriormente válido desde 504, pode portanto ser elevado para 1/2 a partir de 850. Como forma quadrática na janela, a atualização é

\[
R_{\mathrm{P259}}=R_{\mathrm{P258}}+\frac14 P_{|t|>850}
=R_{\mathrm{P258}}+\frac14(I-G_{850}).
\]

## 3. Representação trigonométrica exata de Bessel

Para n=2i+1, a transformada do modo ímpar ortonormal é, à parte de uma fase imaginária comum, \(u_n(t)=(-1)^i\sqrt{2L(2n+1)}\,j_n(Lt)\).

Os polinômios inteiros \(P_n,Q_n\) satisfazem

\[
j_n(z)=\frac{P_n(1/z)\sin(z)+Q_n(1/z)\cos(z)}{z},
\]
\[
P_0=1,\quad Q_0=0,\quad P_1(y)=y,\quad Q_1(y)=-1,
\]
\[
P_{n+1}=(2n+1)yP_n-P_{n-1},\qquad
Q_{n+1}=(2n+1)yQ_n-Q_{n-1}.
\]

Para a testemunha diádica \(v\), definimos as amplitudes \(A,B\), combinações dos polinômios e dos fatores \(\sqrt{2n+1}\), de modo que

\[
F_v(t)=\frac{\sqrt{2L}}{Lt}
\big[A(1/(Lt))\sin(Lt)+B(1/(Lt))\cos(Lt)\big].
\]

O código envolve cada raiz quadrada pelo piso racional de raiz inteira com 320 bits e propaga **racionalmente** o erro \(<2^{-320}\). O sinal não depende de funções trigonométricas avaliadas por float.

## 4. Energia em 23 períodos

Para \(k=0,\ldots,22\), utilize os intervalos de um período inteiro:

\[
I_k=\left[850+\frac{2\pi k}{L},850+\frac{2\pi(k+1)}{L}\right].
\]

A identidade de Machin \(\pi=16\arctan(1/5)-4\arctan(1/239)\), junto da série alternada de arctan, fornece limites racionais certificados para os extremos desses períodos. Todos ficam abaixo de 998,531.

Em cada intervalo um ponto racional \(y_k\) fixa as amplitudes \(A_k,B_k\). Limites de derivada dos polinômios e do erro dos radicais dão \(\epsilon_{A,k},\epsilon_{B,k}\). O auditor verifica exatamente

\[
1024(\epsilon_{A,k}+\epsilon_{B,k})^2 < A_k^2+B_k^2
\]

em **todos os 23 períodos**.

A identidade de energia para seno e cosseno em um período completo é

\[
\int_{I_k}|A_k\sin Lt+B_k\cos Lt|^2dt
=\frac{\pi}{L}(A_k^2+B_k^2).
\]

Aplicar a desigualdade triangular em \(L^2\), com \(\sqrt2<2\), conserva ao menos \((15/16)^2=225/256\) desta energia. Pela normalização angular de Plancherel,

\[
M_{850}(v)=\frac{1}{\pi\|v\|^2}\int_{850}^{\infty}|F_v(t)|^2dt
\ge\sum_{k=0}^{22}\frac{225(A_k^2+B_k^2)}
{128L^2 t_{k+1}^2\|v\|^2}
>9{,}5780\times10^{-24}.
\]

O último limite decorre de uma soma **exata de frações**. O arquivo P259_RESULTADO.json registra um piso racional arredondado para baixo.

## 5. Conclusão e limitação

Tomando o limite **inferior** do Rayleigh antigo e acrescentando \(M_{850}(v)/4\), a auditoria conclui

\[
\boxed{\frac{\langle v,R_{\mathrm{P259}}v\rangle}{\|v\|^2}
>6{,}6443\times10^{-25}>0.}
\]

Portanto a **direção antiga** de 29 coordenadas deixa de ser negativa. Este resultado **não garante positividade do setor ímpar 377×377**, pois outra direção pode ser negativa. Também não certifica extensão do domínio da forma completa de Weil nem constitui prova de RH.

A auditoria foi executada localmente com a entrada P258 autenticada, mais 12 testes normais e 12 com otimização Python. Para reprodução no GitHub, baixe a mesma entrada do Actions e execute o módulo P259. O próximo trabalho (Gate B) é calcular G850 com Arb e efetuar LDL/Schur do setor completo ou obter uma nova testemunha negativa.
