# NUMERIA–RIEMANN ES8-P82 — Auditoria analítica independente

Esta etapa audita a quadratura Gauss–Legendre/Bernstein que foi executada no GitHub Actions nas janelas a=0,30/0,34/0,40 com 24/26/44 modos de Legendre por partes. **Não utiliza gaussianas**. Os dados integrais e artefatos permanecem no pacote ES8-P81 e na revisão de pesquisa.

## Lema uniforme do integrando

Dividir [0,30] em doze trechos, cada qual com meia-largura L=5/4. Para rho=34/25, os semieixos das elipses escaladas são A=1781/1360, B=531/1360. Todos os pontos satisfazem |Im r|<=B<2/5 e |r|<=40881/1360<31. As duas digamas psi(1/4 ± i r/2) têm Re(w)>=149/2720>1/19.

A expansão psi(w)=-gamma+(w-1)*sum_{n>=0} 1/((n+1)(n+w)) implica |psi(w)|<=1+(|w|+1)(1/sigma+2). Logo |psi(1/4+15i)|<103, |[psi(1/4+ir/2)+psi(1/4-ir/2)]/2|<358. O fator digama do integrando é, em módulo, <461.

A identidade integral de Bessel esférica e |P_n(t)|<=1 em [-1,1] provam |j_n(z)|<=cosh(|Im z|). Para cada um dos dois Bessel e o seno/cosseno de fase, |Im argumento|<2/5, portanto |fator|<25/23. Consequentemente, para **todos** os pares das três bases:

    |f_ij(r)| < 461*(25/23)^3 = 7203125/12167 < 600.

## Erro analítico

A soma da cauda de Chebyshev acima do grau 2n-1 é <= 2M rho^(1-2n)/(rho-1). O integral e a regra de Gauss–Legendre têm norma 2L, e a regra é exata até o grau 2n-1. Logo:

    |integral - GL_n| <= 8*L*M*rho^(1-2n)/(rho-1).

Com M=600 e n=110, a soma dos erros nos 12 trechos é **< 1,2×10^-24** por integral. As comparações são racionais exatas e independentes do Arb.

## Identificação da forma (Bombieri §12)

Para G=g*tilde(g), a fórmula explícita expressa o lado da forma quadrática como
T[G] = P[g]-S[g]+(1/2pi) integral_R m(r)|hat(g)(r)|²dr,
m(r)=Re psi(1/4+ir/2)-log(pi).

Usando Parseval e monotonicidade de m para r>=0:
T[G] = P[g]-S[g]+m(30)||g||²-A[g]+R_high[g], R_high[g]>=0.

Para o modo phi_(c,h,p)(x)=sqrt((2p+1)/(2h)) P_p((x-c)/h), sua transformada é sqrt(2h(2p+1)) i^p e^(irc) j_p(hr). Os sinais de seno/cosseno e o fator 2/pi seguem dessa identidade. O termo dos polos é P_ij=U_i V_j+V_i U_j; para a=0,40 o único primo é 2 e acopla os modos dos blocos esquerdo/direito com mesmo grau, por S_ij=log(2)/sqrt(2). Em a=0,30 e 0,34 a contribuição dos primos é nula.

O prolongamento da forma para modos por partes usa variação limitada, decaimento Fourier O(1/r) e aproximação por mollificação na norma logarítmica.

## Evidência e limites

- ES8-P81: 1.641 intervalos executados, 1.328 contidos nas caixas históricas e 313 zeros de paridade, com auditoria adicional sem FLINT.
- ES8-P82: 15 testes de estimativas racionais e SHA256, 208 comparações exatas de coeficientes Legendre–Bessel (graus 0...25).
- Pendente: auditoria de todas as chamadas de funções especiais FLINT e das cotas do operador complementar/Schur. Isto **não prova** a Hipótese de Riemann.

Fontes: https://dlmf.nist.gov/10.54 ; https://flintlib.org/doc/algorithms/ball-arithmetic/integration.html ; https://www.bdim.eu/item?id=RLIN_2000_9_11_3_183_0 .
