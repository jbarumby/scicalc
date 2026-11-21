# P262 — Prova dos restos de integração e limites lógicos

**Fonte:** Konstantopoulos–Patie–Sarkar (2024), Annales de l'Institut Fourier, página impressa 397, seção 4.1, equações (4.15)–(4.18). https://www.numdam.org/item/10.5802/aif.3600.pdf . Também página 395, definição (4.8) de W_phi e J_Psi.

## 1. Transformação exata

F(k) = integral_1^infty (log x)^k x^(-3/4) soma_{m>=1} exp(-pi*m²*x) dx.

Fazendo x=exp(u), obtemos

F(k) = integral_0^infty u^k exp(u/4) soma_{m>=1} exp(-pi*m² exp(u)) du.

Cada termo de m fixo é inteiro em u. A integração em Arb calcula a parte finita m=1,2,3 e 0<=u<=5 com subdivisão em cinco segmentos unitários.

## 2. Resto da soma theta, m>=4

Para u>=0, pi>3 e exp(u)>=1+u. Para m>=4, a razão entre termos sucessivos é <=exp(-27) < 1/2. Portanto

soma_{m>=4} exp(-pi*m² exp(u)) <= 2 exp(-48(1+u)).

A integral do resto, para qualquer inteiro k>=0, é majorada por

E_m(k) = 2*k! / (2^48 * 47^(k+1)),

pois exp(-48)<2^-48, exp(u/4) exp(-48u)<=exp(-47u), e integral_0^infty u^k exp(-47u)du=k!/47^(k+1).

## 3. Resto u>=5

Se u=5+v>=5, então exp(5)>100, pi>3 e exp(v)>=1+v. A série theta inteira satisfaz

theta_0(exp(u)) <= 2 exp(-300) exp(-300v).

Também exp(5/4)<4 e exp(-300)<2^-300. Como 300-1/4>299,

E_u(k) = 8*2^-300 * soma_{j=0}^k C(k,j) * 5^(k-j)*j! / 299^(j+1)

majora toda a cauda u>=5.

A caixa Arb da parte finita é aumentada pelo erro total racional E_m+E_u, sem usar decisões de float.

## 4. Consistência da equação impressa

A Eq. (4.18), tal como publicada, diz

phi(n+1) = - G(2n)/(8(n+1)G(2n+2)),
G(2n)=64n(2n-1)F(2n-2)/F(2n)-1.

Se uma certificação intervalar provar G(2)>0 e G(4)>0, a expressão impressa phi(2) será negativa. Mas uma função de Bernstein, por definição, satisfaz phi(x)>=0 para x>=0. Essa conclusão é apenas um **gate de consistência da fonte impressa**.

A Eq. (4.8), página 395, define W_phi(1)=1 e W_phi(n+1)=produto_{k=1}^n phi(k). Se a Eq. (4.17) definir a_n=1/W_phi(n+1), então a recorrência implica, algebricamente,

phi(n+1) = a_n / a_(n+1).

Antes de formular qualquer alegação de impossibilidade da classe D_P, deve-se verificar o sinal, a indexação e a normalização de Xi/ Xi(0) no artigo e no coeficiente a_n. O status da Hipótese de Riemann **não é afetado** por um potencial problema de transcrição ou de convenção.

## 5. Obrigações abertas

- Reproduzir mais momentos F(2n) com rigor e aumento controlado de precisão.
- Auditar a Eq. (4.16) pelo coeficiente quadrático de Xi calculado diretamente da função completada, e discutir eventual normalização.
- Só testar monotonicidade/concavidade de sequências phi(n) **após** resolver a consistência das definições.
- Em nenhum caso extrapolar verificações finitas para todos os n ou para positividade da forma de Weil.
