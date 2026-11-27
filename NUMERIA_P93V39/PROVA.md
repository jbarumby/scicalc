# NUMERIA–RIEMANN P93-V.39 — localização adaptativa, harmônica e momentos de borda

Escopo: continuidade da P93-V.38; verificamos SHA-256 do ZIP recebido e reproduzimos as novas somas de Schur por inteiros Python. Os recintos espectrais e a fórmula explícita são premissas herdadas da V.38; a varredura IEEE não foi regenerada por MPFR/Arb. Os prefixos são verificados por mpmath.iv sob as premissas de quadratura da V.38. Não é prova global de RH, nem cruzamento efetivo de 16. Nenhuma alteração no GitHub ou Lean.

## 1. Localização adaptativa

Para beta=2, T=4500, a<=1.387, sejam K_i as compressões positivas dos componentes de déficit e suas reflexões, P=P_600, Q=I-P. Para k_i>1, Cauchy–Schwarz no seminorma de K_i e Young dão

    K <= B + Q L Q,
    B = sum_i k_i P K_i P,
    L = sum_i [k_i/(k_i-1)] K_i.

O código hybrid39.py propõe multiplicadores por massa m_i = integral do déficit no componente positivo, e pesos de Schur por iteração numérica. O candidato escolhido dentre a grade testada usa k_i=2 para m_i<1, k_i=14 para m_i>=1. São 1110 e 120 componentes respectivamente. A função objetivo exploratória é sum k_i m_i/l_i^2. Sua razão em relação ao uniforme k=24 é 0.38238944; não é uma comparação de normas exatas de matrizes. Os limites de L são

    0 <= L <= C I, C=34164906632832863/18157142016000000,
    L <= D H, D=27093941/6500000, H=I-G_T.

O primeiro decorre do teste de Schur com kernel majorado na grade 2^-20, pi>333/106, massas e pesos inteiros. O segundo decorre da cota pontual máxima do símbolo amplificado. replay39.py reconstrói os componentes da tabela, verifica a regra dos fatores e repete as 1230 somas de linha em inteiros Python, 1.512.900 termos. Os resultados normal e -O coincidem. A reserva localizada 2-C=0.1183764161, versus 0.0088255179 do uniforme k=14 da V.38.

## 2. Harmônica híbrida e interpolação

No espaço Q, seja H_Q=QHQ, 0<=H_Q<=I_Q. Então QLQ<=CI_Q e QLQ<=D H_Q. Regularize com epsilon I e aplique a ordem inversa a ambas as cotas; tome combinação convexa das inversas e inverta. No limite epsilon->0 obtemos

    QLQ <= F(H_Q), F(x)=(C/lambda)s x/(1+s x),
    s=lambda D/[(1-lambda)C].

Para lambda=91/100, a reserva 2-F(1)=0.02065029399. Como D>C, F(1)>C; a harmônica não melhora a reserva da cauda para a mesma localização. Em compensação, F(x)<C precisamente para x<C/D, o que pode aliviar os modos intermediários. Uma família válida, theta em [0,1], é

    K <= B+Q[(1-theta) C I_Q+theta F(H_Q)]Q.

Não se usa o mínimo pontual entre operadores não comutantes. Para N>=600, a compressão de F(H_Q) é majorada por F(H_{Q,N}), pela desigualdade de compressão da inversa positiva (Schur). Na base ortonormal, escreva J=diag(I_P,(I+sH_{Q,N})/(1+s)). A cabeça suficiente é A0_N-B_N-[(1-theta)C Q_N^loc+theta F(H_{Q,N})], onde Q_N^loc é a projeção sobre o bloco posterior aos 600 modos dentro da cabeça (não confundir com o complemento da cabeça no espaço infinito). A congruência exata é

    J[A0_N-B_N-((1-theta)C+theta C/lambda)Q_N^loc]J
      +theta C/[lambda(1+s)] J Q_N^loc.

Como 0<J<=I, se a congruência >=hI, a cabeça suficiente >=hI. congruence39.py implementa a expressão e testa exatamente a identidade para theta=0,1/2,1 em exemplo racional 3x3. Não é uma montagem física de cabeça de 600+ modos. Para N<=600 a compressão do operador híbrido coincide com a do adaptativo localizado, independentemente de theta.

## 3. Prefixos e majorante de momentos com sinais

A montagem da forma de nove termos usa a autocorrelação racional V.38, GL40 com oito painéis, nós cercados por sinais do polinômio de Legendre, 110 dígitos intervalares. prefix39.py substitui a expansão simbólica por uma fórmula exata de coeficientes: P_n(2z-1)=sum_j (-1)^(n-j) binom(n,j)binom(n+j,j) z^j; integra-se P_n(2(w+u)-1)P_m(2w-1) de w=0 a 1-u e multiplica-se por 2a. Os testes de C(0) são mantidos. A versão antiga de perda usa 16^-j em Young e kdH; a prova suficiente ficou inconclusiva no quarto pivô para k=24. A forma completa de nove termos passou nos oito modos por setor.

Para polinômio f de paridade fixa, integração por partes termina exatamente. Agrupando j pares e ímpares nas contribuições de borda, os fatores trigonométricos são seno e cosseno, e os coeficientes relativos em cada grupo são (-1)^floor(j/2). Até fases comuns de módulo 1 (irrelevantes para a cota), Cauchy–Schwarz em duas coordenadas dá

    |fhat(t)|² <= 4 sum_{r=0,1} |sum_{j congruente r mod 2} (-1)^floor(j/2) f^(j)(a) t^(-j-1)|².

Assim, ao contrário de uma soma de majorantes absolutos das derivadas, são mantidos todos os produtos cruzados com sinais dentro de cada grupo. Para substituir a forma de nove termos pela cabeça localizada/híbrida no prefixo N<=600, basta majorar a perda de alta frequência e a penalidade B. Se W(t)>=Psi9(t)-beta e W(t)>=0 para t>=T, a perda é majorada por uma matriz Gram com momentos

    J_e=(4/3)[integral_T^infty W(t)t^-e dt + sum_cells k_i d_cell integral_l^r t^-e dt], e=j+k+2, e par.

Usamos 1/pi<1/3. O fator de Fourier 4 é da desigualdade anterior. Os termos em j,k com mesma paridade são combinados exatamente por frações; os produtos de derivadas são racionais. A integral de cada célula, l=L/16, r=R/16, d_cell em milionésimos, tem fórmula racional (para e>1)

    k_i d_cell 16^(e-1) [L^(-(e-1))-R^(-(e-1))]/(e-1).

Cada contribuição é encerrada na grade 2^-512 por divisão inteira dirigida; os limites são somados antes de convertê-los a intervalos. A cota W(t)=2 log t+25 é a da V.38. moments39.py usa seus momentos analíticos. Não se comete o erro de majorar entradas separadamente com sinal: encerra-se o momento escalar inteiro e calcula-se a matriz Gram com os produtos assinados. A desigualdade é de formas, válida para coeficientes reais ou complexos. A matriz suficiente é a forma completa menos essa perda; LDL intervalar deslocado por h vezes a massa prova a margem do prefixo, não do espaço inteiro.

O teste uniforme k=24, na largura 6931/5000, passou em oito modos por setor com h=10^-60. O teste adaptativo na largura racional a16 da V.38 passou em oito modos por setor e, no teste N=12 com a cota antiga de W, passou em 12 modos pares e 11 ímpares. A falha no 12º pivô ímpar não prova negatividade do minorante original.

## 4. Cota superior de Binet

Para Re z>0, a derivada da segunda fórmula de Binet fornece

    psi(z)=log z-1/(2z)-2 integral_0^infty u/[(u²+z²)(exp(2pi u)-1)] du.

Para z=1/4+it/2, t>0, |Im z²|=t/4, portanto |u²+z²|>=t/4. Como integral_0^infty u/(exp(2pi u)-1) du=1/24, o módulo do resto é <=1/(3t). Além disso Re[-1/(2z)]<=0 e Re log z<=log(t/2)+1/(8t²). Logo

    Psi9(t)-2 <= log(t/(2pi))+1/(8t²)+1/(3t)+massa-2
                 < log t+8  (t>=4500),

pela massa <9.944 e pi>1. O último orçamento racional é verificado no replay39.py. moments39_binet.py troca W por log t+8; seus momentos de cauda são [ (log T+8)/(e-1)+1/(e-1)^2 ] / T^(e-1). A nova cota fechou 12 pivôs positivos em ambas as paridades na largura racional a16, com deslocamento h=10^-60. O teste anterior havia fechado somente 11 pivôs ímpares. Os prefixos de 5 modos pares e 6 ímpares (ver JSON para as margens maiores) não substituem a certificação de cabeças completas. O argumento não usa a hipótese de Riemann. O teste continua dependente dos recintos de déficit herdados e da fórmula explícita na convenção histórica.

## 5. Schur condicional

Para o complemento da cabeça, a reserva assintótica é 2-[(1-theta)C+theta F(1)]. Pela identidade resolvente, o acoplamento da parte harmônica é <=theta D/(1-lambda)*sqrt(g), onde g é a fuga da banda da V.38. A parte B não cruza a cabeça quando N>=600. Assim b²<=2(50+theta D/(1-lambda))²g+2epsilon0², e o desvio de A0 é <=50g+epsilon0. Testamos a condição de Schur de V.38 com h=10^-60 e m=9h/10, por frações exatas. O pior caso de g é o setor par; o ímpar tem o fator adicional rho^-2. As cotas fatoriais para polos são verificadas. A busca em rho de 1.180 a 1.320, passo .005, dá rho=33/25, N=3309,3310,3312 para theta=0,1/2,1 respectivamente. Não se alega ótimo global. Os números são dimensões suficientes condicionadas à cabeça certificada na largura a16. A transferência de V.38 para 16 conserva o orçamento de perda (o novo minorante continua abaixo do operador auxiliar e a reserva alta é a mesma); não é cruzamento efetivo.

## 6. Reprodução e pendências

Na pasta do pacote, Python 3, NumPy, SymPy, mpmath; BLAS com um thread. A largura racional a16 é 69314718055994530941723212145817656807550013436025525412068000949339362197/50000000000000000000000000000000000000000000000000000000000000000000000000.

    sha256sum -c SHA256SUMS.txt
    python replay39.py
    python -O replay39.py
    OPENBLAS_NUM_THREADS=1 python hybrid39.py
    python congruence39.py
    python schur_adaptive39.py
    python prefix39.py --N 8 --k 24
    python moments39.py --N 8
    python moments39.py --N 12 --adaptive --width <a16>
    python moments39_binet.py --N 12 --adaptive --width <a16>

Os JSONs de matrizes são registros decimais dos intervalos, não certificados de extremos diádicos independentes. Para reprodução, recalcular os intervalos pelos scripts. Próximo gate: montar B em 600 modos e A0/G_T por bandas na largura racional a16, testar cabeças intermediárias e certificar a congruência ou a versão localizada. O operador híbrido ainda não foi fatorado acima de 600 modos. O MPFR/Arb independente, as duas cabeças completas, a extensão efetiva a 16 e a extensão cofinal continuam pendentes.
