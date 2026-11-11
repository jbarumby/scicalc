# NUMERIA–RIEMANN — ES8-P151
## Aproximação admissível explícita por uma única escala, translações bilaterais e Gram de Toeplitz

**Status:** prova matemática construtiva para a **classe-candidata** de funções reais suaves compactamente suportadas com dois momentos nulos. A passagem para o critério clássico de Weil exige confirmação da classe e dos sinais da forma. Nenhuma prova de positividade universal ou de RH é apresentada.

## 1. Família regularizada

Seja `φ ≥ 0` uma função par, `C_c^∞(R)`, suportada em `[-1,1]` e com integral 1. A função bump normalizada da P146 serve de exemplo. Defina, para `d>0`,

- `φ_d(x) = d^(-1) φ(x/d)`;
- `Λ_d(x) = max(1-|x|/d,0)` (chapéu triangular);
- `B_d = Λ_d * φ_d`, uma função par, não negativa, suave, com suporte em `[-2d,2d]`.

Vale `∫B_d = d`, `||B_d||∞ ≤ 1`, `||B_d'||∞ ≤ 1/d`, `||B_d'||_1 ≤ 2` e `B_d > 0` no interior do suporte.

O fator exponencial do chapéu possui expressão exata: `∫Λ_d(x)e^(x/2)dx = 8(cosh(d/2)-1)/d`. Pela convolução, `m_d=[8(cosh(d/2)-1)/d]·∫φ_d(x)e^(x/2)dx`. Para a transformada de Fourier com convenção `e^(-i ξx)`, vale `B̂_d(ξ)=d·sinc²(dξ/2)·φ̂(dξ)`. Essas identidades podem ajudar na próxima recertificação de autocorrelações.

Fixe `d_m=2^(-m)` e `g_(m,j)(x)=B_(d_m)(x-j d_m)`. Escolha

`S_m = span{g_(m,j) : j∈Z, |j d_m|≤m}`.

Este espaço possui `K_m=2m·2^m+1` translações distintas, incluindo **centros positivos e negativos**, e permite `d_m→0`.

## 2. Aproximação C¹ com taxa explícita

Para `f∈C_c²(R)`, seja `I_d f` a interpolação linear global nos pontos `jd`:

`I_d f(x)=∑_j f(jd) Λ_d(x-jd)`.

Ponha `P_d f= (I_d f)*φ_d = ∑_j f(jd) B_d(x-jd)`. A soma é finita se `f` tem suporte compacto. Seja `M1=||f'||∞`, `M2=||f''||∞` e defina `||f||C1=||f||∞+||f'||∞`.

Em cada intervalo da malha, a fórmula do erro de interpolação e o teorema do valor médio dão

`||I_d f-f||∞ ≤ (d²/8)M2`, `||(I_d f)'-f'||∞ ≤ (d/2)M2`

(a segunda desigualdade vale para a derivada definida quase em toda parte). Pela paridade, normalização e suporte de `φ_d`, `∫ |x|φ_d(x)dx ≤ d`. Assim,

`||f*φ_d-f||∞ ≤ d M1`, `||f'*φ_d-f'||∞ ≤ d M2`.

A convolução com `φ_d` não aumenta a norma suprema. Portanto:

**(T1)** `||P_d f-f||C1 ≤ E_d(f) := d M1 + (3d/2+d²/8) M2`.

Se `supp(f)⊆[-L,L]`, então `supp(P_d f)⊆[-L-2d,L+2d]`. Para `m≥max(3,ceil(L))`, o aproximante `P_(d_m) f` pertence a `S_m` e os suportes permanecem numa mesma janela.

## 3. Dois momentos, sem recorrência a corretores externos

Defina `ℓ_±(f)=∫_R f(x)e^(±x/2)dx`. Para o mesmo `d=d_m`,

`ℓ_±(g_(m,j))=m_d e^(±jd/2)`, onde `m_d=∫B_d(x)e^(x/2)dx>0`.

Fixe `d0=1/2`, `u0(x)=B_(d0)(x-1)`, `v0(x)=B_(d0)(x+1)`. Para cada `m≥3`, seja `U_m=P_(d_m)u0` e `V_m=P_(d_m)v0`. Como `1/d_m` é inteiro, ambos pertencem **ao próprio S_m**. Por translação da grade, existe uma função `H_m=P_(d_m) B_(d0)`, par e não negativa, tal que `U_m(x)=H_m(x-1)` e `V_m(x)=H_m(x+1)`.

Consequentemente, definindo `a_m=∫H_m(x)e^(x/2)dx`, a matriz exata de momentos é

`M_m=a_m [[e^(1/2),e^(-1/2)],[e^(-1/2),e^(1/2)]]`.

**Uniformidade explícita:** `a_m ≥∫H_m =d_m Σ_j B_(d0)(j d_m)`. Pelo erro da regra retangular para função de variação limitada,

`|d_m Σ_j B_(d0)(j d_m) -∫ B_(d0)| ≤ d_m ||B_(d0)'||_1 ≤2d_m`.

Como `∫B_(d0)=1/2` e `d_m≤1/8`, concluímos `a_m≥1/4`. Logo `det(M_m)>0` e

`||M_m^-1||∞=1/(2a_m sinh(1/2)) ≤4`.

Pela interpolação convexa e pela convolução normalizada, `||U_m||C1, ||V_m||C1≤3`.

Se `f` já tem os dois momentos nulos, defina o aproximante admissível

`P_m^0 f=P_(d_m) f-[U_m,V_m] M_m^-1 (ℓ_+(P_(d_m)f),ℓ_-(P_(d_m)f))^T`.

Então `P_m^0 f∈S_m` e **exatamente** `ℓ_+(P_m^0 f)=ℓ_-(P_m^0 f)=0`.

Escolha o inteiro `A=max(ceil(L)+1,3)`. Os suportes de `f`, `P_(d_m)f`, `U_m`, `V_m` ficam em `[-A,A]`. Como `|ℓ_±(w)|≤2A e^(A/2)||w||∞`, e `e^(A/2)≤2^A`, a norma da correção é no máximo `48 A·2^A E_(d_m)(f)`.

Portanto:

**(T2)**

`||P_m^0 f-f||C1 ≤ (1+48A·2^A) [2^(-m)M1 + (3·2^(-m)/2+2^(-2m)/8)M2]`.

O lado direito tende a **zero explicitamente**; o suporte é uniformemente controlado em `[-A,A]`. A aproximação é obtida dentro do **núcleo dos dois momentos do espaço de uma única escala**.

## 4. Base local explícita para esse núcleo

Escreva `q_m=e^(d_m/2)>1`. Para cada índice interior `j`, defina o triplo local:

`A_(m,j)=g_(m,j-1) - (q_m+q_m^-1)g_(m,j) + g_(m,j+1)`.

**(T3)** `ℓ_±(A_(m,j))=0` exatamente por cancelamento de polinômios de Laurent.

As `K_m` translações `g_(m,j)` são independentes: a translação com centro mais à esquerda é não nula em um intervalo onde todas as demais ainda são nulas. A aplicação dos dois momentos possui posto 2 (os valores `q_m^j` e `q_m^-j` não são proporcionais). Os `K_m-2` triplos locais são independentes, pois a combinação não nula mais à esquerda apresenta coeficiente exclusivo em seu índice esquerdo. Consequentemente, os triplos locais formam uma base do núcleo, de dimensão

`dim(S_m ∩ ker ℓ_+ ∩ ker ℓ_-) = K_m-2 = 2m·2^m - 1`.

Pelo Teorema T2, essa sequência de núcleos de escala única é suficiente para aproximar qualquer função `C_c²` da classe-candidata com ambos os momentos nulos em norma `C¹` com suporte fixo. Não é preciso que a sequência `S_m` seja aninhada: basta a existência de uma sequência de aproximantes.

## 5. Redução estrutural a matrizes Toeplitz

**Sob a hipótese explícita** de que a forma bilinear real `B_Weil(f,g)` da normalização adotada é simétrica e invariante por translação simultânea, a matriz na família original `g_(m,j)` tem a forma Toeplitz `W_(i,j)=w_m(i-j)`, com `w_m(-k)=w_m(k)`.

Se `c_m=q_m+q_m^-1=2cosh(d_m/2)`, a matriz no núcleo, na base de triplos, também é Toeplitz, com lag

**(T4)** `h_m(k) = w_m(k-2)+w_m(k+2)-2c_m[w_m(k-1)+w_m(k+1)]+(c_m²+2)w_m(k)`.

Para `K_m` centros, a Gram bruta tem `K_m(K_m+1)/2` elementos distintos quando armazenada genericamente, mas, **se a invariância por translação e a simetria forem verificadas para a normalização de Weil**, são necessárias apenas `K_m` avaliações de lags `w_m(0),...,w_m(K_m-1)`. O verificador `p151_toeplitz_exact.py` confere a identidade de congruência em frações racionais, usando valores sintéticos de `w_k` e `q`.

Essa compressão é **algébrica**. Não determina os lags da verdadeira forma de Weil, não implica que a forma seja positiva, e não substitui certificados intervalares das integrais arquimedianas.

## 6. Limite lógico e próximo alvo

A hipótese ainda aberta é demonstrar que as matrizes Toeplitz restritas correspondentes à **forma clássica de Weil corretamente identificada** sejam semidefinidas positivas **para todo nível `m`** (ou provar diretamente positividade em toda a classe candidata). Uma verificação finita, mesmo de todos os elementos de uma matriz grande, não implica esse enunciado universal.

Os resultados aqui são os lemas de aproximação, anulamento exato dos momentos e compressão condicional de Gram. **Não constituem demonstração de RH**.
