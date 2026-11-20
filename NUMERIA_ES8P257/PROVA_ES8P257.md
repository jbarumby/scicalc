# NUMERIA–RIEMANN — ES8-P257: suavização, salto em norma e melhoria do orçamento de borda

**Referências de base:** P256 (`NUMERIA_ES8P256/PROVA_P256.md`, branch `numeria-p256-weighted-von-mangoldt-bridge`), P253 (`NUMERIA_ES8P253/README.md`, branch `numeria-p253-prime7-threshold`) e ponte Guinand–Weil da P251 Gate C. Esta é uma investigação de operadores na janela local `log 7 < 2L < log 8`; não é um resultado cofinal, nem prova RH.

## Teorema A — Uma suavização exata do símbolo

Considere `X=e^{2L}`, `L>0`, `d_q=log q`, `a_q=Λ(q)/√q` para as potências primas `q` e o símbolo duro da P249:

\[
S_L(t)=\sum_{q<X}2a_q\cos(td_q).
\]

Introduza o peso de *taper* `w_L(q)=(1-q/X)_+` e o símbolo regularizado

\[
\widetilde S_L(t)=\sum_{q<X}2a_q(1-q/X)\cos(td_q).
\]

Então, por integração termo a termo de uma soma finita,

\[
\boxed{\widetilde S_L(t)=2\int_0^L e^{-2(L-u)}S_u(t)\,du.}\tag{P257.1}
\]

De fato, um termo `q` surge em `S_u` somente para `u>a_q^*=\tfrac12\log q`, e `2\int_{a_q^*}^Le^{-2(L-u)}du=1-e^{-2(L-a_q^*)}=1-q/e^{2L}`. Os termos com `q=X` têm peso zero. Esta regularização é relacionada à média de von Mangoldt da P256, mas **não é o mesmo peso:** a P256 usa `Λ(q)(X/q−1)`, enquanto aqui o peso completo é `2Λ(q)/√q (1−q/X)`.

## Teorema B — Descontinuidade em norma no limiar do primo 7

Seja `d=log7`, `κ=d/√7`, `log7<2L<log8`, `δ=2L-d>0`, e `H_L=L²([-L,L])`. Prolongue `f∈H_L` por zero. Ponha

\[
(Jf)(x)=f(x-d)\quad (x∈[-L,L]), \qquad K_{7,L}=κ(J+J^*).
\]

A translação `J` é uma isometria parcial da faixa `E_-=[-L,L-d]` sobre `E_+=[-L+d,L]`. Ambas têm comprimento `δ`, são disjuntas porque `δ<d`, e `J²=0`. Na decomposição ortogonal `L²(E_-)\oplus L²(E_+)\oplus L²(E^c)`,

\[
J+J^*=\begin{pmatrix}0&I&0\\I&0&0\\0&0&0\end{pmatrix},
\qquad \boxed{\|K_{7,L}\|=κ},\tag{P257.2}
\]

para **todo** `δ>0`, mesmo arbitrariamente pequeno. Especificamente, para qualquer `u∈L²(E_-)` normalizado, `(u±Ju)/√2` são autovetores de `J+J*` com autovalores `±1`. A parcela `−K_{7,L}` tem direções negativas com quociente `−κ`, assim como direções positivas.

A regularização (P257.1) dá ao primo 7 o peso

\[
w_7=1-7/X=\frac{\varepsilon}{7+\varepsilon},\qquad X=7+\varepsilon,\quad 0<\varepsilon<1.
\]

A translação regularizada desse primo é `\widetilde K_{7,L}=κ w_7(J+J*)`. Logo

\[
\boxed{\|K_{7,L}-\widetilde K_{7,L}\|
=κ\frac7{7+\varepsilon}\longrightarrow\frac{\log7}{\sqrt7}>0}\tag{P257.3}
\]

quando `\varepsilon↓0`. **Não há convergência em norma** dessas parcelas, embora a contribuição suavizada tenda a zero. Esse limite deve ser entendido em cada espaço `H_L` ou após imersão em um espaço fixo maior; a norma permanece grande. No limiar exato `2L=d`, as faixas têm medida zero, e a translação comprimida desaparece. A passagem à direita é descontínua em norma.

Um contraexemplo de 2×2 demonstra por que não se pode transferir positividade por essa suavização: para `β=1/100`, `κ=3/4`, `w=1/100`, a matriz `βI-κw\begin{psmallmatrix}0&1\\1&0\end{psmallmatrix}` é definida positiva (menor autovalor `1/400`), enquanto `βI-κ\begin{psmallmatrix}0&1\\1&0\end{psmallmatrix}` é indefinida (menor autovalor `−74/100`). É um exemplo algébrico, **não** a matriz real da P253.

## Teorema C — Orçamento de borda aprimorado: 2225δ em lugar de 2520δ

**Hipóteses explícitas:** as identidades da fórmula Guinand–Weil da P251 Gate C, o certificado de sinal da P250 para o símbolo de **quatro** potências primas e o limite Binet `Ψ₄(t)−β≥B(t)−B(504)` para `t≥504`, onde `B(t)=log(t/(2π))−1/t` e `β=B(504)−Σ_{q∈{2,3,4,5}}2Λ(q)/√q`. Fixe `Ω=2016` e suponha `log7<2L<log8`, de modo que **7 é o único novo deslocamento de Mangoldt**.

Escreva `R_{4,L}` para o modelo auxiliar de quatro termos, com os polos originais e banda até 504 (mesma normalização da P253). A forma completa de Weil para `f∈C_c^∞([-L,L];ℂ)` satisfaz exatamente

\[
Q_W(f)=R_{4,L}(f)+\frac1{2π}\int_{|t|>504}
  (Ψ_4(t)-β)|\widehat f(t)|²dt
  -κ\langle f,(J+J^*)f\rangle .\tag{P257.4}
\]

Pela P250, o integrando do termo de cauda é não negativo para `|t|≥504`. Para `|t|≥Ω`, por Binet e monotonicidade de `B`, tem-se

\[
Ψ_4(t)-β\ge B(Ω)-B(504)=\log4+\frac1{672}>\frac{11}{8}.\tag{P257.5}
\]

Seja `P_Ω` a projeção de Fourier angular para frequências `|t|≤Ω` em `L²(ℝ)`, com transformada `\hat f(t)=\int f(x)e^{itx}dx` e Plancherel `\|f\|²=(2π)^{-1}\int|\hat f|²`. Ponha `P_>=I-P_Ω` e `E=E_-∪E_+`, com medida `|E|=2δ`. A norma de Hilbert–Schmidt do operador `1_E P_Ω`, de núcleo de Fourier conhecido, dá

\[
\|1_E P_Ω\|²\le\|1_E P_Ω\|_{HS}²=\frac{2δΩ}{π}.\tag{P257.6}
\]

Como `|⟨f,(J+J*)f⟩|≤\|1_E f\|²`, para todo `η>0`, a desigualdade de Young e (P257.6) dão

\[
\begin{aligned}
κ|\langle f,(J+J^*)f\rangle|
&\le κ\|1_E f\|²\\
&\le κ(1+η)\frac{2δΩ}{π}\|f\|²
+κ(1+η^{-1})\|P_>f\|².
\end{aligned}\tag{P257.7}
\]

**Cotas racionais independentes de ponto flutuante.** No `p257_strip_taper.py`, provamos com as séries geométricas de `arctanh`, a identidade de Machin para `π` e um intervalo racional para `√7`:

\[
\frac{73}{100}<κ<\frac{37}{50},\quad π>\frac{157}{50},\quad
 B(Ω)-B(504)>\frac{11}{8}.
\]

A prova usa `log7=2 artanh(3/4)`, `log4=2 artanh(3/5)`, `π=16 arctan(1/5)−4 arctan(1/239)` e `66/25<√7<53/20`, com restos racionais de séries. Os testes só conferem as contas finitas; a identidade de Machin, a representação integral da digama e Plancherel são as entradas matemáticas clássicas declaradas.

Reserve `1/12` da cauda alta. Seja `τ=11/8−1/12=31/24` e tome

\[
\eta=\frac{(37/50)}{\tau-37/50}=\frac{444}{331},\qquad
\frac{37}{50}(1+\eta^{-1})=\tau.
\]

Por (P257.7), o coeficiente de `δ\|f\|²` é estritamente menor que

\[
\frac{2Ω}{157/50}\frac{37}{50}(1+\eta)
=\frac{115617600}{51967}\approx2224{,}827295<2225.\tag{P257.8}
\]

Combinando (P257.4)–(P257.8), **sem pressupor positividade de `R_{4,L}`**, resulta

\[
\boxed{Q_W(f)\ge R_{4,L}(f)-2225\,δ\|f\|²+
\frac1{12}\|P_{|t|>2016}f\|².}\tag{P257.9}
\]

**Corolário estritamente condicional:** se **outro certificado independente** mostrar `R_{4,L}(f)≥m\|f\|²` em alguma largura auxiliar com `2L>log7`, então se poderia atravessar uma faixa com `δ≤m/4450`, garantindo pelo menos `m/2` (respeitando também a disponibilidade da largura auxiliar). A cota P253 usava `2520δ`, ou `δ≤m/5040`. O novo orçamento reduz a perda em `295/2520≈11,7%` e aumenta a tolerância admissível em cerca de 13,3%. **A P253 encontrou direções negativas certificadas para modelos auxiliares testados**. A P257 não cria a hipótese `m>0` que ainda falta: não há nova janela positiva demonstrada.

### O que não se segue

1. Uma soma suavizada não converge em norma para a parcela original do primo 7 no limiar — Teorema B é uma obstrução exata.
2. Positividade de operadores suavizados, de matrizes de células ou de uma banda de Fourier não implica positividade de Weil sem a comparação completa da cauda e dos polos.
3. A cota `2225δ` é um **aperfeiçoamento condicional de forma**, não um certificado de coerção para uma nova largura `L` e muito menos uma demonstração de RH.
4. Não se invoca nem se reutiliza RH presumida, tampouco a formalização de Alpöge, como premissa para essas desigualdades.

## Próxima decisão de pesquisa

Em vez de tentar aproximar o deslocamento original por um taper em norma, procurar (a) uma **reserva espectral efetiva** com bandas múltiplas, adaptada ao vetor com energia nas tiras de borda, ou (b) construir **um novo auxiliar com piso de cauda variável** cuja coercividade possa ser certificada acima de `log7/2`, levando a correção do primo 7 explicitamente. Em paralelo, é preciso revisar a identificação global da forma de Weil com os operadores antes de usar qualquer melhoria da janela como argumento cofinal.
