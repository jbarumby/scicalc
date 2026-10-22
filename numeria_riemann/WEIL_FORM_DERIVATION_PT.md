# NUMERIA–RIEMANN: derivação da forma quadrática de Weil implementada

**Data:** 2026-10-21 (revisão crítica proposta em 2026-10-22). **Escopo:** deduzir, a partir de um teorema citado (fórmula explícita de Guinand–Weil), cada termo da matriz que `weil_interval_audit.py` constrói, incluindo sinais, fatores de normalização e restrições, e mapear cada passo para o código e para um teste. Esta revisão separa explicitamente a identidade computada da forma hermitiana usada no critério global de Weil.

**O que "rigoroso" significa aqui.** Há três camadas, e este documento não as mistura:

1. **Teorema citado** (§1): a fórmula explícita não é reprovada aqui; ela vem da literatura. A leitura exata das hipóteses nas fontes precisa de revisão humana (ver §7).
2. **Deduções** (§2–§4): feitas por cálculo direto, com cada passo escrito. Podem ser conferidas à mão.
3. **Verificações numéricas** (§5): testes independentes (mpmath, zeros reais da zeta). São evidência, não prova.

Nada aqui demonstra a hipótese de Riemann nem a positividade global de Weil.

## 1. O teorema citado (Guinand–Weil)

Seja $h$ inteira (ou analítica numa faixa $|\operatorname{Im}r|\le\tfrac12+\varepsilon$) com decaimento suficiente, e $g(u)=\frac1{2\pi}\int_{\mathbb R}h(r)e^{-iur}dr$ par, contínua e de suporte compacto. Para cada zero não trivial $\rho$, defina $\gamma_\rho=-i(\rho-\tfrac12)$, de modo que $\rho=\tfrac12+i\gamma_\rho$. **Não se assume que $\gamma_\rho$ seja real.** A realidade de todos os $\gamma_\rho$ equivale à Hipótese de Riemann. A soma em (W) deve ser entendida com multiplicidades e sob as hipóteses de convergência da versão exata da fórmula explícita citada; a verificação dessas hipóteses para a função teste usada aqui permanece pendente.
$$\sum_\rho h(\gamma_\rho)=h\!\left(\tfrac i2\right)+h\!\left(-\tfrac i2\right)-g(0)\log\pi+\frac1{2\pi}\int_{\mathbb R}h(r)\operatorname{Re}\psi\!\left(\tfrac14+\tfrac{ir}2\right)dr-2\sum_{n\ge1}\frac{\Lambda(n)}{\sqrt n}\,g(\log n).\tag{W}$$

*Fontes consultadas* (a forma é a mesma, a menos de notação):

- Lema 4 de arXiv:1702.04099 ("Bounding $S_n(t)$ on the Riemann hypothesis") enuncia a fórmula de Guinand–Weil com $\hat h(x)=\int h(t)e^{-2\pi itx}dt$. Como $g(\log n)=\frac1{2\pi}\hat h\!\big(\tfrac{\log n}{2\pi}\big)$ e $g(0)=\frac1{2\pi}\hat h(0)$, os termos de primos e o termo $\log\pi$ coincidem com (W), e o termo de primos $\frac1{2\pi}\sum\frac{\Lambda(n)}{\sqrt n}\big[\hat h(\tfrac{\log n}{2\pi})+\hat h(-\tfrac{\log n}{2\pi})\big]$ é $2\sum\frac{\Lambda(n)}{\sqrt n}g(\log n)$ para $g$ par.
- arXiv:2608.13637, eq. (2.1), traz a mesma fórmula com $\mu(\tau)=\frac1{2\pi}\operatorname{Re}\psi(\tfrac14+\tfrac{i\tau}2)-\frac{\log\pi}{2\pi}$, isto é, o $\log\pi$ já absorvido no peso arquimediano, como em (W) depois do §3.1.
- arXiv:2609.04908 usa a variante $-F(0)\log2\pi-\mathcal A_\infty(F)$ com $\frac1{2\pi}\int H\operatorname{Re}\psi=-F(0)\log2-\mathcal A_\infty(F)$; substituindo, $-F(0)\log2\pi-\mathcal A_\infty=-F(0)\log\pi+\frac1{2\pi}\int H\operatorname{Re}\psi$, de novo (W).
- A prova da fórmula está em Poitou e em Lang (*Algebraic Number Theory*, cap. XVII), conforme citado em arXiv:0903.3105; a formulação original é de Weil (1952).
- O critério de Weil: a forma hermitiana $\sum_\rho\hat f(\gamma_\rho)\overline{\hat g(\bar\gamma_\rho)}$ ser positiva em todas as funções teste compactas é equivalente a RH (Weil 1952; Bombieri 2000, conforme citado em arXiv:2608.13637).

## 2. A função teste do protótipo

Base: $g_k(x)=\sin(w_kx)/\sqrt L$ em $[0,2L]$, $w_k=k\pi/(2L)$, e $h_k(t)=g_k(t+L)$ em $[-L,L]$ (zero fora). Para $c\in\mathbb R^N$ seja
$$f=\sum_kc_kh_k,\qquad F_c(r)=\int f(t)e^{irt}dt=\sum_kc_kF_k(r),\qquad h(r):=F_c(r)F_c(-r),\qquad g_\star(u):=\int f(t+u)f(t)\,dt .$$
(Uso $g_\star$ para a função teste de (W), para não confundir com a base $g_k$.)

- $f$ é real, contínua (anula-se em $t=\pm L$ porque $\sin(0)=\sin(k\pi)=0$) e Lipschitz, com suporte em $[-L,L]$. Logo $g_\star$ é par, de classe $C^1$ e tem suporte em $[-2L,2L]$.
- $F_c$ é inteira e, pelo Lema 1 de `TAIL_PROOF.md`, $F_c(r)=O(r^{-2})$. Como $f$ é real, $F_c(-r)=\overline{F_c(r)}$ para $r$ real, então $h(r)=|F_c(r)|^2\ge0$ é par, inteira, $O(r^{-4})$ em $\mathbb R$, e $g_\star=\frac1{2\pi}\int h\,e^{-iur}dr$ (a transformada de Fourier da autocorrelação é $|F_c|^2$).

Essas propriedades são as hipóteses usuais de (W) (suporte compacto, $g_\star$ contínua e $C^1$, $h$ analítica na faixa com decaimento). A conferência literal contra o enunciado exato de cada fonte fica para revisão humana.

## 3. Dedução termo a termo

### 3.1 Termo arquimediano e o $\log\pi$

Por Fourier, $g_\star(0)=\frac1{2\pi}\int_{\mathbb R}h\,dr$. Então
$$\frac1{2\pi}\int h\operatorname{Re}\psi-g_\star(0)\log\pi=\frac1{2\pi}\int_{\mathbb R}h(r)\,\Omega(r)\,dr,\qquad \Omega(r)=\operatorname{Re}\psi(\tfrac14+\tfrac{ir}2)-\log\pi .$$
Como $h$ e $\Omega$ são pares ($\psi(\bar z)=\overline{\psi(z)}$), isto vale $\frac1\pi\int_0^\infty h\,\Omega\,dr$. E como $h=|F_c|^2=\sum_{i,j}c_ic_j\operatorname{Re}(F_i\overline{F_j})$,
$$\frac1\pi\int_0^\infty h\,\Omega\,dr=\sum_{i,j}c_ic_jG_{ij},\qquad G_{ij}=\frac1\pi\int_0^\infty\operatorname{Re}\!\big(F_i\overline{F_j}\big)\,\Omega\,dr,$$
que é exatamente o que o código integra (`_archimedean_matrix_adaptive`), com a cauda de `TAIL_PROOF.md`. O Lema 1 dá $\operatorname{Re}(F_i\overline{F_j})=\cos(\tfrac{(i-j)\pi}2)A_iA_j$ (`_phase`).

*Consistência da normalização.* Como $\int h_ih_j\,dt=\delta_{ij}$ (Gram, §3.2 com $a=0$), Parseval dá $\frac1\pi\int_0^\infty A_k^2\,dr=\int h_k^2\,dt=1$. É a identidade que faz o $-g_\star(0)\log\pi$ virar exatamente o "$-\log\pi$" de $\Omega$.

### 3.2 Termo de primos

$g_\star(a)=\int f(t+a)f(t)dt$. Com $x=t+L$, para $a\ge0$ e $x\in[0,2L-a]$:
$$g_\star(a)=\sum_{i,j}c_ic_jC_{ij}(a),\qquad C_{ij}(a)=\frac1L\int_0^{2L-a}\sin\!\big(w_i(x+a)\big)\sin(w_jx)\,dx .$$
Por produto-soma, $\sin A\sin B=\tfrac12[\cos(A-B)-\cos(A+B)]$ com $A-B=(w_i-w_j)x+w_ia$ e $A+B=(w_i+w_j)x+w_ia$, então $C_{ij}(a)=\frac1{2L}(\text{first}-\text{second})$ com $\text{first}=\int_0^{2L-a}\cos((w_i-w_j)x+w_ia)dx$ e $\text{second}$ análogo. É exatamente `_overlap` (o fator $1/(2L)$ é o que o erro de 2026-10-20 tinha como $1/(4L)$).

Como $g_\star$ é par, o termo de primos de (W) é $-2\sum_n\frac{\Lambda(n)}{\sqrt n}g_\star(\log n)$, e como $\operatorname{supp}g_\star\subset[-2L,2L]$ só entram os $n$ com $\log n\le2L$ (para $L=0{,}45$: apenas $n=2$). Resulta
$$P_{ij}=-2\sum_{\log n\le2L}\frac{\Lambda(n)}{\sqrt n}\,C^{\rm sym}_{ij}(\log n),\qquad C^{\rm sym}_{ij}=\tfrac12(C_{ij}+C_{ji}),$$
com $\Lambda(n)=\log p$ se $n=p^m$ (`_prime_power_base`). A simetrização não altera a forma quadrática $\sum c_ic_jP_{ij}$. O corte $\log n\le2L$ é decidido em aritmética de bolas, e o código aborta se não consegue separar.

### 3.3 Termos de polo e restrições

Como $f$ é real, a continuação analítica de $h$ é $h(z)=F_c(z)F_c(-z)$, logo $h(\pm\tfrac i2)=F_c(\tfrac i2)F_c(-\tfrac i2)$. Se $F_c(\tfrac i2)=0$ e $F_c(-\tfrac i2)=0$, os dois termos de polo de (W) se anulam. Essas duas condições são **lineares** em $c$:
$$F_k(\mp\tfrac i2)=\int_0^{2L}g_k(x)\,e^{\pm(x-L)/2}dx=\frac{w_k\big(e^{\mp L/2}-(-1)^ke^{\pm L/2}\big)}{\sqrt L\,(w_k^2+\tfrac14)}$$
(integração direta; `plus` e `minus` em `_constraint_nullspace`). O subespaço $\{c:\;F_c(\pm\tfrac i2)=0\}$ é parametrizado resolvendo as duas primeiras coordenadas em função das $N-2$ restantes (sistema $2\times2$ cujo determinante é verificado como separado de zero em bolas), o que define $B$ com $C\,B=0$. O sinal da exponencial é irrelevante porque as **duas** condições são impostas.

Observação de escopo: exigir os dois fatores nulos é **suficiente** para anular o termo de polo, mas não necessário (o produto poderia anular-se com um só fator nulo, o que descreve uma união de hiperplanos, não um subespaço). O protótipo trabalha no subespaço, e é nele que a positividade é testada.

### 3.4 Forma restrita

Para $c=Bc'$, a identidade formal de (W), se aplicável à função teste e com a convergência justificada, dá
$\sum_\rho F_c(\gamma_\rho)F_c(-\gamma_\rho)=c'^{\mathsf T}Ac',\qquad A=B^{\mathsf T}(G+P)B,$
sem termos de polo. É a matriz projetada (`_project`) cuja positividade definida o código certifica.

**Distinção essencial:** se $\gamma_\rho\in\mathbb R$, então $F_c(-\gamma_\rho)=\overline{F_c(\gamma_\rho)}$ e o termo é $|F_c(\gamma_\rho)|^2\ge0$. Se $\gamma_\rho$ é complexo, essa igualdade em geral é falsa; portanto a matriz positiva calculada não pode ser interpretada, sem argumento adicional, como positividade da soma de módulos quadrados sobre todos os zeros. A forma hermitiana geral do critério de Weil envolve a conjugação e a simetria corretas da função teste e dos zeros; a identificação entre essa forma e a forma quadrática construída aqui ainda precisa ser demonstrada. Não usar a equação acima como prova do critério global.

## 4. Sinais e fatores: conferência item a item

| Item | Valor | Origem | Onde é conferido |
|---|---|---|---|
| Peso arquimediano | $\frac1{2\pi}\int_{\mathbb R}\to\frac1\pi\int_0^\infty$ | paridade de $h$ e $\Omega$ | teste ponta a ponta (§5) |
| Termo $\log\pi$ | $-g_\star(0)\log\pi\to-\log\pi$ em $\Omega$ | Fourier: $g_\star(0)=\frac1{2\pi}\int h$ | `test_normalization_parseval…` |
| Termo de primos | $-2\Lambda(n)n^{-1/2}g_\star(\log n)$ (o 2 vem de $\pm\log n$) | (W) com $g_\star$ par | ponta a ponta; sensibilidade em §5 |
| Overlap | $\frac1{2L}(\text{first}-\text{second})$ | produto-soma | `test_zero_shift_overlap_is_identity…`, `test_prime_overlap_matches_direct_quadrature` |
| Fase de $F_i\bar F_j$ | $\cos(\tfrac{(i-j)\pi}2)$ | Lema 1 | `test_tail_proof_lemma1…` |
| Restrições | $F_c(\pm\tfrac i2)=0$ | §3.3 | `test_pole_constraint_vectors_match_direct_integrals` |
| Corte de primos | $\log n\le2L$ | $\operatorname{supp}g_\star\subset[-2L,2L]$ | `_prime_matrix` (separação em bolas) |

## 5. Verificação numérica ponta a ponta (evidência, não prova)

`test_quadratic_form_equals_sum_over_real_zeta_zeros` (`NUMERIA_SLOW=1`) compara, para $c'$ aleatório, $c'^{\mathsf T}Ac'$ mais a cauda arquimediana além de $R$ com $2\sum_{\gamma>0}|F_c(\gamma)|^2$ sobre zeros conhecidos na linha crítica e uma cauda de zeros estimada pela densidade de Riemann–von Mangoldt. Esta é uma **checagem de consistência condicionada aos zeros reais usados**, não uma verificação independente da fórmula para zeros fora da linha crítica, nem uma prova de que o resto da soma de zeros esteja limitado. Os números reportados devem ser vinculados ao log/artifact exato da execução antes de serem tratados como reprodução atual.

*Sensibilidade.* Nesses mesmos vetores o termo de primos vale $-0{,}097$ e $-0{,}069$. Um sinal trocado ou um fator 2 a mais ou a menos deslocaria a comparação em $0{,}07$ a $0{,}19$, cerca de sete ordens de grandeza acima da diferença observada. Por isso o teste distingue os erros de fator e de sinal que já ocorreram nesta linha de trabalho.

Limites: o teste usa os zeros reais (portanto assume RH nos zeros usados), é numérico, e o resto da soma de zeros é estimado, não limitado.

## 6. O que os certificados dizem e não dizem

- O certificado de positividade de $A$ é um fato **sobre esta matriz**, calculado sem usar RH nem a localização de zeros. A interpretação pela fórmula explícita depende das hipóteses e convergência de (W), ainda pendentes de auditoria independente.
- Se RH for assumida, os parâmetros $\gamma_\rho$ são reais e a soma correspondente se torna uma soma de módulos quadrados para esta escolha de $h$. Sem RH, não se pode usar essa identidade de módulos quadrados. Em qualquer caso, esta família é de dimensão finita e tem suporte pequeno ($L=0{,}45$, só $n=2$ entra), enquanto o critério global exige a classe completa de funções teste apropriadas; a positividade desta matriz não implica RH.
- Os relatórios exportam extremos decimais arredondados **para fora** (inferior por baixo, superior por cima) e ponto médio e raio exatos. `verify_certificate.py` refaz o teste final em racionais exatos a partir do JSON, sem Arb e sem ponto flutuante.

## 7. O que continua aberto

1. Conferir literalmente as hipóteses de (W) contra o enunciado primário de Poitou/Lang (ou de Weil) para esta classe de $h$, incluindo a variável complexa $\gamma_\rho$, multiplicidades, simetrias, convergência e eventuais termos de fronteira; aqui só se verificou que $g_\star$ é contínua, $C^1$, par, de suporte compacto e que $h$ é inteira com $O(r^{-4})$ na reta real.
2. Revisão independente deste documento e de `TAIL_PROOF.md`: ambos são novos.
3. A correção da avaliação do Arb é assumida; `verify_certificate.py` cobre só o passo final.
4. Nada aqui diz algo sobre outros $L$, outras bases ou o critério global.

5. Derivar, sem assumir RH, a correspondência exata entre a forma quadrática real $F_c(r)F_c(-r)$ implementada e a forma hermitiana de Weil para a classe de teste escolhida. Até isso ser feito, os resultados são certificados de matriz finita e verificações numéricas restritas, não uma verificação do critério global.
