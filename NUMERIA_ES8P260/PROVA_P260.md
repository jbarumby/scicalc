# NUMERIA–RIEMANN P260 — Demonstração da auditoria arquimediana por frações exatas

## 0. Enquadramento bibliográfico

Ahlfors, *Complex Analysis*, terceira edição, capítulos 4 (resíduos e princípio do argumento) e 5 (§2.2 produtos infinitos, §2.4 Gama, §2.5 Stirling, §4 função zeta). O sumário público foi consultado. O livro integral **não** foi fornecido: a demonstração abaixo é uma derivação independente pela expansão clássica de Euler–Weierstrass, também documentada em [NIST/DLMF](https://dlmf.nist.gov/5.7). Não reproduzimos texto protegido do livro.

## 1. Monotonicidade rigorosa da parte real da digama

A função digama, definida por ψ(z)=Γ′(z)/Γ(z), satisfaz para Re(z)>0:

\[
\psi(z)=-\gamma+\sum_{n\ge0}\left(\frac1{n+1}-\frac1{n+z}\right).
\]

Com a=1/4, b=t/2, x_n=n+a, denotemos D(t)=Re ψ(1/4+it/2). Temos:

\[
D(t)=-\gamma+\sum_{n\ge0}\left(\frac1{n+1}-\frac{x_n}{x_n^2+(t/2)^2}\right).
\]

Para t₂>t₁≥0:

\[
D(t_2)-D(t_1)=\sum_{n\ge0}
\frac{x_n((t_2/2)^2-(t_1/2)^2)}
{(x_n^2+(t_1/2)^2)(x_n^2+(t_2/2)^2)}>0.
\]

Cada parcela é positiva. Para t₁=504 e t₂=T inteiro par, a parcela n é **exatamente**:

\[
\frac{16(4n+1)(T^2-504^2)}
{((4n+1)^2+4\cdot504^2)((4n+1)^2+4T^2)}.
\]

O auditor soma n=0,…,2047, arredondando cada parcela **para baixo** por divisão inteira em unidades de 10⁻³⁰. A cauda omitida é positiva. Além disso, D(t) cresce com |t|, de modo que a cota em T é uniforme para todos |t|≥T.

## 2. Cota independente na frequência 504

Em t=504, b=252. A parcela da série é

\[
\frac1{n+1}-\frac{n+1/4}{(n+1/4)^2+252^2}
=\frac{16\cdot252^2-3(4n+1)}
{(n+1)((4n+1)^2+16\cdot252^2)}.
\]

Somamos os primeiros N=8192 termos com arredondamento inferior por inteiro. Para a cauda:

\[
\frac1{n+1}-\frac{x_n}{x_n^2+252^2}
\ge\frac1{n+1}-\frac1{x_n}
=-\frac{3/4}{(n+1)(n+1/4)}
\ge-\frac{3/4}{(n+1/4)^2}.
\]

A comparação integral de uma série decrescente dá

\[
\sum_{n=N}^{\infty}\frac1{(n+1/4)^2}
\le\frac1{N+1/4}+\frac1{(N+1/4)^2}.
\]

O auditor também usa γ<H₄₀₉₆−log4096. O somatório harmônico é aproximado **para cima** por inteiros, log4096=12log2 é limitado **para baixo** pela série racional de log, e

\[
\log x=2\sum_{k=0}^{K-1}\frac{v^{2k+1}}{2k+1}+R_K,\quad
v=\frac{x-1}{x+1},\quad
0<R_K\le\frac{2v^{2K+1}}{(2K+1)(1-v^2)}
\]

é utilizado para limitar log2, log3, log5, log(355/113).

A identidade de Machin, π=16atan(1/5)−4atan(1/239), e os somatórios alternados dão 3,14159<π<355/113, de modo que logπ<log(355/113). Para q∈{2,3,4,5}, sqrt(q) é limitado inferiormente pelo piso exato de sqrt(q×10⁶⁰)/10³⁰; log p é limitado superiormente, portanto cada amplitude 2log p/sqrt(q) recebe majorante racional.

Conclui-se, independentemente do limite de Binet da P250:

\[
\boxed{
D(504)-\log\pi-\sum_{q\in\{2,3,4,5\}}\frac{2\Lambda(q)}{\sqrt q}
>0.002605374463140167.
}
\]

Cada sinal nesta desigualdade é uma comparação de inteiros / frações exatas.

## 3. Pisos uniformes do símbolo de quatro potências primas

Defina

\[
\Psi_4(t)=D(|t|)-\log\pi-
2\sum_{q\in\{2,3,4,5\}}\frac{\Lambda(q)}{\sqrt q}\cos(t\log q).
\]

Como cos≤1, para todo |t|≥T:

\[
\Psi_4(t)\ge D(504)-\log\pi-\sum_q 2\Lambda(q)/\sqrt q
+ D(T)-D(504).
\]

O software certifica as seguintes cotas **estritas**, expressas abaixo em grade decimal racional conservadora:

| Intervalo | Cota inferior de Ψ₄(t) | Piso escolhido |
|---|---:|---:|
| |t|≥850 | 0,511694434607381370 | 1/2 |
| |t|≥1100 | 0,755782079936074253 | **3/4** |
| |t|≥1500 | 1,037826255173832180 | **1** |
| |t|≥2000 | 1,281505957163487627 | **5/4** |

**Atenção:** a cota obtida por esta técnica em t=504 é apenas positiva; o piso 1/4 na faixa 504–850 continua apoiado no certificado anterior P253. A nova prova, isoladamente, justifica os quatro patamares **a partir de 850**.

## 4. Auditoria dos polos e da equação funcional: limites

Do capítulo de Gama, a recorrência Γ(z+1)=zΓ(z) implica ψ(z+1)=ψ(z)+1/z; logo os polos de ψ em z=0,−1,−2,… têm resíduo −1. A forma clássica da função zeta completada é

\[
\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s),
\quad\xi(s)=\xi(1-s).
\]

Esses fatos fundamentam o tratamento arquimediano. **Não determinam isoladamente os fatores +2 no setor par e −2 no setor ímpar da matriz P251**: aqueles dependem das convenções da fórmula explícita e da base Fourier–Legendre. A auditoria do sinal e da normalização da contribuição de polos inteira permanece uma obrigação independente. Não se atribui a Ahlfors nenhuma demonstração do certificado NUMERIA.

## 5. Consequência matricial condicional

A P259 usa pisos 1/4 a partir de 504 e 1/2 a partir de 850. O novo minorante pode ser reforçado por projeções espectrais comprimidas:

\[
\boxed{
R_{\mathrm{P260}}=R_{\mathrm{P259}}
+\tfrac14(I-G_{1100})
+\tfrac14(I-G_{1500})
+\tfrac14(I-G_{2000}).
}
\]

Por Plancherel, cada parcela é semidefinida positiva. **Não se provou que a nova matriz 377×377 seja SPD**, que a antiga direção negativa de 57 coordenadas tenha sido eliminada nem que a forma completa de Weil tenha positividade além da janela P252. Essas são tarefas do próximo Gate B e exigirão Arb/MPFR e replay de Rayleigh independente.

**Escopo rigoroso:** novos pisos para um símbolo auxiliar, nada sobre a validade global da Hipótese de Riemann.
