# NUMERIA P83 — demonstração elementar independente da cota curta 0,605

**Resultado matemático** (para qualquer **X real >=10^9**):

\[
\boxed{\quad
\sum_{X<n<5X}\frac{\lambda_X(n)^2}{n^2}
\le \frac{121}{200X}=\frac{0.605}{X},\qquad
\lambda_X(n)=\sum_{\substack{d\mid n\\d\le X}}\mu(d)\quad(n>1).
\quad}
\]

Trata-se da mesma cota curta **b3** do Lema 5.6 de Ramaré, utilizada em KLN (2018), Lema 3.4 (3.8), e Chourasiya–Simonič v2, Lema 6. **A prova abaixo não usa o Lema 5.6 como premissa nem a estimativa fina do Lema 3.1 de Ramaré.** Usa apenas identidades elementares de Möbius, a identidade clássica ζ(2)=π²/6, soma por partes, uma cota racional para π e comparação integral. O valor médio de Montgomery–Vaughan/Ramaré **permanece externo**.

## Etapa 1 — a função de contagem de inteiros livres de quadrados

Ponha `Q(x)=sum_{n<=x}mu(n)^2`, `r=sqrt(x)` e `m=floor(r)`. Da identidade

\[
\mu(n)^2=\sum_{d^2\mid n}\mu(d)
\]

segue, por soma finita,

\[
Q(x)=\sum_{d\le r}\mu(d)\left\lfloor\frac{x}{d^2}\right\rfloor.
\]

Como \(\sum_{d\ge1}\mu(d)/d^2=1/\zeta(2)=6/\pi^2=:c\),

\[
|Q(x)-cx|
\le m+x\sum_{d>m}d^{-2}
\le m+\frac{x}{m}.
\]

Para `x>=10^9`, temos `r>30000` e `m>=r-1`. Assim,

\[
m+\frac{x}{m}
\le r+\frac{r^2}{r-1}
=r\left(2+\frac1{r-1}\right)
<\frac{201}{100}\sqrt{x}.
\]

**Esta estimativa de erro com 2,01 é elementar**, bem menos precisa que a cota 0,1333 de Ramaré, mas é suficiente na altura elevada exigida pelo projeto.

## Etapa 2 — primeira parte X<n<=2X

Para `X<n<=2X`, a identidade \(\sum_{d\mid n}\mu(d)=0\) implica
\(\lambda_X(n)=-\mu(n)\). Portanto:

\[
S_1:=\sum_{X<n\le2X}\frac{\lambda_X(n)^2}{n^2}
=\sum_{X<n\le2X}\frac{\mu(n)^2}{n^2}.
\]

Escrevendo `R(t)=Q(t)-ct` e usando soma por partes exata:

\[
S_1=\frac{c}{2X}
+\frac{R(2X)}{4X^2}
-\frac{R(X)}{X^2}
+2\int_X^{2X}\frac{R(t)}{t^3}\,dt.
\]

O limite para `R` dá

\[
XS_1\le\frac3{\pi^2}
+\frac{201}{100\sqrt X}
\left(1+\frac{\sqrt2}{4}+\frac43(1-2^{-3/2})\right)
<
\frac3{\pi^2}
+\frac{201}{100\sqrt X}\frac{65}{24}.
\]

A desigualdade entre parênteses usa somente \(\sqrt2<3/2\) e \(1-2^{-3/2}<1\). As frações do fator `65/24` são verificadas em CI.

## Etapa 3 — segunda parte 2X<n<5X

Para `2X<n<5X`, é possível provar \(|\lambda_X(n)|\le1\) **por casos exaustivos e uniformes**, sem enumerar valores de X:

\[
\lambda_X(n)=-\sum_{\substack{k\mid n\\ k<n/X}}\mu(n/k)
\]

e, já que `n/X<5`, só podem aparecer `k=1,2,3,4`.

- Se \(\mu(n)\ne0\), `n` é livre de quadrados: a parcela `k=4` não existe; se `k=2` ou `3` divide n, \(\mu(n/k)=-\mu(n)\). A soma vale \(-\mu(n)(1-I_2-I_3)\), com `I_2,I_3 in {0,1}`: módulo até 1.
- Se \(\mu(n)=0\) e `4` não divide `n`, não podem ser simultaneamente não nulos \(\mu(n/2)\) e \(\mu(n/3)\); caso contrário `n` seria livre de quadrados. Módulo até 1.
- Se \(\mu(n)=0\) e `4` divide `n`, o termo \(\mu(n/3)\) é zero. Escreva `n=4h`. Se `h` é ímpar, \(\mu(n/2)=\mu(2h)=-\mu(h)=-\mu(n/4)\), de modo que os termos em `k=2,4` se cancelam (se ambos permitidos). Se `h` é par, \(\mu(n/2)=0\), sobrando no máximo o termo \(\mu(n/4)\) de módulo até 1.

No último subcaso, a condição de presença do termo `k=4` garante `n>4X`, que também garante a presença de `k=2` sempre que `2\mid n`, portanto o cancelamento acima é legítimo.

Logo, por monotonicidade de `1/t²`,

\[
S_2:=\sum_{2X<n<5X}\frac{\lambda_X(n)^2}{n^2}
\le\sum_{2X<n<5X}\frac1{n^2}
\le\int_{2X}^{5X}\frac{dt}{t^2}+\frac1{(2X)^2}
=\frac{3}{10X}+\frac1{4X^2}.
\]

## Etapa 4 — fechamento com uma margem racional estrita

Uma cota para \(\pi\) que não depende de aritmética de ponto flutuante segue da identidade de Machin \(\pi=16\arctan(1/5)-4\arctan(1/239)\) e da série alternada de arctan:

\[
\pi>16\left(\frac15-\frac1{375}\right)-\frac4{239}
=\frac{281476}{89625}>\frac{157}{50}.
\]

Como \(\sqrt X>30000\) para `X>=10^9`, obtemos o limite racional **uniforme**:

\[
X(S_1+S_2)
<
\underbrace{\frac{7500}{24649}}_{3/(157/50)^2}
+\underbrace{\frac{201}{100}\frac{65}{24}\frac1{30000}}_{erro\ S_1}
+\underbrace{\frac3{10}+\frac1{4\cdot10^9}}_{S_2}
=\frac{121}{200}
-\frac{161666728553}{295788000000000}.
\]

O termo final é uma **margem positiva de aproximadamente 0,00054656284**. Todas as operações são frações exatas, recalculadas por `b3_proof.py` e verificadas sob Python normal e `-O`.

**Conclusão:** `b3=0.605` está justificada por esta prova elementar para **todo X real >=10^9**. Ramaré Lema 5.6/KLN Lema 3.4 (3.8) deixa de ser uma dependência externa **para esse componente** no NUMERIA. Ainda não provamos independentemente a desigualdade Hilbert/Montgomery–Vaughan (MV), nem o teorema Ramaré–Zúñiga sobre a soma dupla de Möbius.

## Auditoria de origem do teorema de valor médio (dependência ainda aberta)

[Ramaré, *An explicit density estimate for Dirichlet L-series*, Lema 6.5, página impressa 15](https://ramare-olivier.github.io/Maths/DensityEstimate-34.pdf), demonstra que, para coeficientes **reais** `u_n`,

\[
\int_0^T\left|\sum_n u_n n^{it}\right|^2dt
\le\sum_n|u_n|^2[T+\pi m_0(n+1)],\quad
m_0=\sqrt{1+\frac23\sqrt{\frac65}}.
\]

O argumento do artigo usa a desigualdade tipo Hilbert de Montgomery–Vaughan, referida pelo próprio Ramaré como [40]. A redução do fator `2π` para `π` tem explicação transparente: ao integrar os termos cruzados, a expressão sem fase \(\sum_{m\ne n}u_mu_n/(\log m-\log n)\) é nula por antissimetria quando todos os coeficientes são reais; resta apenas a expressão com `exp(iT(log m−log n))`. A separação das frequências satisfaz, para `m≠n`,

\[
|\log m-\log n|\ge\log(1+1/\min(m,n))
\ge\frac1{\min(m,n)+1}.
\]

Isto localiza exatamente o uso de `n+1` e da constante `m0` em KLN Lema 3.6 (3.14). **Ainda não é uma demonstração independente da desigualdade geral de Hilbert**; não substituir a referência MV/KLN por testes finitos.

O integrando do artigo deve ser interpretado com `1/(i\log(m/n))` após a integração; o fator `i` pode ficar implícito ou ausente na tipografia do PDF, mas não altera a antissimetria nem o majorante final.

## Originais e encadeamento

- [Ramaré, Lema 5.6 (b3), Lema 6.5 (MV)](https://ramare-olivier.github.io/Maths/DensityEstimate-34.pdf).
- [Kadiri–Lumley–Ng, Lema 3.4 (3.8), Lema 3.6 (3.14)](https://www.cs.uleth.ca/~kadiri/articles/Explicit-zero-density-for-the-Riemann-zeta-function-JMAA-May2018.pdf).
- [Chourasiya–Simonič v2, Lema 6](https://arxiv.org/html/2507.15184v2).
- [NUMERIA P82 — cadeia condicional F_X](https://github.com/jbarumby/scicalc/pull/82).

**Próximo gate P84:** auditar a desigualdade ponderada de Hilbert com \(\lambda_n=\log n\) e \(m_0\), usando fonte primária Montgomery–Vaughan para fechar a outra metade da estimativa F_X. Depois o único insumo não demonstrado deste ramo será o majorante de Ramaré–Zúñiga para `S_ε(X)` (P78/P81).

**Não prova RH global**, nem aumenta automaticamente o intervalo finito de Suzuki `18<=t<=79`.
