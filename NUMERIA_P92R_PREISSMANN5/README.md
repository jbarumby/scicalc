# NUMERIA–RIEMANN P92-R — Reconstrução independente do Lema 5 quantitativo e passagem candidata ao Lema 6 geral

**Estado de auditoria:** o Lema 5 quantitativo abaixo tem uma demonstração analítica autocontida; a passagem ao Lema 6 também é deduzida por um argumento analítico explícito e deve receber **revisão matemática independente**, pois resultados computacionais e CI não substituem a avaliação de uma prova geral. Não há reivindicação de prova de RH ou de novidade histórica: os Lemas 5–6 já são teoremas publicados de Preissmann (1984). Esta é uma **reconstrução independente** no projeto, sujeita a revisão.

**Fonte de comparação:** E. Preissmann, *Sur une inégalité de Montgomery–Vaughan*, *L'Enseignement Mathématique* **30** (1984), pp. 95–113; Lemas 5 e 6, pp. 104–110. PDF do usuário “Reconstrução Detalhada: A Passagem do Lema 5 para o Lema 6”, 5 páginas, auditado na P92-Q. Os três reparos naquele PDF estão registrados em `NUMERIA_P92Q_LATTICE_THINNING/LEMA5_TO_LEMA6_UPLOAD_AUDIT.md`. Esta nota NÃO depende de usar suas fórmulas danificadas.

## 1. O lema fundamental de Bellman para funções convexas

Seja `f:(0,+∞)→[0,+∞)` **convexa e decrescente**, com suporte compacto, e `w,h>0`. Defina

\[
U_0(w,h)=h\sum_{j=0}^{\infty}f(w+jh),\qquad U_1(w,h)=U_0(w,h)-hf(w)=h\sum_{j=1}^{\infty}f(w+jh).
\]

**Lema R1.** Para `w>0` fixo, `U_0(w,h)` é não decrescente em `h>0` e `U_1(w,h)` é não crescente. Em particular, para qualquer `b,h>0`,

\[
\boxed{\min(h,b)f(w)+U_0(w+b,b)\le U_0(w,h).}\tag{R1}
\]

**Demonstração independente.** Primeiro considere uma função-rampa `f_T(u)=(T-u)_+`, com `L=(T-w)_+`. Se `L>0` e `n=\lceil L/h\rceil`, então exatamente

\[
U_0=hnL-\frac{h^2n(n-1)}2,\qquad U_1=h(n-1)L-\frac{h^2n(n-1)}2.
\]

Em todo intervalo em que `n` é constante, as derivadas são

\[
\partial_h U_0=n[L-(n-1)h]>0,\qquad
\partial_h U_1=(n-1)[L-nh]\le0.
\]

As expressões são contínuas nas fronteiras de mudança de `n`. Se `L=0`, ambas são nulas. Logo as duas monotonias valem para cada rampa. Uma função convexa, decrescente, não negativa e de suporte compacto é uma combinação integral não negativa de rampas: `f(u)=\int (T-u)_+\,d\mu(T)`, onde `\mu` é a medida de segunda derivada distribucional; não é preciso que `f` seja `C³`. Pelo teorema de Tonelli, as monotonias passam para `f`.

Se `b<=h`, então

\[
\min(h,b)f(w)+U_0(w+b,b)=U_0(w,b)\le U_0(w,h).
\]

Se `b>=h`, o lado esquerdo é `hf(w)+U_1(w,b)\le hf(w)+U_1(w,h)=U_0(w,h)`. Isto prova (R1).

A mesma prova vale para funções convexas decrescentes expressáveis como misturas positivas de rampas com séries finitas; em particular **`f(u)=u^{-2}`**, pois

\[
u^{-2}=\int_{u}^{\infty}\frac{6(T-u)}{T^4}\,dT.
\]

## 2. Comparação unilateral para espaçamentos reais arbitrários

Seja `0=z_0<z_1<⋯<z_N`, com `g_j=z_j-z_{j-1}>0`, `a=g_1=z_1`, e

\[
\delta_j=\min(g_j,g_{j+1})\quad(1\le j<N),\qquad\delta_N=g_N.
\]

**Lema R2 (generalização independente do controle unilateral).** Para toda `f` da classe R1,

\[
\boxed{\sum_{j=1}^{N}\delta_j f(z_j)\le a\sum_{k=1}^{\infty} f(ka).}\tag{R2}
\]

**Prova:** começando pelo último termo, `\delta_N f(z_N)=g_N f(z_N)\le U_0(z_N,g_N)`; substitua na soma do penúltimo termo e aplique (R1) com `w=z_{N-1},h=g_{N-1},b=g_N`. Continue retroativamente. Ao chegar a `z_1=a`, resta a majorante `U_0(a,a)=a\sum_{k\ge1}f(ka)`. Todas as desigualdades são válidas para quaisquer **valores reais positivos** de `g_j`, sem procurar pontos críticos do funcional.

Para `f(u)=u^{-2}`, a representação integral anterior e Tonelli permitem usar (R2), obtendo

\[
\boxed{\sum_{j=1}^{N}\frac{\delta_j}{z_j^2}\le\frac{\zeta(2)}a.}\tag{R2b}
\]

## 3. Lema 5 quantitativo, com erro explícito mais forte

**Teorema R3.** Com as hipóteses de R2, seja `X=z_N` e `K=\lfloor X/a\rfloor`. Então

\[
\boxed{\sum_{j=1}^{N}\frac{\delta_j}{z_j}\le H_K+\frac{X-Ka}{X}<H_K+\frac aX,\qquad H_K=\sum_{k=1}^{K}\frac1k.}\tag{R3}
\]

A primeira desigualdade **pode ser igualdade**, e a segunda é estrita (porque `0<=X-Ka<a`).

**Prova integral e sem erro assintótico implícito.** Defina para `u>0`

\[
f_X(u)=\max\left(\frac1u-\frac1X,0\right).
\]

Ela é convexa, decrescente, não negativa e de suporte compacto quando restrita a `u>=a>0`. Sua segunda derivada como medida no intervalo positivo é `2u^{-3}du` para `u<X`, mais um átomo de massa `X^{-2}` em `u=X`. Portanto cumpre R1. Pela comparação unilateral,

\[
\sum_{j=1}^N\delta_j\left(\frac1{z_j}-\frac1X\right)
\le a\sum_{k=1}^{K}\left(\frac1{ka}-\frac1X\right)
=H_K-\frac{Ka}{X}.
\]

Adicionando `(\sum_j\delta_j)/X` e observando `\delta_j\le g_j` e `\sum g_j=X`, obtemos exatamente R3. **Não há constante oculta:** o erro de borda é `(X-Ka)/X<a/X`, mais forte que a expressão `O(a/X)` com constante `e^5` presente no roteiro do usuário.

O coeficiente `1` no erro simples `a/X` não pode ser reduzido universalmente: escolha `g_1=a, g_2=a-\varepsilon` e faça `\varepsilon\downarrow0` (com `X/a=2-\varepsilon/a` não inteiro). A diferença da soma para `H_{\lfloor X/a\rfloor}` tende a `1/2`, assim como `a/X`.

## 4. Comparação rigorosa das caudas do núcleo do Lema 6

Seja `h>0` e `\phi_h(u)=[u^2(u+h)^2]^{-1}` para `u>0`. A função é positiva, decrescente e estritamente convexa; pode-se aplicar R2 à cauda de pontos à esquerda ou à direita de um par de âncoras. Seja `g` a primeira lacuna da cauda e `d` a distância ao vizinho mais próximo na âncora; `d<=g`. Por R2,

\[
\sum_{j\ \mathrm{cauda}}\delta_j\phi_h(\mathrm{dist}_j)
\le g\sum_{k\ge1}\phi_h(kg)
\le d\sum_{k\ge1}\phi_h(kd),\tag{R4}
\]

pois para cada `k>=1` a função `v\mapsto v\phi_h(kv)=1/[k^2 v(kv+h)^2]` é decrescente. Portanto as duas caudas de **qualquer configuração finita** são majoradas pelas duas caudas aritméticas com passos `\delta_t` e `\delta_s`. Os pontos entre as âncoras e as distâncias `\delta_t,\delta_s` dos extremos permanecem inalterados. A comparação usa somente R2, não os Lemas 5–6 publicados.

## 5. Compensador de dois índices — demonstração sem notação O

Fixe a extensão bilateral `(...,y_{t-1},y_t,...,y_s,y_{s+1},...)` assim: `t<s`, `y_t=0`, `y_s=h>0`, pontos internos arbitrários. Para `r<t`, `y_r=(r-t)L` com `L=\delta_t`; para `r>s`, `y_r=h+(r-s)R` com `R=\delta_s`. Os espaçamentos extremos internos satisfazem `a_t=y_{t+1}-y_t\ge L` e `a_s=y_s-y_{s-1}\ge R`. Todos os `\delta` interiores e os dos extremos são os da configuração original.

Escolha um inteiro `A>=1` e ponha

\[
Z=h+AL,\qquad J=\left\lfloor\frac ZR\right\rfloor,\qquad Z'=h+JR.
\]

Trunque os índices de `r=t-A` até `r=s+J`, e defina (sempre excluindo `t,s`)

\[
\Sigma_t^{A}=\sum_r\frac{\delta_r}{y_r-y_t},\qquad
\Sigma_s^{A}=\sum_r\frac{\delta_r}{y_r-y_s}.
\]

A soma unilateral **à esquerda de s**, incluindo o nó t antes de excluí-lo, satisfaz pelo Teorema R3

\[
P_s=\sum_{r=t-A}^{s-1}\frac{\delta_r}{y_s-y_r}
\le H_{\lfloor Z/a_s\rfloor}+\frac{a_s}{Z}.
\]

A soma unilateral **à direita de t**, incluindo o nó s, satisfaz

\[
P_t=\sum_{r=t+1}^{s+J}\frac{\delta_r}{y_r-y_t}
\le H_{\lfloor Z'/a_t\rfloor}+\frac{a_t}{Z'}.
\]

As partes aritméticas das caudas opostas são exatamente `H_J` e `H_A`. Ao remover os termos t e s nas duas somas, obtemos a **identidade exata**

\[
\Sigma_t^A-\Sigma_s^A
=P_t+P_s-H_A-H_J-\frac{L+R}{h}.
\]

Como `a_s>=R`, tem-se `\lfloor Z/a_s\rfloor<=J`; logo `H_{\lfloor Z/a_s\rfloor}<=H_J`. Escreva `K'=\lfloor Z'/a_t\rfloor`. Para `K'<=A` também `H_{K'}-H_A<=0`; para `K'>A`,

\[
H_{K'}-H_A\le\frac{K'-A}{A}.
\]

Além disso, `a_t>=L` e `JR<=h+AL` dão

\[
K'-A\le \frac{h+JR}{L}-A\le\frac{2h}{L}.
\]

Assim, nos dois casos:

\[
\boxed{
\Sigma_t^A-\Sigma_s^A
\le-\frac{L+R}{h}
+\frac{a_s}{Z}+\frac{a_t}{Z'}+\frac{2h}{AL}.}\tag{R5}
\]

A expressão de erro é **racional quando todos os espaçamentos são racionais**, é não negativa, e tende a `0` quando `A→∞`. O número de pontos da cauda direita também cresce (`J→∞`). As diferenças individuais

\[
\delta_r\left(\frac1{y_r-y_t}-\frac1{y_r-y_s}\right)
\]

são absolutamente somáveis, pois nas duas caudas equidistantes decaem como `O(r^{-2})`. Portanto R5 implica a **compensação infinita de dois índices**

\[
\boxed{\Sigma_t-\Sigma_s\le-\frac{\delta_t+\delta_s}{h}.}\tag{R6}
\]

O argumento inclui a parte positiva necessária quando `K'<=A`, reparando o erro da eq. (21) no PDF recebido. Não depende da alegação falsa de que as lacunas pontuais aumentam após prolongar as caudas.

## 6. Consequência: reconstrução candidata do Lema 6 afiado para todos os pares

A identidade algébrica **exata**, para `r,s,t` distintos e `h=y_s-y_t>0`, é

\[
\frac1{(y_r-y_s)^2(y_r-y_t)^2}
=\frac1{h^2}\left(\frac1{(y_r-y_s)^2}+\frac1{(y_r-y_t)^2}\right)
+\frac2{h^3}\left(\frac1{y_r-y_t}-\frac1{y_r-y_s}\right).
\]

A soma sobre `r\notin\{s,t\}` é absolutamente convergente (a combinação do último termo converge). Por R6 o termo misto é limitado por `-2(\delta_s+\delta_t)/h^4`. Por R2b nas duas direções da âncora e `\zeta(2)=\pi^2/6`, para `i=s,t`,

\[
\sum_{r\ne i}\frac{\delta_r}{(y_r-y_i)^2}\le\frac{2\zeta(2)}{\delta_i}=\frac{\pi^2}{3\delta_i}.
\]

Excluindo a outra âncora dos somatórios, retira-se também `\delta_t/h^2` do lado de `s` e `\delta_s/h^2` do lado de `t`. Somando,

\[
\boxed{\sum_{r\ne s,t}\frac{\delta_r}{(y_r-y_s)^2(y_r-y_t)^2}
\le\frac{\pi^2}{3h^2}\left(\frac1{\delta_s}+\frac1{\delta_t}\right)
-\frac{3(\delta_s+\delta_t)}{h^4}.}\tag{R7}
\]

A comparação de caudas R4 transfere R7 da configuração prolongada `y` a **qualquer configuração finita original** `x`, preservando as distâncias dos extremos e todos os pontos internos. Todos os passos R1–R7 são analíticos e não requerem assumir os Lemas 5 ou 6 publicados. **Classificação:** reconstrução autônoma candidata à validação matemática independente; os testes e a CI não são a prova do enunciado geral, que depende da verificação das justificações analíticas acima. Não se deve declarar este texto uma prova formalizada em Lean ou uma publicação revisada por pares.

## 7. Regressões racionais e limites do projeto

`certificate.py` usa apenas `fractions.Fraction`, verificando 222 configurações do Lema 5, 400 desigualdades Bellman para rampas, 12 compensadores de truncamentos assimétricos e 60 configurações Lema 6 completo usando um limite **inferior** exato para pi pela fórmula de Machin. `test_p92r.py` acrescenta 10 testes unitários com casos-limite, simetria/reflexão e invariância afim. A execução independente exploratória conferiu ainda 350 pares internos/exteriores com lacunas arbitrárias racionais.

Reprodução na raiz do repo:

    python -m unittest -v NUMERIA_P92R_PREISSMANN5.test_p92r
    python -O -m unittest -v NUMERIA_P92R_PREISSMANN5.test_p92r
    python -m NUMERIA_P92R_PREISSMANN5.certificate --out p92r-certificate.json

**Nota sobre as duas implementações:** no repositório GitHub a execução compacta de revisão usa `python -m unittest -v NUMERIA_P92R_PREISSMANN5.test_p92r_github` e `python -m NUMERIA_P92R_PREISSMANN5.verify --out p92r-review.json`. Os comandos `test_p92r`/`certificate` acima correspondem ao pacote local completo entregue ao usuário. O teorema e os argumentos analíticos são os mesmos; os hashes dos JSONs diferem por terem esquemas e apresentações distintos.

**Pendências científicas:** revisão externa da passagem global R1–R7, comparação contra o enunciado original sem OCR, rigor formal opcional em Lean/Isabelle/Coq. Não inferir nem declarar a Hipótese de Riemann demonstrada. Não realizar merge automático.