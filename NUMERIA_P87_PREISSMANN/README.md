# NUMERIA P87 — Lema 1 de Preissmann reconstituído; Lema 6 em auditoria separada

**Novo resultado analítico desta etapa:** a desigualdade de energia diagonal correta
\[
\boxed{\qquad Q_i=\sum_{j\ne i}\frac{\delta_i\delta_j}{(x_i-x_j)^2}\le\frac{\pi^2}{3},\quad
\delta_i=\min_{j\ne i}|x_i-x_j| \qquad}
\]
é demonstrada para **todo N finito ≥2 e qualquer sequência de pontos reais distintos**. A cota `Q_i<=4` da P86 era válida, porém insuficiente para a constante `m0`; **agora a desigualdade mais forte está demonstrada** por uma prova analítica geral, não por enumeração de matrizes.

## Referências e relação com a fonte original

- [Wijit Yangjit, *On the Montgomery–Vaughan weighted generalization of Hilbert's inequality*, arXiv:2203.14950v1 (2022), §2.2, Lemas 2.2–2.5, pp. 4–7](https://arxiv.org/pdf/2203.14950). O autor explica que o resultado é consequência do Lema 1 de Preissmann e fornece uma **demonstração alternativa pelo método de Shan**. Reconstituímos integralmente as três etapas de Shan abaixo em vez de tomar o Lema 1 como hipótese.
- [E. Preissmann, *Sur une inégalité de Montgomery et Vaughan*, Enseign. Math. 30 (1984), 95–113](https://www.e-periodica.ch/cntmng?pid=ens-001:1984:30::29). Fonte histórica do Lema 1 e da desigualdade cruzada do Lema 6. O acesso direto ao PDF original apresentou restrição de site; para a **fórmula exata do Lema 6** e a combinação espectral, utilizamos Yangjit §3.2, Lema 3.2 e Teorema 1.2.
- [Yangjit, §3.1–3.2, pp. 8–9](https://arxiv.org/pdf/2203.14950), passagem algébrica da desigualdade com três denominadores à constante de Preissmann.
- Fase-base P86: [PR #86](https://github.com/jbarumby/scicalc/pull/86); P84/P85 permanecem preservadas.

## 1. Teorema unidirecional de Shan — versão geral e prova completa

Seja `f:[1,+∞)->[0,+∞)` **decrescente, convexa** e com `sum_(k>=1)f(k)<∞`. Para uma sequência positiva `a=(a_1,a_2,...)` com **a_1>=1** e para qualquer `N>=1`, defina
\[
s_n=\sum_{k=1}^n a_k,\quad
F_N(a)=\sum_{n=1}^N\min(a_n,a_{n+1})f(s_n).
\]
**Proposição:** \(F_N(a)\le\sum_{j=1}^{\infty}f(j)\) para toda sequência admissível.

### 1a. Lema de interpolação e contagem por inteiros

Se `b_n>=1` para todos `n`, escreva `t_n=sum_(k<=n)b_k`, `m_n=floor(t_n)`, `theta_n=t_n-m_n`, com `t_0=m_0=theta_0=0`. Pela convexidade:
\[
f(t_n)\le (1-\theta_n)f(m_n)+\theta_n f(m_n+1).
\]
Como `f` decresce e `b_n>=1`,
\[
b_nf(t_n)\le (b_n-\theta_n)f(m_n)+\theta_n f(m_n+1).
\]
Usando `b_n-\theta_n=(m_n-m_{n-1}-1)+(1-\theta_{n-1})`, `m_n>=m_{n-1}+1` e a monotonicidade de `f`, a primeira parcela fica limitada por
\[
\sum_{j=m_{n-1}+1}^{m_n}f(j)-\theta_{n-1}f(m_{n-1}+1).
\]
A soma telescópica em `n` fornece a desigualdade **finita e explícita**
\[
\boxed{\sum_{n=1}^N b_nf(t_n)
 \le\sum_{j=1}^{\lfloor t_N\rfloor}f(j)
      +\{t_N\}f(\lfloor t_N\rfloor+1)
 \le\sum_{j=1}^{\infty}f(j).}
\]
A última passagem exige `f>=0`, razão de essa hipótese constar no enunciado.

### 1b. Lema de eliminação dos descensos

Se `a_{ν-1}>a_ν` para algum `ν>=2`, aumente `a_ν` em
`0<ε<=a_{ν-1}-a_ν`, sem mudar os demais termos, obtendo `b`. O peso `min(a_(ν-1),a_ν)` aumenta em `ε`; o peso seguinte não diminui; os argumentos `s_n` para `n>=ν` aumentam em `ε`. Pela não negatividade de `f`, e por convexidade,
\[
\begin{aligned}
F_N(b)-F_N(a)
&\ge \varepsilon f(s_{\nu-1})
+\sum_{n=\nu}^N\min(a_n,a_{n+1})[f(s_n+\varepsilon)-f(s_n)]\\
&\ge \varepsilon f(s_{\nu-1})
+\varepsilon\sum_{n=\nu}^N\frac{\min(a_n,a_{n+1})}{a_n}
   [f(s_n)-f(s_{n-1})]\\
&\ge \varepsilon f(s_{\nu-1})
+\varepsilon\sum_{n=\nu}^N[f(s_n)-f(s_{n-1})]\\
&=\varepsilon f(s_N)\ge0.
\end{aligned}
\]
A terceira linha usa o fato de o incremento de `f` ser **não positivo** e `min(a_n,a_(n+1))/a_n<=1`. Se `ν=N+1`, apenas o último peso aumenta, e o argumento é ainda mais simples.

Eliminando sucessivamente todos os descensos em `a_2,...,a_(N+1)`, chegamos ao envelope **não decrescente**
\(\overline a_n=\max_{1\le k\le n}a_k\), com \(F_N(a)\le F_N(\overline a)\).

### 1c. Fechamento do teorema

Como o envelope é não decrescente e \(\overline a_1=a_1\ge1\), vale
\(\min(\overline a_n,\overline a_{n+1})=\overline a_n\) e todos os termos são ao menos 1. Aplicando 1a:
\[
F_N(a)\le F_N(\overline a)
=\sum_{n=1}^N\overline a_n f\!\left(\sum_{k=1}^n\overline a_k\right)
\le\sum_{j=1}^{\infty}f(j).
\]
Isso é uma prova geral completa, para quaisquer entradas **reais** que satisfaçam as hipóteses, independente de testes computacionais.

## 2. Da desigualdade de Shan ao Lema 1 de Preissmann

Ordene \(x_1<\cdots<x_N\) e fixe `i`. Separe a contribuição à direita, `j>i`. Escreva `h_n=x_(i+n)-x_(i+n-1)`, `a_n=h_n/delta_i`. Temos `a_1>=1`, e
\[
\frac{\delta_{i+n}}{\delta_i}=\min(a_n,a_{n+1})
\]
para nós interiores. Para o último nó finito, adote `a_(N-i+1)=a_(N-i)`, a **extensão virtual** que codifica corretamente o espaçamento de borda. Além disso,
\[
\frac{x_{i+n}-x_i}{\delta_i}=\sum_{k=1}^na_k.
\]
Portanto, com `f(t)=t^{-2}` (positiva, decrescente, convexa),
\[
\sum_{j>i}\frac{\delta_i\delta_j}{(x_i-x_j)^2}
=F_{N-i}(a)\le\sum_{k\ge1}\frac1{k^2}=\frac{\pi^2}{6}.
\]
Simetricamente, a contribuição à esquerda também é no máximo \(\pi^2/6\); a soma conclui
\[
\boxed{\forall N\ge2,\ \forall x_i\in\mathbb R\text{ distintos}:\ Q_i\le\pi^2/3.}
\]

A constante é **uniformemente ótima**: para `N=2M+1` com pontos igualmente espaçados, o ponto central satisfaz
\[
Q_i=2\sum_{k=1}^M k^{-2}\longrightarrow\pi^2/3
\]
quando `M→∞`. Esse argumento também explica por que a cota 4 da P86 não era ótima.

## 3. Lema 6 de Preissmann: enunciado exato e estado da auditoria

A fórmula original foi localizada independentemente da digitalização primária por Yangjit, Lema 3.2 (equação 3.6, p. 9). Para `l!=m`:
\[
\sum_{k\notin\{l,m\}}
\frac{\delta_k}{(x_k-x_l)^2(x_k-x_m)^2}
\le
\frac{\pi^2(\delta_l+\delta_m)}
     {3\delta_l\delta_m(x_l-x_m)^2}
-\frac{3(\delta_l+\delta_m)}{(x_l-x_m)^4}.
\tag{P87-L6}
\]
A P87 **ainda não prova (P87-L6) para todos os N e espaçamentos**. Em vez de mascarar essa lacuna:
- `lemma6_bridge.py` verifica um conjunto explícito de instâncias finitas **com `fractions.Fraction`**, substituindo `π` por um limite racional **estritamente inferior** obtido com 12 termos da identidade de Machin. É um teste mais exigente que o enunciado para cada instância, não prova universal.
- O certificado descreve cada par testado, os valores racionais exatos e a margem positiva, e mantém `Preissmann_Lemma6_all_N_proved_here=false`.

## 4. Combinação espectral de Preissmann — prova condicional, álgebra exata

No esquema das matrizes da P86, a identidade dos quadrados da matriz antissimétrica e o Lema 1 produzem
\[
\|B\|_2^2\le\frac{\pi^2}{3}+2\frac UV,
\quad
U=\sum_{m\ne n}\frac{\delta_m^{3/2}\delta_n^{1/2}t_mt_n}{(x_m-x_n)^2},
\quad V=\sum_n t_n^2,
\]
com `t_n>=0`. A estimativa de Cauchy e de Shan com `σ=4` fornecem um termo \(S\le\pi^4 V/45\). **Se** (P87-L6) é admitido, o outro termo satisfaz \(T\le(2\pi^2/3)U\). Assim,
\[
U^2\le V(S+T)\le V\left(\frac{\pi^4}{45}V+
                                  \frac{2\pi^2}{3}U\right).
\]
Para \(V>0\), escrevendo \(u=U/(\pi^2V)\), temos
\(u^2\le1/45+(2/3)u\), cuja raiz positiva é
\[
u\le\frac{1+\sqrt{6/5}}3.
\]
Inserindo o resultado:
\[
\|B\|_2^2\le\frac{\pi^2}{3}
+2\pi^2\frac{1+\sqrt{6/5}}3
=\boxed{\pi^2\left(1+\frac23\sqrt{\frac65}\right)}.
\]
Logo \(\|B\|_2\le\pi m_0\), exatamente a constante usada em `F_X`, **condicionalmente à prova geral do Lema 6**. A P87 verifica a identidade da equação quadrática **no anel `Q[sqrt(6/5)]`**, sem qualquer `float`.

**O que foi fechado:** Lema 1 em todas as dimensões por uma prova matemática completa no README; todas as transformações da combinação espectral *assumindo* Lema 6. **O que permanece externo:** o Lema 6 uniforme, demonstrado em Preissmann (1984) mas não refeito independentemente aqui.

## 5. Novos uploads S3: deduplicação, não recertificação repetida

Os cinco arquivos de entrada desta rodada foram conferidos por SHA-256 e correspondem aos originais já auditados na P84:
- `p78_s3_certificado_oficial (2).json`: `6c0be960a3824310aa73b04bbb2f5bf1eb54a688b845ab8391bb542b48fdef1a` — igual ao `(1)`.
- `p78_s3_adocao_table (2).json`: `35a8bd405b7a368587a842318448f8af3513b406e78d77882c8a39eb0d065fb2`.
- `p78_s3_doc1_certificado (2).json`: `093dd70b4fc6ddf0385094cd112699034f2fad2c36a3c51758c7f7c044652137`.
- `p78_s3_doc2_certificado (1).json`: `ac356d2997e93945e8ce0751dd9e2f00054b45b9284a800d3ed690b57d8c6846`.
- `p78_s3_nota_adocao_repo_certified (1).pdf`: `c916aa2ff5dc59ff19c0e89177b8bf9335d82322fd91e1a3be7f9355d4d1bb89`.

O rótulo `REPO_CERTIFIED` no S3 original não substitui os certificados Arb da P84/P85. Não reexecutamos o mesmo kit nem sobrescrevemos a TABLE antiga; os sete coeficientes S2 robustos já estão em versão separada, com quadratura intervalar, SHA e cobertura finita condicional Suzuki até `t=79`.

## 6. Reproduzir P87

```bash
python -m unittest -v NUMERIA_P87_PREISSMANN.test_p87
python -O -m unittest -v NUMERIA_P87_PREISSMANN.test_p87
python -m NUMERIA_P87_PREISSMANN.shan_spacing --out p87-shan.json
python -m NUMERIA_P87_PREISSMANN.lemma6_bridge --out p87-lemma6.json
```

A CI executa e revalida os exemplos com frações exatas, a contagem dos pares verificados e as identidades de \(m_0\), gerando hashes SHA-256. **Não confundir a aprovação desses testes com uma prova automática da demonstração geral escrita no README.**

## Próxima prioridade — P88

Reconstruir *de fato* o **Lema 6 em dimensão arbitrária**, em especial sua parcela negativa
\(-3(\delta_l+\delta_m)/(x_l-x_m)^4\), por decomposição dos intervalos à esquerda, entre e à direita dos dois índices. Só depois declarar a desigualdade de Preissmann inteiramente rederivada e retirar sua dependência de \(F_X\).

Ainda abertos: Ramaré–Zúñiga para \(S_\varepsilon(X)\) e a faixa \(3×10^{10}\le X<3×10^{12}\). **Não há prova global da Hipótese de Riemann.**
