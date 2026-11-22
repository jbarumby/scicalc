# NUMERIA P85 — Hilbert/Montgomery–Vaughan e limiares b4* rigorosos

**Base:** PR #84 (Arb 192/256/384 bits, TABLE S2 robusta, Suzuki até t=79 condicional). Este PR **não mescla** anteriores e **não** reivindica uma prova global da Hipótese de Riemann.

## Parte A — origem do fator pi*m0 no limite de F_X

**Fontes primárias:**

- [Olivier Ramaré, *An explicit density estimate for Dirichlet L-series*, Lema 6.5, p. 15](https://ramare-olivier.github.io/Maths/DensityEstimate-34.pdf). A constante publicada é `m0=sqrt(1+(2/3)*sqrt(6/5))`. Se os coeficientes são reais, o peso `2*pi*m0*(n+1)` melhora para `pi*m0*(n+1)`.
- [E. Preissmann, *Sur une inégalité de Montgomery–Vaughan* (1984), pp. 95–113](https://doi.org/10.5169/seals-53823), desigualdade ponderada de Hilbert, **ainda usada como teorema externo**.
- [H. Montgomery–R. Vaughan, *Hilbert’s inequality* (1974)](https://doi.org/10.1112/jlms/s2-8.1.73).
- [E. Carneiro–F. Littmann, *Monotone extremal functions and the weighted Hilbert’s inequality* (2023)](https://arxiv.org/html/2302.14658), oferece prova alternativa do limite geral com constante `2*pi`, **insuficiente** para substituir a constante `pi*m0≈4.13` na P84.

### Proposição — transferência elementar do teorema ponderado

Seja `N≥2`, `T≥0`, `u_n∈R` e `lambda_n=log n`. A expansão da integral finita é exata:

```
I_N(T) = int_0^T |sum_(n≤N) u_n exp(it log n)|^2 dt
       = T sum_(n≤N)u_n²
         + sum_(m!=n) u_m*u_n (exp(iT(log m-log n))-1)/(i(log m-log n)).
```

**Atenção à tipografia da fonte:** a expressão display do PDF de Ramaré p. 15 omite visualmente o `i` no denominador depois da integração. Pela derivação elementar `int_0^T exp(itd)dt=(exp(iTd)-1)/(id)`, esse `i` é **obrigatório**.

O termo inferior `sum_(m!=n)u_m*u_n/(log m-log n)` é **exatamente zero**, agrupando `(m,n)` e `(n,m)`, pois `u_m u_n` é simétrico e o denominador antissimétrico. `hilbert_source.py` testa quatro configurações independentes com `fractions.Fraction`; a demonstração geral dispensa os testes.

Se `delta_n` é a separação mínima das frequências usadas, como `m` e `n` são inteiros distintos,

```
|log m - log n| >= log(1+1/n)
                = integral_0^(1/n) dt/(1+t)
                >= 1/(n+1).
```

Assim `1/delta_n≤n+1`. **Importando** a desigualdade ponderada de Preissmann, que tem constante `pi*m0`, e escolhendo `a_n=u_n exp(iT log n)`:

```
|sum_(m!=n) a_m*conj(a_n)/(log m-log n)|
    ≤ pi*m0 * sum_n |a_n|²/delta_n
    ≤ pi*m0 * sum_n(n+1)u_n².
```

A diferença de apenas um termo de fronteira (porque o inicial cancelou) dá

```
I_N(T) ≤ sum_(n≤N)u_n² [T + pi*m0*(n+1)].
```

**Aplicação a F_X:** para `sigma>1`, tome `u_n=lambda_X(n)n^(-sigma)`, com

```
lambda_X(n) = sum_(d|n,d≤X) mu(d) - 1_(n=1), 
|lambda_X(n)| ≤ floor X.
```

A convergência absoluta e uniforme do polinômio limite decorre de `sum n^(-sigma)<∞`; além disso `sum (n+1)|u_n|²<∞` porque `sum n^(1-2sigma)<∞`. O limite `N→∞` é legítimo por convergência uniforme no compacto `[0,T]` e convergência das séries de majorantes. Portanto,

```
F_X(sigma,0,T,1/2) ≤ pi*m0 * sum_n lambda_X(n)²/n^(2sigma-1)
                   +(T+pi*m0)*sum_n lambda_X(n)²/n^(2sigma).
```

A transição do teorema ponderado para `F_X` está **fechada com prova elementar**, e a norma constante satisfaz rigorosamente `263/200 < m0 < 329/250`, verificado com frações exatas.

**Limite explícito do alcance:** a **desigualdade ponderada de Preissmann** *não foi provada independentemente* nesta P85; apenas sua aplicação foi reconstruída. Sem essa entrada, não afirmar que o ramo de Hilbert está provado a partir de primeiros princípios.

## Parte B — limiares b4* com Arb, sem ajustar float

Na P84 ficou pendente o valor do limiar `b4*` para o `B1` original e para a TABLE S2 robusta.

Nas sete linhas afetadas `i=10,…,16`, o cálculo P79 satisfaz, **na região K3>1**,

```
K3(b4) = alpha*b4 + beta, alpha>0, beta>0,
e = (2*sigma2-1)/(2-sigma2) = 2i/(48-i) >0,
B1_majorant(b4) = D * (alpha*b4+beta)^e, D>0.
```

**Prova da fatoração:** os fatores da fórmula P77 `max(1,K3^(-3delta/L))` e `max(1,K3^(-3d/L))` valem 1 quando K3>1; o `max(K3^e(sigma1),K3^e(sigma2))` se reduz a `K3^e(sigma2)`, pois `e(sigma)=(2sigma-1)/(2-sigma)` é estritamente crescente, `e'(sigma)=3/(2-sigma)^2>0`. Todos os outros fatores podem ser agrupados em `D`, independente de b4. Assim, a derivada é

```
dB1/db4 = D * e * alpha * K3^(e-1) > 0.
```

Para um alvo impresso `P`, a raiz da **expressão majorante implementada** é

```
b4*(P) = ((P/D)^(1/e)-beta)/alpha.
```

`b4star_arb.py` reconstrói `D,alpha,beta` diretamente das fórmulas P77/P79 usando quadratura de K1. Ele forma a raiz como intervalo Arb em `ctx.workprec` 192/256/384, exporta endpoints `man_exp` **exatos** e fornece um intervalo racional com extremos de 10 casas decimais, ambos verificados por sinais intervalares estritos:

```
B1(raiz_inferior).upper() < P
B1(raiz_superior).lower() > P.
```

A exigência extra nas sete linhas é a separação **rigorosa**, para cada precisão:

```
b4*_old.upper < b4_public < b4*_S2robust.lower.
```

Com isso, **a falha dos antigos B1 e o sucesso da TABLE robusta** ficam demonstrados também na linguagem de limiares, sem aceitar `b4star_upper` float do DOC1/S3 como certificado.

**Cuidado:** a raiz é da expressão \(B_1\) *majorante escolhido pelo algoritmo P77*, com `A0`, `K1` e as fontes fixadas, não um limiar universal sobre a zeta.

## Reprodução

```bash
python -m pip install python-flint==0.9.0 mpmath==1.3.0
python -m unittest -v NUMERIA_P85_HILBERT_B4STAR.test_p85
python -O -m unittest -v NUMERIA_P85_HILBERT_B4STAR.test_p85
python -m NUMERIA_P85_HILBERT_B4STAR.hilbert_source --out p85-hilbert-source.json
for bits in 192 256 384; do
  python -m NUMERIA_P85_HILBERT_B4STAR.b4star_arb --bits "$bits" --out "p85-b4star-\${bits}.json"
done
```

Após executar as três precisões, a CI exige compatibilidade das 14 raízes em cada faixa, margens positivas, hashes de scripts e insumos. **Apenas uma CI aprovada** permite declarar a recertificação computacional da Parte B.

## O que fica aberto

- Prova independente da **desigualdade ponderada geral de Preissmann**. A prova por Fourier de Carneiro–Littmann dá constante `2pi`, mais fraca que a exigida.
- Prova independente das estimativas de Ramaré–Zúñiga para `S_epsilon(X)` e cobertura pública de `3e10≤X<3e12`.
- O argumento original sobre monotonicidade global de `G(epsilon)` e os 28 testes zeta do DOC2 S3 não são reconstruídos aqui.
- A verificação Suzuki P84 até `t=79` permanece finita e condicional; não há prova global de RH.

**Fontes matemáticas utilizadas:** Ramaré Lema 6.5, Preissmann (1984), Montgomery–Vaughan (1974), Carneiro–Littmann (2023), e os próprios certificados P77/P79/P84 do repositório.
