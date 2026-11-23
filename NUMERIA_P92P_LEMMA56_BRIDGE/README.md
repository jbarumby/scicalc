# NUMERIA–RIEMANN P92-P — Ponte exata entre as energias, potenciais e o Lema 6 de Preissmann

**Status:** identidade e equivalência algébricas gerais **demonstradas** para toda configuração **finita**, um critério suficiente analítico para certos pares, e certificações racionais finitas. **Não** é uma prova independente dos Lemas 5 e 6 completos, nem da Hipótese de Riemann.

**Referências:** E. Preissmann, *Sur une inégalité de Montgomery–Vaughan*, *L'Enseignement Mathématique* 30 (1984), 95–113, [texto original](https://www.e-periodica.ch/cntmng?pid=ens-001%3A1984%3A30%3A%3A31), especialmente pp. 104–110; NUMERIA [P92, PR #92](https://github.com/jbarumby/scicalc/pull/92) (prova para separação de até 4 índices); [P88_UPLOAD_AUDIT](https://github.com/jbarumby/scicalc/blob/numeria-p92-preissmann-four-gaps-exact/NUMERIA_P92_FOUR_GAPS/P88_UPLOAD_AUDIT.md) (erro de sinal no comutador).

## 1. Lema de compensação exata P92-P.1 (todas as configurações finitas)

Sejam `x_0<...<x_{N-1}`, N>=2, reais, e `delta_i=min_{j!=i}|x_i-x_j|`. Para `0<=s<t<N`, ponha `h=x_t-x_s>0` e defina

\[
E_i=\sum_{r\ne i}\frac{\delta_r}{(x_i-x_r)^2},\qquad
T_i=\sum_{r\ne i}\frac{\delta_r}{x_i-x_r},\qquad d=\delta_s+\delta_t.
\]

**Identidade EXATA, sem qualquer hipótese espectral:**

\[
\boxed{h^2\sum_{r\ne s,t}\frac{\delta_r}{(x_r-x_s)^2(x_r-x_t)^2}
+\frac{3d}{h^2}
=E_s+E_t-\frac{2(T_s-T_t)}h.}\tag{P1}
\]

### Demonstração em três passos, sem desigualdades

Para os índices r diferentes de s,t, escreva

\[
a_r=\frac1{x_r-x_s},\qquad b_r=\frac1{x_r-x_t},\qquad
I=\sum_{r\ne s,t}\frac{\delta_r}{(x_r-x_s)(x_r-x_t)}.
\]

(1) A identidade `a_r-b_r=-h/[(x_r-x_s)(x_r-x_t)]` dá

\[
h^2\sum_{r\ne s,t}\frac{\delta_r}{(x_r-x_s)^2(x_r-x_t)^2}
=\sum_{r\ne s,t}\delta_r(a_r-b_r)^2
=E_s+E_t-\frac{d}{h^2}-2I.
\]

(2) Separar nos potenciais T_s,T_t os termos em r=t e r=s respectivamente dá

\[
T_s-T_t=hI-\frac{d}{h},\qquad
I=\frac{T_s-T_t}{h}+\frac{d}{h^2}.
\]

(3) Substituição de (2) em (1) produz exatamente P1. **O coeficiente 3 aparece sem Young, sem perda de sinal e sem aproximações:** uma unidade provém de subtrair os dois termos de energia correspondentes ao outro extremo, e duas unidades da identidade de frações parciais.

Em configurações infinitas, P1 exige cuidado adicional com a convergência de T_s,T_t; **não extrapolar** a fórmula finita para séries principais sem justificação de convergência.

## 2. Ponte exata para o Lema 6: novo alvo variacional

A forma usada no projeto é

\[
\sum_{r\ne s,t}\frac{\delta_r}{(x_r-x_s)^2(x_r-x_t)^2}
\le\frac{\pi^2}{3h^2}\left(\frac1{\delta_s}+\frac1{\delta_t}\right)
-\frac{3d}{h^4}.\tag{L6}
\]

Multiplicando por h² e usando P1, temos a **equivalência matemática**

\[
\boxed{E_s+E_t-\frac{2(T_s-T_t)}h
\le\frac{\pi^2}{3}\left(\frac1{\delta_s}+\frac1{\delta_t}\right).}\tag{P2}
\]

Se `G_i=(pi²/3)/delta_i-E_i` é o déficit de energia, a mesma desigualdade é

\[
\boxed{G_s+G_t+\frac{2(T_s-T_t)}h\ge0.}\tag{P3}
\]

**P1, P2 e P3 não são uma prova geral de L6; P2/P3 são reformulações exatamente equivalentes.** No entanto, P3 isola rigorosamente o problema de compensação *entre dois índices* deixado aberto nas investigações P87–P92. A origem do termo negativo afiado fica explícita.

### Corolário P92-P.2 — critério suficiente para separações arbitrárias

A conhecida desigualdade de energia (Lema 1 de Preissmann; tratada na NUMERIA P87) estabelece

\[
E_i\le\frac{\pi^2}{3\delta_i},\quad\text{logo}\quad G_i\ge0.
\]

Assim, **se `T_s>=T_t`, o Lema 6 vale para esse par, sem restrição de `t-s`**. Isto é uma consequência rigorosa da P1 + lema de energia; não depende de importar os Lemas 5/6 completos.

O critério é estritamente suficiente e não necessário: para uma configuração uniforme finita, frequentemente `T_s<T_t` e ainda assim L6 é verdadeiro. Não reivindicamos monotonicidade dos potenciais para todos os pares.

**Testemunho racional além do escopo P92:**

\[
(x_0,\ldots,x_6)=(0,7,8,16,17,32,33),\quad (s,t)=(0,5),
\]

com seis lacunas consecutivas `(7,1,8,1,15,1)`, cinco passos entre os índices. Uma conta inteiramente com Fraction dá

\[
\boxed{T_0-T_5=\frac{31333}{261800}>0.}
\]

Logo o critério energético demonstra L6 para este par de separação cinco, **sem afirmar o teorema para todas as configurações de separação cinco**. Outro testemunho de separação seis é

\[
(0,8000,8300,32300,32500,32510,33510,33515),\quad(s,t)=(0,6),
\]

com `T_0-T_6>47/100>0`, verificado exatamente.

O caso de uma rede uniforme bilateral satura formalmente a compensação no limite (energias `pi²/3` e potenciais simétricos nulos): a perda do termo `-3d/h⁴` destruiria a afiação deste balanço.

## 3. Formulação operatorial genuinamente equivalente via complemento de Schur

Para `r not in {s,t}` defina o vetor de coeficientes **racionais quando os pontos são racionais**

\[
y_r=\frac1{x_r-x_s}-\frac1{x_r-x_t},\quad
P=\operatorname{diag}(\delta_r^{-1}),\quad
A_\kappa=\kappa\left(\frac1{\delta_s}+\frac1{\delta_t}\right)-\frac{3d}{h^2}.
\]

A **matriz explícita** (real simétrica) é

\[
\boxed{\mathcal S_{s,t}(\kappa)=
\begin{pmatrix} P&y\\ y^\top&A_\kappa\end{pmatrix}.}\tag{S1}
\]

Como `P` é positiva definida, seu complemento de Schur é precisamente

\[
A_\kappa-y^\top P^{-1}y
=A_\kappa-\sum_{r\ne s,t}\delta_r y_r^2
=\kappa\left(\frac1{\delta_s}+\frac1{\delta_t}\right)
-E_s-E_t+\frac{2(T_s-T_t)}h.\tag{S2}
\]

Logo

\[
\boxed{\mathcal S_{s,t}(\pi^2/3)\succeq0
\iff \text{L6 é verdadeiro para esse par}.}\tag{S3}
\]

Isto é **uma equivalência exata**, não uma prova universal da positividade. O operador não é abstrato ou arbitrário como na proposta de decompor `H=H0+H0_perp`: seus blocos foram construídos diretamente a partir do núcleo do Lema 6, incluindo o termo corretivo e sem erro de sinal.

### Certificação racional suficiente

Se `kappa_-<pi²/3` é um número racional inferior rigoroso e se a matriz racional `S_st(kappa_-)` é **positiva definida** por LDLᵀ exata, então L6 está rigorosamente certificado para essa configuração racional.

O código usa a identidade de Machin

\[
\frac\pi4=4\arctan(1/5)-\arctan(1/239)
\]

junto com restos racionais da série alternada de arctan para obter `pi_lower<pi<pi_upper` e `kappa_-=pi_lower²/3`. Nunca se utiliza o arredondamento de `float` como majorante rigoroso. Todos os pivôs LDLᵀ usam `fractions.Fraction`.

**A matriz Schur oferece um certificado exato caso a caso e uma formulação para buscar novas desigualdades por blocos. Não converte automaticamente certificados finitos em prova de todos os pares.**

## 4. Por que o texto do Gemini não basta

A decomposição `u=u0+u_perp`, `T=T0+R`, seguida de Cauchy–Schwarz/Young, é formalmente válida sob hipóteses de limitação/autoadjunção. Entretanto:

1. O *Lema 5 de 1984* é um resultado sobre extremalização de `T(x1,...,xN)=sum(delta_k/x_k²)` com restrições nas lacunas, não um teorema abstrato `\langle T0u,u\rangle<=C0||u||²`.
2. Sem cotas explícitas para a compressão ao complemento ortogonal e os acoplamentos, não há estimativa afiada para o operador completo.
3. Young geralmente perde a compensação do coeficiente `-3`; P1 e S2 a preservam exatamente.

O objeto correto a demonstrar no passo seguinte é a positividade de S1 para **todos os pares e configurações**, ou equivalentemente a desigualdade de potenciais P2 por uma comparação variacional do tipo Lema 5. A formulação não supõe essa positividade de antemão.

## 5. Reprodutibilidade e alcance

Do diretório-raiz do repositório:

```bash
python -m unittest -v NUMERIA_P92P_LEMMA56_BRIDGE.test_p92p
python -O -m unittest -v NUMERIA_P92P_LEMMA56_BRIDGE.test_p92p
python -m NUMERIA_P92P_LEMMA56_BRIDGE.certificate --out p92p-certificate.json
```

- **63 configurações com frações exatas**: todas verificam P1, P2/P3 e S2 por duas construções independentes; 
- **27 fatorações LDLᵀ densas exatas**, com pivôs positivos;
- testes incluem um certificado que **rejeita uma cota falsa** (pivô final negativo), invariância afim e testemunhos com separação 5 e 6;
- `pi_lower` é racional e demonstrado estritamente inferior a π;
- os campos do JSON declaram expressamente que Lemas 5/6 gerais e RH permanecem não demonstrados de forma independente.

### Próxima barreira real

Provar uniformemente, para todas as configurações e s<t,

\[
G_s+G_t\ge\frac2h\max(0,T_t-T_s).
\]

O caso `T_s>=T_t` está resolvido pelo lema de energia. O **caso `T_s<T_t`**, especialmente com muitos pontos interiores, requer o mecanismo variacional que o Lema 5 original forneceu. Em vez de ampliar apenas verificações de norma Delta, estudar diretamente este funcional de **dois potenciais + dois déficits de energia**, ou produzir uma família de regiões dos espaçamentos coberta por certificados algébricos.

**Escopo científico:** sem prova independente geral de L5/L6, Preissmann global ou RH. Manter PR como *draft* e não realizar merge automático.