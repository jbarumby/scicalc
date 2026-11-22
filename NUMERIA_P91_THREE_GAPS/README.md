# NUMERIA–RIEMANN P91 — Lema 6 afiado para dois pontos intermediários

**Base:** [PR #90](https://github.com/jbarumby/scicalc/pull/90), branch numeria-p90-preissmann-interior-variational.

## Teorema novo (escopo delimitado)

Para qualquer configuração finita de pontos reais distintos e quaisquer índices a < b com **b-a <= 3**, o Lema 6 afiado de Preissmann vale com a constante original:

\[
\sum_{k\ne a,b}\frac{\delta_k}{(x_k-x_a)^2(x_k-x_b)^2}
\le \frac{\pi^2}{3h^2}\left(\frac1{\delta_a}+\frac1{\delta_b}\right)
-\frac{3(\delta_a+\delta_b)}{h^4},
\quad h=x_b-x_a.
\]

P90 demonstrou os casos b-a=1,2. A **contribuição de P91 é b-a=3**, isto é, **três lacunas e dois pontos internos**. Pontos exteriores são arbitrários em número e espaçamento. Para essa classe não é necessário importar o Lema 5 ou o Lema 6 completo de Preissmann.

## 1. Redução ao orçamento interior

Normalize h=1. Sejam a,b,c>0, a+b+c=1, as três lacunas entre os extremos.

Defina

\[
Z=\zeta(2)=\pi^2/6,\quad
S(s)=\sum_{j\ge1}\frac1{j^2(1+js)^2},\quad
D(s)=\frac{2Z-S(s)}s-3s.
\]

P89 controla as caudas exteriores e P90 prova que D é decrescente. Basta demonstrar, **para todo a,b,c>0 com soma 1**:

\[
\boxed{\mathcal M(a,b,c):=
\frac{\min(a,b)}{a^2(1-a)^2}+
\frac{\min(b,c)}{c^2(1-c)^2}
\ \le D(a)+D(c).}
\tag{P91-GATE}
\]

De fato, as distâncias normalizadas dos extremos são no máximo a,c; logo a monotonicidade de D e as caudas de P89 implicam o Lema 6 afiado.

## 2. Minorante racional obtido por convexidade

De P90, e por frações parciais:
\[
S(0)=Z,\quad S(1/3)=2Z-31/12,\quad
S(1/2)=2Z-11/4,\quad S(1)=2Z-3.
\]
As cordas de S convexo entre nós sucessivos são majorantes. Como Z>3/2 e o coeficiente de Z no intervalo inicial é (1/s-3)>=0, segue:

\[
D(s)\ge d(s)=
\begin{cases}
3/(2s)+13/4-3s,&0<s\le1/3 \quad(d1),\\
9/(4s)+1-3s,&1/3\le s\le1/2 \quad(d2),\\
5/(2s)+1/2-3s,&1/2\le s\le1 \quad(d3).
\end{cases}
\]

As expressões concordam nas fronteiras. Em particular, **d(1/3)=D(1/3)=27/4**. Portanto, a desigualdade suficiente mais forte a provar é inteiramente racional: M(a,b,c) <= d(a)+d(c).

## 3. Certificado universal de nove triângulos

Por reflexão (a,b,c) ↔ (c,b,a), reduza ao domínio

\[
\Omega=\{(a,c):0\le a\le c,\ a+c\le 1\},
\]

com área exata 1/4. As três ordens possíveis são abc (a<=b<=c), acb (a<=c<=b), bac (b<=a<=c). Intersectando essas regiões com as faixas d1,d2,d3 e dividindo um quadrilátero em dois triângulos, resulta uma partição exata em **nove triângulos**:

| Triângulo | Ordem | Envelopes (a,c) | Vértices (a,c) |
|:--|:--|:--|:--|
| 1 | abc | d1,d2 | (1/3,1/3), (1/4,1/2), (0,1/2) |
| 2 | abc | d1,d3 | (0,1/2), (1/4,1/2), (0,1) |
| 3 | acb | d1,d1 | (0,0), (1/3,1/3), (0,1/3) |
| 4 | acb | d1,d2 | (0,1/3), (1/3,1/3), (0,1/2) |
| 5 | bac | d1,d2 | (1/3,1/3), (1/3,1/2), (1/4,1/2) |
| 6 | bac | d1,d3 | (1/4,1/2), (1/3,1/2), (1/3,2/3) |
| 7 | bac | d1,d3 | (1/4,1/2), (1/3,2/3), (0,1) |
| 8 | bac | d2,d2 | (1/3,1/3), (1/2,1/2), (1/3,1/2) |
| 9 | bac | d2,d3 | (1/3,1/2), (1/2,1/2), (1/3,2/3) |

Os interiores não se sobrepõem porque as condições de ordem e as faixas são disjuntas, exceto nas fronteiras; a soma exata das áreas é 1/4, cobrindo o domínio todo.

Num triângulo de vértices (a_i,c_i), tome u,v,w>=0 e H=u+v+w>0:

\[
A=\sum_i a_i u_i,\quad C=\sum_i c_i u_i,\quad
a=A/H,\quad c=C/H,\quad B=H-A-C,\quad L=H-A,\quad R=H-C.
\]

Escreva d(s)=alpha/s+beta-3s em cada faixa e defina

\[
P(a,c)=4a^2c^2(1-a)^2(1-c)^2
\,[d(a)+d(c)-\mathcal M(a,1-a-c,c)].
\]

O polinômio homogêneo exato de grau 9 H^9 P(A/H,C/H) é:

\[
\begin{aligned}
&4\alpha_a AC^2L^2R^2H^2 +4\alpha_c A^2CL^2R^2H^2\\
&+4(\beta_a+\beta_c)A^2C^2L^2R^2H
-12(A^3C^2+A^2C^3)L^2R^2\\
&-4C^2R^2\min(A,B)H^4-4A^2L^2\min(B,C)H^4.
\end{aligned}
\]

Em cada região a ordem fixa os mínimos como formas lineares. O arquivo [certificate.py](certificate.py) constrói essas expressões por multiplicações esparsas, usando **apenas fractions.Fraction**, e registra **cada coeficiente com expoentes**. São respectivamente 53, 42, 41, 44, 54, 55, 49, 54, 55 coeficientes não nulos, **todos positivos**: 447 no total. Os outros coeficientes são zero. Cada polinômio é não negativo para u,v,w>=0. O fator de P é positivo para lacunas estritamente positivas: assim M<=d(a)+d(c)<=D(a)+D(c) para **todos os reais admissíveis**.

A checagem programática de coeficientes é um certificado algébrico finito, não uma extrapolação de amostras. Os testes em pontos interiores apenas conferem se a expansão corresponde à fórmula original.

## 4. Reprodução e fronteira honesta

Executar a partir da raiz do repositório:

    python -m unittest -v NUMERIA_P91_THREE_GAPS.test_p91
    python -O -m unittest -v NUMERIA_P91_THREE_GAPS.test_p91
    python -m NUMERIA_P91_THREE_GAPS.certificate --out p91-three-gaps.json

O workflow P91 executa a auditoria, cria o JSON integral e o manifesto SHA-256. A igualdade de rede uniforme ocorre em a=b=c=1/3 e é preservada exatamente.

**Não demonstrados:** Lema 5 geral, Lema 6 quando b-a>=4, constante espectral geral sem Preissmann externo, Ramaré–Zúñiga no intervalo de X pendente e RH global. A próxima fronteira P92 é quatro lacunas internas (três pontos interiores) e exige nova redução variacional ou uma prova indutiva; não presumir extensão automática. Preservar PR em draft, sem merge automático.
