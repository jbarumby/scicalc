# P92-F.2 — Critério espectral refinado para qualquer quantidade de escalas

Esta nota é uma extensão da [prova P92-F.1](README.md). A cota `||Delta|| <= (m-1)*zeta(2)` descartava as razões reais entre as escalas; retendo-as, obtém-se um teste mais eficiente **mesmo com m>3**.

## Teorema analítico universal baseado no conjunto dos valores de delta

Seja S o conjunto FINITO de valores positivos assumidos pelas distâncias `delta_i=min_{j!=i}|x_i-x_j|` de uma configuração arbitrária de pontos distintos reais. Defina, para cada `a in S`,

\[
F_S(a)=\sum_{b\in S,\ b<a}\big(1-\sqrt{b/a}\big)^2\big(1+b/a\big)
+\zeta(2)\sum_{b\in S,\ b>a}\big(1-\sqrt{a/b}\big)^2.
\]

Então a prova dos casos `b<a` e `b>a` da P92-F.1, **antes de remover os fatores dependentes das razões**, dá

\[
\boxed{\|\Delta\|_2\le\max_{a\in S}F_S(a).}\tag{F8}
\]

A afirmação é válida para *toda* configuração cujos `delta_i` pertencem a S, **independentemente da ordem, da quantidade e da posição dos pontos**, e também vale quando só um subconjunto de S ocorre (omitindo termos positivos só diminui a soma, e o máximo pode ser tomado sobre S). É invariável pela dilatação `S -> cS` para qualquer c>0. É uma prova analítica que não depende de cálculo numérico nem de forma em blocos.

## Certificado racional para cinco valores distintos de delta

O código [spectrum_gate.py](spectrum_gate.py) usa `Fraction` para espectros `S={s_i^2}` com raízes `s_i` racionais e usa a cota **superior** 

\[
\zeta(2)=\pi^2/6<5/3 \quad (\text{pois }\pi<22/7\text{ implica }\pi^2<10).
\]

Considere, para qualquer `c>0`, o conjunto

\[
S=c\{1,4,9,16,25\}.
\]

O máximo exato da estimativa racional de F8, verificado por cinco linhas de aritmética de frações, é:

\[
\boxed{\|\Delta\|_2\le\frac{6829}{2160}\approx3.16157.}
\]

Esse valor satisfaz o **limite residual da P88** de maneira rigorosa, por cotas racionais inferiores `pi>25/8` e `sqrt(6/5)>109/100`:

\[
\frac{6829}{2160}
<\frac{(25/8)^2}{3}\,\frac{109}{100}
<\frac{\pi^2}{3}\sqrt{6/5}.
\]

Logo a estimativa residual P88 vale para **qualquer quantidade de pontos e qualquer arranjo** em que os `delta_i` pertençam a um múltiplo comum desse conjunto de **cinco** escalas. Por exemplo, folgas consecutivas `1,4,4,9,9,16,16,25,25` realizam os cinco valores de delta, mas a prova NÃO depende desse exemplo.

O mesmo verificador encontra certificados para **seis** escalas próximas, como `S={100,121,144,169,196,225}` (ou seus múltiplos positivos), novamente sem limitar a quantidade de pontos. Ele também marca conjuntos de escalas extremamente separadas como **inconclusivos pelo teste**, o que não é um contraexemplo ao limite procurado.

## Uma condição simples que também permite m>3

Suponha que `h=min S`, `H=max S`. Ponha `epsilon=1-sqrt(h/H)` em [0,1). Em cada parcela do somatório F8, `|1-sqrt(min(a,b)/max(a,b))|<=epsilon`; e `1+b/a<=2` no caso b<a. Como `zeta(2)<2`, todas as parcelas são <=`2 epsilon^2`. Assim

\[
\boxed{\|\Delta\|_2\le2(m-1)\left(1-\sqrt{h/H}\right)^2.}
\tag{F9}
\]

Essa estimativa é independente do N e do arranjo, e permite **m arbitrário**, desde que o contraste máximo de escalas seja suficientemente pequeno em função de m. É grosseira, mas analítica; para certificados mais fortes, preferir F8 com as razões individuais.

## Reproduzir e auditar

    python -m unittest -v NUMERIA_P92F_THREE_SCALES.test_spectrum
    python -O -m unittest -v NUMERIA_P92F_THREE_SCALES.test_spectrum
    python -m NUMERIA_P92F_THREE_SCALES.spectrum_gate --out p92f-spectrum.json

O JSON registra as parcelas exatas por escala, o máximo racional, a comparação rigorosa contra o alvo P88 e campos explícitos que NÃO atribuem validade geral à estimativa residual para todas as configurações.

**Limites:** F8, F9 não fornecem uma constante residual uniforme quando o número de valores distintos e seu contraste crescem sem controle. Lema 5/6 completos, caso geral de Preissmann, P93 geométrica e RH continuam em aberto. A dedução da norma de B usa o elo de energia da P87 com a correção analítica já documentada, além da identidade Rayleigh com sinal de comutador reparado na P88.
