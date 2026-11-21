# NUMERIA–RIEMANN · Extensão escalar da P259

A partir da [P259](../NUMERIA_ES8P259/README.md), o certificado [JSON](SUZUKI_PREFIX_100000000.json) cobre todas as **5.762.859 potências primas** $q\le 10^8$ em Arb a 192 bits. São **5.684.125 prefixos adicionais**. A cota quadrática da seção seguinte também é estritamente positiva em cada prefixo. A parte comum até $10^6$ reproduz o SHA-256 dos limites diádicos da P259, entrada por entrada.

| Teste | Limite inferior uniforme (decimal arredondado para baixo) | Prefixo limitante da cota registrada |
| --- | ---: | ---: |
| Margem tangente $\Delta_q$, faixa toda | 2,69846 | 2 |
| Cota quadrática $Q_q$, faixa toda | 1,79065 | 2 |
| Margem tangente $\Delta_q$, prefixos novos | 40,84264 | 1.000.003 |
| Cota quadrática $Q_q$, prefixos novos | 40,83789 | 1.000.033 |
| Função completa $\Psi$, intervalo novo | 0,02149858 | 34.186.367 |

## Margem aritmética e cota quadrática

Sejam $A_X=\sum_{n\le X}\Lambda(n)/\sqrt n$, $B_X=\sum_{n\le X}\Lambda(n)\log n/\sqrt n$ e $P(t)=\sum_n \Lambda(n)(t-\log n)_+/\sqrt n$. Escreva $F(X)=4\sqrt X-P(\log X)$. Para $X\ge2$, a mínima distância entre a linha do prefixo $A_Xt-B_X$ e $4e^{t/2}$ é

$$
\Delta_X=B_X-2A_X\bigl(\log(A_X/2)-1\bigr),\qquad
t_* = 2\log(A_X/2).
$$

A minimização vale para todo $t\in\mathbb R$. A inclinação $A_X$ é positiva e a segunda derivada de $4e^{t/2}$ é positiva. O ponto crítico $t_*$ é assim o mínimo global. Em $t\ge0$, a linha ativa para $P(t)$ é a do último $q\le e^t$, com a linha vazia antes de $\log 2$. Portanto o certificado implica $P(t)<4e^{t/2}$ em todo $0\le t\le\log(10^8)$, inclusive entre potências primas.

Com $u=A_X/(2\sqrt X)$ e $h(u)=u\log u-u+1$, a substituição dá a identidade exata

$$
\Delta_X=F(X)-4\sqrt X\,h(u),\qquad
0\le h(u)\le (u-1)^2.
$$

A primeira desigualdade segue da convexidade e a segunda de $\log u\le u-1$. Logo

$$
Q_X:=F(X)-\frac{(A_X-2\sqrt X)^2}{\sqrt X}\le\Delta_X.
$$

O teste de $Q_q$ usa a soma ponderada, uma raiz quadrada e o desvio do primeiro momento. Ele evita o logaritmo do primeiro momento, que aparece no teste de $\Delta_q$. As duas expressões foram avaliadas separadamente. A desigualdade $Q_X\le\Delta_X$ é analítica e a positividade de $Q_q$ vem dos limites inferiores calculados em Arb.

## Função completa de Suzuki na faixa nova

Na fórmula explícita de Suzuki, o termo de índice zero da série de Lerch cancela exatamente $4e^{-t/2}$. Para $t>0$, sobram

$$
\begin{aligned}
\Psi(t)&=4e^{t/2}-ct+C-P(t)-\varepsilon(t),\\
c&=\tfrac12(\log\pi-\Gamma'(1/4)/\Gamma(1/4)),&
C&=\tfrac14\zeta(2,1/4)-8,\\
\varepsilon(t)&=4\sum_{k\ge1}\frac{e^{-(4k+1)t/2}}{(4k+1)^2}.
\end{aligned}
$$

Para $t\ge\log S$, com $S=10^6$, a cauda satisfaz $0<\varepsilon(t)\le4/[25(S^2-1)\sqrt S]$, pela série geométrica de razão $S^{-2}$ com denominador mínimo $5^2$. O certificador substitui $A_q,B_q$ por $A_q+c,B_q+C$ na minimização anterior e deduz essa cota de cauda, menor que $1{,}601\cdot10^{-16}$. Isso dá um único limite inferior para toda a célula em que o prefixo é ativo. Inclui-se o último prefixo até $S$ para cobrir o ponto de junção e o trecho até a primeira potência prima nova. Os **5.684.126 prefixos ativos** cobrem todo o intervalo novo. Em conjunto com a parte inicial da P259, obtém-se **$\Psi(t)>0$ para $0<t\le\log(10^8)$**, e $\Psi(0)=0$. [Suzuki 2023, fórmula (1.1)](https://doi.org/10.1112/jlms.12785) · [Certificado da parte inicial](../NUMERIA_ES8P259/data/SUZUKI_FINITE.json)

## O que falta para a cauda infinita

Defina $e(X)=A_X-2\sqrt X$. Integrar a soma finita termo a termo dá $P(\log X)=\int_1^X A_y\,dy/y$. A cota quadrática vira

$$
Q_X=4-\int_1^X e(y)\,\frac{dy}{y}-\frac{e(X)^2}{\sqrt X}.
$$

Assim, seria suficiente provar nos prefixos grandes uma cota unilateral $\int_1^q e(y)\,dy/y\le 4-\kappa\log q$ com alguma constante $\kappa>0$, junto com $e(q)=o(q^{1/4}\sqrt{\log q})$. A expressão exata de $Q_q$ permite outras estimativas conjuntas. O teorema dos números primos, por soma parcial, dá $e(X)=o(\sqrt X)$ e $F(X)=o(\sqrt X)$. Essas escalas ainda deixam indeterminado o sinal da diferença e da média ponderada.

O sinal eventual de $Q_q$ oferece uma formulação equivalente de RH. Se $Q_q>0$ para todos os prefixos suficientemente grandes, a minimização acima implica $P(\log X)\le4\sqrt X$ para todo $X$ suficientemente grande. O Teorema 1 de Suzuki então implica RH. Na direção inversa, **sob RH**, a fórmula explícita (46) de Suzuki dá $F(X)=c\log X+O(1)$, com $c>0$. A estimativa clássica sob RH $\sum_{n\le X}\Lambda(n)=X+O(\sqrt X\log^2 X)$ dá por soma parcial $e(X)=O(\log^3 X)$. Assim o termo quadrático é $O(\log^6 X/\sqrt X)=o(1)$, e $Q_q/\log q\to c$. A implicação inversa é condicional e serve para dimensionar a margem esperada. A certificação finita foi obtida incondicionalmente. [Suzuki 2025, Teorema 1 e equações (46)–(47)](https://doi.org/10.1007/s11139-025-01238-9)

## Comparação bibliográfica

Mittermeier relata aritmética MPFR dirigida para a função completa de Suzuki até $X=10^{10}$ e estuda a cauda por uma identidade de balanço. Sua continuação pelo lado dos zeros e o preprint de Pedro Caceres relatam faixas finitas da ordem de $2\cdot10^{29}$. Esses alcances externos foram consultados nas fontes originais. A reprodução deste repositório usa primos até $10^8$ e inclui o teste da cota quadrática $Q_q$ em cada prefixo. [Mittermeier, verificação por células](https://zenodo.org/records/23377174) · [Mittermeier, balanço](https://zenodo.org/records/23377241) · [Mittermeier, verificação pelo lado dos zeros](https://zenodo.org/records/23544568) · [Caceres, preprint Hexagon](https://hexagonmath.org/2610.00661)

## Execução e auditoria

Da raiz do repositório, em Python 3.12

```sh
python -m pip install -r NUMERIA_ES8P249/requirements.txt
python -m unittest NUMERIA_SUZUKI_PREFIX.test_prefix NUMERIA_ES8P259.test_p259 -v
python -O -m unittest NUMERIA_SUZUKI_PREFIX.test_prefix NUMERIA_ES8P259.test_p259 -v
python -O -m NUMERIA_SUZUKI_PREFIX.audit_sieve
sha256sum -c NUMERIA_SUZUKI_PREFIX/SHA256SUMS
python -O -m NUMERIA_SUZUKI_PREFIX.certify --block-size 65521 \
    --verify NUMERIA_SUZUKI_PREFIX/SUZUKI_PREFIX_100000000.json
```

Para gerar outro certificado, use `python -m NUMERIA_SUZUKI_PREFIX.certify --limit 100000000 --out /tmp/SUZUKI_PREFIX.json`. Para reproduzir a parte inicial, use `python -m NUMERIA_ES8P259.suzuki --limit 1000000 --out /tmp/SUZUKI_FINITE.json`.

O crivo segmentado conserva um bloco, os primos até $\sqrt X$ e as potências de expoente pelo menos dois. O auditor usa um crivo independente para ímpares e divisão por tentativa para a base das potências próprias, e compara todos os pares $(q,p)$. Pequenos prefixos têm ainda uma checagem independente com frações e séries para logaritmos. A avaliação da cota da função completa é comparada com a avaliação direta da série de Lerch usada pela P259. Os limites são diádicos exatos `[m, e]`, com valor $m2^e$. O SHA-256 inclui as linhas da margem tangente, com JSON compacto `[q,[m,e]]`, um `\n` por linha e `m` como texto decimal. O replay recompõe todas as somas e os limites. O escopo de todas as verificações computacionais é a faixa registrada no certificado.
