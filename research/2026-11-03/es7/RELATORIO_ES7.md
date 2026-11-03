# NUMERIA–RIEMANN — ES7: certificação da parcela arquimediana gaussiana

**Data:** 03/11/2026  
**Status:** prova computacional de *enclosures* racionais sob lemas analíticos clássicos explícitos. Sem Lean; não é prova de RH nem positividade universal de Weil.

## 1. Continuação de ES6

A ES6 implementou, para `h_a(t)=exp(-a t²)`, a parcela arquimediana da fórmula explícita de Guinand–Weil, mas avaliou a integral por `mpmath` sem certificação da quadratura. Nesta etapa certificamos, **somente**, para `a>0`,

\[
A(a)=\frac{1}{\pi}\int_0^\infty e^{-a t^2}\,\Omega(t)\,dt,\qquad
\Omega(t)=\operatorname{Re}\psi\left(\frac14+\frac{it}{2}\right)-\log\pi,
\]

além da derivada real

\[
A'(a)=-\frac{1}{\pi}\int_0^\infty t^2e^{-a t^2}\Omega(t)\,dt.
\]

A fórmula global da ES6 tem fator `1/(2π)` na integral em toda a reta; por paridade, essa convenção coincide com `1/π` na semirreta. **Isso não prova que as entradas de matrizes V24 usem a mesma função-teste.**

## 2. Série positiva de digama e valor em zero

Pela expansão em frações parciais da digama (DLMF 5.7.6),

\[
\Omega(t)=\Omega(0)+\sum_{n\ge0}b_n(t),\quad
b_n(t)=\frac{t^2}{q_n(4q_n^2+t^2)},\quad q_n=n+\tfrac14,
\]

\[
\Omega(0)=-\gamma-\frac\pi2-3\log2-\log\pi.
\]

Cada `b_n(t)` é não negativo, cresce com `t>=0` e diminui com `n`. Após `m` termos, o teste integral fornece

\[
\frac12\log\left(1+\frac{t^2}{4(m+1/4)^2}\right)
\le\sum_{n=m}^\infty b_n(t)
\le\frac12\log\left(1+\frac{t^2}{4(m-3/4)^2}\right).
\]

Esses logaritmos são cercados por séries racionais de `atanh` com resto geométrico explícito. `π` é cercado pela fórmula de Machin, `log 2` e `log π` pela série racional e `γ` por um intervalo de Young aplicado a `H_4096 - log 4096`:

\[
\frac{1}{2(4097)} < H_{4096}-\log4096-\gamma <\frac{1}{8192}.
\]

Os arredondamentos dirigidos operam em inteiros com escala `10^28`; os resultados exatos estão preservados em `RESULTADOS_ES7.json`, não apenas nas apresentações decimais.

## 3. Erro de integração e diferenciação

Da série acima, usando `\sum_{n\ge0}q_n^{-3}<=64+8=72`, obtemos as cotas globais

\[
|\Omega(t)|\le 7+18t^2,\qquad
|\Omega'(t)|\le36t,\qquad
|\Omega''(t)|\le36.
\]

O último majorante usa `|1-3u|/(1+u)^3<=1` para `u>=0` em cada segunda derivada `b''_n(t)`. Para `f_a(t)=e^{-at^2}\Omega(t)` e `y=at^2`,

\[
|f_a''(t)|\le e^{-y}\{36+180y+72y^2+28ay+14a\}.
\]

A quadratura de ponto médio composta em `[0,T]` é cercada pelo teorema de erro, com majorantes *por subintervalo* e soma dos erros em aritmética racional dirigida. Para a derivada paramétrica, o integrando é `-t² f_a(t)` e usamos

\[
|(t^2f_a)''|\le2|f_a|+4t|f_a'|+t^2|f_a''|,
\qquad |f_a'|\le e^{-at^2}\{36t+2at(7+18t^2)\}.
\]

Com `I_0<=exp(-aT²)/(2aT)` e integração por partes, a cauda é majorada por

\[
\frac{e^{-aT^2}}{\pi}\left[\frac7{2aT}
+18\left(\frac{T}{2a}+\frac{1}{4a^2T}\right)\right].
\]

Para a derivada, aplicamos as correspondentes cotas de `I_2` e `I_4`. Em ambos os casos o programa usa aritmética racional para os limites e para o exponencial negativo (Taylor alternado em intervalo reduzido e quadraturas sucessivas).

## 4. Resultados certificados

Parâmetros `T=48`, `M=8192` subintervalos, 192 parcelas de digama e 28 casas decimais de escala fixa.

| a | Certificado para A(a) | Certificado para A'(a) |
|:--|:--|:--|
| 1/50 = 0,02 | **[-1,968639397083; -1,967605299413]** | **[-6,727200145263; -6,587641747029]** |
| 1/25 = 0,04 | **[-1,950194337532; -1,949460944506]** | **[3,798673129248; 3,848183534206]** |

Os extremos decimais foram arredondados **para fora**, e o JSON fornece os limites inteiros originais. São intervalos relativamente largos porque os majorantes de segunda derivada foram deliberadamente conservadores. A obtenção de um intervalo mais estreito exige refinar a quadratura ou melhorar os majorantes, não simplesmente aumentar dígitos de ponto flutuante.

O cálculo independente por `mpmath` da ES6, **não certificado**, fornece aproximadamente:

| a | A(a) numérico | A'(a) numérico |
|:--|:--|:--|
| 0,02 | -1,96812233304665481 | -6,65742109602382 |
| 0,04 | -1,94982762864192756 | 3,82342821773342 |

Todos estão contidos nos respectivos intervalos. Essa concordância é um teste independente da implementação, **não** a origem dos certificados racionais.

## 5. Cotejo com V24: somente parcialmente aprovado

A V24 usa formalmente

\[
A_{ij}^{V24}=\frac1\pi\int_0^\infty \Omega(r)F_i(r)\overline{F_j(r)}\,dr.
\]

A ES7 usa o mesmo kernel `Ω` e normalização `1/π`, mas o peso escalar **gaussiano** `e^{-a r²}`. Já a V24 utiliza produtos de transformadas de Fourier de modos originalmente compactamente suportados, com suavização especial e posterior projeção matricial.

**Gate: `NOT_ESTABLISHED`.** O fato de o kernel ser igual não autoriza transportar intervalos de ES7 a `A_ij`, aos autovalores da matriz projetada, à estabilidade LDL ou à positividade de Weil. Também não certificamos a cauda de zeros nem o acoplamento de termos de polos e primos em entradas individuais. `es7_v24_gate.py` confere apenas âncoras documentais e guarda hashes dos arquivos examinados, **não** constitui prova de equivalência.

## 6. Próximos gates analíticos

A ES8 deverá trocar o peso gaussiano pelo produto `F_i\overline{F_j}` das funções suavizadas, demonstrar a fórmula de traço compatível com o teste admissível, verificar fatores aritméticos e polos contra fonte independente, transportar limites pela projeção `B* (G+P) B`, e certificar autovalores ou pivôs **sem supor RH**. Só então se poderá discutir a positividade das matrizes finitas; mesmo ela não implica positividade para todos os testes.

## Referências

- NIST DLMF, §5.7.6, série de digama: https://dlmf.nist.gov/5.7.E6
- NIST DLMF, §5.4 e §5.9, valores especiais e representação integral de digama: https://dlmf.nist.gov/5.4 e https://dlmf.nist.gov/5.9
- Young/Cesàro para limites da constante de Euler, conforme: https://link.springer.com/article/10.1186/1029-242X-2013-222
- Carneiro–Chirre–Milinovich, fórmula explícita de Guinand–Weil (convenção da ES6): https://doi.org/10.1007/s00209-021-02820-9
- V24: https://github.com/jbarumby/scicalc/blob/main/v24/README_V24.md
- ES6: https://github.com/jbarumby/scicalc/commit/5c059bbaf1cac7a37e7be5fcd166020efca83957

**Conclusão:** parcela arquimediana e sua derivada agora possuem enclosures numéricos auditáveis; a equivalência matricial e a certificação global permanecem abertas. RH segue sem demonstração.
