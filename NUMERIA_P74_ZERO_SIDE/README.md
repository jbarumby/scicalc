# NUMERIA–RIEMANN · P74 · Certificação pelo lado dos zeros

A P74 parte do [PR #73](https://github.com/jbarumby/scicalc/pull/73), commit `e4c5cf326f8d7e05795ab5457a37deee6beb7a72`. O cálculo usa a identidade incondicional de Suzuki, a verificação finita publicada de Platt–Trudgian e cotas publicadas de contagem e densidade. As conclusões desta etapa dizem respeito a intervalos finitos. [Suzuki](https://doi.org/10.1112/jlms.12785) · [Platt–Trudgian](https://doi.org/10.1112/blms.12460) · [Bellotti–Wong](https://doi.org/10.1090/mcom/4133) · [Kadiri–Lumley–Ng](https://doi.org/10.1016/j.jmaa.2018.04.071)

## Auditoria da P73

O cálculo racional até $t=46{,}6$ reproduziu os bytes publicados a 192 bits. São **565 células**, 566 avaliações de extremos, profundidade 12 e limite inferior estrito de **$4{,}1984667188770\cdot10^{-5}$**. A 256 bits, as 20 caixas de Arb ficam estritamente contidas nas originais e a mesma partição conserva a positividade. Uma avaliação separada a 320 bits verificou sinais opostos de Hardy $Z$ nos extremos diádicos das 20 caixas originais. A contagem exata $N(78)=20$ e a disjunção das caixas confirmam a completude e a simplicidade nessa faixa. [Auditoria JSON](p73_replay.json) · [README da P73](../NUMERIA_P73_ADAPTIVE/README.md) · [Plano da P73](../NUMERIA_P73_ADAPTIVE/CODEX_NEXT_STEPS.md) · [Isolamento de zeros no FLINT](https://flintlib.org/doc/acb_dirichlet/zeta-zeros.html)

## Experimentos em Arb

A cota por faixas certificou sucessivamente **$t=47{,}48,\ 50,\ 55,\ 60$** com 20, 50, 100, 200 e 500 zeros, tanto a 192 como a 256 bits. Cada execução recompõe as caixas FLINT, confere a ordem e a parte real, verifica sinais de Hardy $Z$ nos extremos ao dobro da precisão e conta todos os zeros até um corte que separa o último retido do próximo. O JSON registra o SHA-256 das caixas e os limites diádicos. Os tempos são racionais exatos, incluindo $47{,}48=1187/25$. [Certificados](experiments.json)

| Zeros | Células adaptativas até 60, em ambas as precisões | Limite inferior em ambas as precisões, arredondado para baixo | Cota sem faixas até 47,48, a 192 bits |
| ---: | ---: | ---: | --- |
| 20 | 2.690 | 0,00697180 | Inconclusiva |
| 50 | 2.690 | 0,01033589 | Inconclusiva |
| 100 | 2.690 | 0,01022401 | Inconclusiva |
| 200 | 2.690 | 0,00635827 | Inconclusiva |
| 500 | 2.841 | 0,00000117533 | Positiva, com cota 0,0000961132 |

Com 500 zeros e a grade uniforme de passo $1/200$, a mesma janela $[18,60]$ fecha em **8.400 células**, com cota **0,01609245**. A grade adaptativa usa 2.841 células e aceita uma margem menor. As margens medem a qualidade da cota de certificação. Adicionar zeros aumenta a soma pontual, enquanto aumenta também a perda de curvatura estimada para uma célula de largura fixa.

Uma rodada adicional com **500 zeros a 256 bits** cobriu $[18,67{,}5]$ em **3.457 células**, com limite inferior **$1{,}17533\cdot10^{-6}$**. Isso reproduz o extremo temporal relatado por Mittermeier com uma implementação independente e um prefixo menor. Como Arb verifica $e^{18}<10^8$, o certificado se junta à P71, que cobre $(0,\log(10^8)]$. A união certifica $\Psi(t)>0$ para $0<t\le67{,}5$ e a fórmula dá $\Psi(0)=0$. [P71](../NUMERIA_SUZUKI_PREFIX/README.md) · [Mittermeier](https://zenodo.org/records/23544568)

Em $t=50$, a cota sem faixas é maior que $0{,}1065459$. A soma positiva retida, mesmo incluindo todos os zeros da reta, é menor que $0{,}092533$ pelo limite $4\bigl(\sum_{j=1}^{500}\gamma_j^{-2}+\overline W(812)\bigr)$. A contagem é $N(812)=500$. Assim, refinar a malha ou aumentar apenas o número de zeros é insuficiente com essa cota. A diferença negativa diz respeito ao limite inferior escolhido e deixa o sinal de $\Psi$ em aberto por esse método. [Certificado da obstrução](experiments.json)

## Fórmula e cotas

Escreva $z=(\rho-1/2)/i$, com multiplicidade. A soma de Suzuki é absolutamente convergente para cada $t$ e utiliza $(1-\cos(zt))/z^2$. Na reta crítica, os pares de ordenadas $\pm\gamma$ dão

$$
S_n(t)=2\sum_{j=1}^n\frac{1-\cos(\gamma_jt)}{\gamma_j^2}
      =4\sum_{j=1}^n\frac{\sin^2(\gamma_jt/2)}{\gamma_j^2}.
$$

As contribuições omitidas dos zeros sobre a reta são não negativas. Pelo Teorema 1 de Platt–Trudgian, os zeros potencialmente fora da reta têm $|\gamma|>H=3\cdot10^{12}$. Nesse domínio, $\Re(1/z^2)>0$, $|z|^2\ge\gamma^2$ e $|\cos(zt)|\le\cosh(|\Re\rho-1/2|t)$ para $t\ge0$. Logo cada termo fora da reta perde no máximo $\gamma^{-2}\cosh(|\Re\rho-1/2|t)$. A contagem de ordenadas positivas é multiplicada por dois para cobrir os sinais de $\gamma$. [Suzuki, equação (1.3)](https://arxiv.org/html/2206.03682v4) · [Platt–Trudgian, Teorema 1](https://arxiv.org/abs/2004.09765)

Para o peso $W(H)=\sum_{\gamma>H}\gamma^{-2}$, escreva $N(u)=M(u)+R(u)$, com $M(u)=u\log(u/(2\pi e))/(2\pi)$. Bellotti–Wong dá $|R(u)|\le E(u)=a\log u+b\log\log u+c$ para $u\ge e$, com os decimais publicados $(a,b,c)=(0.10076,0.24460,8.08344)$. A soma por partes conserva o termo de fronteira $-N(H)/H^2$. Como $E(u)\le E(H)+(a+b/\log H)\log(u/H)$,

$$
W(H)\le\overline W(H)=\frac{\log(H/(2\pi))+1}{2\pi H}
      +\frac{2E(H)+(a+b/\log H)/2}{H^2}.
$$

A cota sem faixas é $2\overline W(H)\cosh(t/2)$. Melhorar apenas o erro de contagem altera a correção de ordem $\log H/H^2$. A auditoria também verificou a cota mais larga da P73, $5/H+38/H^2$ para $W(H)$, e seu fator $3+e^{t/2}\ge2(1+\cosh(t/2))$. [Bellotti–Wong](https://doi.org/10.1090/mcom/4133) · [P73](../NUMERIA_P73_ADAPTIVE/README.md)

Para a cota por faixas, usamos 16 limiares da Tabela 1 de Kadiri–Lumley–Ng, com parâmetro $k=1$ e altura $H$ superior ao mínimo exigido de $3{,}0610046\cdot10^{10}$. Para cada $\sigma_i$, integramos a cota de $N(\sigma_i,u)$ contra $2u^{-3}\,du$ em $[H,\infty)$ para limitar $D_i=\sum_{\gamma>H,\,\Re\rho>\sigma_i}\gamma^{-2}$. A gama incompleta avalia a integral da potência de $\log u$ em Arb. Também vale $D_i\le W(H)/2$ pela reflexão $\beta\leftrightarrow1-\beta$. Tomamos limites superiores decrescentes $\overline D_i$ e escrevemos $\theta_i=\sigma_i-1/2$, com último extremo $\theta_{17}=1/2$. Para $t\ge0$ uma cota é

$$
U(t)=2\overline W(H)\cosh(\theta_1t)
 +4\sum_{i=1}^{16}\overline D_i
       \bigl(\cosh(\theta_{i+1}t)-\cosh(\theta_it)\bigr).
$$

Os incrementos são positivos. Quando a parte real coincide com um limiar estrito, os incrementos anteriores já cobrem o peso nesse ponto. O fator quatro combina a simetria da ordenada com a reflexão em torno da reta crítica. A implementação toma o menor limite entre essa cota e a cota sem faixas. Ambas são crescentes em $t$. [Kadiri–Lumley–Ng, Teorema 1.1 e Tabela 1](https://www.cs.uleth.ca/~kadiri/articles/Explicit-zero-density-for-the-Riemann-zeta-function-JMAA-May2018.pdf) · [Mittermeier, método por faixas](https://zenodo.org/records/23544568)

Em uma célula $[a,b]$, $S_n''\le 2n$. Se $\ell_a$ e $\ell_b$ são limites inferiores nos extremos, a parábola

$$
p(u)=(1-u)\ell_a+u\ell_b-n(b-a)^2u(1-u),\qquad
u=(t-a)/(b-a)\in[0,1]
$$

limita $S_n(t)$ por baixo. O verificador minimiza essa parábola e subtrai a cota da cauda em $b$. Com $C=n(b-a)^2$ e $d=|\ell_b-\ell_a|$, a perda além do menor extremo é $(C-d)^2/(4C)$ se $d<C$, e zero se $d\ge C$. Esse limite recupera a perda uniforme $C/4$ quando os extremos têm a mesma altura. Arb limita todas as avaliações e decisões de subdivisão. O percurso determinístico avança da esquerda para a direita, com extremos racionais exatos.

## Busca com orçamento fixo

`discovery.py` compara a malha uniforme, a subdivisão adaptativa e três sequências pseudoaleatórias com sementes 73, 74 e 101. Em $[46,48]$, usa o indicador de 20 zeros com a cauda sem faixas. Com 512 avaliações, a largura entre as amostras que cercam a primeira transição de sinal observada foi $4{,}88\cdot10^{-4}$ na adaptativa, $3{,}91\cdot10^{-3}$ na uniforme e $2{,}73\cdot10^{-3}$ a $7{,}92\cdot10^{-3}$ nas sementes aleatórias. Com 256 avaliações, as larguras foram respectivamente $9{,}77\cdot10^{-4}$, $7{,}84\cdot10^{-3}$ e $8{,}03\cdot10^{-3}$ a $1{,}37\cdot10^{-2}$. São diagnósticos do indicador amostrado. A prova é executada em Arb sem usar esses pontos. [Dados da busca](discovery.json)

Keating–Snaith calcula momentos dos polinômios característicos de matrizes CUE. As comparações de Odlyzko com espaçamentos GUE fornecem evidência empírica. Esses modelos podem sugerir prioridades de exploração. A garantia local depende das fases dos zeros isolados e das cotas da cauda. O ganho de um critério de matrizes aleatórias diante dos três controles medidos fica em aberto. [Keating–Snaith](https://doi.org/10.1007/s002200000261) · [Odlyzko](https://doi.org/10.1090/S0025-5718-1987-0866115-0)

## Limitação e próxima etapa

Para um prefixo fixo, a aproximação simultânea das razões $\gamma_j/\gamma_1$ faz a soma retida retornar arbitrariamente perto de zero em tempos crescentes. A cota absoluta por densidade cresce exponencialmente em $t$. Uma etapa analítica seguinte é preservar informação sobre as fases na contribuição acima de $H$, sobretudo perto das quase recorrências. O Teorema 1.7 de Suzuki exige uma condição global de sinal que a certificação finita deixa em aberto. [Suzuki](https://doi.org/10.1112/jlms.12785)

## Reprodução

A partir da raiz do repositório, em Python 3.12

```sh
python -m pip install python-flint==0.9.0 mpmath==1.3.0
sha256sum -c NUMERIA_P74_ZERO_SIDE/SHA256SUMS
python -m NUMERIA_P74_ZERO_SIDE.audit_p73 --out /tmp/p73-replay.json
cmp /tmp/p73-replay.json NUMERIA_P74_ZERO_SIDE/p73_replay.json
python -m NUMERIA_P74_ZERO_SIDE.experiments --verify NUMERIA_P74_ZERO_SIDE/experiments.json
python -O -m NUMERIA_P74_ZERO_SIDE.experiments --verify NUMERIA_P74_ZERO_SIDE/experiments.json
python -m unittest NUMERIA_P74_ZERO_SIDE.test_p74 -v
python -O -m unittest NUMERIA_P74_ZERO_SIDE.test_p74 -v
python -m NUMERIA_P74_ZERO_SIDE.discovery --out /tmp/p74-discovery.json
```

Os limites diádicos são pares `[m,e]` com valor $m2^e$. O SHA-256 da cobertura inclui uma linha JSON compacta por célula, com extremo esquerdo, extremo direito e cota inferior. O verificador refaz as avaliações direcionadas e compara todos os campos. Um valor negativo em `endpoint_envelope_upper_dyadic` mostra onde o limite inferior escolhido perde a capacidade de certificar o sinal.
