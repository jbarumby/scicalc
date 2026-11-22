# P92 — Auditoria dos uploads P88

**Escopo:** auditoria dos seis novos arquivos P88 fornecidos em 22/11/2026. Esta nota NÃO amplia a prova P92 nem comprova o Lema 5/6 geral.

## Arquivos íntegros
Todos os cinco hashes anunciados no manifesto SHA256SUMS conferem:
- p88_verificacao.py: 890779c7ddf80ad18a7bbc5e9a959a4e334e859a9f9c33d809b1aae479077680
- p88_execucao.log: 42917b2efec911615893505e9541ee570e3ec93fae0c893d1841cf88114f2410
- p88_ledger.json: 72d292248927fb55c6d93e385a9b1c72aeb2cee97dcf8a8463e91b6605cdf0df
- p88_nota_demonstracao.tex: 446252ab24dd3005ccd97b7dc25b413ec2559faf9b1ca79d6f4c1f98edae39b0
- p88_nota_demonstracao.pdf: 492791d581f17aede0ce5fab69b190f3988dfe09aa1cbccf084d523ca3f9dcc6

O ledger registra 124 gates, zero falhas, mas marca expressamente L3_cota_Delta_analitica=false.

## 1. ERRO DE SINAL NO TEOREMA P88.1

Definições originais:
\[
B_{ij}=\sqrt{\delta_i\delta_j}/(x_i-x_j),\quad
T_i=\sum_{k\ne i}\delta_k/(x_i-x_k),\quad
R_{ij}=\sqrt{\delta_i\delta_j}(T_i-T_j)/(x_i-x_j),
\]
\[
C_{ij}=\sqrt{\delta_i\delta_j}(\delta_i+\delta_j)/[2(x_i-x_j)^2].
\]

A identidade matricial correta é
\[
\boxed{-B^2=\operatorname{diag}(Q)-R+2C.}
\]
A nota e o comentário de código escrevem +R, portanto a identidade para VETORES ARBITRÁRIOS está errada.

**Contraexemplo exato:** pontos x=(0,1,2), delta=(1,1,1): entrada (0,1) tem (-B²)01=1/2, R01=3/2, C01=1. Vale -R01+2C01=1/2; com +R01 resulta 7/2, errado.

O cancelamento do comutador R=K Lambda-Lambda K, K=-B, Lambda=diag(T), continua válido no autovetor u de B, pois <u,Ru>=0. Logo o RESULTADO de Rayleigh
\[
\|B\|_2^2=\sum_iQ_i|u_i|^2+2\langle u,Cu\rangle
\le\max_iQ_i+2\|C\|_2
\]
continua correto com o sinal reparado. Isso **NÃO PROVA** a alegação operadorial lambda_max(E)<=2||C|| para todo vetor: a anulação de <u,Ru> é específica dos autovetores de B.

## 2. T2 útil, residual Delta ainda aberto

A identidade C=D+Delta é válida, com D_ij=delta_i delta_j/(xi-xj)^2, Delta_ij=sqrt(delta_i delta_j)(sqrt(delta_i)-sqrt(delta_j))²/[2(xi-xj)²]. Matrizes simétricas de entradas não negativas dão
\[
\|C\|_2\le\|D\|_2+\|\Delta\|_2
\le\max_iQ_i+\|\Delta\|_2\le\pi^2/3+\|\Delta\|_2.
\]
Mas a cota universal ||Delta||<=(pi²/3)sqrt(6/5) NÃO foi provada analiticamente. B4 usa autovalores float64 em amostras finitas, e B6 transforma um float por arredondamento em uma fração sem demonstrar que tal fração majoriza o autovalor verdadeiro. A comparação B6 também usa uma cota SUPERIOR para o alvo. Não promover a certificação de exemplos finitos a teorema universal.

## 3. Correções adicionais verificadas

- B3 usa uma cota superior de pi. Recalculei as três famílias com Fraction e pi_lower²/3: passam, com margens decimais indicativas 0.13343201125092316 (uniforme_30), 0.917870263695595 (diadico), 1.3116727098347754 (geometrico_15). O script original requer alteração para expor as margens racionais exatas.
- A função families_float retorna **22** famílias predefinidas, não 25; existem mais 60 buscas aleatórias.
- O campo dobrada_gamma05_r70 do ledger ficou com valor 0.22604845816999916, pertencente à família CLUSTER após sobrescrita da variável r2; o log do teste dobrado relata r(70) aproximadamente 8.182.
- Para a família x_k=2^k-1, delta_k=2^(k-1) quando k>=1, NÃO 2^(k-2). O argumento do Lema L4 pode ser reparado: para k>=4, Delta_0k >= 2^(-k/2-9/2), e para gamma>=1/2,
\[
(\Delta p)_0/p_0\ge\sum_{k=4}^K 2^{(\gamma-1/2)k-\gamma-9/2},
\]
que diverge. Para gamma<1/2, a família âncora+cluster dá a obstrução complementar. A conclusão delimitada é: nenhum peso fixo p_j=delta_j^gamma produz um teste Schur–Collatz uniforme nessa família; outros métodos permanecem possíveis.

## 4. Consequência para a sequência NUMERIA

- A P88 oferece uma rota SPECTRAL paralela que vale estudar, após corrigir o sinal, a linguagem sobre lambda_max(E) e os testes de pi.
- O objetivo nela é provar uma estimativa UNIVERSAL sobre ||Delta||_2, possivelmente com decomposição multiescala e normas intervalares, não extrapolar autovalores float.
- **P93** continua com o objetivo diferente: cinco lacunas interiores do Lema 6 afiado. A P88 não fecha esse orçamento nem muda a demonstração existente P92 para b-a<=4.
- Não há prova nova da RH, nem do Lema 5/6 completos ou da constante espectral geral independente do resultado externo.

**Status:** nota de auditoria, sem modificação dos certificadores matemáticos P92, PR permanece draft.
