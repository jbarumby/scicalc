# NUMERIA–RIEMANN · P93-V.40 · Minorantes para a passagem de 16

A V.40 importa os arquivos da V.39 sem alterar seus bytes. O ZIP recebido tem SHA-256 `574ad27dfdd137e7938548ef3bfe8cae9bf57a6be386dbc5454b45aa90d84443` e as 35 entradas do manifesto conferem. O PR #109 está aberto em rascunho, no commit `cbc67f2ec6388aca71d72402680b6eb586b4825e`. Ele trata da etapa P92-T em Lean. A base matemática desta implementação é o ZIP, importado numa branch separada a partir da `main`. [V.39](../NUMERIA_P93V39/PROVA.md) · [PR #109](https://github.com/jbarumby/scicalc/pull/109)

## Resultados executados

A montagem física certificou os prefixos pedidos de **25 e 50 modos por paridade**, com deslocamento $h=10^{-60}$. A rodada de 100 modos, com GL140 e 768 bits, reaproveitou a mesma montagem para provar **65 modos pares e 75 ímpares**.

Na dimensão seguinte de cada setor, a rodada verificou vetores diádicos com quociente de Rayleigh nos recintos abaixo. Os recintos decimais da tabela foram arredondados para fora a partir do JSON.

| Setor | Último prefixo provado nesta rodada | Suporte do vetor | Quociente de Rayleigh da cabeça sem deslocamento |
| --- | ---: | ---: | ---: |
| Par | 65 | 66 | $(6{,}55608,\ 6{,}55610)\cdot10^{-61}$ |
| Ímpar | 75 | 76 | $(7{,}74850,\ 7{,}74852)\cdot10^{-61}$ |

Esses intervalos ficam inteiramente entre zero e $h$. Logo a desigualdade pedida com $hI$ falha em todas as cabeças maiores que contêm esses vetores, incluindo os tamanhos 100 a 3309. As cinco variantes harmônicas têm a mesma forma nesses vetores porque seu suporte fica nos primeiros 600 modos. Esse teste resolve o sinal da forma deslocada na direção testemunha. A positividade da matriz sem deslocamento em todas as direções maiores permanece em aberto. [Certificado dos 100 modos e testemunhas](prefix_replay_100.json)

Uma segunda avaliação por Wronskians recompõe $v^\top Bv$ para cada testemunha. Somando-a ao recorte anterior, os quocientes da cabeça sem penalidade $A_0$ são menores que $6{,}5610\cdot10^{-61}$ e $7{,}7571\cdot10^{-61}$, respectivamente. Assim, reduzir ou redistribuir apenas uma penalidade positiva não recupera $h=10^{-60}$. Com a margem tentativa $10^{-64}$, a transferência ainda teria saldo para a largura do ZIP, e a cota 18 exige 3316, 3317, 3318, 3318 e 3320 modos na ordem dos cinco $\theta$. A positividade dessas cabeças com margem menor continua como condição. [Origem da obstrução](obstruction_origin.json) · [Orçamento com margem menor](transfer_relaxed.json)

As duas testemunhas normalizadas têm transformadas de ordem 1 em algumas frequências comuns, mas os módulos amostrados nos primeiros 30 zeros de zeta retornados pelo FLINT ficam abaixo de $5{,}05\cdot10^{-33}$ no setor par e de $2{,}92\cdot10^{-35}$ no ímpar. As massas espectrais além de 4500 ficam abaixo de $4{,}65\cdot10^{-63}$ e $3{,}95\cdot10^{-63}$. Esse diagnóstico identifica direções quase ortogonais às primeiras amostras e muito concentradas até o corte. A prova da obstrução usa os quocientes de Rayleigh verificados diretamente no integrando. [Diagnóstico espectral](witness_spectral.json)

Os programas originais reproduziram os seis JSONs das matrizes e dos resultados dos dois perfis de 12 modos, com os SHA-256 do manifesto. Com Binet, as duas paridades fecharam os 12 pivôs históricos. O ensaio anterior, sem Binet, parou no décimo segundo pivô ímpar. O auditor usa uma cópia temporária e mantém os certificados originais. A montagem da V.40 calcula diretamente uma forma espectral de cabeça em base ortonormal. Ela é diferente da forma suficiente por momentos usada para o replay histórico. [Replay dos 12 modos](historical_replay.json)

## Déficit espectral e Binet

Para $\beta=2$ e $T=4500$, [o certificado espectral](spectral_replay.json) valida em Arb as **9.833 células com déficit tabulado** e a união das lacunas, coberta por **681.817 folhas de déficit zero**. O recinto vale no intervalo inteiro $4500\le |t|\le10^6$. Acima de $10^6$, a cota inferior de Binet e a soma das amplitudes dos nove termos garantem $\Psi_9(t)>2$. Usamos os extremos diádicos do replay. A análise inteira preserva a convenção $\Psi_9=\Re\psi(1/4+it/2)-\log\pi-2\sum_q\Lambda(q)q^{-1/2}\cos(t\log q)$.

Em Binet, $|\Im(1/4+it/2)^2|=|t|/4$ e $\int_0^\infty u/(e^{2\pi u}-1)\,du=1/24$. Assim, o módulo da correção integral da digama é no máximo $1/(3|t|)$. A soma das amplitudes primas é menor que $9{,}944$, o que confirma a cota histórica $\Psi_9(t)-2<\log|t|+8$. Separar a integral em $u\le |t|/4$ e sua cauda também permite $\log|t|+6{,}1059$. O teste de Schur recompõe exatamente os multiplicadores $2$ para 1.110 componentes e $14$ para 120 componentes. A reserva localizada é **$0{,}118376416138\ldots$**. [Auditoria Binet e reserva](binet_replay.json) · [DLMF, fórmula de Binet](https://dlmf.nist.gov/5.9.E15)

## Montagem, sinais e controle dos erros

As duas penalidades adaptativas completas, de ordem 600, foram montadas a 384 bits, com os multiplicadores e todos os termos de Wronskian. O maior raio de entrada após a montagem fica abaixo de $1{,}722\cdot10^{-102}$ em ambas as paridades. O certificado registra SHA-256 da sequência canônica de extremos diádicos de cada triângulo inferior. Ele verifica a reconstrução de $B$ e suas cotas de arredondamento. A positividade de $B$ decorre da representação como soma ponderada de Gram. [Replay das duas penalidades](penalty_replay_600.json)

Para $\phi_n(x)=\sqrt{(2n+1)/(2a)}P_n(x/a)$, a transformada é $i^n\sqrt{2a(2n+1)}j_n(at)$. Para dois modos de mesma paridade, o produto carrega a fase real $(-1)^{(n-m)/2}$. A cabeça $A_0$ integra o símbolo $\Psi_9-\beta$ na faixa $[0,T]$, adiciona $\beta I$ e conserva o termo dos polos $2(-1)^p vv^\top$. A concentração $G_T$ e as bandas da penalidade $B$ usam primitivas de Bessel. Se $n\ne m$,

$$
\big[x^2(j_n'j_m-j_nj_m')\big]'=\big(n(n+1)-m(m+1)\big)j_n j_m.
$$

Nas diagonais, $I_0(x)=\operatorname{Si}(2x)-\sin^2(x)/x$ e $(2n+1)I_n=(2n-1)I_{n-1}-x(j_{n-1}^2+j_n^2)$. Os extremos repetidos acumulam pesos com sinal antes da multiplicação em bloco. As bandas da penalidade têm **11.060 extremos ponderados**. Os testes comparam as primitivas diagonais e os Wronskians a uma integral direta em Arb.

Para $A_0$, são 345 painéis com abscissas e pesos encerrados em Arb. A quadratura usa GL100 a 384 bits e GL140 quando a precisão solicitada é de pelo menos 512 bits. Nas elipses de Bernstein de parâmetro 4 dos painéis de largura $1/4$ em $[0,16]$, vale $|\Im z|\le15/64$ e os argumentos da digama têm parte real $\ge17/128$. Binet dá $|\psi|<13$. Nos painéis de largura até 16, vale $|\Im z|\le15$ e $7\le\Re z\le4510$. Deslocar cada argumento da digama em 8 dá parte real $\ge3/4$ e permite a cota $|\psi|<14$. Como $|\widehat\phi_n(z)|\le \sqrt{2a}\,e^{a|\Im z|}$ por Cauchy–Schwarz, os módulos dos integrandos são menores que $10^4$ e $10^{41}$ nessas duas famílias. Aproximação de Chebyshev até grau $2m-1$ e positividade dos pesos de GL de ordem $m$ dão por painel a cota $(16/3)|I|\max_{E_4}|f|\cdot 4^{-2m}$. Ela é somada e incluída em cada entrada.

A tentativa inicial de 200 modos com recorrência intervalar ingênua perdeu o recinto no pivô 1. A perda pôde ser reproduzida no cálculo de Bessel próximo da transição de ordem. Agora a recorrência parte do ponto médio exato. Uma norma de energia controla o acúmulo de arredondamentos na região oscilatória e frações contínuas tratam a região restante. Ao fim acrescentamos metade do raio do argumento, por $|j'_n(x)|\le\frac12$, a partir da representação integral. Isso elimina o alargamento provocado por propagar a incerteza do nó em todas as etapas.

A LDL intervalar direta ainda pode alargar os recintos nos pivôs 29 e 30. O código usa a fatoração do ponto médio apenas para propor uma matriz triangular diádica $R$ com diagonal não nula. A decisão refaz $R^\top(A-hI)R$ em Arb e exige uma margem de Gershgorin estritamente positiva. Uma única congruência triangular preserva todos os prefixos, pois suas primeiras $N$ colunas são suportadas nas primeiras $N$ coordenadas. Em caso de falha, o código distingue um recinto contendo zero de um vetor diádico com numerador de Rayleigh estritamente negativo, para a forma deslocada e para a forma sem deslocamento.

## Desigualdades de operadores e transferência

Com $P=P_{600}$, $Q=I-P$ e $K_i\ge0$, Cauchy–Schwarz e Young no seminorma de $K_i$ dão

$$
K_i\ \le\ k_iP K_iP+\frac{k_i}{k_i-1}QK_iQ.
$$

O replay de Schur fornece $L\le CI$ e $L\le D(I-G_T)$. Regularização e inversão da ordem dão $QLQ\le F(H_Q)$, com $F(x)=\frac C\lambda[1-(1+sx)^{-1}]$. A compressão até $N$ majora a compressão de $F(H_Q)$ por $F(H_{Q,N})$, pelo complemento de Schur da inversa. A congruência usa $J=\mathrm{diag}(I,(I+sH_{Q,N})/(1+s))$ e evita inverter uma matriz intervalar. A família usa $\theta=0,1/4,1/2,3/4,1$. Até 600 modos as cinco cabeças coincidem. Os testes independentes comparam a identidade de congruência a uma inversa racional exata.

A transferência aperta de 50 para **18** a cota de $|\Psi_9-\beta|$ na faixa baixa. A monotonicidade de $\Re\psi(1/4+it/2)$ para $t\ge0$ segue termo a termo da série da digama. A fuga da banda usa o resto de Chebyshev de $e^{iatx}$ em $E_\rho$, com $\rho=33/25$. Para um modo ortonormal, ortogonalidade ao truncamento de grau menor dá a cota em módulo $\sqrt{2a}\,2e^{a|t|(\rho-\rho^{-1})/2}\rho^{-n}/(1-\rho^{-1})$. Integrando e somando a progressão de cada paridade obtém-se a cota de traço $g$ usada pelo teste. Os polos têm recintos fatoriais pelo resto de Taylor de $e^{x/2}$.

Para as cinco escolhas de $\theta$ na ordem acima, a condição de Schur com a cota 18 requer cabeças de **3307, 3308, 3309, 3310 e 3311** modos por paridade. Com a cota anterior 50, o mesmo cálculo dá 3309, 3310, 3310, 3311 e 3312. São dimensões suficientes para a passagem ao complemento, dadas as hipóteses matriciais. [Cálculo dirigido](transfer.json)

A largura racional do ZIP satisfaz $2a-\log16>0$, cerca de $1{,}22114\cdot10^{-75}$, e $2a<\log17$. O único termo novo é o de $16=2^4$, de peso $\kappa=\log2/4$. O produto das translações vive em duas faixas de borda de comprimento $\delta=2a-\log16$. Para $\Omega=4\cdot10^6$, a reserva espectral acima de $\Omega$ supera 1 sobre o auxiliar. Cauchy–Schwarz nas bordas dá perda no máximo $\kappa\|f\|^2_{\mathrm{bordas}}$. Aplicar Young com $1/4$ à decomposição de Fourier absorve a parte alta, pois $5\kappa<1$. A norma suprema da parte baixa é limitada em quadrado por $\Omega\|f\|^2/\pi$, o que dá perda restante $\Omega\delta/6$. Ela é menor que $8{,}141\cdot10^{-70}$. Assim, se as duas cabeças requeridas fecharem, a margem transferida será maior que $8{,}99999999\cdot10^{-61}$. O cruzamento efetivo ainda depende da cabeça completa.

### Teste com margem menor e duas faixas

Com $h=8\cdot10^{-70}$ a rodada a 768 bits certificou 91 modos pares e 100 ímpares. A primeira direção par abaixo do alvo tem suporte em 92 modos e quociente em $(5{,}05603,5{,}05604)\cdot10^{-70}$. Esse valor também está abaixo de $\Omega\delta/6$ para o orçamento anterior, mesmo se a condição de Schur retiver toda a margem disponível da cabeça. [Replay com margem reduzida](prefix_strip_floor_100.json)

O operador de projeção na faixa $[-\Omega,\Omega]$ tem núcleo $\sin(\Omega(x-y))/(\pi(x-y))$. Restringindo a união das duas faixas de borda, cada bloco diagonal tem linhas com integral absoluta no máximo $\delta\Omega/\pi$. O outro bloco contribui no máximo $\delta/[\pi(\log16-\delta)]$, pela distância entre as faixas. O teste de Schur dá a soma dessas duas cotas como limite para a norma da compressão. Com $\Omega=1{,}6\cdot10^6$ e parâmetro de Young $8/15$, o custo da parte alta fica abaixo da reserva dirigida de Binet e a perda uniforme da parte baixa cai para menos de $1{,}652483\cdot10^{-70}$. As dimensões suficientes para a transferência com $h=8\cdot10^{-70}$ e a cota baixa 18 tornam-se 3326, 3327, 3328, 3329 e 3330, sob a condição matricial. [Orçamento dirigido das duas faixas](transfer_separated_strips.json)

Um controle mais fino da concentração conserva $\rho=33/25$. Escreva $d_\rho=\rho-\rho^{-1}$, $Z=aTd_\rho$ e $n_0=2N+p$. A fórmula de Cauchy para os coeficientes de Chebyshev dá $|c_n|\le 2\rho^{-n}I_0(atd_\rho/2)$. Cosseno e seno preservam as respectivas paridades. Projetar ortogonalmente em Legendre tem erro em $L^2$ limitado pelo da série de Chebyshev truncada. O erro quadrático da série na norma ponderada de Chebyshev é $(\pi/2)\sum_{n\ge n_0}|c_n|^2$, e essa norma majora a norma sem peso. Na definição integral de $I_0$, $\cos u\le 1-2u^2/\pi^2$ para $0\le u\le\pi$. Isso dá $I_0(v)\le\min(e^v,e^v\sqrt{\pi}/(2\sqrt{2v}))$. Integrando em duas partes, com divisão em $T/2$, resulta a cota de traço

$$
g\le \frac{2\rho^{-2n_0}}{d_\rho(1-\rho^{-4})}
\left[e^{Z/2}-1+\frac{\pi}{2Z}(e^Z-e^{Z/2})\right].
$$

Com as duas faixas, margem condicional $h=8\cdot10^{-70}$, retenção $h/4$ no complemento e a cota baixa 18, os cinco limiares caem para **3305, 3306, 3307, 3308 e 3309**. Se as hipóteses matriciais valessem, sobraria mais de $3{,}47\cdot10^{-71}$ após o termo de 16. [Recintos e teste de Schur refinados](transfer_chebyshev.json)

Uma rodada adicional com $h=10^{-70}$ provou 96 modos pares e 100 ímpares e achou, no suporte par de 97 modos, um quociente rigoroso em $(8{,}61969,8{,}61970)\cdot10^{-71}$. Esse limite superior fica abaixo até da estimativa uniforme refinada das faixas para a largura do ZIP. Ele exclui os alvos de cabeça necessários por **essa estimativa**, que cobra a mesma perda em todas as direções. Na rodada **sem deslocamento**, a congruência certificou as duas matrizes de **200 modos**. Os resultados apontam para uma transferência que acompanhe o custo de borda nas direções quase nulas, ou para uma janela mais próxima de $\log16/2$ depois de obter coercividade para uma cabeça suficiente. [Replay de 100 modos](prefix_two_strip_floor_100.json) · [Replay de 200 modos](prefix_unshifted_200.json)

### Custo de borda dentro da cabeça

Para um polinômio de um setor, as duas bordas têm o mesmo supremo de $|f|$. Seja $u_n=\phi_n(a)$ e seja $M_N^2=\sum_{n\in \mathrm{setor}}(2n+1)n^4/(2a^3)$, com a soma restrita à cabeça. Pela desigualdade de Markov, $\|f'\|_\infty\le M_N\|f\|_2$. Para $x$ em uma das duas faixas de comprimento $\delta$, isso implica $|f(x)|\le|f(a)|+\delta M_N\|f\|_2$. O módulo da contribuição nova é no máximo $\kappa\|f\|_{\mathrm{bordas}}^2$. Aplicar $(x+y)^2\le 2x^2+2y^2$ dá a majorante matricial $4\kappa\delta\,uu^\top+4\kappa\delta^3M_N^2 I$. Sua parcela dominante tem posto um em cada paridade. Um teste independente integra diretamente a translação de polinômios de grau até três com Gauss–Legendre exato no grau do produto.

Para a testemunha par de 97 modos, a avaliação direta com coeficientes exatos e raio de derivada dá custo absoluto de 16 menor que $4{,}366\cdot10^{-144}$ após normalizar $\|f\|_2=1$. A folga em uma direção isolada indica onde a transferência uniforme perde precisão. A passagem no espaço inteiro exige controlar também os termos cruzados com o complemento. [Diagnóstico de borda](witness_spectral_two_strip.json)

Subtrair a majorante de posto um e seu resto escalar das cabeças da V.40 a 768 bits certificou **a forma de Weil com o termo de 16 em 100 modos de cada paridade**, na largura racional estritamente acima de $\log16/2$. A forma de Weil é positiva nesses subespaços porque a cabeça penalizada é um limite inferior da forma completa ali e a congruência verificou a positividade estrita. [Certificado com a contribuição de 16](prefix_prime16_100.json)

Para estender o tratamento da borda ao espaço inteiro, tome $g=P_{[-\Omega,\Omega]} f$ na decomposição de Young. Por Cauchy–Schwarz em Fourier, $\|g'\|_\infty^2\le\Omega^3\|f\|_2^2/(3\pi)$. Integrar $2|g(\pm a)|^2+2x^2\|g'\|_\infty^2$ em cada faixa, com $0\le x\le\delta$, dá uma alternativa à perda escalar anterior. Ela cobra $2\kappa(1+\eta)\delta$ vezes a soma dos quadrados das duas avaliações de $g$ e um resto escalar menor que $2{,}804\cdot10^{-208}$ na largura do ZIP. Essas avaliações são funcionais limitados, representados pelo núcleo sinc centrado em $\pm a$. Falta incluir a correção de posto dois e seus termos cruzados no teste de Schur da cabeça completa. [Constantes dirigidas da alternativa](transfer_chebyshev.json)

## Material complementar

O JSON A1/A4 da V.38 foi incorporado byte a byte. O auditor reconstruiu 25 perfis, 50 constantes de Schur, seis produtos de Perron com pesos positivos em aritmética inteira, 212 sinais e 120 entradas da grade de dimensões. Para A4 ($T=4500$, $k=20$, $\lambda=9/10$), ambas as condições de transferência fecham em $N=3350$ na grade testada. A aritmética dirigida também recupera os valores estritamente positivos que foram impressos como zero por subfluxo de `float64`: 20 ocorrências de $g$ e 16 de $b^2$. A faixa do déficit de 1600 a 4500 foi verificada em Arb e complementa a validação da V.40 de 4500 até $10^6$. A passagem de A4 continua condicionada à positividade das cabeças completas. [Entrada original](v38_inputs/P93V38_A1_A4_resultados.json) · [Replay V.38](v38_audit.json)

O post indicado reproduz imagens da primeira versão do guia de Rubinstein. A auditoria relaciona a notação das imagens ao PDF integral, reproduz cálculos finitos da exposição e registra as hipóteses usadas no argumento de Alpöge. [Auditoria do guia](RUBINSTEIN.md)

## Reprodução

Da raiz do repositório, em Python 3.12

```sh
python -m pip install -r NUMERIA_P93V40/requirements.txt
sha256sum -c NUMERIA_P93V40/SHA256SUMS
python -m unittest discover -s NUMERIA_P93V40 -t . -v
python -O -m unittest discover -s NUMERIA_P93V40 -t . -v
python -m NUMERIA_P93V40.historical_replay --verify NUMERIA_P93V40/historical_replay.json
python -O -m NUMERIA_P93V40.spectral_replay --verify NUMERIA_P93V40/spectral_replay.json
python -m NUMERIA_P93V40.digamma_binet --precision 256
python -m NUMERIA_P93V40.transfer
python -O -m NUMERIA_P93V40.prefix --verify NUMERIA_P93V40/prefix_replay_50.json
python -O -m NUMERIA_P93V40.prefix --verify NUMERIA_P93V40/prefix_replay_100.json
python -O -m NUMERIA_P93V40.prefix --verify NUMERIA_P93V40/prefix_strip_floor_100.json
python -O -m NUMERIA_P93V40.prefix --verify NUMERIA_P93V40/prefix_two_strip_floor_100.json
python -O -m NUMERIA_P93V40.prefix --verify NUMERIA_P93V40/prefix_unshifted_200.json
python -O -m NUMERIA_P93V40.prefix --verify NUMERIA_P93V40/prefix_prime16_100.json
python -m NUMERIA_P93V40.origin --verify NUMERIA_P93V40/obstruction_origin.json
python -m NUMERIA_P93V40.penalty --verify NUMERIA_P93V40/penalty_replay_600.json
python -m NUMERIA_P93V40.transfer --margin 8e-70 --separated-strips
python -m NUMERIA_P93V40.transfer --margin 8e-70 --separated-strips --chebyshev-tail --retain 1/4
python -m NUMERIA_P93V40.witness_spectral NUMERIA_P93V40/prefix_replay_100.json
python -O -m NUMERIA_P93V40.v38_audit --verify NUMERIA_P93V40/v38_audit.json
python -O -m NUMERIA_P93V40.rubinstein_replay --verify NUMERIA_P93V40/rubinstein_replay.json
```

A certificação é para a forma finita especificada. Ela deixa em aberto o cruzamento completo e a Hipótese de Riemann global.
