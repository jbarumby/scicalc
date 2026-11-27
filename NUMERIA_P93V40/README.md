# NUMERIA–RIEMANN · P93-V.40 · Minorantes para a passagem de 16

A V.40 importa os arquivos da V.39 sem alterar seus bytes. O ZIP recebido tem SHA-256 `574ad27dfdd137e7938548ef3bfe8cae9bf57a6be386dbc5454b45aa90d84443` e as 35 entradas do manifesto conferem. O PR #109 está aberto em rascunho, no commit `cbc67f2ec6388aca71d72402680b6eb586b4825e`. Ele trata da etapa P92-T em Lean. A base matemática desta implementação é o ZIP, importado numa branch separada a partir da `main`. [V.39](../NUMERIA_P93V39/PROVA.md) · [PR #109](https://github.com/jbarumby/scicalc/pull/109)

## Resultados executados

A montagem física por bandas certificou **12, 25 e 50 modos por paridade**, com deslocamento de $10^{-60}$. Para 50 modos, a LDL direta ficou inconclusiva, e a congruência triangular verificou ambos os setores. Os cinco valores de $\theta$ têm a mesma compressão até 600 modos. As dimensões acima da maior cabeça certificada continuam em aberto. [Prefixos em Arb](prefix_replay_50.json)

Os programas originais reproduziram os seis JSONs das matrizes e dos resultados dos dois perfis de 12 modos, com os SHA-256 do manifesto. Com Binet, as duas paridades fecharam os 12 pivôs históricos. O ensaio anterior, sem Binet, parou no décimo segundo pivô ímpar. O auditor usa uma cópia temporária e mantém os certificados originais. A montagem da V.40 calcula diretamente uma forma espectral de cabeça em base ortonormal. Ela é diferente da forma suficiente por momentos usada para o replay histórico. [Replay dos 12 modos](historical_replay.json)

## Déficit espectral e Binet

Para $\beta=2$ e $T=4500$, [o certificado espectral](spectral_replay.json) valida em Arb as **9.833 células com déficit tabulado** e a união das lacunas, coberta por **681.817 folhas de déficit zero**. O recinto vale no intervalo inteiro $4500\le |t|\le10^6$. Acima de $10^6$, a cota inferior de Binet e a soma das amplitudes dos nove termos garantem $\Psi_9(t)>2$. Usamos os extremos diádicos do replay. A análise inteira preserva a convenção $\Psi_9=\Re\psi(1/4+it/2)-\log\pi-2\sum_q\Lambda(q)q^{-1/2}\cos(t\log q)$.

Em Binet, $|\Im(1/4+it/2)^2|=|t|/4$ e $\int_0^\infty u/(e^{2\pi u}-1)\,du=1/24$. Assim, o módulo da correção integral da digama é no máximo $1/(3|t|)$. A soma das amplitudes primas é menor que $9{,}944$, o que confirma a cota histórica $\Psi_9(t)-2<\log|t|+8$. Separar a integral em $u\le |t|/4$ e sua cauda também permite $\log|t|+6{,}1059$. O teste de Schur recompõe exatamente os multiplicadores $2$ para 1.110 componentes e $14$ para 120 componentes. A reserva localizada é **$0{,}118376416138\ldots$**. [Auditoria Binet e reserva](binet_replay.json) · [DLMF, fórmula de Binet](https://dlmf.nist.gov/5.9.E15)

## Montagem, sinais e controle dos erros

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
python -m NUMERIA_P93V40.prefix --out /tmp/v40-head.json --modes 50 --bits 384
```

A certificação é para a forma finita especificada. Ela deixa em aberto o cruzamento completo e a Hipótese de Riemann global.
