# Correção de borda e complemento de Schur

Esta rodada usa a mesma meia-largura racional $a$ e o mesmo déficit espectral da P93-V.40. Com $\delta=2a-\log16>0$, $\Omega=1\,600\,000$, $\eta=8/15$ e $\kappa=\log(2)/4$, a estimativa de Binet absorve o custo de alta frequência do termo de 16 e deixa reserva superior a $0{,}0056685$. Na parte de baixa frequência, seja $g=P_{[-\Omega,\Omega]}f$. A desigualdade $\|g'\|_\infty^2\le\Omega^3\|f\|_2^2/(3\pi)$ limita o resto nas duas faixas de borda. Por paridade, os módulos das avaliações de $g$ em $a$ e $-a$ são iguais. A forma completa de Weil é limitada inferiormente pelo minorante auxiliar menos

$$
\gamma \lvert\langle e,f\rangle\rvert^2+\epsilon\|f\|^2,\qquad
e(x)=\frac{\sin(\Omega(a-x))}{\pi(a-x)},\quad
\gamma=4\kappa(1+\eta)\delta,\quad
\epsilon=\frac{4\kappa(1+\eta)\delta^3\Omega^3}{9\pi}.
$$

Usamos extremos superiores diádicos de $\gamma$ e $\epsilon$. O resto escalar é menor que $2{,}804\cdot10^{-208}$. O núcleo da borda esquerda é o reflexo de $e$. A correção ocupa até duas dimensões no espaço inteiro e uma em cada setor.

## Eliminação dirigida no espaço finito

Escreva o minorante após a dedução de $\epsilon I$ na forma de blocos $\bigl(\begin{smallmatrix}A&C\\C^*&D\end{smallmatrix}\bigr)$, e decomponha a projeção do núcleo $e$ da respectiva paridade em $(u,v)$. Após verificar $D\succ0$ por congruência, uma resolução intervalar encontra recintos para $D^{-1}C^*$ e $D^{-1}v$. O bloco da cauda corrigido é positivo se $r=1-\gamma v^*D^{-1}v>0$. Nesse caso o complemento exato é

$$
S=A-CD^{-1}C^*-
  \frac{\gamma}{r}(u-CD^{-1}v)(u-CD^{-1}v)^*.
$$

A positividade de $S$ por congruência, combinada com a condição sobre a cauda, prova a positividade da restrição finita com o termo de 16. A verificação usa o vetor de borda acoplado $u-CD^{-1}v$. Também converte a margem de Gershgorin transformada em um limite inferior na base ortonormal, dividindo-a por uma cota superior da norma de Frobenius ao quadrado da congruência. Os testes comparam a fórmula a uma eliminação racional independente e verificam a rejeição dos casos em que o denominador é singular, negativo ou inconclusivo.

A montagem de 200 modos por paridade foi feita em Arb a 768 bits. A cabeça usa os modos 1 a 100 e o complemento finito, os modos 101 a 200. A carga $\gamma v^*D^{-1}v$ ficou abaixo de $3{,}221\cdot10^{-72}$ nos dois setores. A eliminação corrigida também fornece um piso para a forma finita antes de qualquer transformação. As margens abaixo foram arredondadas para baixo.

| Paridade | Piso do complemento de Schur | Piso da forma finita com o termo de 16 |
| --- | ---: | ---: |
| Par | $2{,}379\cdot10^{-75}$ | $2{,}285\cdot10^{-76}$ |
| Ímpar | $3{,}506\cdot10^{-71}$ | $3{,}360\cdot10^{-72}$ |

O piso da forma finita usa $D-\gamma vv^*\succeq rD$ e a triangular de eliminação com bloco $(D-\gamma vv^*)^{-1}(C-\gamma uv^*)^*$. Seu módulo em norma de Frobenius é verificado por um novo produto com os vetores da mesma resolução intervalar. [Certificado dos 200 modos](edge_schur_finite_200.json)

## Complemento infinito

O teste usa o minorante localizado com $\theta=0$. A localização efetiva fica dentro da cabeça. Assim os termos localizados de déficit e a penalidade escalar no complemento da localização têm bloco cruzado nulo. A concentração residual calculada em Arb por elipses de Bernstein dá um piso $D\succeq dI$ para a cauda. Para a norma quadrada da parte caudal do núcleo de borda, $b\ge\|v\|^2$, usamos as normas analíticas dos núcleos sinc projetados por paridade e subtraímos a soma dirigida dos quadrados dos coeficientes da cabeça.

Se $d-\gamma b>0$ e $CC^*\preceq W$, a fórmula aplicada à substituição $D\mapsto dI$ e Young com parâmetro $s=1/16$ dão uma condição suficiente de positividade na cabeça

$$
A-\left(\frac1d+\frac{\gamma(1+1/s)b}{d(d-\gamma b)}\right)W
  -\frac{\gamma d(1+s)}{d-\gamma b}uu^*\succ0.
$$

A opção escalar usa $W=\|C\|_{\rm cota}^2 I$. Na opção direcional, para cada banda de frequência $j$, a norma do bloco cruzado é limitada na ordem de operadores por $C_jC_j^*\preceq 18^2 g_jG_{P,j}$. Aqui $g_j$ é a concentração na banda no complemento e $G_{P,j}$ é a Gram da banda restrita à cabeça. Aplicar Cauchy–Schwarz à soma dos blocos, incluindo o polo, dá $W=8(\sum_j18^2 g_j G_{P,j}+b_{\rm polo}^2 I)$ para sete bandas. Cada $G_{P,j}$ usa as primitivas de Wronskian. Os $g_j$ usam a cota integrada de Chebyshev com início no grau $2N+p$, onde $N$ é o corte por paridade e $p$ é a paridade.

O piso da cauda corrigida excede $0{,}118$ em 3.309 modos por paridade. A compressão da condição suficiente aos primeiros 200 modos foi montada e verificada em Arb a 768 bits, com pisos superiores a $2{,}379\cdot10^{-75}$ no setor par e $3{,}506\cdot10^{-71}$ no ímpar. Isso certifica a condição para a submatriz principal de 200 modos. A prova no espaço inteiro ainda exige os modos 201 a 3.309 e seu acoplamento aos 200 primeiros. Com corte em 400 modos a estimativa da cauda permanece inconclusiva. O teste usa o minorante na janela local especificada sem tomar a Hipótese de Riemann como entrada. [Orçamento do complemento](edge_schur_complement.json) · [Compressões verificadas](edge_schur_finite_200.json)

O orçamento escalar par em 3.309 modos exigiria uma margem superior a $4{,}5218\cdot10^{-72}$. Uma testemunha diádica de 103 modos tem quociente do minorante base no recinto $(3{,}95683,3{,}95684)\cdot10^{-72}$. Para essa mesma direção normalizada, o custo do acoplamento por bandas é inferior a $4{,}441\cdot10^{-146}$ e a correção de borda custa menos que $6{,}640\cdot10^{-147}$. O teste direcional fecha nessa direção. [Recintos da testemunha](edge_schur_directions.json)

## Reprodução

```sh
python -m pip install -r NUMERIA_P93V40/requirements.txt
sha256sum -c NUMERIA_P93V40/SHA256SUMS
python -O -m NUMERIA_P93V40.edge_schur --verify NUMERIA_P93V40/edge_schur_complement.json
python -O -m NUMERIA_P93V40.edge_schur --verify NUMERIA_P93V40/edge_schur_finite_200.json
python -O -m NUMERIA_P93V40.edge_schur --verify NUMERIA_P93V40/edge_schur_directions.json
```

O replay da testemunha usa os recintos diádicos da montagem finita. O SHA-256 do certificado de origem é conferido ao reproduzir o JSON direcional. A reprodução da montagem finita recalcula independentemente esse certificado de origem a partir das tabelas e quadraturas.
