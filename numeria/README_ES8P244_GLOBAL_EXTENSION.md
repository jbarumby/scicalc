# P244: auditoria do bloco finito e uma rota cofinal para a positividade de Weil

**Autor:** Manus AI
**Data:** 18 de novembro de 2026

A reprodução independente refeita nesta sessão confirma que os dois blocos P237 de dimensão 377 foram processados por C/MPFR a 1024 bits e que os rastros gerados coincidem byte a byte com os anexados. A auditoria racional exata aceita 754 pivôs positivos e a transferência de Weyl dá margem inferior `1/10^18` para todas as matrizes dentro das caixas P237 fornecidas. Isso é um resultado sobre dois blocos finitos. Não prova, por si só, a positividade da forma de Weil para todos os suportes nem a Hipótese de Riemann.

## O que a auditoria P244 verificou

Todos os hashes listados no manifesto do ZIP recebido passaram. Reexecutei `p244_analytic_gates.py`; o JSON produzido coincidiu byte a byte com `P244_ANALYTIC_CHECK.json`. Recompilei o programa C/MPFR, exportei de novo os centros diádicos do `P237_FINITE_LDL_2048.json.gz` versionado e regenerei as duas fatorações/inversas. Os rastros par e ímpar coincidiram byte a byte com os originais, e `P244_REPLAY_CERT.json` local também coincidiu byte a byte com o certificado recebido. Os hashes dos programas MPFR e verificador racional do arquivo P242 coincidem com os hashes fixados no wrapper recebido.

Há uma vulnerabilidade real no auditor P242 isolado. Sua função `load_mpfr_trace` aceita qualquer número positivo no campo `inverse_frobenius_squared_upper`; o valor é então usado como se fosse um limite superior demonstrado para a norma de Frobenius do inverso. Repeti o ataque ponta a ponta com os dados P237 e os rastros anexados: mantive todos os intervalos dos 754 pivôs e troquei apenas o campo `ν` por `1`. O auditor exato P242 aceitou os dois setores e imprimiu “CERTIFIED”. O wrapper P244 rejeitou os mesmos rastros por divergência de SHA-256. Portanto, o certificado original não foi refutado: o problema é que o verificador P242, sozinho, não autentica esse campo. O replay P244 fecha especificamente essa substituição ao reexecutar o código C e comparar os rastros fixados.

O wrapper P244 recebido verifica hashes para o C e o verificador Python, mas não fixa diretamente o hash de `export_centers.py`. A cópia integrada em `numeria/p244/p244_replay_bound.py` acrescenta pins para esse exportador e para o certificado P242. A versão integrada também distingue “rastro existente conferido” de “MPFR reexecutado nesta chamada”. Ela ainda não fixa compilador, bibliotecas, sistema operacional ou flags completas. Os hashes de saída e a reprodução observada são evidência forte para esta execução, não uma atestação criptográfica de todo o ambiente.

Nesta execução, o sandbox não tinha `gmp.h`; a compilação usou um cabeçalho mínimo temporário e um symlink para a biblioteca GMP de runtime já instalada. Nenhum arquivo recebido foi alterado. A reprodução completa passou, incluindo o teste 2×2, mas uma repetição independente de longo prazo deve usar os cabeçalhos oficiais de desenvolvimento GMP/MPFR e registrar as versões do compilador e das bibliotecas.

Os portões racionais P244 conferem consequências exatas dos limites escritos: por exemplo, `M=1384/75<20`, o orçamento GL60 por linha abaixo de `4×10⁻⁴⁷` em cada paridade e o piso de alta frequência `335/169344` em `T=504`. Esses cálculos dependem de lemas analíticos clássicos, da fórmula explícita e da redução de Zhu. O código não prova sozinho a identidade de Binet, a extensão holomorfa uniforme de cada integrando nos painéis, a identidade Legendre–Bessel, a exatidão de Gauss–Legendre ou a correspondência completa da quadratura P237 com a forma original. Tampouco refaz as 1.008 integrais por uma segunda implementação intervalar. A bola positiva de `β*` para `L=161/200` é importada de P237/P239; P244 não a reconstrói independentemente.

## Corolário explícito de cauda para qualquer janela finita

Há, contudo, uma generalização analítica simples e útil do *gate* de alta frequência. Ela torna explícito que esse gate pode ser fechado para cada `L<∞`, embora não feche o bloco de baixas frequências. Defina

\[
A_L=\sum_{\log n<2L}\frac{2\Lambda(n)}{\sqrt n},\qquad
T_L=\max\{504,\,4\pi e^{A_L}\},\qquad
\beta_L=\log\frac{T_L}{2\pi}-\frac1{T_L}-A_L.
\]

A soma é finita e seus termos são não negativos. Como `T_L≥4πe^{A_L}`, temos

\[
\beta_L\ge \log 2-\frac1{T_L}
>\frac12-\frac1{504}=\frac{251}{504}>0.
\]

A estimativa necessária para a cauda decorre da representação de Binet da digama [2]. Para `t≥504`, ponha `z=1/4+it/2`. Então

\[
z^2=\frac1{16}-\frac{t^2}{4}+i\frac t4,
\qquad
\psi(z)=\Log z-\frac1{2z}-2\int_0^\infty
\frac{u\,du}{(u^2+z^2)(e^{2\pi u}-1)}.
\]

Na parte `0≤u≤t/4`,

\[
|u^2+z^2|\ge |\Re(u^2+z^2)|
\ge \frac{3t^2}{16}-\frac1{16}\ge \frac{t^2}{8}.
\]

Usando `∫₀∞ u/(e^{2πu}−1)du=1/24`, o módulo da contribuição dessa parte do resto de Binet é no máximo `2/(3t²)`. Na parte `u≥t/4`, temos `|u²+z²|≥t/4`. Como `π>3` e `u≥126`, vale `1/(e^{2πu}−1)≤2e^{-6u}`, de modo que essa contribuição é limitada por

\[
\frac{16}{t}\int_{t/4}^{\infty}u e^{-6u}\,du
=\left(\frac23+\frac4{9t}\right)e^{-3t/2}
\le \left(\frac23+\frac4{9t}\right)\frac{16}{9t^3}
<\frac1{3t^2}.
\]

A penúltima desigualdade usa `e^x≥x³/6`; a última equivale a `27t²−96t−64>0`, válida para `t≥504`. Portanto o resto integral tem módulo estritamente menor que `1/t²`. Além disso, `Re(Log z)=log|z|≥log(t/2)` e `Re(1/(2z))=2/(1+4t²)≤1/(2t²)`. Concluímos

\[
\Re\psi\!\left(\frac14+\frac{it}{2}\right)
>\log\frac t2-\frac{3}{2t^2},\qquad t\ge504.
\]

Com o símbolo de Zhu,

\[
\Psi_L(t)=\Re\psi\!\left(\frac14+\frac{it}{2}\right)-\log\pi-
\sum_{\log n<2L}\frac{2\Lambda(n)}{\sqrt n}\cos(t\log n),
\]

segue, para `t≥T_L`, que

\[
\Psi_L(t)-\beta_L>
\log\frac{t}{T_L}+\frac1{T_L}-\frac{3}{2t^2}
+\sum_{\log n<2L}\frac{2\Lambda(n)}{\sqrt n}
\bigl(1-\cos(t\log n)\bigr)
\ge \frac1{T_L}-\frac{3}{2T_L^2}>0.
\]

O último passo usa `Λ(n)≥0`, `t≥T_L≥504` e a definição de `β_L`. Assim, para cada janela finita, existe uma escolha explícita de corte com `β_L>0` e símbolo de alta frequência positivo. Isso fornece a condição de sinal usada na redução *one-stroke*. O resultado é uma dedução neste dossiê, não uma reivindicação de novidade: Zhu já parametriza o mesmo corte por `A_L` e reduz cada janela a um bloco finito mais caudas quantitativas [1]. O nosso `T_L` é apenas uma escolha suficiente e simples, não ótima.

O custo permanece importante. `T_L` contém `e^{A_L}`; o tamanho do bloco finito requerido pela redução depende de `L T_L`. Fechar a cauda de frequência para cada `L` não demonstra que os blocos finitos resultantes sejam positivos, nem que o acoplamento com os modos descartados caia abaixo de sua menor margem.

## Critério cofinal por complemento de Schur

A passagem correta de janelas para positividade global pode ser escrita como um alvo quantitativo. Seja `L_j→∞` uma sequência cofinal de suportes e suponha que cada espaço de testes `H_{L_j}` se inclua no seguinte por extensão zero, preservando a forma `Q`. Em cada janela, decomponha `H_{L_j}=K_j⊕T_j`. Suponha que o bloco principal satisfaça `A_j≥m_j I`, a cauda satisfaça `D_j≥η_j I` com `η_j>0`, o acoplamento tenha norma `||B_j||≤b_j`, e o erro total de identificação/aritmética seja limitado por `ε_j I`. Então

\[
Q_j(x,y)\ge (m_j-\varepsilon_j)\|x\|^2
-2b_j\|x\|\|y\|+\eta_j\|y\|^2.
\]

Completando o quadrado,

\[
Q_j(x,y)\ge
\eta_j\left(\|y\|-\frac{b_j}{\eta_j}\|x\|\right)^2
+\left(m_j-\varepsilon_j-\frac{b_j^2}{\eta_j}\right)\|x\|^2.
\]

Portanto, se `m_j−ε_j>b_j²/η_j` para todo `j`, cada janela cofinal é positiva. Para qualquer teste de suporte compacto, escolhe-se uma janela `L_j` que contenha seu suporte e usa-se a compatibilidade por extensão zero; daí `Q(f)≥0` em toda a união dos domínios. Se `Q` se estender continuamente ao fecho pertinente, a desigualdade passa ao fecho. Não é necessária uma margem uniforme positiva: a quantidade `m_j−ε_j−b_j²/η_j` pode tender a zero, desde que permaneça não negativa em cada etapa.

Esse enunciado é uma aplicação direta do complemento de Schur, não uma alegação de novo princípio. Zhu já prova uma redução de janela fixa com cauda e acoplamento explícitos [1]. O critério de Weil, por sua vez, exige positividade em toda a classe pertinente de testes, não apenas numa janela finita [3]. O passo ainda ausente para este projeto é demonstrar, para uma sequência `L_j→∞`, a identificação exata dos blocos com as compressões da forma de Weil, uma cota coerciva válida para cada `D_j`, controles `b_j` e `ε_j` que satisfaçam a desigualdade, e a compatibilidade do domínio.

## O que não se pode concluir

P244 confirma uma cadeia computacional finita do centro às saídas MPFR e às desigualdades racionais de Weyl. A matriz de entrada, contudo, veio da montagem P237; o replay não reavalia as integrais que a originaram. Os próprios artefatos registram que a ponte completa Weil–Zhu e a prova de RH não foram executadas. O acordo float64 de doze entradas é somente diagnóstico e não certifica arredondamento nem quadratura contínua.

A revisão bibliográfica também impede chamar o mecanismo geral de inédito. Zhu apresenta redução de cada janela compacta a um bloco finito mais limites de cauda e acoplamento [1]. Um preprint recente de Liu relata certificados de positividade em janelas fixas até meia-largura `17/16`, declara que não prova RH e diz que a reprodução externa ainda não foi concluída [5] [6]. Isso não valida independentemente o preprint, mas torna imprudente reivindicar prioridade para o método de Schur ou para certificados fixos de janela. O que pode ser relatado com rigor é: a execução P242/P244 foi reproduzida aqui; a vulnerabilidade do verificador P242 foi demonstrada; e a família explícita de cortes acima é uma forma conveniente de organizar o gate de alta frequência. A contribuição substantiva para positividade global só existirá quando os parâmetros do critério cofinal forem verificados para todas as janelas de uma exaustão.

## Referências

[1]: https://arxiv.org/html/2608.24827v2 "Xuefeng Zhu, Weil positivity in compact windows: a finite reduction, certified two-sided bounds, and a Landau–Widom decay law, version 2"
[2]: https://dlmf.nist.gov/5.9 "NIST Digital Library of Mathematical Functions, Section 5.9, integral representations for the psi function"
[3]: https://eudml.org/doc/252338 "Enrico Bombieri, Remarks on Weil's quadratic functional in the theory of prime numbers, I"
[4]: https://people.math.ethz.ch/~hiptmair/Seminars/RAP_22/TRE13.pdf "Lloyd N. Trefethen, Approximation Theory and Approximation Practice, Theorem 8.1"
[5]: https://www.alphaxiv.org/abs/2609.weil-positivity-riemann-zeta-bounds "Vincent Liu, Certified Weil Positivity Beyond the Unit Window, preprint listing"
[6]: https://github.com/luciferyu666/certified-weil-positivity "Frozen reproduction repository for the Liu preprint; external reproduction status as stated by its README"
