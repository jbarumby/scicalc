# NUMERIA ES8-P261 — relatório matemático de Lee–Yang, Fisher, Laguerre–Pólya e Bost–Connes

## 1. O que afirma o artigo primário

Konstantopoulos, Patie e Sarkar, *A new class of solutions to the van Dantzig problem, the Lee–Yang property, and the Riemann hypothesis*, Annales de l'Institut Fourier 74 (2024), 377–421, DOI 10.5802/aif.3600.

**Convenção de notação:** neste relatório, a zeta completada é \(\xi_{\rm R}(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)\), e a função na variável de Fourier é \(\Xi(t)=\xi_{\rm R}(1/2+it)\). O artigo denomina a última função \(\xi\) de Landau. Essa distinção impede confundir as variáveis complexas.

Na Seção 2.1, eq. (2.7), os autores dão uma densidade de Fourier positiva \(\Phi(x)\) e a representação \(\Xi(t)=\int_{\mathbb R}e^{itx}\Phi(x)\,dx\). Consequentemente \(\Xi(t)/\Xi(0)\) é **incondicionalmente** uma função característica positiva-definida. A classe \(\mathcal D_L=\mathbb{LP}_e\cap\mathbb P_+\), Definição 2.2, exige **adicionalmente** pertencer à classe de Laguerre–Pólya.

- **Teorema 2.3:** \(\Xi/\Xi(0)\in\mathcal D_L\) se e somente se RH.
- **Corolário 2.4:** a função \(t\mapsto\Xi(0)/\Xi(it)\) é função característica associada a densidade de frequência de Pólya simétrica se e somente se RH.
- **Seção 4.1, equações (4.15)–(4.18):** a possibilidade de uma representação \(\Xi=\mathcal J_\Psi\) requer a identificação de uma função de Bernstein a partir de razões de momentos theta. **Não é declarado que essa identificação esteja provada.**

Uma equivalência com RH é um **gate de pesquisa**, não uma solução independente.

## 2. Lee–Yang: prova de caso finito com frações exatas

Escolhemos a cadeia aberta de Ising de quatro spins \(\sigma_j\in\{-1,1\}\), acoplamento \(K\), \(a=e^{-2K}>0\) e fugacidade \(z\). Após remover o fator comum \(e^{3K}\):

\[
P_4(z,a)=\sum_\sigma a^{\#\{j:\sigma_j\ne\sigma_{j+1}\}}z^{\#\{j:\sigma_j=-1\}}.
\]

Enumeração explícita com pesos **racionais** fornece:

\[
P_4(z,a)=z^4+A z^3+B z^2+A z+1,
\quad A=2a+2a^2,\quad B=2a+2a^2+2a^3.
\]

Dividir por \(z^2\) e fazer \(y=z+z^{-1}\) produz a quadrática

\[
Q(y)=y^2+A y+(B-2).
\]

As raízes \(z\) estão todas no círculo unitário **se e somente se** as duas raízes de \(Q\) forem reais e pertencerem ao intervalo \([-2,2]\). O auditor verifica exatamente essa condição por discriminante, coordenada do vértice e sinais nos extremos.

Para \(a=1/2\), \(A=3/2\), \(B=7/4\), e
\[
Q(y)=y^2+\frac32y-\frac14,\qquad
y_{\pm}=\frac{-3\pm\sqrt{13}}4.
\]
Como \(3<\sqrt{13}<4\), temos \(-7/4<y_-<-3/2\) e \(0<y_+<1/4\). Todas as raízes de \(P_4\) estão no círculo. Isso é um **caso-modelo** que não prova o teorema de Lee–Yang geral, já estabelecido na literatura.

Para \(a=2\) (antiferromagnetismo: \(K<0\)), temos \(A=12\), \(B=28\),
\[
Q(y)=y^2+12y+26,\qquad y_\pm=-6\pm\sqrt{10}.
\]
Ambas as raízes de \(Q\) são menores que \(-2\); nenhuma corresponde a \(z\) no círculo unitário. **Todos os coeficientes de \(P_4\) ainda são positivos**. Logo, nem positividade dos coeficientes nem simetria recíproca garantem Lee–Yang sem hipóteses ferromagnéticas.

## 3. Fisher: a variável complexa é outra

Na ausência de campo \(z=1\), a soma sobre a cadeia se fatoriza escolhendo o primeiro spin (duas opções) e, independentemente em cada uma das \(N-1\) arestas, manter o sinal (peso 1) ou invertê-lo (peso \(a\)). Portanto:

\[
\boxed{P_N(1,a)=2(1+a)^{N-1}}.
\]

A única raiz como função de temperatura complexificada \(a\) é \(-1\), multiplicidade \(N-1\). **Lee–Yang e Fisher são estudos de zeros em variáveis diferentes**: campo/fugacidade \(z\) versus temperatura/acoplamento \(a\). Nenhuma lei do círculo universal para Fisher é obtida da lei Lee–Yang.

## 4. Contraexemplo exato à transferência automática para RH

Considere a medida simétrica
\[
\mu=\frac34\delta_0+\frac18\delta_{-1}+\frac18\delta_{+1}.
\]
Sua função característica é
\[
\varphi(t)=\frac34+\frac14\cos t,
\]
e é positiva-definida **para todos os conjuntos de pontos** pelo teorema elementar
\[
\sum_{j,k}\overline c_j c_k\varphi(t_j-t_k)
=\int\left|\sum_jc_j e^{it_jx}\right|^2d\mu(x)\ge0.
\]

Em \(t_j=j\pi/2\), \(j=0,1,2,3\), a matriz de Toeplitz obtida em frações exatas tem **15 menores principais não negativos**: menores de ordem 1 valem 1; os de ordem 2 valem \(7/16\) ou \(3/4\); os de ordem 3 valem \(3/16\); o determinante de ordem 4 é zero.

Entretanto, \(\varphi\) **não** está na classe Laguerre–Pólya. Para uma função real inteira em LP, uma condição necessária é a desigualdade de Laguerre
\[
\mathcal L_1[\varphi](t)=\varphi'(t)^2-\varphi(t)\varphi''(t)\ge0\qquad (t\in\mathbb R).
\]
No nosso exemplo, \(\mathcal L_1[\varphi](\pi)=-1/8<0\), enquanto \(\mathcal L_1[\varphi](0)=1/4>0\).

Equivalentemente, com \(z=e^{it}\), a equação dos zeros é
\[
8z\varphi(t)=z^2+6z+1=0,\qquad z=-3\pm2\sqrt2.
\]
Ambas as raízes têm módulo diferente de 1; existem zeros de \(\varphi\) com \(t\) não real.

**Implicação para NUMERIA:** \(\Xi/\Xi(0)\) já define kernels positivos de Toeplitz. Sem uma nova hipótese de LP (forte como RH), não se obtém localização dos zeros. E a positividade desses kernels não demonstra a positividade da forma de Weil P259, que possui um multiplicador e termos de polos distintos.

## 5. Bost–Connes: campo físico e faixa crítica

Bost–Connes (1995) constrói um sistema \(C^*\)-dinâmico cuja função de partição para o regime físico de inversa de temperatura real \(\beta>1\) é \(Z(\beta)=\zeta(\beta)=\sum_{n\ge1}n^{-\beta}\). A expressão da série vale analiticamente em \(\operatorname{Re}s>1\). **Não** implica um teorema de Lee–Yang para a continuação meromorfa no plano inteiro.

Mesmo a truncagem elementar \(Z_2(s)=1+2^{-s}\) tem zeros em
\[
s=-\frac{(2k+1)\pi i}{\log2},\quad k\in\mathbb Z,
\]
todos de parte real zero, **não** na reta crítica. Portanto o nome “função de partição” por si só não induz zeros na reta \(1/2\), e aproximações de volume finito não transportam automaticamente localização dos zeros.

## 6. Que critérios poderiam contribuir para RH?

A estratégia teoricamente relevante é fornecer uma **propriedade nova e independente** de \(\Xi\), por exemplo uma demonstração global de \(\Xi/\Xi(0)\in\mathcal D_L\), ou um certificado universal de uma cadeia de desigualdades necessárias que, com um teorema adicional de suficiência, force todos os zeros reais. Verificar apenas alguns determinantes finitos, sinais de derivadas em pontos amostrados, ou zeros numericamente encontrados não basta.

O artigo KPS, eqs. (4.15)–(4.18), sugere um problema de **momentos theta / funções de Bernstein**. Um próximo Gate B útil será:
1. reproduzir rigorosamente a normalização dos momentos theta \(F(2n)\) e dos coeficientes de Taylor de \(\Xi\);
2. deduzir e auditar, com intervalos dirigidos, as sequências candidatas de Bernstein no índice \(n\);
3. separar condições necessárias finitas de alguma condição comprovadamente suficiente **para todos os índices**;
4. exigir a normalização exata e um operador explícito antes de qualquer suposta transferência para a forma completa de Weil.

**Nada nesta P261 altera** os certificados locais de Weil da P252, a nova testemunha negativa do minorante P259, a auditoria arquimediana P260 ou o status da RH. O Gate A apresenta exclusivamente identidades exatas e contraexemplos para impedir saltos lógicos.


## 7. Um teste de Laguerre que é automaticamente positivo para \(\Xi\)

A representação positiva do artigo KPS mostra que a função \(\varphi_\Xi(t)=\Xi(t)/\Xi(0)\) é a transformada característica de uma distribuição simétrica com densidade \(\Phi(x)/\Xi(0)\), cujo segundo momento é finito pela decadência rápida da densidade theta. Logo, **sem RH**,

\[
\varphi_\Xi(0)=1,\qquad \varphi_\Xi'(0)=0,\qquad
\varphi_\Xi''(0)=-\mathbb E[X^2]<0,
\]

e, consequentemente,

\[
\mathcal L_1[\varphi_\Xi](0)
=\varphi_\Xi'(0)^2-\varphi_\Xi(0)\varphi_\Xi''(0)
=\mathbb E[X^2]>0.
\]

Isso é uma **obstrução metodológica**: conferir apenas a primeira desigualdade de Laguerre na origem, ou o primeiro menor de matriz de momentos, não acrescenta informação à RH. O contraexemplo exato da Seção 4 também possui \(\mathcal L_1(0)>0\), mas falha em \(t=\pi\).

## 8. Gate B quantitativo baseado na Equação (4.18) de KPS

A Seção 4.1 do artigo apresenta, como condição de uma possível representação em uma classe de processos de Lévy/Bernstein, para os índices em que as razões são definidas,

\[
\phi(n+1)=-\frac{G(2n)}{8(n+1)G(2n+2)},\quad
G(2n)=64n(2n-1)\frac{F(2n-2)}{F(2n)}-1,
\]

com \(F(k)=\int_1^\infty(\log x)^k x^{-3/4}\theta_0(x)dx\) e \(\theta_0(x)=\sum_{m\ge1}e^{-\pi m^2 x}\).

Uma função de Bernstein é não negativa, crescente e côncava, logo suas amostras devem ter diferenças de primeira ordem não negativas e de segunda ordem não positivas. **Próxima verificação realmente informativa:** obter enclosures intervalares rigorosos dos momentos \(F(0),F(2),\ldots,F(2N)\) com resíduos das integrais e das somas theta; construir os quocientes sem dividir por intervalos contendo zero; e testar essas **condições necessárias** com aritmética racional/Arb.

Uma violação rigorosa poderia **excluir essa representação candidata de KPS**, mas **não refutaria a RH**, pois a classe candidata \(\mathcal D_P\) não foi demonstrada necessária para a RH. E um número finito de diferenças com os sinais esperados não demonstra que \(\phi\) é uma função de Bernstein em todos os argumentos.
