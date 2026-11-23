# NUMERIA–RIEMANN P92-Q — Lema das tendas e Lema 6 para submalhas

Base: PR #105 (P92-P), P89 (cotas das caudas exteriores), P90 (orçamento interior e identidade harmônica). Esta etapa desenvolve uma nova comparação variacional que permite buracos internos e pontos isolados. NÃO prova os Lemas 5/6 gerais nem RH.

## Teorema Q1: compactação convexa em malha esparsa

Seja M>=2 inteiro e J={j_1<...<j_p} subconjunto arbitrário de {1,...,M-1}. Acrescente os extremos 0 e M, e ponha
\[
d_i=\min(j_i-j_{i-1},j_{i+1}-j_i),\quad j_0=0,\quad j_{p+1}=M.
\]
Para toda função convexa NÃO NEGATIVA f:(0,M)->[0,infinito) finita nos inteiros interiores,
\[
\boxed{\sum_{i=1}^{p}d_i f(j_i)\le\sum_{k=1}^{M-1}f(k).}\tag{Q1}
\]

**Prova analítica completa:** cada ponto interior j_i define a tenda
\[
w_i(k)=\max\{0,1-|k-j_i|/d_i\},\quad 1\le k<M.
\]
Como d_i é inteiro, a soma aritmética simétrica dá massa e centro exatos
\[
\sum_k w_i(k)=d_i,\qquad\sum_k k\,w_i(k)=j_i d_i.
\]
O suporte da tenda fica entre os vizinhos de j_i. Portanto, para todo k, apenas uma ou as duas tendas de nós adjacentes podem ser positivas. Se j_i<k<j_{i+1} com g=j_{i+1}-j_i, como d_i,d_{i+1}<=g:
\[
w_i(k)+w_{i+1}(k)\le[1-(k-j_i)/g]+[1-(j_{i+1}-k)/g]=1.
\]
Se apenas uma tenda está ativa, vale igualmente a cota 1. Em todo caso,
\[
0\le\sum_i w_i(k)\le1.
\]
Pela desigualdade de Jensen e a identidade do baricentro:
\[
d_i f(j_i)\le\sum_k w_i(k)f(k).
\]
Somando sobre i e usando a capacidade das tendas com f>=0, obtém-se Q1. A prova vale para **todos** M,J, sem enumerá-los.

## Corolário Q2: Lema 6 afiado para submalhas interiores, exterior arbitrário

Considere quaisquer pontos reais distintos ordenados x_0<...<x_(N-1), e s<t. Suponha que haja a>0, M inteiro positivo, tais que
\[
h=x_t-x_s=Ma,\quad
\{x_k:s<k<t\}\subseteq\{x_s+ja:1\le j<M\},
\quad\delta_s\le a,\ \delta_t\le a.
\]
Os pontos exteriores a (x_s,x_t) são **totalmente arbitrários**. Para cada ponto interior x_s+j_i a, sua distância ao vizinho mais próximo é a d_i, porque nenhum exterior ao intervalo pode estar mais perto do ponto interior que um dos extremos.

A função
\[
f(u)=1/[u^2(M-u)^2]
\]
é positiva e estritamente convexa no intervalo (0,M) (a segunda derivada de log f é 2/u²+2/(M-u)²>0 e f>0). Portanto Q1 dá
\[
h^3\sum_{s<k<t}\frac{\delta_k}{(x_k-x_s)^2(x_k-x_t)^2}
\le M^3\sum_{j=1}^{M-1}\frac1{j^2(M-j)^2}=2D(1/M).\tag{Q2}
\]
A **igualdade final** é a identidade harmônica já demonstrada na P90, válida para todo M:
\[
M^3\sum_{j=1}^{M-1}\frac1{j^2(M-j)^2}
=2M H_{M-1}^{(2)}+4H_{M-1}=2D(1/M).
\]

Como \(\delta_s/h,\delta_t/h\le1/M\) e D é decrescente (P90),
\[
2D(1/M)\le D(\delta_s/h)+D(\delta_t/h).
\]
Aplicando o orçamento interior P90-BUDGET e o teorema das caudas P89, segue o Lema 6 **na classe Q2**, sem importar os Lemas 5/6 gerais:
\[
\boxed{
\sum_{k\ne s,t}\frac{\delta_k}{(x_k-x_s)^2(x_k-x_t)^2}
\le \frac{\pi^2}{3h^2}\left(\frac1{\delta_s}+\frac1{\delta_t}\right)
-\frac{3(\delta_s+\delta_t)}{h^4}.}\tag{Q3}
\]
Este resultado permite qualquer quantidade de pontos interiores e lacunas de vários tamanhos (inteiros em unidades a), inclusive pontos isolados com \(\delta_k>a\). O caso M=1 sem ponto interior já está contido no P90.

**Exemplo de 100 unidades:** x_s=0, x_t=100; nós internos em 1,2,5,19,20,46,73,74,99; exteriores em -1 e 101. Há buracos de comprimento 26 e 27 e pontos isolados, portanto os espaçamentos internos NÃO são iguais. A Q1 ainda majoriza a soma pelo modelo uniforme de 100 lacunas.

## Relação com a ponte de potenciais P92-P

A igualdade P92-P para os déficits de energia,
\[
G_s+G_t+\frac{2(T_s-T_t)}h
=\frac{\pi^2}{3}\left(\frac1{\delta_s}+\frac1{\delta_t}\right)
-\frac{3(\delta_s+\delta_t)}{h^2}
-h^2\sum_{r\ne s,t}\frac{\delta_r}{(x_r-x_s)^2(x_r-x_t)^2},
\]
mostra que Q3 fornece uma classe nova em que a compensação é rigorosamente não negativa, sem exigir T_s>=T_t. O complemento de Schur explícito da P92-P é PSD nessa classe.

## Verificação e fronteiras

O arquivo certificate.py verifica Q1 por tendas racionais e identidade harmônica para muitas submalhas, inclusive com buracos e nós isolados. Calcula também a desigualdade completa Q3 para exemplos com pontos exteriores, usando uma cota **inferior** de pi obtida por Machin. Nenhum teste finito é apresentado como demonstração da generalidade de Q1/Q3: as provas são as contas analíticas acima, com hipóteses explícitas.

Reproduzir da raiz:
- python -m unittest -v NUMERIA_P92Q_LATTICE_THINNING.test_p92q
- python -O -m unittest -v NUMERIA_P92Q_LATTICE_THINNING.test_p92q
- python -m NUMERIA_P92Q_LATTICE_THINNING.certificate --out p92q-certificate.json

Limites não alcançados: posições interiores reais arbitrárias que não sejam de uma mesma submalha; remoção da hipótese de distâncias dos extremos <=a; lemas 5/6 gerais, constante de Preissmann geral por demonstração independente e RH. Próximo alvo P92-R: transportar esse argumento a lacunas reais arbitrárias, com um majorante de ocupação e centro adequado. Não realizar merge automático.
