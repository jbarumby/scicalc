# NUMERIA ES8-P229 — limite independente de erro GL60/Bernstein

## Resultado

Derivei um majorante uniforme simples, independente do arquivo P193 ausente, para todos os integrandos do bloco de 200 modos, para \(T\in\{150,200\}\) e \(0\le n,m\le399\). Em cada elipse de Bernstein mapeada com \(\rho=3\), vale \(|f_{nm}(t)|<100.000\). Aplicando o teorema de Gauss–Legendre para funções analíticas de Trefethen, obtém-se o erro composto exato

\[
E_{T}=\frac{4T\,100.000}{15(3^{120}-1)}.
\]

O denominador \(3^{120}-1\), em vez de \(3^{120}\), torna a cota ligeiramente maior e mais simples de comparar racionalmente. O arquivo `p229_quadrature_bound.py` calcula as duas cotas com frações exatas e registra a comparação P193 em `P229_BOUND.json`.

Para \(T=200\), a cota é \(2,967892467610\times10^{-51}\); para \(T=150\), \(2,225919350707\times10^{-51}\). Elas são aproximadamente 10.674 vezes menores que os `E_GL_total` guardados nos certificados P193. O P228 foi atualizado para usar esta cota nova no lado P227; a referência P193 continua expandida pelo seu próprio erro certificado.

## Derivação do majorante \(M=100.000\)

O integrando para cada entrada é

\[
f_{nm}(t)=\frac{\Psi_L(t)-\beta_T}{\pi}\,u_n(t)u_m(t),\qquad
u_n(t)=(-1)^{\lfloor n/2\rfloor}\sqrt{2L(2n+1)}\,j_n(Lt),\quad L=4/5.
\]

Para prolongar analiticamente o termo digama que é real na reta, uso

\[
\Psi_L^{\mathrm{hol}}(t)=\frac{\psi(1/4+it/2)+\psi(1/4-it/2)}{2}-\log\pi-\sum_{q=2,3,4}a_q\cos(t\log q),
\]

com \(a_2=2\log 2/\sqrt2\), \(a_3=2\log3/\sqrt3\), \(a_4=\log2\). Para real \(t\), isto coincide com a expressão usada na integral.

Cada painel real tem comprimento \(1/2\), é centrado em \(c_k=(2k+1)/4\), e usa \(t=c_k+x/4\). A elipse \(E_3\) no plano \(x\) tem semieixos \(5/3\) e \(4/3\). Logo, em todos os painéis, \(-1/6\le\Re t\le T+1/6\) e \(|\Im t|\le1/3\). Para \(z_\pm=1/4\pm it/2\), resulta \(\Re z_\pm\ge1/12\) e \(|z_\pm|<101\); portanto, a elipse não contém os polos do digama.

Para limitar o digama, desloca-se \(z\) por 200 unidades:

\[
\psi(z)=\psi(z+200)-\sum_{k=0}^{199}\frac1{z+k},\qquad
\psi(w)=-\gamma+\sum_{k=0}^{\infty}\left(\frac1{k+1}-\frac1{k+w}\right).
\]

Com \(\Re z\ge1/12\), \(|z|<101\), \(\Re w>200\) e \(|w-1|<300\), a série deslocada satisfaz

\[
|\psi(w)|\le\gamma+300\frac{H_{199}}{199}<1+\frac{2100}{199}<12.
\]

A soma de retorno satisfaz \(\sum_{k=0}^{199}|z+k|^{-1}\le12+H_{199}<19\). Assim \(|\psi(z_\pm)|<31\).

Também valem \(|\log\pi|<2\), \(a_2+a_3+a_4<7\) e \(|\cos(t\log q)|<2\) nesta elipse. Para a última desigualdade, \(|\Im(t\log q)|<(\log4)/3<2/3\) e \(\log2>2/3\), esta última obtida da série \(\log2=2(1/3+1/(3^3\cdot3)+\cdots)\); logo \(e^{2/3}<2\) e \(\cosh(2/3)<e^{2/3}<2\). Além disso, para \(T=150,200\), \(1<T/(2\pi)<34\), de modo que \(|\log(T/(2\pi))|<33\) pela desigualdade \(\log x<x-1\) para \(x>1\). Com \(1/T<1\) e \(a_2+a_3+a_4<7\), temos \(|\beta_T|<41\). Portanto, \(|\Psi_L^{\mathrm{hol}}(t)-\beta_T|<31+2+14+41=88<100\).

A representação integral

\[
j_n(z)=\frac1{2i^n}\int_{-1}^{1}e^{izx}P_n(x)\,dx
\]

e \(|P_n(x)|\le1\) em \([-1,1]\) dão \(|j_n(z)|\le e^{|\Im z|}\). Para os graus usados, \(n\le399\), temos \(\sqrt{2L(2n+1)}<36\), enquanto \(|\Im(Lt)|\le4/15<1/3\) e \(e^{1/3}<3/2\). Logo \(|u_n(t)|<54\). Como \(\pi>3\), concluímos

\[
|f_{nm}(t)|<\frac{100\cdot54^2}{3}=97.200<100.000.
\]

O mesmo limite vale para as 40.200 entradas superiores, nos dois setores.

## Aplicação do teorema de quadratura

Trefethen, Teorema 4.5, dá para a regra de Gauss de ordem \(r\), em \([-1,1]\), quando \(f\) é analítica e limitada por \(M\) na elipse \(E_\rho\):

\[
|I-Q_r|\le\frac{64M}{15(\rho^2-1)\rho^{2r}}.
\]

Com \(r=60\), \(\rho=3\), e o fator de escala \(1/4\) de cada painel, a cota por painel é no máximo \(2M/[15(3^{120}-1)]\). Há \(2T\) painéis; somando as cotas obtemos a fórmula P229. A regra GL60 integra exatamente polinômios até grau 119; a hipótese restante é a analiticidade e o majorante já demonstrados acima. A forma publicada do teorema é indicada em [1].

## Achado de proveniência P193

Os certificados P193 do ZIP P220 contêm `majorant_M_upper`, `E_GL_total` e o campo `majorant_file: P193_majorant_exact.json`. Esse arquivo, contudo, não está no ZIP P220 e não foi encontrado nas branches/objetos Git acessíveis. A igualdade

\[
E_{\mathrm{P193}}=\frac{9600\,T\,M_{\mathrm{P193}}}{3^{120}-1}
\]

é reproduzível exatamente a partir dos números do certificado, mas isso apenas reconstrói a fórmula numérica; não valida a derivação do majorante P193. O P229 contorna essa lacuna com uma nova cota que não usa \(M_{\mathrm{P193}}\).

Também não endosso aqui a fórmula mais estreita \(2T\cdot120/(3^{120}-1)\) usada no `p220_merge.py` para as novas caixas P220: ela não decorre do majorante P229 e exigiria uma justificativa própria. O limite deste relatório é para os enclosures P227, cuja reconstrução P228 foi atualizada.

## Limites do resultado

A prova do majorante acima é analítica, explícita e deliberadamente grosseira, mas ainda não foi formalizada num assistente de provas. O resultado final continua limitado à quadratura finita e depende também da exatidão dos enclosures Arb exportados pela P227, da fórmula de redução do operador e da avaliação intervalar dos termos não quadratura. Não verifica as caudas/acoplamentos P207, não prova positividade global de Weil e não prova a Hipótese de Riemann.

### Referências

[1]: https://people.maths.ox.ac.uk/trefethen/publication/PDF/2008_127.pdf "L. N. Trefethen, Is Gauss Quadrature Better than Clenshaw–Curtis, SIAM Review 50, Theorem 4.5, 2008"
