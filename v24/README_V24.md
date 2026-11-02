# NUMERIA V24 — erro global de suavização em um modelo da forma de Weil

**Status: envelope de erro explícito e racional, condicionado à identificação e aos lemas analíticos descritos. Não é certificado de RH.**

A V23 provou apenas uma cota uniforme \(\|\widehat{\phi}-\widehat{\chi_\delta\phi}\|_\infty=O(\delta^2)\). Não se pode multiplicá-la por um peso arquimediano que cresce com \(r\) e integrá-la em \([0,\infty)\) sem controle adicional. Esta V24 introduz **uma segunda estimativa de decaimento para a transformada de Fourier** das funções de teste e propaga a diferença aos três tipos de termos do **modelo**: arquimediano, primos e polos.

## 1. Classes de funções e corte

\(\phi_k(x)=\sin(k\pi(x+L)/(2L))/\sqrt L\), estendida por zero fora de \([-L,L]\), com \(L=9/20\), \(k=1,\ldots,10\). Para \(k=1\) corresponde a \(\cos(\pi x/(2L))/\sqrt L\) e para \(k=2\) a \(-\sin(\pi x/L)/\sqrt L\), sinais sem efeito nas cotas absolutas.

\(\chi_\delta(x)=\eta((L-|x|)/\delta)\), com o corte padrão \(C^\infty\) da V23, \(\eta(t)=\exp(-1/t)/[\exp(-1/t)+\exp(-1/(1-t))]\) para \(0<t<1\), estendido por 0 e 1. A função \(\chi_\delta\phi_k\) pertence a \(C_c^\infty\); \(0<\delta<L\).

Defina \(\omega_k=k\pi/(2L)\), \(a_k=\omega_k/\sqrt L\), \(c_k=\omega_k^2/\sqrt L\) e
\[
e_k=2\delta^2 a_k,\quad B=2\sqrt L.
\]
Por valor médio e suportes nas bordas,
\[
\|\phi_k-\chi_\delta\phi_k\|_1\le e_k,\qquad
|\widehat\phi_k(r)-\widehat{\chi_\delta\phi_k}(r)|\le e_k.
\]
A cota \(e_k\) é conservadora: não usa a constante ótima.

## 2. Novo lema de derivadas do corte e cauda de Fourier

Em \(0<t<1/2\), defina \(h=-1/t+1/(1-t)\), \(\eta=(1+e^{-h})^{-1}\). Então \(h'\le t^{-2}+4\), \(|h''|\le2t^{-3}+16\) e \(\eta(1-\eta)\le e^{2-1/t}\). Para \(x=1/t\ge2\), maximizar \(e^{2-x}x^p\) para \(p=2,3,4\), usando \(e>8/3\), fornece
\[
|\eta'|\le8,\qquad|\eta''|\le128.
\]
A simetria \(\eta(1-t)=1-\eta(t)\) estende ambas a \([0,1]\).
Essas desigualdades são **lemas clássicos por cálculo**, documentados; ainda não estão formalizados em Lean.

Duas integrações por partes, com \(\phi_k(\pm L)=0\), dão para \(|r|\ge1\):
\[
|\widehat\phi_k(r)|\le D_k/r^2,\quad
D_k=2a_k+2Lc_k.
\]
Para a função suavizada, \(2\int|\chi'_\delta||\phi'_k|\le 4(8)a_k=32a_k\) e
\(\int|\chi''_\delta||\phi_k|\le2(128)a_k=256a_k\). Assim,
\[
|\widehat{\chi_\delta\phi_k}(r)|\le
D'_k/r^2,\quad D'_k=2Lc_k+288a_k.
\]
O importante é que **\(D'_k\) é independente de \(\delta\)**. O salto na derivada da função original aparece explicitamente no termo \(2a_k\) de \(D_k\).

## 3. Peso arquimediano de r=0 ao infinito

Adotamos \(\Omega(r)=\operatorname{Re}\psi(1/4+ir/2)-\log\pi\), a convenção arquimediana reproduzida na V17/V23. A identidade de série real positiva e a comparação integral para \(u_j=j+1/4\) fornecem
\[
|\Omega(r)|\le|\Omega(0)|+4+\tfrac12\log(1+4r^2).
\]
A V23 fornece \(-6<\Omega(0)<0\). Como \(\tfrac12\log(1+4r^2)\le\log(1+2r)\),
\[
|\Omega(r)|\le12\quad(0\le r\le1),\qquad
|\Omega(r)|\le10+3\sqrt r\quad(r\ge1).
\]
A integral ponderada explícita
\[
\int_1^\infty(10+3\sqrt r)\,r^{-2}\,dr=16
\]
é convergente. Logo, para o **modelo**
\(A_{ij}=\pi^{-1}\int_0^\infty \Omega(r)F_i(r)\overline{F_j(r)}\,dr\),
\[
|\Delta A_{ij}|\le
\frac{12B(e_i+e_j)+16(e_iD'_j+e_jD'_i)}{\pi}.
\]
O código usa majorantes racionais externos para \(\pi,\sqrt L\).

## 4. Polos e primos com suporte compacto

Na variável logarítmica, defina \(M_i(s)=\int\phi_i(x)e^{(s-1/2)x}\,dx\). Para \(s=0,1\),
\[
|M_i(s)|\le Bq,\quad|\Delta M_i(s)|\le e_iq,\quad q=(1-L/2)^{-1},
\]
pela desigualdade \(e^{L/2}\le(1-L/2)^{-1}\). A **soma de dois termos de polos** no modelo tem erro absoluto
\[
|\Delta P_{ij}|\le 2Bq^2(e_i+e_j).
\]
É obrigação de cotejo adicional verificar que os coeficientes e os cancelamentos projetados no E7/V17 usam precisamente esta identificação.

Para a convolução logarítmica \(k_{ij}(t)\), há
\[
|\Delta k_{ij}(t)|\le(e_i+e_j)/\sqrt L.
\]
Como \(\log2<2L=0,9<\log3\), só \(n=2\) contribui na soma de primos e potências com \(|t|\le2L\). Os dois termos \(t=\pm\log2\) dão, com \(\Lambda(2)/\sqrt2=\log2/\sqrt2\),
\[
|\Delta {\rm Prime}_{ij}|\le 2\frac{\log2}{\sqrt L}(e_i+e_j).
\]
O fator \(1/\sqrt2\le1\) foi descartado para ampliar a cota.

**A identidade completa com a implementação de todas as fases E7/V16 não foi demonstrada.** São envelopes rigorosos do modelo explicitamente definido, **condicionais** à fórmula explícita e às suas normalizações.

## 5. Certificação finita de diferenças, sem extrapolação

São testados \(k=1,\dots,10\), \(\delta\in\{10^{-3},10^{-4},10^{-5},10^{-6}\}\). Geram-se envelopes racionais para cada parcela, termos diagonais e máximo da soma por linha da matriz de erros \(10\times10\). Para uma perturbação hermitiana, \(\|E\|_2\le\|E\|_\infty\). Uma aplicação de estabilidade da positividade precisa **além disso** demonstrar que a matriz original/projetada de Weil tem autovalor mínimo maior que essa cota, e provar que a mesma matriz está sendo perturbada. A V24 **não** declara essa prova nem a positividade geral.

## Reprodutibilidade

    python v24/AUDITAR_V24_CONTINUIDADE.py --out v24/RESULTADO_V24_CONTINUIDADE.json
    python v24/TESTAR_V24.py

O teste fixa por SHA-256 os dois scripts matemáticos da V23, verifica simetria, monotonia \(O(\delta^2)\), entradas inválidas e reprodução byte a byte. Um manifesto de hashes acompanha a evidência versionada após sua geração no GitHub.

**Próximas obrigações:** cotejar literalmente projeções e normalização das fórmulas G+P da V16, controlar as aproximações em todas as dimensões, conferir a regra de erro sob polo e excluir circularidade entre suposições de RH e positividade de Weil. A Hipótese de Riemann permanece sem demonstração.
