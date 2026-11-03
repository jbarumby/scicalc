# NUMERIA–RIEMANN ES6 — fórmula explícita global gaussiana, diferenciais e auditoria aritmética

**Data:** 03/11/2026  
**Escopo:** estender a ES5 da derivada logarítmica em `Re(s)>1` à **fórmula global clássica de Guinand–Weil**, usando uma função-teste analítica e diferenciação paramétrica, sem pressupor RH. A reprodução global permanece exploratória; apenas o termo primo recebeu certificado aritmético rigoroso nesta etapa.

## 1. Formulação e convenção verificadas

O teorema de Guinand–Weil na Proposição 5 de Carneiro–Chirre–Milinovich usa

\[
\widehat h(y)=\int_{-\infty}^\infty h(t)e^{-2\pi ity}\,dt.
\]

Para `h_a(t)=exp(-a t²)` e `g_a(u)=\widehat h_a(u/(2\pi))/(2\pi)=exp(-u²/(4a))/(2\sqrt{\pi a})`, define-se `z_ρ=(ρ-1/2)/i` para **todo** zero não trivial `ρ`, sem assumir que `z_ρ` seja real. A fórmula resulta em

\[
\boxed{\sum_{\rho}e^{-az_\rho^2}
=2e^{a/4}+\frac{1}{2\pi}\int_{-\infty}^{\infty}e^{-at^2}\Omega(t)dt
-\frac{1}{\sqrt{\pi a}}\sum_{n=2}^{\infty}\frac{\Lambda(n)}{\sqrt n}e^{-(\log n)^2/(4a)}}
\]

com `Ω(t)=Re ψ(1/4+it/2)-log π`. Esse teorema **não** é novo e não equivale a ter construído um operador Hilbert–Pólya. A integral arquimediana é calculada numericamente no intervalo `[-T,T]`, usando a paridade. A gaussiana é admissível no espaço holomorfo da fórmula; **não** pertence automaticamente à família de convoluções compactamente suportadas usada nas matrizes V24.

## 2. Uso de diferenciais — nova identidade de teste

A diferenciação em `a` de ambos os lados (justificada pelo decaimento gaussiano e pelo controle de zeros da zeta na faixa crítica) fornece:

\[
\frac{d}{da}\sum_{\rho}e^{-a z_\rho^2}=-\sum_{\rho}z_\rho^2e^{-az_\rho^2}
\]

\[
=\frac12 e^{a/4}-\frac{1}{2\pi}\int_{\mathbb R}t^2 e^{-at^2}\Omega(t)dt
+\frac1{\sqrt{\pi a}}\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}e^{-(\log n)^2/(4a)}\left(\frac1{2a}-\frac{(\log n)^2}{4a^2}\right).
\]

A segunda identidade fornece uma checagem independente mais sensível aos sinais e à normalização, mas **não** prova que somente existam zeros críticos.

## 3. Certificação racional do termo aritmético

Para `a=1/50`, `N=64`, a quantidade positiva

\[
P(a)=\frac1{\sqrt{\pi a}}\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}e^{-(\log n)^2/(4a)}
\]

recebe intervalo racional de arredondamento **para fora**. Implementação:

1. `log n`: série racional de `atanh`, com resto explícito;
2. `π`: fórmula de Machin com série alternada e resto limitado;
3. `sqrt`, inversas e produtos: limites externos por inteiros e `isqrt`;
4. `exp(-x)`: Taylor alternado de um argumento dividido por potência de dois e quadraturas sucessivas com arredondamento externo;
5. cauda `n>N`: usar `0<=Λ(n)<=log n` e, para `N>=64`, a monotonia de `log x/sqrt x * exp[-(log x)^2/(4a)]`.

Após `y=log x`, a integral que limita a cauda é menor que

\[
4\sqrt{a/\pi}\,e^{-(\log N)^2/(8a)}.
\]

Para `N>=64`, `log N>4`, `π>3` e `a=1/50`, resulta `4 sqrt(a/pi)<1/3` e o limite é menor que `exp(-100)/3<10^{-40}`. As desigualdades usadas também são verificadas no código com aritmética racional. O resultado foi encerrado em inteiros na escala `S=10^40`:

\[
\boxed{0.0048201252095083525923710170549320614490\;\le\;P(1/50)\;\le\;0.0048201252095083525923710170549320614627.}
\]

**O intervalo se refere somente à parcela positiva dos primos.** Na fórmula global, ela entra com sinal negativo. Não foi feita uma certificação de igual rigor da integral de digama nem da completude do lado espectral.

## 4. Reproduções numéricas não certificadas

Usamos `mpmath 1.3.0`, `55 dps`, `N=128`, `T=128`, `J=24` zeros críticos conhecidos, pareados por conjugação. Os termos arquimedianos foram computados com quadratura numérica de digama; as somas sobre zeros são apenas parciais.

| a | Polos | Arquimediano | Primos | RHS | 24 pares críticos, soma parcial |
|---:|---:|---:|---:|---:|---:|
| 0.02 | 2.010025041719 | -1.968122333047 | -0.004820125210 | 0.037082583463 | 0.037082583463 |
| 0.04 | 2.020100334168 | -1.949827628642 | -0.069596090972 | 0.000676614554 | 0.000676614554 |

Os resíduos das quantidades calculadas (ponto flutuante, não majorantes de erros totais) foram `−1.85×10^-56` e `1.17×10^-56`. Para a **derivada em a**, foram `4.08×10^-55` e `1.73×10^-55`. Essas magnitudes **não** são certificados de correção: faltam prova independente de quadratura, delimitação da cauda de zeros e prova de ausência de zeros fora da reta; nenhuma dessas propriedades pode ser inferida de coincidência finita.

## 5. Controle adversarial ES4

Calculamos também a soma `2 Σ_{j<=J} exp(-a E_j²)` com o operador diagonal ES4, cujo espectro reproduz a lei assintótica de contagem, mas não os zeros da zeta. Com os mesmos 24 termos:

| a | Soma ES4 | Soma dos 24 zeros críticos | ES4 menos zeros |
|---:|---:|---:|---:|
| 0.02 | 0.003465056409 | 0.037082583463 | -0.033617527053 |
| 0.04 | 0.000005853313 | 0.000676614554 | -0.000670761241 |

Assim, acertar a contagem média dos níveis **não** reproduz a fórmula explícita nem os dados aritméticos. Nenhuma identificação operador–zeta foi estabelecida.

## 6. O que a ES6 permite transportar à V24

Já existe compatibilidade **formal** do kernel arquimediano `Ω(t)` e do peso `Λ(n)/sqrt(n)` com os objetos descritos em `v24/README_V24.md`. Essa coincidência, por si só, não fixa fatores de Fourier, sinais após projeções, polos nem a matriz da forma de Weil. A ES6 usa gaussianas, não funções de suporte compacto; não se devem transportar cotas numéricas ES6 aos pivôs LDL de E7/V24.

**Pendências prioritárias para ES7:** (i) quadratura intervalar de digama com limitação de cauda e arredondamento; (ii) verificação de termos da fórmula explícita com funções compactamente suportadas suavizadas e suas transformadas; (iii) cotejo entrada a entrada com a matriz projetada e as convenções de E7/V24; (iv) sem extrapolação de exemplos finitos para todos os testes. Um esquema de intervalos rigorosos na quadratura poderá utilizar Arb/FLINT, mas não foi instalado nem executado nesta ES6.

## 7. Fontes e reprodutibilidade

- Carneiro, Chirre e Milinovich, **Proposition 5**, *Bounding the log-derivative of the zeta-function*, Mathematische Zeitschrift, https://doi.org/10.1007/s00209-021-02820-9 (fórmula global e convenção de Fourier).
- Paul Garrett, *Guinand's explicit formula*, University of Minnesota, https://www-users.cse.umn.edu/~garrett/m/complex/notes_2020-21/12g_guinand_explicit_fml.pdf
- V24: https://github.com/jbarumby/scicalc/blob/main/v24/README_V24.md
- ES5: https://github.com/jbarumby/scicalc/commit/2484a37170a077089e8c7626a11edb7d92c5ee23

**Conclusão:** ES6 executou uma reprodução auditável da fórmula de Weil com derivada paramétrica, certificou apenas a soma global dos primos em `a=1/50` e estabeleceu outro teste adversarial para ES4. Nenhuma demonstração nova de RH, de positividade universal da forma de Weil ou de equivalência de operadores foi alcançada.
