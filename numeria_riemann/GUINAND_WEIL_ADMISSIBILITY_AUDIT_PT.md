# NUMERIA–RIEMANN — auditoria de admissibilidade da fórmula de Guinand–Weil

**Data:** 2026-10-22  
**Branch:** `numeria-riemann/weil-hermitian-audit`  
**Base preservada:** `6b1da033352ccb02731bc1dc9b81525eb207c80a`  
**Escopo:** auditoria do item 2 — regularidade, decaimento, convergência e correspondência da fórmula explícita com a forma quadrática implementada. Não altera o algoritmo nem aumenta a dimensão da matriz.

## Resultado executivo

A análise indica que a família finita do protótipo parece satisfazer as condições analíticas usuais da versão de Guinand–Weil usada no projeto, desde que a prova de decaimento uniforme abaixo seja incorporada e a versão exata do teorema seja citada de forma estável. A soma sobre zeros converge absolutamente sob a estimativa local padrão de contagem de zeros. A normalização do termo de primos, do peso arquimediano e dos polos coincide algebricamente com uma forma publicada da fórmula.

Isso **não** valida a positividade global de Weil nem prova RH. A positividade intervalar certifica somente a matriz finita projetada para os parâmetros escolhidos. Também não constitui uma revisão formal independente por especialista.

## 1. Classe de funções teste do protótipo

Para um vetor real finito (c), o protótipo define
[
f(t)=sum_{k=1}^{N}c_k h_k(t),qquad
F_c(z)=int_{-L}^{L}f(t)e^{izt},dt,
]
com (f) uma combinação finita de senos em ([-L,L]), estendida por zero fora do intervalo. Cada seno se anula nos extremos, logo (f(pm L)=0). A função é contínua, suave por partes, Lipschitz e tem suporte compacto. A derivada é suave por partes e limitada em cada lado dos extremos.

Defina
[
H_c(z)=F_c(z)F_c(-z),qquad
g_c(u)=int_{mathbb R}f(t+u)f(t),dt.
]
Então (H_c) é inteira e par; (g_c) é par, tem suporte em ([-2L,2L]) e é (C^1). Com a convenção de Fourier usada no código,
[
g_c(u)=rac{1}{2pi}int_{mathbb R}H_c(r)e^{-iur},dr.
]

## 2. Decaimento uniforme em faixas horizontais

Fixe (a>0) e tome (z=x+iy), com (|y|le a). O fator (e^{izt}) e suas derivadas em (t) são uniformemente limitados por constantes dependentes de (a,L) e da função (f), mas não de (x).

Uma integração por partes dá
[
F_c(z)=left[\frac{f(t)e^{izt}}{iz}\right]_{-L}^{L}
-\frac{1}{iz}\int_{-L}^{L}f'(t)e^{izt},dt.
]
O termo de fronteira é zero porque (f(pm L)=0). Integrando por partes uma segunda vez em cada intervalo suave (ou usando a derivada lateral nos extremos),
[
F_c(z)=O_a((1+|x|)^{-2})
]
uniformemente para (|y|le a). Os termos de fronteira da segunda integração envolvem apenas os valores laterais limitados de (f'); a integral restante envolve (f''), que é integrável por partes. Para (|x|) limitado, a cota é obtida pela continuidade e compacidade.

Consequentemente, em qualquer faixa fixa (|\operatorname{Im}z|le 1/2+\varepsilon),
[
H_c(z)=F_c(z)F_c(-z)=O_\varepsilon((1+|\operatorname{Re}z|)^{-4}).
]
Isso satisfaz a hipótese usual (O((1+|\operatorname{Re}z|)^{-1-\delta})) para, por exemplo, (delta=2). O decaimento necessário não é apenas pontual na reta real: a estimativa é uniforme na faixa relevante.

## 3. Convergência absoluta da soma sobre zeros

Escreva (\gamma_\rho=(\rho-1/2)/i). Como (0<\operatorname{Re}\rho<1), temos (|\operatorname{Im}\gamma_\rho|<1/2). Portanto os argumentos (\gamma_\rho) pertencem a uma faixa horizontal coberta pela estimativa acima.

A estimativa local padrão de contagem de zeros,
[
N(T+1)-N(T)=O(\log(T+2)),
]
implica que a soma dos módulos dos termos em blocos de altura (T) é
[
O(\log(T+2)(1+T)^{-4}).
]
A série (sum_{Tge1}\log(T+2)/T^4) converge. Assim, a soma (sum_\rho H_c(\gamma_\rho)), com multiplicidades, converge absolutamente. Isso também permite agrupar os zeros por blocos de altura sem ambiguidade de ordenação para esta classe de teste.

**Limite da conclusão:** isso prova convergência da soma, não a localização dos zeros nem positividade global.

## 4. Correspondência com a forma hermitiana

Como (f) é real,
[
F_c(-z)=\overline{F_c(\bar z)}.
]
Logo
[
H_c(\gamma_\rho)
=F_c(\gamma_\rho)\overline{F_c(\overline{\gamma_\rho})}.
]
Essa é a forma diagonal hermitiana relevante; quando (\gamma_\rho) é real, reduz-se a (|F_c(\gamma_\rho)|^2). Quando (\gamma_\rho) é complexo, não se deve substituir o produto por um módulo quadrado no mesmo argumento.

A identidade algébrica é válida sem assumir RH. O que ela não faz é transformar positividade numa subfamília finita em positividade para todas as funções teste.

## 5. Conferência da normalização da fórmula

A forma usada no documento de derivação é
[
\sum_\rho H_c(\gamma_\rho)=H_c(i/2)+H_c(-i/2)
+\frac{1}{2\pi}\int_{\mathbb R}H_c(r)
\left(\operatorname{Re}\psi(1/4+ir/2)-\log\pi\right)dr
-2\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}g_c(\log n).
]

A conferência com a versão publicada que escreve o termo arquimediano como integral de (\operatorname{Re}\Gamma'/\Gamma) mais (-\log\pi\,\widehat H(0)/(2\pi)) dá:

- **Termo arquimediano:** (g_c(0)=(2\pi)^{-1}\int H_c(r)dr), então o termo (-g_c(0)\log\pi) é absorvido no peso (\Omega(r)=\operatorname{Re}\psi(1/4+ir/2)-\log\pi).
- **Primos:** para a convenção (g_c(u)=(2\pi)^{-1}\int H_c(r)e^{-iur}dr), a transformada de Fourier normalizada em (2\pi) avaliada em (\log n/(2\pi)) é (2\pi g_c(\log n)). O coeficiente (-1/\pi) na versão publicada converte-se em (-2g_c(\log n)), como no código.
- **Polos:** como (H_c) é par, (H_c(i/2)=H_c(-i/2)); a soma dos dois termos é (2H_c(i/2)). As condições lineares (F_c(i/2)=F_c(-i/2)=0) são suficientes para anular ambos. São suficientes, não necessárias para anular o produto: exigir ambos define um subespaço conveniente, mais restrito que o conjunto em que apenas um dos fatores zera.
- **Suporte/primos:** (g_c) tem suporte em ([-2L,2L]). Para (L=0{,}45), (2L=0{,}9), de modo que somente (n=2) satisfaz (log n\le0{,}9); (n=3) já fica fora. O corte implementado por aritmética de bolas é coerente com o suporte.

A derivação por convenções de Fourier deve continuar citando uma única versão de referência com todos os fatores explícitos. A concordância numérica é um teste cruzado, não substitui essa derivação.

## 6. Matriz e escopo do certificado

O código forma (W=G+P) e depois (A=B^{\mathsf T}WB). A certificação intervalar da positividade definida de (A) estabelece positividade para toda matriz real simétrica dentro do enclosure certificado e, portanto, para a matriz exata representada pelo procedimento, se todas as hipóteses do certificado forem satisfeitas.

Há três níveis que devem permanecer separados:

1. **Analítica da fórmula:** a família de teste parece admissível pelas estimativas acima; convém obter uma revisão independente formal.
2. **Computação finita:** os enclosures e a positividade da matriz são resultados computacionais rigorosos sobre os parâmetros certificados.
3. **RH:** não segue da positividade de uma única matriz nem de uma família finita de dimensão (N-2). O critério global requer positividade para a classe completa de funções teste admissíveis.

## 7. Pontos a resolver antes de ampliar a pesquisa

1. Fazer uma revisão independente da prova de cauda em `TAIL_PROOF.md`, incluindo as identidades da digama, constantes numéricas e propagação da cota pela projeção intervalar.
2. Confirmar no log final da CI que o teste lento e todas as auditorias intervalares concluíram com sucesso.
3. Verificar que o certificado intervalar cobre a matriz exata e a transformação (B) com os mesmos critérios de inclusão de bolas usados no relatório.
4. Manter explícita a distinção entre a positividade desta família finita e a positividade global do critério de Weil.
5. Não aumentar (N), (L) ou (R) como substituto de uma prova; qualquer ampliação futura deve responder a uma pergunta matemática concreta.

## Referências externas para a revisão

- E. Bombieri, *Remarks on Weil’s quadratic functional in the theory of prime numbers, I*, seção 2, fórmula explícita de Guinand–Weil: https://www.bdim.eu/item?fmt=pdf&id=RLIN_2000_9_11_3_183_0
- *Bounding the log-derivative of the zeta-function*, Proposition 5 (enuncia a fórmula e suas hipóteses analíticas): https://doi.org/10.1007/s00209-021-02820-9
- Bombieri–Lagarias, *Complements to Li’s Criterion for the Riemann Hypothesis*: https://doi.org/10.1006/jnth.1999.2392
- Formulação de convenções de Fourier e zeros possivelmente complexos, discussão matemática: https://mathoverflow.net/questions/62816/the-guinand-weil-explicit-formula-without-entire-function-theory

*Este documento é uma auditoria interna de pesquisa, não uma prova de RH nem um parecer formal independente.*
