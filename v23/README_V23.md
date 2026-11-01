# NUMERIA V23 — majorantes independentes e suavização

Status: onze majorantes de segunda derivada B3 reproduzidos por fonte independente, com hipóteses clássicas explícitas. Não é prova de RH.

## 1. Verificador de derivadas

AUDITAR_V23_DERIVADAS.py lê o ZIP E7 pelo SHA-256 original e reconstrói todos os 11 majorantes de |f_k''(r)| (G11 e G22) com aritmética Fraction. Não importa os programas B3 históricos. Recalcula π por Machin, log por atanh com resto racional, γ por Euler–Maclaurin alternada, Ω via série real de digama e cauda integral, e Ω' e Ω'' via diferenciação e desigualdades globais. Trata os polos removíveis da transformação seno/cosseno usando representação integral; nas bandas afastadas do polo diferencia o quociente explícito.

Os principais insumos analíticos permanecem: identidades de digama, limite alternado de Euler–Maclaurin, representação de Fourier dos modos, critério de convergência para diferenciar a série e comparação integral. São hipóteses documentadas, não inteiramente provadas em Lean.

A regra do produto usada é:
  |(A² Ω/π)''| <= ((2a1²+2a0a2)|Ω| +4a0a1|Ω'|+a0²|Ω''|)/π.

Os 11 valores independentes ficam abaixo dos majorantes B3 exportados. A correção de arredondamento de 1/10^45 no erro G22 [96,160] da V22 continua necessária.

## 2. Suavização C∞ de funções-teste

Os modos φ1(x)=cos(πx/(2L))/sqrt L e φ2(x)=sin(πx/L)/sqrt L são nulos nas bordas x=±L, mas suas extensões por zero não são C∞.

Defina b(t)=exp(-1/t) para t>0 e b(t)=0 para t<=0. Defina η(t)=b(t)/(b(t)+b(1-t)) em 0<t<1, com η=0 em t<=0 e η=1 em t>=1. Para 0<δ<L, χδ(x)=η((L−|x|)/δ) é um corte suave, igual a 1 se |x|<=L−δ e zero se |x|>=L. Não apresenta singularidade em x=0 porque vale 1 em uma vizinhança desse ponto.

A aproximação satisfaz a cota conservadora:
  ||φk−χδ φk||_1 <= 2 δ² ||φk'||∞
  ||φk'||∞ <= kπ/(2L sqrt L), k=1,2.
Pela desigualdade Fourier, ||hatφk−hat(χδφk)||∞ <= kπ δ²/(L sqrt L).
O segundo programa gera limites racionais para δ=1/100 e 1/1000.

**Limite:** erro Fourier uniforme NÃO é cota global da forma de Weil com peso Ω(r), que exige decaimento, integração por partes e controle dos termos de polos/primos. Não foi provado que L=0,45 esteja dentro do suporte suficiente do Teorema 3 de Burnol.

## 3. Teste adversarial de simetrias

O polinômio exato P(s)=[(s−1/2)²+17/16]²−(s−1/2)²/4 satisfaz P(s)=P(1−s) e tem zeros em 1/4±i, 3/4±i. É contraexemplo à inferência de que a equação funcional e a conjugação sozinhas forçam zeros na reta crítica.

## Rodar

    python v23/AUDITAR_V23_DERIVADAS.py --e7 NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip --out v23/RESULTADO_V23_DERIVADAS.json
    python v23/AUDITAR_V23_SUAVIZACAO.py --out v23/RESULTADO_V23_SUAVIZACAO.json
    python v23/TESTAR_V23.py

O workflow reproduz os JSONs duas vezes e testa entradas inválidas, adulteração da origem e igualdade byte a byte.

## Próximas obrigações

Formalizar os lemas analíticos em Lean ou submetê-los à revisão externa; transportar as cotas de suavização à integral arquimediana inteira; repetir 15.754 enclausuramentos B5 e os pivôs LDL nas funções-teste admissíveis; demonstrar positividade universal. Nenhum desses passos foi substituído pela V23.

**A Hipótese de Riemann permanece aberta.**
