# NUMERIA–RIEMANN — ES8-P88: oito translações, Gram e alvo racional

## Arquivos examinados
- Avaliar avanços para resolver a hipótese de Riemann.zip (8 itens, varredura h=0,1…1,0);
- matriz_weil_8_centros_h040.json;
- certifica_matriz_intervalar.py (aceita entry_intervals rigorosos, ausentes do JSON enviado).

## Descobertas reproduzidas
- Família suave b_h(x)=exp(1−1/(1−(x/h)²)) para |x|<h, zero fora.
- Centros log(1), log(2), log(3), log(5), log(7), log(11), log(13), log(17).
- h=0,4: matriz 8×8 reproduzida bit a bit pelo programa original; λmín(M_float64)=4,7524966561986686×10⁻⁵.
- Segunda implementação por Gauss–Legendre fixo sem cache quantizado, ordens 160×240 e 240×320: discrepância máxima ≈3,601×10⁻¹² na matriz em relação ao original; os dois cálculos fixos diferem no máximo ≈5,44×10⁻¹³. Isto NÃO é cota rigorosa.
- Com Gram das oito translações, menor λ generalizado (M,G) ≈9,022228952484222×10⁻⁵ para h=0,4.
- Dez valores h=0,1…1,0 reproduzidos, positivos numericamente. κ₂(G) cresce de ≈1,002 para ≈580,745.
- A coluna original λmín(M)/||b||² não é λmín(M,G) quando há sobreposição entre translações.
- O suporte de uma combinação de oito centros é [−0,4, log(17)+0,4], comprimento ≈3,6332; NÃO equivale à janela [-0,4,0,4] das fases de Legendre.

## Teorema condicional exato (novo)
Seja A₀ a matriz de h=0,4 arredondada a 14 casas decimais. Verificou-se por Bareiss inteiro e Sylvester exato (8 menores principais) que

A₀ − (47/10⁶)I ≻ 0.

**SE** uma integração independente rigorosa demonstrar |M_ij−(A₀)_ij|≤10⁻⁶ para as 36 entradas independentes, então

||M−A₀||₂ ≤ ||M−A₀||_F ≤8·10⁻⁶,

e M ≽ (39/10⁶)I ≻0. Além disso, ||G||₂≤tr G≤32/5, logo λmín(M,G)≥39/6.400.000, condicionalmente.

A demonstração algébrica é rigorosa, mas **os limites sobre as integrais ainda NÃO foram obtidos**. O certificador fornecido recebe INPUT_ERROR quando se passa o JSON de float64, que não inclui entry_intervals.

## Execução ES8-P88
- 14/14 testes locais aprovados.
- Reexecução integral da varredura com Gram; calculadora fixa independente; script de certificado racional de tolerância.
- Pacote original, relatório, resultados JSON e manifesto SHA-256 preservados no arquivo da conversa NUMERIA_ES8P88_OITO_CENTROS_AUDITORIA_2026-11-07.zip.
- Sem novas inclusões Arb; sem merge na main; sem prova da Hipótese de Riemann.

## Porta ES8-P89
Certificar rigorosamente as 36 entradas independentes de M a erro absoluto ≤10⁻⁶ usando aritmética intervalar, incluindo as bordas não analíticas da bump e a regularização da digama perto t=0. Só então certificar a positividade **desta** matriz finita.

Referência primária: Bombieri (2000), Remarks on Weil's quadratic functional in the theory of prime numbers, I, §12.
