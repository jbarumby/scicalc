# NUMERIA ES8-P83 — Auditoria adicional
## Execução verificada
- GitHub Actions: https://github.com/jbarumby/scicalc/actions/runs/42422317716
- python-flint 0.9.0, 256 bits, **75/75 verificações**: 50 Bessel (j/i), 5 incrementos pares da digama, 20 identidades para nós/pesos de Gauss–Legendre.
- As referências Bessel/digama usam séries de **Fraction** com limites de resto, não chamadas ao FLINT. O teste avalia interseção de intervalos.
- A comparação de pesos usa a equação de Legendre e o valor 2 para a soma dos pesos.
- Teste adicional externo de 21 verificações locais sem FLINT e segundo cálculo racional de caudas/Schur.
## Resultados condicionais
- a=0,30 (24 modos): forma limitada inferiormente por 0,0063447, condicional à identificação de Weil e às margens finitas da ES8-P81.
- a=0,34 (26 modos): 0,00135877, nas mesmas condições.
- a=0,40 (44 modos): 0,00014279, nas mesmas condições.
O lema do majorante adaptativo decorre da identidade de Rodrigues
j_p(z) = z^p/(2^(p+1)p!) ∫_{-1}^1 (1−t²)^p exp(i z t) dt,
logo |j_p(z)| ≤ |z|^p cosh(|Im z|)/(2p+1)!!.
As caudas do operador são limitadas independentemente por primeiro termo e série geométrica, preservando os blocos de Legendre por partes e o primo 2 na janela a=0,40.
## Limites
Os testes de componentes não auditam formalmente FLINT nem demonstram a equivalência universal com a forma de Weil. Não provam a Hipótese de Riemann. O PR #17 permanece em rascunho.
