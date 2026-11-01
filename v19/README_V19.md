# NUMERIA–RIEMANN — E8-K V19: replay das cotas arquimedianas

**Status científico:** P5_P6_CONSTANTS_AND_TAIL_ENCLOSURES_REPLAYED_CONDITIONALLY.

V19 reexecuta a aritmética racional exata das Proposições 5 e 6 de `fontes/TAIL_PROOF.md` contido no ZIP E7 original (SHA fixado). A prova `referencias/TAIL_PROOF_ENVIADO_2026-11-01.md` é uma versão distinta: não contém a Proposição 6; ambos os hashes são conferidos separadamente.

## Resultados verificados

- 30 pares independentes de mesma paridade (i <= j <= 10), e 25 pares de paridade diferente (identicamente zero por fase).
- P5 (conservadora): |T_ij| <= (16*w_i*w_j/(pi*L))*(8+log(3*R)+1/3)/(3*R^3).
- P6 (do E7): |T_ij| <= [4*w_i*w_j/(pi*L*(1-eps_i)*(1-eps_j))]*(log R+3.14+1/3)/(3*R^3), eps_k=(w_k/R)^2.
- Rejeita qualquer cálculo se R < N*pi/L (testado: L=0.45, R=160, N=30 é inválido); N=20 nesse R satisfaz a pré-condição, o que NÃO significa positividade da matriz.
- Coteja as cotas históricas G11/G22 da Proposição 6 com erro absoluto <1e-22, sem alegar igualdade byte a byte.
- Pi: série alternada de Machin; log: série de atanh e cauda racional; gamma: limites clássicos harmônico-log **assumidos**; cosh: soma parcial inferior; decisões com Fraction, sem float.

## Uso

Executar na raiz do repositório:

    python v19/AUDITAR_E8K_V19_CAUDA.py --e7 NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip --submitted-proof referencias/TAIL_PROOF_ENVIADO_2026-11-01.md --out v19/RESULTADO_E8K_V19.json
    python v19/TESTAR_E8K_V19.py

O workflow audita, testa adulterações de origem e gera duas vezes o JSON para comparação byte a byte. O resultado é publicado como artefato GitHub Actions; não modifica o manifesto V14 da raiz.

## Limites e pendências

Os lemas analíticos (identidades da digamma, integral de Binet e as desigualdades da constante de Euler) não foram formalizados integralmente em Lean na V19. A V19 também não reexecuta as 1.600 células de quadratura da Fase 4, as cotas complexas de Cauchy, as matrizes projetadas ou a fórmula explícita global de Weil. Nenhuma positividade em dimensão arbitrária nem RH foram demonstradas.
