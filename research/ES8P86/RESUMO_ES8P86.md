# NUMERIA–RIEMANN ES8-P86 — auditoria ponta a ponta da fórmula explícita

**Data:** 2026-11-07. **Resultado:** consistência numérica corroborada em três representações, não certificação universal de Weil/RH.

## Representações independentes

Para g real de suporte compacto em [-a,a], G(t)=∫g(x+t)g(x)dx, N=G(0),
U=∫e^(x/2)g(x)dx, V=∫e^(-x/2)g(x)dx, com m(r)=Re ψ(1/4+ir/2)−logπ:

    Q(g)=2UV−2(log2)/sqrt(2) G(log2) 1_{2a>log2}
           +(1/2π)∫_{ℝ}m(r)|ĝ(r)|²dr.

A identidade integral da digama fornece **outra expressão**, no espaço real:

    (1/2π)∫_{ℝ}m(r)|ĝ(r)|²dr
     = −(γ+logπ)N
       + ∫_0^∞ [e^(−t)N−e^(−t/4)G(t/2)]/(1−e^(−t))dt.

A parte t>4a se integra exatamente como −N log(1−e^(−4a)).

## Testes locais concluídos

- 5 funções genuinamente C_c^∞ feitas com exp(1−1/(1−(x/h)²)) em |x|<h, zero fora: **maior diferença** entre integral real e espectral <7×10^(-12) em ponto flutuante. Nas janelas 0,30, 0,34 e 0,40, e em dois testes de suporte duplo separados por log2.
- 5 funções polinomiais (1−(x/h)²)^4_+; **atenção: estas são C³, NÃO C_c^∞**: diferenças <3×10^(-12).
- Para h=9/200 e dois suportes em ±log2/2: identidade exata G_±(log2)=±8192/303875, logo o termo do primo 2 é ±[2log2/sqrt2]*(8192/303875).
- Terceira verificação: soma de 350 zeros produzidos por mpmath.zetazero, aplicada às funções polinomiais; nas funções centrais a diferença do valor aritmético é ~10^(-15), enquanto nos dois suportes duplos há cauda remanescente ~4×10^(-10).
- **20/20 testes de regressão passaram**, com resultados em JSON e scripts no pacote ES8-P86 entregue na conversa.

## Posição em relação à literatura

O Teorema 12 de Bombieri (2000), §12, já garante positividade para C_c^∞ em intervalos de comprimento menor que log2. Consequentemente os suportes [-0,30,0,30] e [-0,34,0,34] pertencem ao regime conhecido. A janela [-0,40,0,40] ultrapassa o limiar e requer tratar o primo 2.

**Limitações:** as verificações usam SciPy/mpmath em ponto flutuante; não são bolas Arb nem comprovam uma cota rigorosa da cauda de todos os zeros. As observações não estabelecem RH. A identificação da forma de Weil ainda requer prova completa de domínio, normalização e avaliação intervalar independente.

Fonte: Enrico Bombieri, *Remarks on Weil's quadratic functional in the theory of prime numbers, I* (2000), https://www.bdim.eu/item?id=RLIN_2000_9_11_3_183_0.

**Próxima porta ES8-P87:** certificar em Arb as duas representações para exemplos C_c^∞, dar erro explícito da cauda de Fourier, preservar as três bases históricas e o tratamento exato do primo 2.
