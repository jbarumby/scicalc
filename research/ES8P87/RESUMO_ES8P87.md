# NUMERIA–RIEMANN — ES8-P87 | Auditoria de translações log 2 / log 3

Dados fornecidos: matriz 3×3 para b(x)=exp(1−1/(1−(x/h)^2)) no interior de |x|<h, zero fora, com h=0,4, centros 0, log 2 e log 3.

## Achados

- Reexecução original reproduziu todas as matrizes do JSON fornecido.
- Segunda implementação em NumPy/SciPy, usando Gauss–Legendre de ordem fixa **sem cache por arredondamento de s**: maior diferença na matriz total da fonte ≈1,91×10⁻¹². Testes de resolução 160×192 e 256×384 nós diferiram na matriz em menos de 10⁻¹⁴.
- Conferência do termo arquimediano por integração no espaço de Fourier com digama: quatro deslocamentos, diferenças <1,4×10⁻¹³.
- Os índices de primos ativos são 2,3,4,5 (o 4 é potência prima); nenhum n>6 alcança o suporte e Λ(6)=0.
- Autovalores da matriz ordinária de coeficientes: 0,0002715057167; 0,0019226064407; 0,0074665927192 (exploratórios).
- Como as translações não são ortonormais, o problema coerente com a norma L² usa **M c=λ G c**, com G_ij=∫b(x−x_i)b(x−x_j)dx. Autovalores generalizados numéricos: **0,0006373873032; 0,0063514379487; 0,0168190967774**.
- Suporte global das combinações está em [−0,4, log3+0,4] (comprimento ≈1,8986), **não é** a mesma janela [-0,4,0,4] das fases anteriores.
- **14 testes locais passaram**.

## Certificado racional condicional — meta para Arb

A matriz M₀, obtida arredondando a matriz fornecida a 12 casas, satisfaz **M₀ − (1/4000)I ≻ 0** por Sylvester com frações exatas.

**Se**, por uma auditoria rigorosa ainda não executada, cada coeficiente verdadeiro M_ij for incluído em M₀_ij ± 1/100000, então

||M−M₀||₂ ≤ 3/100000,
λ_min(M) > 11/50000 = 0,00022.

Como 0≤b≤1 em suporte de comprimento 4/5, ||G||₂≤12/5. Logo λ_min(M,G)>11/120000≈9,17×10⁻⁵ **condicionalmente**.

O arquivo original utiliza cache por grade 10⁻¹² para a correlação R(s); perto t=0, a expressão arquimediana sofre cancelamento e deve ser reescrita de modo estável. A nova quadratura float64 corrigiu isso, mas **não** fornece erro rigoroso.

## Próxima etapa ES8-P88

Construir inclusões Arb das seis entradas independentes e de G, separando interiores analíticos e bordas da bump. Não aplicar diretamente à bump C∞ os contornos de Bernstein sem controlar as singularidades em ±h.

Fonte da fórmula explícita: Bombieri (2000), §12, https://www.bdim.eu/item?id=RLIN_2000_9_11_3_183_0.

**Status:** exploração computacional reproduzida e teste racional condicional; não há certificado intervalar de positividade e não é prova da Hipótese de Riemann.
