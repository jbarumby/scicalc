# NUMERIA–RIEMANN — V21: traço térmico, regularização espectral e Casimir

Status do auditor: TOY_LAPLACIAN_EXACT_RATIONAL_AUDITS_PASSED_OPERATOR_BRIDGE_OPEN.

Este é um modelo **auxiliar**, não um operador que realize Hilbert–Pólya, e não certifica a fórmula explícita de Weil nem a Hipótese de Riemann.

## Modelo e resultados matemáticos

Seja H=-d²/dx² no espaço L²(0,π), domínio H²(0,π) ∩ H¹₀(0,π), com condições de Dirichlet. Os autovalores são λ_n=n², n≥1; o operador é positivo e autoadjunto.

Para t>0, a função de partição térmica (traço de calor) é Θ(t)=Tr(exp(-tH))=Σ_{n≥1} exp(-tn²). Na região de convergência Re(s)>1/2, a zeta espectral é ζ_H(s)=Σ_{n≥1} n^(-2s)=ζ(2s). A transformada de Mellin satisfaz Γ(s)ζ_H(s)=∫₀^∞ t^(s-1)Θ(t)dt, para Re(s)>1/2. A justificativa usa o teorema clássico do laplaciano de Dirichlet, integral gama e intercâmbio válido na região de convergência; esses passos não estão formalizados em Lean.

O módulo AUDITAR_E8K_V21_ESPECTRAL.py verifica com Fraction (sem float nas decisões):

1. Três envoltórias de ζ_H(s) para s=1,2,3, com cutoff N=32 e cotas por comparação integral. As identidades ζ(2)=π²/6 e ζ(4)=π⁴/90 são referências matemáticas externas, não resultados novos.
2. Vinte e quatro igualdades exatas da integral gama por modo, em n=1...8, s=1,2,3.
3. Três envoltórias de Θ(t), para t=1/4,1/2,1, com N=6 e N=10 e largura estritamente menor no corte maior. Exponenciais têm cotas racionais por Taylor alternado, com redução de argumento. Para N≥2, a cauda satisfaz:
   exp[-t(N+1)²] ≤ Σ_{n>N}exp(-tn²) ≤ exp[-t(N+1)²] / (1-exp[-t(2N+3)]).
4. A soma de Abel S(x)=Σ n exp(-nx)=exp(x)/(exp(x)-1)² possui série Laurent S(x)=x^-2 -1/12 +O(x²). O coeficiente -1/12 é recomposto por álgebra de séries racionais, confrontado com -B₂/2=-1/12 pela recorrência de Bernoulli.
5. Em um campo escalar unidimensional Dirichlet, E_Casimir=(πℏc/2a)ζ(-1)=-πℏc/(24a). O código certifica o fator **racional** -1/24; a identificação com a continuação analítica e as escolhas de renormalização/condições de contorno são teoremas e convenções externas.

## Execução

Na raiz do repositório:

    python v21/AUDITAR_E8K_V21_ESPECTRAL.py --out v21/RESULTADO_E8K_V21.json
    python v21/TESTAR_E8K_V21_ESPECTRAL.py

O workflow auditoria-e8k-v21-espectral.yml gera duas vezes o JSON e compara byte a byte, roda mutações e publica a testemunha com hashes SHA-256. O resultado será fixado no GitHub junto com um manifesto SHA256SUMS quando o pipeline confirmar os bytes.

## Fontes e separação dos níveis de evidência

- Bourgade, Paul; Keating, Jonathan P. (2013), Quantum Chaos, Random Matrix Theory, and the Riemann ζ-function, Chaos: Poincaré Seminar 2010, p.125–168; DOI https://doi.org/10.1007/978-3-0348-0697-8_4. A cópia enviada contém **somente as duas primeiras páginas**. A conexão com a fórmula de Weil é uma analogia de fórmulas de traço, não construção do operador Hilbert–Pólya.
- Shen, Hui; Liu, Shuqi; Quan, Ying; Li, Yingying; Wang, Tongbiao; Yu, Tianbao; Tao, Xixia (2025), The progress of the Casimir effect: from theory to application, Nanotechnology 36, 482002; DOI https://doi.org/10.1088/1361-6528/ae1b0e. Revisão física de Casimir que motiva a análise; o exemplo 1D é matemática-padrão reconstruída independentemente, não um teorema extraído integralmente desse artigo.
- R. S. Laugesen, Spectral Theory of Partial Differential Equations, capítulo 15, Spectral zeta function, https://publish.illinois.edu/rlaugesen/files/2023/07/595Lectures.pdf — fonte matemática pública para o modelo do laplaciano do intervalo e transformada de Mellin.

## Pendências

Não construímos um operador com os zeros não triviais de Riemann como espectro. Não demonstramos uma identidade do seu traço com os termos de primos, polos e parcela arquimediana da fórmula de Weil. A V21 não reexecuta as cotas de segunda derivada, 1600 subintervalos de quadratura da Fase 4, nem sua projeção LDL. O resultado negativo de posto ≤1 na V20 para o kernel Gibbs direto ainda impede uma identificação imediata com Weil.
