# NUMERIA E8-K V20 — Gibbs, Euler, Yang–Lee e Gram

**Status científico:** FINITE_GIBBS_EULER_YANG_LEE_GRAM_AUDITED_WEIL_BRIDGE_OPEN.

V20 é uma linha auxiliar independente, não altera a certificação condicional da fórmula explícita. Todos os cálculos decisórios implementados nos módulos usam a aritmética racional exata de Python (Fraction); valores decimais no JSON são arredondados para fora.

## Executar

Na raiz:

    python v20/AUDITAR_E8K_V20.py --out v20/RESULTADO_E8K_V20.json
    python v20/TESTAR_E8K_V20.py
    (cd v20 && sha256sum -c SHA256SUMS_V20.txt)

GitHub Actions gera o JSON duas vezes, compara bytes e publica o resultado como artefato. Não copia PDFs protegidos para o repositório.

## Módulo 1 — Gibbs–Euler

Na região real beta>1:

    Z(beta) = sum n^-beta = zeta(beta) = product_p (1-p^-beta)^-1.
    -d/d(beta) log Z(beta) = sum_n Lambda(n)/n^beta
                           = sum_p log(p)/(p^beta-1).

Testamos beta=2,3,4: série finita + cauda integral contra produto de Euler finito + majorante por união dos múltiplos de primos omitidos. Para beta=3,4 verificamos a expansão von Mangoldt de duas formas com logaritmos intervalares racionais. São 99 identidades geométricas exatas para potências de primos. Euler e a integral são teoremas clássicos tomados como premissas: não afirmamos formalização Lean.

## Módulo 2 — Yang–Lee finito

Para um modelo auxiliar de um único modo de bósons de energia zero, limitando sua ocupação a K, Xi_K(z)=1+z+...+z^K. A identidade (z-1)Xi_K(z)=z^(K+1)-1 é verificada coeficiente a coeficiente em inteiros. Seus zeros são as raízes da unidade de ordem K+1 exceto z=1, e ficam fora do eixo real positivo em K finito. A distância da raiz mais próxima a 1 é menor que 44/[7(K+1)]. Esse *pinching* não implica, sozinho, transição termodinâmica, nem identifica zeros Yang–Lee com zeros da zeta. NÃO usamos o polinômio interagente de X.-Z. Wang.

## Módulo 3 — Gram racional e a obstrução de posto

Para f_i(n)=n^-i e pesos de Gibbs normalizados w_n=n^-beta/sum_{m<=M}m^-beta, seis matrizes são reconstituídas por LDL racional com todos pivôs positivos. A prova Gram é c^T K c = sum w_n (sum c_i f_i(n))^2 >=0.

**Resultado negativo relevante para a ponte com Weil:** as funções seno do NUMERIA com suporte x em [0,2L], L=9/20, avaliadas diretamente em x=log n, dão log 2 < 2L=0,9 < log 3. Só n=1,2 estão no suporte, mas f_k(log 1)=f_k(0)=sin 0=0. Apenas n=2 pode contribuir. O Gram é um produto externo de posto <=1, e B^T K B também tem posto <=1. Logo essa IDENTIFICAÇÃO DIRETA não pode reproduzir a positividade definida de uma matriz Weil projetada de dimensão >=2. Outras possíveis correspondências via operadores/transformações continuam abertas, mas exigem teoremas novos.

## Referências e revisão crítica

- M. Howard Lee, “Polylogarithms and Riemann’s ζ Function”, Phys. Rev. E 56 (1997), 3909; DOI https://doi.org/10.1103/PhysRevE.56.3909.
- X.-Z. Wang, “Critical temperature of Bose–Einstein condensation of a dilute Bose gas”, Physica A 341 (2004), 433–443; DOI https://doi.org/10.1016/j.physa.2004.04.125. O coeficiente 4,92506 deve ser comparado criticamente com valores de outros modelos/estudos, não assumido como consenso para gases homogêneos.
- B. Kerr, O. Klurman, J. Thorner, “Zeros of L-functions and large partial sums of Dirichlet coefficients”, Mathematische Annalen (2026), DOI https://doi.org/10.1007/s00208-026-03542-1; a equação (1.8) generaliza os coeficientes de von Mangoldt em derivadas logarítmicas.
- M. J. Bowick, “Quantum Statistical Mechanics and Multiplicative Number Theory”, Thermal Field Theories (1991), p.315–323. Seu conteúdo integral não foi verificado.
- Detalhes dos hashes dos arquivos recebidos: REFERENCIAS_V20.md.

## Pendências científicas

V20 NÃO executa os majorantes de derivada da Fase 4, não certifica as 1.600 células de quadratura nem as matrizes LDL Weil projetadas, não formaliza os lemas clássicos em Lean, não demonstra uma ponte Gibbs–Weil positiva de dimensão geral, e **não prova a Hipótese de Riemann**.
