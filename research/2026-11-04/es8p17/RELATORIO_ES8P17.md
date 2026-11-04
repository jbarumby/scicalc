# NUMERIA–RIEMANN — ES8-P17 | Revisão estrutural e orçamento de cauda após projeção

**Data:** 04/11/2026. **Status:** auditoria racional exata; teoremas Lean de aritmética condicionais; identidade integral E7 e identificação com Weil ainda pendentes.

## 1. Documento enviado e referências

O arquivo `deepseek_latex_20261104_dd6431.tex` propõe trabalhar com formas fechadas/Dirichlet, cita positividade em janela limitada, e destaca a necessidade de controlar `B_N`. **Apoiado:** limite finito ≠ universal, a forma de Galerkin exige identidade e domínio, o projetor amplifica erro, e a necessidade de uma estrutura para todos os suportes. **Não provado pelo arquivo:** que as matrizes históricas sejam compressões de uma única forma fechada; uma margem inferior uniforme; que a extensão de autoadjuntos identifique zeros; a convergência Hilbert–Schmidt. O artigo de M. Suzuki, *Weil's quadratic form via the screw function* (arXiv:2606.09096v3, setembro/2026), trata do operador de Friedrichs e apresenta como **conjectura** uma identificação espectral no limite. O trabalho de Xuefeng Zhu, arXiv:2608.24827v2, reivindica uma certificação de janela `[-0.8,0.8]` e análise de deterioração da margem; não implica RH. Essas alegações externas não foram reprovidas numericamente nesta etapa.

**Correção quantitativa essencial:** o documento chama a soma dos quadrados dos erros de 8×8 de `≈1,2×10^-11`; porém, o certificado ES8-P16 registra exatamente `12382243/200000000000 = 0.000061911215`, ou seja, `6,1911215×10^-5`. A norma de Frobenius seria a raiz quadrada desse número, cerca de `0,0078684`, não `10^-11`. É menor que `10^-2`, mas três dimensões não provam convergência na dimensão infinita.

Outra distinção: uma forma **fechada e limitada inferiormente não precisa ser positiva**. KLMN assegura fechamento/semilimitação sob hipóteses de relativa limitação, não a desigualdade `Q≥0` automaticamente. Para suporte compacto, termos de primos fora de `log(p^k)≤2L` **desaparecem** em vez de gerar uma cauda numérica de primos a truncar.

## 2. Projeção exata da família L=9/20

O código original (`es8p_projection.py`, função `pole_coefficient_exact`) usa, para cada modo `k≥3`, pivô `p=1` se k ímpar e `p=2` se k par, `X=π²` e `L=9/20`:

`B[p-1,k-3]=-(k/p)(p² X+L²)/(k² X+L²)`, `B[k-1,k-3]=1`.

As outras entradas são zero. Separando as duas linhas de pivô, `B_N=[C_N;I_{N-2}]` e

`||B_N x||²=||x||²+|Σ_{k ímpar} c_k x_k|²+|Σ_{k par} c_k x_k|²`, onde `c_k=|(B_N)_{p,k-3}|`.

A fórmula é **exata** e os suportes das duas somas são disjuntos. Pelo Cauchy–Schwarz, vale

`||B_N||² = 1+max(Σ_{k≥3,ímpar}c_k², Σ_{k≥4,par}c_k²)`.

Como `X=π²≥9` e `k>p`, a fração satisfaz

`c_k ≤ [p²+L²/9]/(pk)`.

Assim `c_k ≤ (409/400)/k` para k ímpar, e `c_k ≤ (1609/800)/k` para k par. Usando

`Σ_{k≥3} 1/k² ≤ Σ_{k≥3}1/[k(k-1)] = 1/2`

e, para `k=2m,m≥2`,

`Σ_{m≥2}1/(2m)² ≤ (1/4)Σ_{m≥2}1/[m(m-1)]=1/4`, obtemos

`Σ_odd c_k² ≤ 167281/320000 < 1`;

`Σ_even c_k² ≤ 2588881/2560000`.

Logo para todo N≥3,

`||B_N||² ≤ 5148881/2560000 ≈2.01128 < 21/10`.

A sequência também é crescente por inclusão de novos termos e possui limite finito; portanto esta **projeção concreta** se estende continuamente como operador de `ℓ²` para `ℓ²`, sob a identificação usual entre coordenadas. Isso NÃO demonstra que a forma infinita de Weil esteja definida ou seja positiva.

## 3. Cauda projetada com N arbitrário

A ES8-P15 registrou, sob a hipótese integral analítica E7 e com `R(N)=200N²`, a cota da soma por linha da cauda bruta:

`δ_raw(N) ≤ 16/(3645N)+2/(30375N³) ≤ 406/91125` para N≥1.

Para a matriz de cauda simétrica real `E_N`, `||E_N||op ≤ δ_raw(N)`. A projeção satisfaz

`||B_Nᵀ E_N B_N||op ≤ ||B_N||²op ||E_N||op ≤ (21/10)(406/91125)=1421/151875 < 1/100`, com margem exata `391/607500`.

A validade de `R(N)≥max(1,2w_max)` também resulta de `π<4`: `2w_max=Nπ/L < (80/9)N ≤200N²` para N≥1.

**Esta conclusão resolve somente o orçamento da CAUDA ARQUIMEDIANA projetada, não o orçamento de erro total após a projeção**. Os erros de digama, interpolação, normalização e avaliação dos termos aritméticos exigem análise adicional. A margem `1/20` da matriz central continua não estabelecida para todo N.

## 4. Evidências, testes e próximo obstáculo

O código `audit_es8p17.py` verifica todas as constantes, reconstrói `BᵀB` diretamente em frações racionais para N até 20, verifica a identidade de Gram por paridade, os coeficientes no intervalo racional `π²∈[9,10]`, a monotonicidade em X, a cota uniforme e os testes adversariais. Foram 13 testes aprovados localmente.

Em Lean, `ES8P17.lean` prova as identidades racionais e a transferência do limite sob **hipóteses explícitas sobre o erro bruto e a norma de projeção**. Não formaliza ainda a desigualdade de Cauchy–Schwarz para a construção específica `B_N`, nem a cota E7 no integral analítico. O código será considerado compilado apenas após execução do GitHub Actions.

**ES8-P18 recomendada:** definir uma candidata global à forma quadrática em domínio comum; provar o vínculo Galerkin para todos os N; formalizar a soma infinita dos coeficientes de B; investigar o menor autovalor das matrizes centrais sem supor uma margem `1/20` uniforme. Verificar separadamente a rota Suzuki/Friedrichs e as restrições do suporte do NUMERIA.
