# NUMERIA–RIEMANN ES8-P242 — reprodução independente das duas caixas finitas

Data: 2026-11-18. **Escopo preciso:** reproduzir, sem `python-flint/Arb`, a positividade das DUAS caixas simétricas intervalares `377×377` **já montadas** em `P237_FINITE_LDL_2048.json.gz`, incluindo a transferência centro–Weyl. O presente cálculo não reavalia as integrais, os 1008 painéis GL60, a fórmula explícita de Weil, ou a redução analítica de Zhu. Não há alegação sobre RH.

## 1. Entradas auditadas

- Fonte: `P237_FINITE_LDL_2048.json.gz`, SHA-256 `2695768a3cca11237dd63cb40c0b6c34427b390329e5ff2ec089e31111c7898a`.
- Fonte de comparação: `P237_CENTER_WEYL_2048.json.gz`, SHA-256 `272caafee875c2d5301782ed56bcb4d895df98ee0551a2e878d99577f44156d3`.
- Auditoria-sombra: `P237_SHADOW_FULL_AUDIT.json`, SHA-256 `e51a28fe0180c1a6215be23b0f268d570bed3c748210dcca0a0d9ff3f59bea26`.
- Matrizes agregadas de quadratura: `P237_FULL_C_even.json.gz` / `odd.json.gz`, SHA-256 registrados em `p242_verify_exact.py` e conferidos no replay.

Cada entrada armazenada tem uma codificação de bola diádica `[ [m,e], [r,k] ]`, equivalente ao centro exato `m·2^e` e ao raio `r·2^k ≥0`. O script `export_centers.py` extrai **somente centros**, inverte a indexação triangular de linhas para colunas e vincula os dados à fonte por SHA-256.

## 2. LDL completamente independente

O código `p242_mpfr_interval_ldl.c` (GCC + MPFR + GMP) implementa um LDL por atualizações do complemento de Schur. Ele lê inteiros decimais e expoentes binários; as conversões e operações são calculadas em **precisão 1024 bits com arredondamento dirigido** (`MPFR_RNDD` para o limite inferior, `MPFR_RNDU` para o superior).

Escrevendo `C` para a matriz real simétrica **formada pelos centros diádicos de P237**, cada pivô exato `d_k` está contido no intervalo independente `[d_k^-,d_k^+]`. A cada passo, a atualização `S_{ji} ← S_{ji} - S_{ik}S_{jk}/d_k` é encerrada intervalarmente e os fatores `L_{ik}` são guardados para a etapa posterior. Verificou-se **`d_k^- > 0` para todos os 377 índices em cada paridade**. Assim, pelo critério LDL, ambos os centros `C_even`, `C_odd` são positivos definidos.

Todas as 754 caixas de pivôs independentes têm interseção não vazia com as caixas publicadas pelo Arb, mas essa comparação não foi usada para decidir a positividade. Foram preservados os intervalos MPFR em forma de quatro inteiros/expoentes por pivô, sem arredondamento decimal para decidir sinais.

## 3. Inversão independente e passo de Weyl

Após a fatoração, o mesmo código C calcula uma inclusão para cada coluna de `C^{-1}` resolvendo, com arredondamento dirigido, os sistemas `L y=e_j`, `D z=y`, `L^T x=z`. A soma dos quadrados dos **limites superiores dos módulos** das entradas produz `ν ≥ ||C^{-1}||_F²`.

O auditor Python **não usa os raios nem as cotas de inversa do certificado centro–Weyl anterior**. Em vez disso, ele soma com `fractions.Fraction` os raios diádicos originais de cada uma das 142.506 entradas triangulares **do arquivo matriz**, repetindo as contribuições simétricas à soma de cada linha. Seja `δ = max_i Σ_j r_ij`. Para qualquer matriz simétrica `A` incluída nessa caixa, `A=C+E`, temos `||E||₂ ≤ ||E||∞ ≤ δ`.

As desigualdades exatas verificadas são

`δ² ν < 1` e `(δ + 1/10^18)² ν < 1`.

De fato, `C` SPD implica `λ_min(C) = 1/||C^{-1}||₂ ≥ 1/||C^{-1}||_F ≥ 1/√ν`. Por Weyl, `λ_min(A) ≥ 1/√ν − δ > 1/10^18`. Ambas as matrizes intervalares P237 satisfazem essa implicação **sob a correção do algoritmo MPFR e da codificação de entrada**.

## 4. Números de apresentação (não usados para decisões)

| Setor | Pivôs positivos independentes | δ, soma de raios exata (aprox.) | ν, limite superior da inversa (aprox.) | δ²ν (aprox.) |
|---|---:|---:|---:|---:|
| Par | 377/377 | 3.661198923268398 × 10⁻⁴⁷ | 4.356188589587664 × 10³⁴ | 5.839199655884685 × 10⁻⁵⁹ |
| Ímpar | 377/377 | 3.657549602819543 × 10⁻⁴⁷ | 3.684275778979701 × 10²⁸ | 4.928702223359698 × 10⁻⁶⁵ |

Os valores de prova são as linhas dyádicas de `P242_PIVOTS_{even,odd}_1024.txt`, a montagem original sujeita aos hashes e as comparações com `Fraction` executadas por `p242_verify_exact.py`. Todas as comparações são **exatas**, não baseadas em `float64`.

## 5. Distinções obrigatórias

**Comprovado computacionalmente, condicionado à implementação das operações MPFR:** cada matriz real simétrica nas duas caixas P237 montadas é SPD com margem `>10⁻¹⁸`, sem reutilização dos pivôs ou inversas calculados pelo Arb.

**Não reproduzido nesta P242:** o integrando Bessel/digama e cada uma das 1008 quadraturas, validade externa das cotas analíticas GL60, a identificação da forma de Weil com a matriz, a passagem ao operador Legendre–Bessel infinito (P239), e a igualdade ou desigualdade correta para o núcleo reduzido de Zhu (P240–P241). As afirmações de suporte fixo dependem dessas entradas analíticas separadas e ainda requerem revisão independente.

**Não provado:** positividade para suportes arbitrários e Hipótese de Riemann. O artigo de Xuefeng Zhu (arXiv:2608.24827v2) fornece a redução de janela em seus próprios termos; citá-lo não verifica automaticamente a configuração `L=0.805` do NUMERIA nem a integral que deu origem à matriz.

## 6. Reprodutibilidade

```
# GNU/Linux, python3, gcc, libgmp-dev e biblioteca MPFR 4.x (runtime libmpfr.so.6)
gcc -O3 -Wall p242_mpfr_interval_ldl.c -o p242_mpfr_ldl -l:libmpfr.so.6 -lgmp
python3 export_centers.py --source /CAMINHO/P237_FINITE_LDL_2048.json.gz --out-dir .
./p242_mpfr_ldl P242_CENTERS_even.txt P242_PIVOTS_even_1024.txt 1024 377
./p242_mpfr_ldl P242_CENTERS_odd.txt P242_PIVOTS_odd_1024.txt 1024 377
python3 p242_verify_exact.py --files /CAMINHO/DOS/ARQUIVOS_P237 --traces . --out P242_AUDIT_CERT.json
python3 -m unittest -v test_p242_mpfr.py
python3 -O -m unittest -q test_p242_mpfr.py
```

A versão C usa declaração manual da ABI pública MPFR porque os cabeçalhos `mpfr.h` não estavam instalados neste ambiente. A variante preferível para manutenção em Linux é compilar contra o `mpfr.h` do pacote de desenvolvimento, sem modificar o algoritmo. O código C e os dados de saída precisam ser revistos por terceiros antes de serem apresentados como um certificado matemático independente de auditoria externa.
