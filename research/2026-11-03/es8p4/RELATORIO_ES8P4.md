# NUMERIA–ES8-P4 — Reconciliação documental, controles algébricos e triagem de artigos

**Data:** 03/11/2026  
**Status matemático:** `DOCUMENTARY_AND_FINITE_V17_GATE_PASS`, `ORIGINAL_E7_V16_ENTRY_IDENTITY_NOT_ESTABLISHED`, `NO_RH_PROOF`.

## 1. Objetivo e resultado científico

A ES8-P3 construiu de novo uma família finita de matrizes reais simetrizadas a partir de (i) um bloco arquimediano com núcleo digama e (ii) o único termo de primo com suporte em `2L=0,9`, seguida da projeção de anulação dos polos. A próxima obrigação era **provar que esses objetos coincidem com os artefatos originais E7/V16**.

**Resultado desta fase:** a inspeção do GitHub confirmou a identidade dos **objetos Git** históricos, e o programa reexecutou as propriedades algébricas e intervalares das matrizes reconstruídas. **Não houve abertura dos bytes do ZIP original E7 nem cotejo das matrizes internas E7/V16, portanto a equivalência literal continua não demonstrada.** O teste estrito falha propositadamente quando esse requisito não é atendido.

## 2. Evidência criptográfica verificável no GitHub

Dois commits examinados pela API de árvores Git:

- `fe6bdb7463172ea7107cd050410fe0f94196de29` — registro histórico auditado de V13–V15.
- `2053129b03cffbce302b142af73df37e993f1613` — versão publicada da ES8-P3.

Em ambos os commits os caminhos `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip` e `numeria_e8k/v15/E7_original.bin` apontam para o **mesmo blob Git SHA-1** `5b6abb59aa09aa53b088d672a3c4e2c240e5ac4d`, tamanho **1.641.486 bytes**. O hash SHA-256 histórico declarado no parecer do E7 é `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7`. Este SHA-256 não foi **recalculado nesta fase** porque o conector não forneceu os bytes do ZIP.

A igualdade dos IDs Git prova que as referências apontam para o mesmo objeto armazenado, não a veracidade das relações do manuscrito ou o conteúdo do ZIP. Os metadados foram consultados diretamente no GitHub e registrados em `GITHUB_TREE_E7_SNAPSHOT.json`. O manifesto da ES8-P4 confere os arquivos de **esta** etapa, não os 32 componentes históricos E7.

## 3. Reconciliação analítica das convenções

As definições documentadas da V17 exigem:

\[
 F_k(r)=\int_{-L}^L f_k(x)e^{irx}\,dx,\qquad
 H_{ij}(r)=\frac{F_i(r)F_j(-r)+F_j(r)F_i(-r)}{2}.
\]

O termo arquimediano real é `1/π ∫_0^∞ Ω(r) H_ij(r) dr`. Como `f_i` são reais, `F_i(-r)=conj(F_i(r))` e `H_ij=Re(F_i conj F_j)`; sinais de paridade oposta se cancelam. O termo primo é `−2 log(2)/sqrt(2) C_ij(log 2)` na convenção V17, com o suporte `log 2<9/10<log 3` já auditado.

Para `N=6,8,10` e `k=3,...,N`, as colunas da projeção têm a forma

\[
 b_k=f_k-\frac{k}{p}\frac{p^2\pi^2+L^2}{k^2\pi^2+L^2}f_p,
 \quad p=1\ (k\ \text{ímpar}),\quad p=2\ (k\ \text{par}).
\]

O programa testa com números racionais **exatos** os coeficientes de `X=π²` e as partes constantes dos numeradores resultantes de `M_k(s)+c_{pk}M_p(s)`, verificando zeros exatos em ambos os polos `s=0,1` para todos os modos. A reconstrução por intervalos racionais verifica `LDL(M_N−(1/20)I)>0` para `N=6,8,10`, paridades opostas nulas e erro de suavização `δ=10⁻⁶` inferior a `1/20`. Esses são testes dos códigos reproduzidos da ES8-P3, não replay de dados E7/V16.

**Nota sobre a métrica:** a matriz `B` não possui colunas ortonormais; a margem `1/20` para a matriz expressa na base `B` é uma afirmação em termos da norma dos coeficientes `c∈R^(N−2)`. Para afirmá-la na norma Hilbert da função `Bc` seria necessário comparar `Bᵀ B` ou obter um limite inferior generalizado `M_N−λ Bᵀ B`.

## 4. Referências novas: CSVs 14 e 15

Foram lidos **31 registros**, com **28 trabalhos distintos** por DOI/título (3 ocorrências duplicadas). Dois são particularmente úteis:

1. **Borwein, Bradley e Crandall (2000)**, *Computational strategies for the Riemann zeta function*, DOI `10.1016/S0377-0427(00)00336-8`: inventário de métodos de avaliação numérica da zeta e suas regiões de validade. Pode apoiar auditorias independentes de `mpmath` e controles da função ξ; **não** resolve o cotejo dos coeficientes da matriz E7.
2. **Chaubey, Malik, Robles e Zaharescu (2018)**, *Zeros of normalized combinations of ξ⁽ᵏ⁾(s) on Re(s)=1/2*, DOI `10.1016/j.jmaa.2017.12.045`: estuda proporções de zeros na linha crítica para certas combinações de derivadas de ξ sob condições no artigo. É útil para a rota diferenciada/ES8-M, **não** demonstra que todos os zeros de ζ estão na reta crítica.

Como pesquisas indiretas foram destacados um artigo sobre **séries de Dirichlet e séries de potências**, um sobre **Mellin e Fourier nas sequências de Appell** e artigos sobre a **função de Möbius generalizada**. Os demais tratam predominantemente de polinômios Apostol–Euler/Bernoulli, Hurwitz–Lerch, q-extensões e números combinatórios; não há relação direta demonstrada com a igualdade das matrizes V17/E7.

A avaliação usa os **títulos e resumos entregues nos CSVs**, confirmando em fonte editorial o escopo dos dois trabalhos prioritários. Não são alegações de leitura integral ou verificação dos lemas de todos os 28 artigos.

## 5. Estado dos gates e rota seguinte

| Gate | Resultado ES8-P4 |
|---|---|
| Hash Git blob E7 consistente entre caminhos e commits | **PASSOU** |
| SHA de origem copiada ES8-P2 igual ao registro ES8-P3 | **PASSOU** |
| Fórmulas de polo e paridade V17 reproduzidas | **PASSOU** |
| LDL racional do modelo V17 para N=6,8,10 e δ=10⁻⁶ | **PASSOU condicionalmente** |
| ZIP E7 original aberto e SHA-256 recalculado nesta execução | **NÃO EXECUTADO** — bytes indisponíveis |
| Entradas E7/V16 comparadas uma a uma às de ES8-P3 | **NÃO ESTABELECIDO** |
| Positividade universal de Weil / RH | **NÃO PROVADA** |

**Próxima obrigação específica:** recuperar o ZIP original e os relatórios/fonte V16 em bytes comparáveis, verificar SHA-256 e manifesto interno, extrair os registros de entradas projetadas e mapear cada normalização e sinal em relação a `es8p3_arch_projected.py` / `es8p2_prime_cert.py`. Caso a representação histórica use outras bases ou janelas, documentar explicitamente a transformação de base antes da comparação de intervalos ou de sinais. A ausência dessas condições impede qualquer promoção do resultado finito para `E7_HISTORICAL_MATCH_CERTIFIED`.

## Fontes

- ES8-P3 original: https://github.com/jbarumby/scicalc/commit/2053129b03cffbce302b142af73df37e993f1613
- E7 histórico: https://github.com/jbarumby/scicalc/commit/fe6bdb7463172ea7107cd050410fe0f94196de29
- V17 normalizações: https://github.com/jbarumby/scicalc/blob/main/v17/RESULTADO_E8K_V17.json
- V24 envelope: https://github.com/jbarumby/scicalc/blob/main/v24/README_V24.md
- Borwein–Bradley–Crandall: https://doi.org/10.1016/S0377-0427(00)00336-8
- Chaubey–Malik–Robles–Zaharescu: https://doi.org/10.1016/j.jmaa.2017.12.045

**Resultado:** reconciliação de objetos Git e do modelo matemático V17; **sem revalidação do interior do ZIP E7**, sem prova da Hipótese de Riemann.
