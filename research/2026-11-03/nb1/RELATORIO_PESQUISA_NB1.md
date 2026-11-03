# NUMERIA–RIEMANN — Pesquisa e auditoria experimental NB1

**Data:** 03/11/2026. **Status:** certificado finito de *erro para coeficientes racionais fixos*; sem certificado de otimalidade, de limite para dimensões arbitrárias ou de RH.

## 1. Achados documentais verificados

1. **Werner Ehm (2024), *On certain Gram matrices and their associated series*.** https://arxiv.org/abs/2405.06349. Fórmulas de matrizes de Gram no critério de Báez–Duarte em termos de séries de tipo Müntz; formulações q=1 e q=2. Prioridade: reproduzir entradas G^(1), F^(1) e cotejá-las com a integração analítica intervalar deste módulo. Os coeficientes de Bettin–Conrey–Farmer aparecem como benchmark na introdução de Ehm, com hipóteses adicionais para certas afirmações assintóticas.
2. **Sébastien Darses e Erwan Hillion (2021), *On probabilistic generalizations of the Nyman–Beurling criterion for the zeta function*.** https://cml.centre-mersenne.org/item/10.5802/cml.71.pdf. No espaço H = L²(0,∞), com ρ_θ(x) = {θ/x}, a função indicadora de (0,1) pertence ao fecho do span de ρ_(1/n) se, e somente se, RH. Ver Teoremas 1.1 e 1.3.
3. **Xuefeng Zhu (2026), *Weil positivity in compact windows: a finite reduction, certified two-sided bounds, and a Landau-Widom decay law*.** https://arxiv.org/abs/2608.24827. Preprint v2 de 02/09/2026; declara positividade incondicional em uma janela de suporte até L=0,8, com cota extremamente pequena. **Não auditado de forma independente neste trabalho.** Prioridade de cotejo com sinais, projeções, regularidade e constantes do E8K/V24, sem importar os resultados como teoremas verificados pelo NUMERIA.
4. **Kim et al. (2026), *A Numerical Realization of Suzuki's Weil-Quadratic-Form Operator*.** https://arxiv.org/abs/2607.24830. Modelagem por elementos finitos P1 e extrapolação de Richardson; o artigo ressalta que não fornece demonstração de RH e diferencia o operador de forma de Weil de um operador com zeros como autovalores. Serve como benchmark espectral ES2, não como equivalência automática entre operadores.
5. **FLINT, *Zeros of the Riemann zeta function*.** https://flintlib.org/doc/acb_dirichlet/index.html. Dispõe de isolamento rigoroso de zeros pela Hardy Z e contagem por método de Turing. Deve funcionar como conferência independente da NL2, não como substituto da positividade global.

## 2. Experimento executado: aproximações finitas na base de Báez–Duarte

Sejam χ = 1_(0,1), ρ_n(x) = {1/(n x)} para n >= 1, c_n racionais fixos, e

    R_N(x) = χ(x) - sum_{n=1}^N c_n ρ_n(x),
    E_N(c) = int_0^infinity |R_N(x)|² dx.

A escolha θ_n=1/n é a base da versão forte do critério de Báez–Duarte; o limite E_N^* = inf_c E_N(c) -> 0 é equivalente a RH. **O presente cálculo apenas limita E_N(c) para três c finitos.**

### 2.1 Separação exata do domínio

Fixamos ε = 1/10000.

* Se x > 1, tem-se ρ_n(x) = 1/(n x). Por isso, o erro integrado é exatamente (sum_n c_n/n)².
* Se 0 < x < ε, como 0 <= ρ_n(x) < 1, temos 0 <= integral |R_N(x)|² dx <= ε(1 + sum_n |c_n|)².
* Para ε <= x <= 1, particionamos em 9999 intervalos abertos (1/(k+1), 1/k), k=1,...,9999. Dentro de cada intervalo, floor(1/(n x)) = floor(k/n), para n inteiro.

Pondo

    A_k = 1 + sum_n c_n floor(k/n)
    B = -sum_n c_n/n

obtemos em cada célula R_N(x) = A_k + B/x e a identidade exata

    integral_{1/(k+1)}^{1/k} |R_N(x)|² dx
      = A_k²(1/k - 1/(k+1)) + 2 A_k B log((k+1)/k) + B².

Assim, o único termo não racional é o logaritmo. Ele é englobado por séries racionais com resto explícito:

    log(r) = 2 * sum_{j=0}^{m} [z^(2j+1)/(2j+1)] + R_m,
    z = (r-1)/(r+1),
    0 <= R_m <= 2 z^(2m+3) / ((2m+3)(1-z²)),  1 < r <= 2.

O algoritmo soma racionalmente com arredondamento **para fora** a 26 casas decimais; o único uso de ponto flutuante é na determinação heurística de coeficientes candidatos, que são então arredondados para números racionais **antes** da certificação.

### 2.2 Resultados

| N | Limite inferior E_N(c) | Limite superior E_N(c) | Estimativa (float, não certificada) da condição da matriz |
|---:|---:|---:|---:|
| 6 | 0,0359686542 | 0,0382842437 | ~117,4 |
| 8 | 0,0241455736 | 0,0273231717 | ~200,7 |
| 10 | 0,0237561967 | 0,0272368780 | ~329,4 |

**Interpretação:** os intervalos acima englobam o erro **dos coeficientes escolhidos**, e não o mínimo E_N^*. Como exemplo, a linha N=10 prova apenas que existe uma escolha racional finita com E_10(c) <= 0,0272368780. A observação de erros pequenos ou decrescentes em um conjunto finito não prova que E_N^* -> 0.

Os valores racionais integrais e os coeficientes estão no arquivo `RESULTADOS_NB1.json`. As quantidades decimais desta tabela são arredondadas para facilitar a leitura; a evidência exata encontra-se em JSON.

### 2.3 Testes e limitações

Comando:

    python nb1_gram_audit.py
    python test_nb1_gram_audit.py

Cinco testes passaram: representação por células, avaliação independente pontual do logaritmo com mpmath, invariância por preenchimento de coeficientes com zeros, caso trivial χ com coeficientes nulos e consistência dos limites gerados.

Dependências: Python 3, NumPy, mpmath (somente nos testes). **python-flint/Arb não esteve disponível no ambiente desta execução**; por isso, o módulo utiliza frações exatas de Python, não declara integração com Arb já realizada e não executou CI no GitHub.

As provas elementares dos envelopes constam deste relatório, e o programa as implementa; para publicação formal, ainda é necessária revisão de fonte independente, testes adversariais adicionais e, idealmente, formalização mecânica.

## 3. Relação com NL2 e E8K/V24

* **NB2 (próxima etapa):** reconstruir rigorosamente as matrizes G_N e os vetores b_N, e certificar os **minimizadores** das equações lineares G_N c=b_N com intervalos (ou Krawczyk para a solução do sistema); então comparar os limites para N=6,8,10,20,30. Não substituir a exigência de controle uniforme quando N -> infinito.
* **ES2:** reproduzir a estrutura de matrizes/operadores de Suzuki e os testes de elementos finitos de Kim et al., controlando distinção entre operadores e hipóteses condicionais.
* **Weil V24:** comparar o preprint de Zhu com as normalizações e a projeção de V24 antes de adotar seu limite de suporte. A positividade de Gram em NB1 decorre da definição de produto interno; **não implica** por si só a positividade da forma de Weil.

## 4. Proveniência

Este pacote foi produzido localmente, sem gravar arquivos no repositório GitHub. O programa e o JSON são disponibilizados para repetição e inspeção independente. Nenhuma demonstração da Hipótese de Riemann foi obtida.
