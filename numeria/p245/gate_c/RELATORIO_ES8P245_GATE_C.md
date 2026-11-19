# NUMERIA–RIEMANN — ES8-P245 Gate C

**Data:** 18/11/2026. **Janela:** meia-largura `L=161/200=0,805`, `T♯=504`; 1008 painéis de largura 1/2, cada qual com regra GL60 de 60 nós. **Escopo:** reconstituição independente de *somas quadraturadas selecionadas*, em aritmética racional intervalar, e auditoria separada do orçamento analítico GL60. A certificação integral do operador de Weil e a RH **não** seguem deste trabalho.

## 1. Linhagem das fontes

O pacote P245 Gate B foi reexecutado e seu certificado de 48 amostras coincidiu byte a byte com o anterior (`SHA-256` do certificado: `65f8c25d15eb0c4975fc2b432bafb0fb0e8bb8e269243a6d733cd94f5ea98b53`). Reutilizam-se, sem modificações, `p245b_certified_gl60.py` e `P245B_ROOT_BRACKETS.json` (as 60 raízes são novamente validadas por sinais exatos de `P_60`).

Foram recuperados os dois blocos P236 de 25 painéis, com SHA-256 correspondente ao manifesto de origem P237:

- par, `P236_CHUNK25_even.json.gz`: `6cfc2d3080389751b08d62dbdcdabdbd16a8be64ea5b23993dac0f27f61383f5`;
- ímpar, `P236_CHUNK25_odd.json.gz`: `76a7c147bf805fca4b67f65cbfd01464cf95134b30e990d3b5f24226a4de1f3d`.

Cada bloco contém 71.253 intervalos triangulares das somas GL60 de painéis 0–24. Os dados originais são **referência para inclusão**, não participam da recomputação independente das novas somas.

## 2. Replay intervalar de painéis 0–24

`p245c_chunk25.py` utiliza 60 nós/ pesos isolados em frações racionais, funções elementares em ponto fixo dirigido de 256 bits, digama por Binet deslocada, Bessel de graus baixos por fórmulas trigonométricas e Bessel de graus 4–33 por série com resto geométrico racional explícito.

Para os graus `n=752,753`, usa cotas bilaterais conservadoras por `|j_n(x)|<=x^n/(2n+1)!!`, em vez de atribuir um sinal não demonstrado. O erro da série de baixa frequência foi estendido de `x<=10` para `x<=11`, suficiente nos painéis 0–24; as razões de resto são majoradas com `x²<=121` e somente frações exatas.

Foram computados 13 pares de índice por paridade:

`(0,0), (0,1), (0,2), (0,3), (1,1), (1,2), (2,2), (3,3), (0,16), (16,16), (0,376), (1,376), (376,376)`.

São **25 × 13 × 2 = 650 avaliações intervalares de entradas de painéis**. Após somar as 25 caixas por paridade, o auditor comparou cada um dos **26 intervalos de soma** com a entrada correspondente do bloco original P236; o intervalo independente **contém a totalidade da caixa Arb** em todos os 26 casos. Não se demonstrou que as outras 71.240 entradas por setor foram reproduzidas.

**Reprodução determinística:** `P245C_CERT.json` e o replay inteiro são idênticos byte a byte; SHA-256 `8077ec76d8750538dea420285917d237394a83d6c66706a303892e6bb8a8b18d`.

## 3. Painéis estratificados no restante da janela

`p245c_stratified.py` calcula independentemente os painéis **100, 250, 503, 750, 875 e 1007**. Para graus pequenos usa recorrência esférica de Bessel iniciada em `j0`,`j1`, agora que `Lt>30`, e para graus 752/753 mantém a cota racional por magnitude. São **6 × 13 × 2 = 156 avaliações intervalares adicionais**.

Nos painéis 503 e 1007, os arquivos P236 originais disponibilizam **52 intervalos ao todo** para comparação, todos rigorosamente incluídos nas novas caixas; as fontes têm hashes SHA-256 conferidos. Para painéis 100, 250, 750 e 875 há apenas o resultado intervalar independente, **sem comparação externa Arb individual neste pacote**.

## 4. Orçamento analítico global GL60 — cálculo racional independente

`p245c_gl60_error_budget.py` recompõe o orçamento dos 754 graus, independentemente do certificado anterior. Assume duas entradas matemáticas externas: (i) extensão holomorfa do integrando à elipse de Bernstein `E_3`, com `|Ψ_L(t)-β|<50`; (ii) cota complexa esférica de Bessel `|j_n(Lt)| <= min(2,2*406^n/(2n+1)!!)` nessa elipse. Não formaliza a demonstração desses lemas; a origem P237 apresenta argumentos para ambos, ainda sujeitos à revisão matemática externa.

A desigualdade clássica `|c_k| <= 2M 3^{-k}` para coeficientes de Chebyshev limita a cauda desde `k=120` por `3M/3^120`. A quadratura GL60 é exata até grau 119; usando pesos positivos somando 2, a diferença integral–quadratura de cada painel é limitada por `12M/3^120`. Como a mudança de variáveis produz `1/(4π)` e `π>3`, `M<(50/12) a_i a_j`, logo:

`|I_ij - GL60_ij| < (1008*50/3^120)*a_i*a_j`, onde

`a_n = ceil(sqrt(2L(2n+1))) * min(2, 2*406^n/(2n+1)!!)`,

para cada grau `n=2i+paridade`, `i=0,...,376`. O código rederiva os tetos das raízes por aritmética inteira e verifica em `Fraction` as somas de 377 amplitudes e a pior soma de linha para cada setor. Resultados de apresentação:

| Setor | Majorante racional da pior soma de linha (aproximação) |
|---|---:|
| Par | `3.661198905758e-47` |
| Ímpar | `3.657549585347e-47` |

Ambos são estritamente menores que `4e-47` **sob as hipóteses analíticas acima**. Este é um replay do *orçamento de erro*, não uma avaliação intervalar de todas as integrais nem uma prova independente das hipóteses de holomorfia.

## 5. Testes e ameaças

`python -m unittest -q test_p245c.py` → **14/14 OK**. `python -O -m unittest -q test_p245c.py` → **14/14 OK**. São verificadas a proveniência das fontes, a inclusão por frações exatas das 26 somas e das 52 amostras estratificadas, o orçamento racional e a rejeição de adulteração de raiz e de arquivo-fonte. O replay estratificado e o orçamento também reproduziram seus respectivos certificados byte a byte.

O programa prova inclusões **se** sua implementação de operações intervalares e seus lemas clássicos estiverem corretos. Não é formalização integral em Lean. As 60 raízes são propostas históricas, mas sua propriedade de isolar todas as raízes de P60 é verificada novamente por frações exatas.

## 6. Escopo final e próxima prioridade

**Confirmado nesta etapa:** 650 entradas de painéis 0–24 computadas e agregadas em 26 somas incluídas nas caixas Arb P236; 156 amostras estratificadas, das quais 52 comparadas a Arb; 754 amplitudes de erro e duas cotas por linha, sob os majorantes analíticos declarados; reprodução byte a byte; 14 testes adversariais.

**Não confirmado:** as outras entradas de cada painel; todos os painéis e toda a matriz P237 por esta segunda implementação; revalidação externa da holomorfia e da fórmula explícita Weil–Zhu; positividade para suportes arbitrários; hipótese de Riemann. A positividade finita P242 e o orçamento infinito P239 continuam condicionados à fidelidade da montagem integral e das identidades importadas.

**Gate D sugerido:** certificar blocos P237 subsequentes em lotes (por exemplo 25–50, 500–525, 975–1000) com segunda implementação de funções especiais e intervalos dirigidos, e revisar externamente os lemas complexos que legitimam o erro GL60 antes de reivindicar certificação independente integral da matriz.

## Reprodução

```bash
cd NUMERIA_ES8P245_GATE_C
python3 p245c_chunk25.py --chunks . --roots P245B_ROOT_BRACKETS.json --out /tmp/P245C_CERT_REPLAY.json
python3 p245c_stratified.py --roots P245B_ROOT_BRACKETS.json --pilots pilot_sources --out /tmp/P245C_STRATIFIED_REPLAY.json
python3 p245c_gl60_error_budget.py --out /tmp/P245C_GL60_ERROR_BUDGET_REPLAY.json
python3 -m unittest -q test_p245c.py
python3 -O -m unittest -q test_p245c.py
```

Referências matemáticas: [DLMF, função esférica de Bessel](https://dlmf.nist.gov/10.54); [Trefethen, *Approximation Theory and Approximation Practice*](https://people.maths.ox.ac.uk/trefethen/trefethen_sample.pdf); [FLINT, erro GL por elipses de Bernstein](https://flintlib.org/doc/algorithms/ball-arithmetic/integration.html).
