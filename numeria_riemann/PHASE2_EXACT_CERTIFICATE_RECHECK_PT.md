# NUMERIA–RIEMANN — fase 2: rechecagem independente dos certificados

**Data:** 2026-10-22  
**Branch:** `numeria-riemann/weil-hermitian-audit`  
**Escopo:** verificar os JSONs exportados sem Arb, usando aritmética racional exata para o passo final de positividade; auditar arredondamento externo; comparar relatórios de 192 e 256 bits e o método legado com o adaptativo.

## 1. Integridade dos pacotes recebidos

Os três pacotes ZIP fornecidos pelo usuário contêm os mesmos 11 caminhos e conteúdo idêntico em todos os arquivos internos (SHA-256 por arquivo igual). Os hashes dos arquivos ZIP inteiros diferem, portanto os ZIPs são embalagens diferentes dos mesmos artefatos, não três execuções independentes.

## 2. Verificação racional exata da positividade

Foi executado `verify_certificate.py` sobre os seis relatórios (N=6, 8, 10 a 192 e 256 bits). O verificador lê os limites decimais exportados como racionais exatos, forma o ponto médio (M), a matriz de raios (ho), calcula (delta=\max_i\sum_j\rho_{ij}) e aplica LDL exato a (M-\delta I).

| N | Precisão | Dimensão restrita | δ (norma de raio) | Menor pivô de (M-\delta I) | Resultado |
|---:|---:|---:|---:|---:|---|
| 6 | 192 bits | 4 | 0.012675985562964343 | 0.2189160458563405 | PASS |
| 6 | 256 bits | 4 | 0.012675985577516258 | 0.2189160458417886 | PASS |
| 8 | 192 bits | 6 | 0.02792852713173488 | 0.20366350428756996 | PASS |
| 8 | 256 bits | 6 | 0.02792852710263105 | 0.2036635043166738 | PASS |
| 10 | 192 bits | 8 | 0.05178402640740387 | 0.17980800501190097 | PASS |
| 10 | 256 bits | 8 | 0.051784026392851956 | 0.1798080050264529 | PASS |

Em todos os seis relatórios: nenhum pivô falhou, o verificador exato concorda com o status declarado e o código de saída foi 0.

## 3. Auditoria do arredondamento externo

Executado `audit_outward_rounding.py` nos seis relatórios:

| N | Precisão | Enclosures verificados | Violações |
|---:|---:|---:|---:|
| 6 | 192 bits | 20 | 0 |
| 6 | 256 bits | 20 | 0 |
| 8 | 192 bits | 42 | 0 |
| 8 | 256 bits | 42 | 0 |
| 10 | 192 bits | 72 | 0 |
| 10 | 256 bits | 72 | 0 |
| **Total** | | **268** | **0** |

Os limites decimais exportados contêm os extremos reconstruídos de `mid_exact` e `rad_exact`, com a ordem correta, em todos os casos auditados.

## 4. Consistência entre precisões e método legado

O verificador de intervalos confirmou interseção de todos os pares de envoltórios da matriz projetada entre 192 e 256 bits para N=6, 8 e 10.

Também foi comparado o relatório legado de N=6 com o relatório adaptativo de N=6. Os envoltórios correspondentes da matriz projetada se intersectam; o teste racional de positividade do relatório legado também passou. Isso é compatibilidade entre os certificados exportados, não independência da derivação matemática ou da quadratura subjacente.

## 5. O que foi e não foi verificado

**Verificado nesta fase**
- A positividade definida de cada matriz finita exportada é implicada pelos respectivos envoltórios, segundo o teste de Rump refeito em aritmética racional exata.
- Os 268 envoltórios auditados foram arredondados para fora sem violações.
- Os relatórios 192/256 bits são intervalarmente compatíveis para N=6, 8 e 10.
- O relatório legado N=6 é compatível com o adaptativo N=6 nos envoltórios comparados.

**Não verificado por estes scripts**
- A quadratura Arb original não foi recalculada independentemente a partir do código-fonte.
- A prova analítica das cotas de cauda e a validade da fórmula explícita continuam dependentes da revisão matemática.
- Não se estabelece positividade global da forma de Weil, nem a Hipótese de Riemann.

## 6. Próxima etapa recomendada

Prosseguir com uma revisão independente da cadeia analítica que produz os envoltórios: (i) fórmula de Guinand–Weil e convenções de Fourier; (ii) admissibilidade da função teste; (iii) cota da cauda arquimediana; (iv) propagação das incertezas pela projeção (B^{\mathsf T}WB); e (v) relação entre os dados exportados e a matriz efetivamente certificada.

Não aumentar a dimensão da matriz nesta etapa. A evidência computacional finita está reforçada; o próximo gargalo é a validação analítica e a reprodução independente do cálculo que gera os envoltórios.
