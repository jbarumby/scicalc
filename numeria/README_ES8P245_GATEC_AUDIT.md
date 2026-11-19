# ES8-P245 Gate C — auditoria independente e alcance

**Resultado da revisão:** hashes do pacote válidos; 14 testes passaram em Python normal e `python -O`; os três replays completos foram byte a byte idênticos aos certificados anexados. O Gate C é uma auditoria parcial útil de somas GL60, não uma nova prova global.

## Reprodução independente

O pacote completo recebido está preservado sem modificação em [`numeria/p245/gate_c/`](p245/gate_c/). O manifesto `SHA256SUMS_P245C.txt` passou. O código Gate B e `P245B_ROOT_BRACKETS.json` coincidem byte a byte com a versão do Gate B já arquivada.

Reexecutei:

1. `p245c_chunk25.py`: 650 avaliações intervalares selecionadas em 25 painéis (13 pares por setor) e inclusão das 26 somas acumuladas — 13 por paridade — nos intervalos Arb do bloco P236 correspondente. JSON reproduzido byte a byte; SHA-256 `8077ec76d8750538dea420285917d237394a83d6c66706a303892e6bb8a8b18d`.
2. `p245c_stratified.py`: 156 avaliações nos painéis 100, 250, 503, 750, 875 e 1007. Nas 52 entradas dos painéis 503 e 1007 com piloto Arb disponível, as caixas Arb foram incluídas nas novas caixas racionais. Nos outros 104 casos não há comparação Arb externa no pacote. JSON reproduzido byte a byte; SHA-256 `d2a2696efc4c374e687fc21390f773d0a431dd777ca491ce7058728531803638`.
3. `p245c_gl60_error_budget.py`: orçamento de erro composto para 1008 painéis e 754 graus, rederivado por frações exatas. JSON reproduzido byte a byte; SHA-256 `84632341c00c635256d9d43fe01e052975486c41a07188445dae6fad1bc5ae55`.

Os 14 testes foram executados em ambos os modos, com resultado `OK`. Eles incluem verificações de hashes, brackets das raízes, inclusões exatas, orçamento racional e rejeição de alguns arquivos adulterados. A suíte também checa os indicadores de escopo dos certificados.

## Leitura matemática

O replay dos primeiros 25 painéis avalia os pares `(0,0), (0,1), (0,2), (0,3), (1,1), (1,2), (2,2), (3,3), (0,16), (16,16), (0,376), (1,376), (376,376)`. As entradas extremas de grau 752/753 são limitadas por majorantes de Bessel centrados em zero; não se lhes atribui um sinal que o cálculo não estabeleça. As outras amostras usam intervalos dirigidos de 256 bits, raízes e pesos GL60 racionais, Binet para digama e séries/recorrências de Bessel com restos controlados.

O orçamento racional dá, como apresentação, limites superiores para a pior soma de linha de aproximadamente `3,661198905758 × 10^-47` no setor par e `3,657549585347 × 10^-47` no ímpar, ambos menores que `4 × 10^-47`. A aritmética das amplitudes e somas de linha é refeita exatamente. A implicação para o erro integral, contudo, depende de hipóteses analíticas explícitas: holomorfia do integrando na elipse de Bernstein `E_3`, `|Ψ_L(t)-β|<50` nela, e `|j_n(Lt)| ≤ min(2, 2·406^n/(2n+1)!!)` na região complexa correspondente.

A geometria da elipse é compatível com essas cotas: para os painéis cobertos, `|Im t|≤1/3`, o argumento digama `z=1/4+it/2` satisfaz `Re z≥1/12` e `|z|<253`. A estimativa publicada no material P237 usa a expansão em frações parciais da digama, limites racionais para os termos de primos e `cosh(2/3)≤9/7`, obtendo `|Ψ_L-β|<274/7<50`, sob identidades clássicas nomeadas. Minha checagem confirma a consistência dessas reduções algébricas e constantes, mas não formaliza as identidades complexas importadas. Além disso, o script de majorante P237 que as verifica não está atualmente incluído no branch versionado; a derivação precisa ser incorporada e revisada de forma durável antes de tratar essa parte como um certificado analítico independente.

## Limite da conclusão

O Gate C não recalcula independentemente todas as 71.253 entradas triangulares de cada setor, todos os 1008 painéis, nem todas as integrais. As comparações dos painéis 0–24 são feitas após somar 25 painéis para os 13 pares selecionados. A amostragem dos outros painéis é esparsa; não se infere dela a validade de todos os painéis ou de todas as entradas.

Em particular, os certificados não verificam a identidade geral entre a redução de Zhu e a forma quadrática de Weil, não demonstram positividade para suportes arbitrários e não provam a Hipótese de Riemann. O orçamento GL60 é uma ponte condicional para controlar o erro da quadratura; ele não substitui a matriz P237 completa nem remove as hipóteses analíticas importadas. Portanto, o resultado útil é uma confirmação computacional independente, mas estratificada e parcial, mais uma reexecução racional do orçamento.

## Reprodução

A partir de `numeria/p245/gate_c`:

```bash
python3 p245c_chunk25.py --chunks . --roots P245B_ROOT_BRACKETS.json --out /tmp/P245C_CERT_REPLAY.json
cmp /tmp/P245C_CERT_REPLAY.json P245C_CERT.json
python3 p245c_stratified.py --roots P245B_ROOT_BRACKETS.json --pilots pilot_sources --out /tmp/P245C_STRATIFIED_REPLAY.json
cmp /tmp/P245C_STRATIFIED_REPLAY.json P245C_STRATIFIED_CERT.json
python3 p245c_gl60_error_budget.py --out /tmp/P245C_GL60_ERROR_BUDGET_REPLAY.json
cmp /tmp/P245C_GL60_ERROR_BUDGET_REPLAY.json P245C_GL60_ERROR_BUDGET.json
python3 -m unittest -q test_p245c.py
python3 -O -m unittest -q test_p245c.py
```

**Prioridade concreta:** versionar e revisar a prova do majorante complexo `E_3`, formalizar ou justificar de modo independente as identidades especiais importadas e, então, reproduzir blocos adicionais de painéis com comparação integral do triângulo superior. Até lá, P245C fortalece a confiança em amostras e no orçamento condicional, não a conclusão global.
