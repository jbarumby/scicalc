# ES8-P245 Gate B — auditoria independente de amostras GL60

**Estado da auditoria:** replay reproduzido localmente, certificado byte a byte idêntico e 12 testes aprovados nos modos normal e `python -O`. **Escopo:** 48 entradas selecionadas da soma GL60, não a matriz integral completa.

## Resultado reproduzido

O pacote original está preservado, sem alterações, em [`p245/package/`](p245/package/); seu `SHA256SUMS_P245B.txt` foi validado antes e depois da cópia. A execução local do certificador, lendo os 60 brackets racionais fornecidos e sem usar `mpmath` para a certificação, produziu um JSON idêntico a `P245B_CERT.json`. O SHA-256 do replay é:

```text
65f8c25d15eb0c4975fc2b432bafb0fb0e8bb8e269243a6d733cd94f5ea98b53
```

O programa verificou exatamente mudanças de sinal de `P_60` nos 60 intervalos racionais disjuntos de largura `2^-224`; a contagem igual ao grau do polinômio estabelece que todos os seus zeros foram isolados. Os pesos GL60 foram encerrados por aritmética intervalar dirigida, e a soma dos pesos contém 2.

Foram recalculadas **48 somas GL60 de painel**: oito pares de índices `(0,0), (0,1), (1,1), (0,16), (16,16), (0,376), (1,376), (376,376)`, para os painéis `0`, `503` e `1007`, em ambas as paridades. Em cada caso, o intervalo racional de ponto fixo produzido independentemente contém ambos os extremos exatos do intervalo Arb do piloto P236, cujos seis arquivos estão vinculados a SHA-256 no certificador. Os 12 testes passaram com Python normal e novamente com `python -O`, incluindo rejeição de brackets e piloto adulterados.

A cota intervalar do beta no certificado corresponde a aproximadamente `0,0012165183470747105186604451`. Para entradas que envolvem o modo de índice 376, o pacote também dá limites racionais estritos de magnitude para cada um dos três painéis. A tabela mostra os expoentes `k` para `|C_(0,376)|`, `|C_(1,376)|` e `|C_(376,376)|` menores que `10^-k` (os dois primeiros compartilham o mesmo expoente):

| Painel | Setor | `k` para as entradas cruzadas | `k` para a diagonal alta |
|---:|---|---:|---:|
| 0 | par | 2360 | 4722 |
| 0 | ímpar | 2364 | 4730 |
| 503 | par | 328 | 658 |
| 503 | ímpar | 329 | 660 |
| 1007 | par | 102 | 206 |
| 1007 | ímpar | 103 | 207 |

Essas são cotas para as **somas GL60 dos painéis indicados**, não para o erro de quadratura da integral nem para as entradas da matriz após todos os painéis.

## Checagem contra a convenção P237

Comparei a expressão do certificador P245 diretamente com o gerador P237 `p237_377_chunk.py`. A transformação é coerente: o nó do painel é `t=(2p+1+x_k)/4`; o argumento de Bessel é `L t`; a identidade `sqrt(pi/(2x)) J_(n+1/2)(x)=j_n(x)` converte a avaliação Arb em Bessel esférica; e a base usa `(-1)^i sqrt(2L(2n+1)) j_n(Lt)`, com `n=2i+paridade`. O peso `w_k/4 · (Psi_L(t)-beta)/pi` coincide com `w_k(Psi_L(t)-beta)/(4pi)` usado pelo P245. Os sinais manuais nos índices amostrados 0, 1, 16 e 376 correspondem a `(-1)^i`.

Assim, o P245 fornece uma checagem independente em aritmética de inteiros/Fractions para esses 48 valores do funcional **já discretizado por GL60**. A digama é avaliada pelo método Binet com deslocamento e resto explícito; as raízes, pesos, logaritmos, trigonometria e funções de Bessel também recebem inclusões intervalares. A comparação com os pilotos Arb é uma checagem cruzada, e não a fonte das caixas independentes.

## Limites que permanecem

O Gate B **não** certifica as outras entradas dos três painéis, os 1005 painéis restantes por setor, as 142.506 entradas triangulares completas, nem a diferença entre a integral verdadeira e a quadratura GL60. Os principais majorantes de Bessel para os graus altos são deliberadamente centrados em zero; as cotas racionais de magnitude mais fortes valem somente para as somas GL60 e os painéis enumerados.

O P245 também não verifica por si só a identidade geral entre a redução de Zhu e a forma quadrática de Weil, não prova positividade global, não substitui o certificado finito P237/P239 e não demonstra a Hipótese de Riemann. Seu avanço concreto é estreito, porém útil: três painéis distribuídos pelo intervalo e modos baixos, intermediários e extremos agora têm uma segunda avaliação rigorosa, baseada em implementação e aritmética diferentes das amostras Arb.

## Reprodução

A partir do diretório `numeria/p245/package`:

```bash
python3 p245b_certified_gl60.py \
  --pilot-dir pilot_sources \
  --roots P245B_ROOT_BRACKETS.json \
  --out /tmp/P245B_REPLAY.json
cmp /tmp/P245B_REPLAY.json P245B_CERT.json
python3 -m unittest -v test_p245b_certified_gl60
python3 -O -m unittest -q test_p245b_certified_gl60
```

O SHA do replay listado acima é idêntico ao do certificado anexado. O relatório e o código originais, incluindo suas ressalvas, estão preservados em `numeria/p245/package/`.
