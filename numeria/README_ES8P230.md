# NUMERIA ES8-P230 — certificado SPD dos blocos finitos de 200 modos

## Resultado

Certifiquei positividade definida para as duas matrizes **finitas** reduzidas de 200×200 reconstruídas a partir dos blocos P227, condicionada à validade de suas caixas Arb, da cota GL60 P229 e dos termos intervalares adicionados. O limite inferior obtido para o menor autovalor é maior que \(5\times10^{-21}\) em ambos os setores. Isto é um avanço sobre a P228: além de comparar entradas, agora se prova uma margem espectral finita por uma âncora SPD verificada independentemente e uma cota de perturbação exata.

Os resultados reproduzíveis estão em `P230_FINITE_SPD.json`; a implementação está em `p230_finite_spd.py`.

| Setor | Pivôs LDL positivos | \(\|B-C\|_\infty\) majorada | Limite inferior uniforme para \(\lambda_{\min}(B)\) |
|---|---:|---:|---:|
| Par | 200/200 | \(5,93578494\times10^{-49}\) | \(>5\times10^{-21}\) |
| Ímpar | 200/200 | \(4,45183870\times10^{-49}\) | \(>5\times10^{-21}\) |

Aqui \(C\) é a matriz racional dos centros P193 usada somente como âncora numérica e \(B\) é a matriz finita alvo contida nas caixas P227 recompostas pela P228, incluindo \(\beta I\), polos e o erro quadratural P229. P230 recalcula SPD de \(C\) por LDL/inversa Arb; não usa os limites espectrais nem a cota `E_GL_total` da P193 para a desigualdade de Weyl. Os arquivos da matriz central são fixados pelos hashes SHA-256 guardados no JSON.

## Certificação da matriz central

Cada entrada de \(C\) é um racional diádico armazenado como mantissa e expoente; a simetria foi checada exatamente. Os centros foram convertidos, a 512 bits, em bolas Arb que contêm os racionais exatos. Uma eliminação simétrica por complementos de Schur (LDL intervalar) verificou que os **200 pivôs** em cada setor têm limite inferior estritamente positivo. Pelo critério de Sylvester, \(C\) é definida positiva.

Em seguida, `arb_mat.inv()` calculou enclosures da inversa. Somando os quadrados dos limites superiores dos módulos de todas as entradas, a aritmética Arb verificou

\[
\|C^{-1}\|_F^2 < 10^{40}.
\]

Como \(C\) é SPD, \(\lambda_{\min}(C)=1/\|C^{-1}\|_2\ge1/\|C^{-1}\|_F>10^{-20}\). Os enclosures numéricos do somatório quadrático são cerca de \(9,46838941\times10^{33}\) (par) e \(1,20273122\times10^{28}\) (ímpar), muito abaixo do teto \(10^{40}\). Esta margem \(10^{-20}\) é deliberadamente bem mais fraca que a margem espectral aproximada informada pelos certificados P193; P230 não depende desses valores de autovalor nem dos lower bounds P193.

## Perturbação racional e desigualdade de Weyl

A P228 foi ampliada para calcular, com frações exatas, a maior soma de linha dos desvios máximos entre o centro de cada caixa P227 e a entrada correspondente de \(C\). Para cada entrada \((i,j)\), toma-se

\[
e_{ij}=\max(|\ell_{ij}-C_{ij}|,|u_{ij}-C_{ij}|),
\]

 onde \([\ell_{ij},u_{ij}]\) é a caixa racional P227 completa. O alvo \(B\) é real simétrico pela própria fórmula: o integrando bilinear contém \(u_n(t)u_m(t)\), os termos de polo são simétricos em \(n,m\), e \(\beta_T\) só multiplica a diagonal. A P228 armazena o triângulo superior e P230 o reflete para a conta de linhas. A soma de cada linha inclui a entrada fora da diagonal também na posição simétrica. O máximo \(R=\max_i\sum_j e_{ij}\) é, portanto, um majorante exato de \(\|B-C\|_\infty\). Para matrizes reais simétricas, \(\|B-C\|_2\le\|B-C\|_\infty\); então Weyl implica

\[
\lambda_{\min}(B)\ge\lambda_{\min}(C)-\|B-C\|_2>10^{-20}-R.
\]

Nos dois setores, \(R<5\times10^{-21}\), logo \(\lambda_{\min}(B)>5\times10^{-21}>0\). `P230_FINITE_SPD.json` guarda as cotas exatas como frações, incluindo \(10^{-20}-R\), e também o limite uniforme mais simples.

## Reprodução

Com os 14 lotes Arb P227, os arquivos de entrada P193 do ZIP P220 e `python-flint` 0.9.0:

```bash
python numeria/p229_quadrature_bound.py \
  --p193 /tmp/p220/NUMERIA_ES8P220/inputs \
  --out numeria/P229_BOUND.json
python numeria/p230_finite_spd.py \
  --chunks /tmp/p227-todos \
  --p193 /tmp/p220/NUMERIA_ES8P220/inputs \
  --out numeria/P230_FINITE_SPD.json
python3 -m unittest -v \
  numeria/test_p228_recompose_audit.py \
  numeria/test_p230_finite_spd.py
```

P230 chama internamente a rotina P228, que revalida cobertura e metadados de todos os lotes, hashes das matrizes P193 e a derivação das caixas antes de computar a norma de linha.

## Escopo e limite

Este é um certificado de positividade definida **dos dois blocos finitos de 200 modos**. A prova numérica intervalar foi executada em Arb a 512 bits, e a recomposição/perturbação usa frações exatas. Ela não foi formalizada em Lean, Coq ou outro assistente de provas e depende de que o modelo do operador, as caixas P227, a derivação analítica P229 e os termos de polo/\(\beta\) representem corretamente a forma quadrática pretendida. A matriz central P193 é um ponto racional de referência cujo hash e simetria são verificados; suas cotas espectrais e de integração próprias não são pressupostos do certificado P230.

Mais importante: P230 não verifica os limites de cauda e acoplamento P207 nem os demais passos da redução global de Weil. Portanto, não é uma prova de positividade global de Weil e não prova a Hipótese de Riemann.
