# NUMERIA ES8-P228 — recomposição racional independente da P227

## Objetivo

Recompor, fora dos runners de CI, os 14 artefatos Arb 0.8 da P227 (8 lotes pares e 6 ímpares) e cotejar a matriz completa de 200 modos com as matrizes de referência P193 fornecidas no pacote P220. O script `p228_recompose_audit.py` converte cada centro/raio Arb em frações diádicas exatas, soma os lotes, acrescenta \(\beta I\), os termos de polo com enclosures racionais e a cota GL60 independente P229; faz então os testes de interseção entrada a entrada, também com frações exatas. Para \(\beta\), usa `mpmath.iv` e exporta os extremos intervalares como frações.

O script valida a cobertura sem lacunas nem duplicatas, os metadados Arb 0.8/1024 bits, os 20.100 elementos do triângulo superior por lote e os hashes gzip e descomprimido das duas matrizes P193.

## Resultado

Artefatos dos workflows GitHub Actions `43996977950` (13 lotes) e `43996977896` (primeiro lote) foram baixados. A auditoria processou 42.000 nós nos 700 painéis e recompôs os dois blocos de 200×200:

| Setor | Entradas triangulares cotejadas | Interseções com P193 | Desacordos | Maior diferença entre centros | \(\|B-C\|_\infty\) majorada |
|---|---:|---:|---:|---:|---:|
| Par | 20.100 | 20.100 | 0 | \(2,68172322\times10^{-75}\) | \(5,93578494\times10^{-49}\) |
| Ímpar | 20.100 | 20.100 | 0 | \(4,30963997\times10^{-76}\) | \(4,45183870\times10^{-49}\) |

O P227 usa agora os erros GL60 P229, de \(2,96789247\times10^{-51}\) (par) e \(2,22591935\times10^{-51}\) (ímpar); as referências P193 continuam ampliadas por seus próprios `E_GL_total`, cerca de 10.674 vezes maiores. As diferenças entre centros são da ordem de \(10^{-75}\), muito menores que ambas as cotas. Além da interseção das **40.200 entradas finitas**, P228 calcula com frações exatas a soma máxima de desvios de cada linha em relação à matriz central P193; a cota de perturbação resultante é usada por P230.

Com P230, os 200 pivôs LDL intervalares da matriz central são estritamente positivos em cada setor, e uma cota Arb da norma de Frobenius da inversa prova \(\lambda_{\min}(C)>10^{-20}\). A norma de linha P228 combinada com Weyl transfere isso para cada matriz finita P227: \(\lambda_{\min}(B)>5\times10^{-21}\) nos setores par e ímpar. Os detalhes e o JSON reprodutível estão em [README_ES8P230.md](README_ES8P230.md) e `P230_FINITE_SPD.json`.

## Como reproduzir

Requer Python 3 e `mpmath` instalado. Baixe os artefatos usando o GitHub CLI autenticado e descompacte o pacote P220 do usuário, que fornece as referências P193:

```bash
gh run download 43996977950 --repo jbarumby/scicalc --dir /tmp/p227-rest
gh run download 43996977896 --repo jbarumby/scicalc --dir /tmp/p227-primeiro
mkdir -p /tmp/p227-todos
cp -a /tmp/p227-rest/* /tmp/p227-todos/
cp -a /tmp/p227-primeiro/* /tmp/p227-todos/
unzip NUMERIA_ES8P220.zip -d /tmp/p220
python numeria/p229_quadrature_bound.py \
  --p193 /tmp/p220/NUMERIA_ES8P220/inputs \
  --out numeria/P229_BOUND.json
python numeria/p228_recompose_audit.py \
  --chunks /tmp/p227-todos \
  --p193 /tmp/p220/NUMERIA_ES8P220/inputs \
  --out numeria/P228_AUDIT.json
python numeria/p230_finite_spd.py \
  --chunks /tmp/p227-todos \
  --p193 /tmp/p220/NUMERIA_ES8P220/inputs \
  --out numeria/P230_FINITE_SPD.json
python3 -m unittest -v \
  numeria/test_p228_recompose_audit.py \
  numeria/test_p230_finite_spd.py
```

`P228_AUDIT.json` conserva os hashes dos 14 arquivos-fonte e os extremos racionais usados na verificação. Os artefatos de CI e os grandes arquivos P193 não foram copiados para o Git.

## Escopo e ressalvas

Este é um resultado de **matriz finita**. A positividade P230 é certificada para os blocos de 200 modos, mas depende da validade das caixas Arb P227, do limite analítico P229 e de que os termos recompostos representem a forma reduzida pretendida. A interseção P228, isoladamente, não prova positividade nem verifica os limites P207 de cauda/acoplamento. Mesmo com a margem finita P230, não se provou positividade global de Weil nem a Hipótese de Riemann.

O próximo gate prioritário é verificar a ponte de cauda/acoplamento P207 e a redução global de Weil; a margem espectral finita de 200 modos já recebeu o certificado separado P230.
