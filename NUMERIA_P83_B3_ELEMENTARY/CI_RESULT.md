# NUMERIA P83 — resultado da prova elementar de b₃=0,605

**Execução verificada:** [GitHub Actions #45036247458](https://github.com/jbarumby/scicalc/actions/runs/45036247458), job `independent-b3` **completed/success**, 22/11/2026.

**Artefato:** [numeria-p83-independent-b3-proof](https://github.com/jbarumby/scicalc/actions/runs/45036247458/artifacts/14759083078), ID `14759083078`, contendo `p83-b3-proof.json` e `p83-sha256sums.txt`.

SHA-256 do JSON: `3d91250c2bdd0319d331be4f07a3d0451700bb6ae8f28723458a68a70167fc8f`.

## Resultado matemático

Obtivemos **uma demonstração analítica alternativa**, independente dos resultados específicos de Ramaré sobre a soma curta e a contagem fina de quadrados livres, de que **para todo X real ≥ 10⁹**:

```
sum_(X<n<5X) lambda_X(n)^2/n^2 <= 0.605/X
```

A prova geral está no [README.md](README.md), e não se baseia na enumeração de valores de X. A fórmula de contagem de livres de quadrados, o controle do resto por `2.01*sqrt(x)`, a soma por partes e o limite `|lambda_X(n)|≤1` para `2X<n<5X` conduzem ao majorante racional:

```
X * sum_(X<n<5X) lambda_X(n)^2/n^2
  < 178790073271447/295788000000000
  ≈ 0.604453437162586
  < 121/200 = 0.605.
```

A **folga estrita exata** é

```
161666728553/295788000000000
≈ 0.000546562837414.
```

A cota `π>157/50` também recebe uma justificativa matemática independente via identidade de Machin e o primeiro par de termos da série alternada de arctan.

## Testes de reprodução

- **8/8 testes** aprovados em Python normal e **8/8** sob `python -O`.
- Frações exatas confirmaram o orçamento global e a folga estrita.
- Testes finitos conferiram valores de Möbius, condições adversariais e o limite `|lambda|≤1` em `2X<n<5X`. Eles **não são** usados como prova para X grande; esta vem das identidades gerais no README.
- Os artefatos e fontes receberam manifesto SHA-256.

## Implicações e limitações

A P82 tinha três premissas analíticas importadas:
1. Cota curta `b3=0.605` (Ramaré Lema 5.6 / KLN Lema 3.4 (3.8)): **agora suprida por demonstração elementar independente** na P83 para o mesmo domínio relevante.
2. Desigualdade de valor médio `T+πm0(n+1)` (Ramaré Lema 6.5 / KLN Lema 3.6): **ainda importada**, embora a antissimetria para coeficientes reais e a separação `log((n+1)/n)≥1/(n+1)` estejam explicadas.
3. Soma dupla de Möbius `S_ε(X)` (Ramaré–Zúñiga Alterman v1/v3): **ainda importada** nas faixas restritas `X≥3×10¹²`.

Os novos dados **não** alteram os coeficientes B1 da P79, não ampliam a janela finita de Suzuki até t=79 e não representam uma prova da Hipótese de Riemann global.

**Próxima prioridade P84:** demonstrar/reconstituir criticamente o teorema ponderado de Hilbert que fornece a constante `m0=sqrt(1+(2/3)*sqrt(6/5))`, usando fonte primária Montgomery–Vaughan e tomando cuidado com hipóteses para séries infinitas.
