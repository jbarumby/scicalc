# P84 — RESULTADO FINAL: S3 robusto com Arb, adoção versionada e downstream

**GitHub Actions:** [run 45038163289](https://github.com/jbarumby/scicalc/actions/runs/45038163289) — **completed/success**.
**Artefato:** [14759814197](https://github.com/jbarumby/scicalc/actions/runs/45038163289/artifacts/14759814197) — `numeria-p84-s3-robust-interval-grade-and-downstream`.

## 1. Documentos 1 e 2: novo certificado reproduzível no padrão P77

Os certificados enviados com `schema=REPO_STANDARD` são **fontes legadas não verificadas**, não evidência de que o projeto executou Arb. Foram preservados byte a byte:
- DOC1 original SHA-256 `093dd70b4fc6ddf0385094cd112699034f2fad2c36a3c51758c7f7c044652137`
- DOC2 original SHA-256 `ac356d2997e93945e8ce0751dd9e2f00054b45b9284a800d3ed690b57d8c6846`.

A CI P84 realmente reconstruiu as **16 linhas em cada uma de 192/256/384 bits** (`48` verificações de linha) utilizando `python-flint==0.9.0`, `ctx.workprec` e a quadratura `acb.integral` da P79/P77. O script calcula e armazena:
- intervalos dirigidos de integral `K1` e sua cauda, `K2` analítico, `K3` e `B1_upper`;
- extremos **diádicos exatos** na forma `[mantissa_string, exponent_int]`;
- margem exata como `fractions.Fraction`, com `B1_robusto - B1_upper >= 1/1000` nas sete linhas S2;
- positividade de `dK3/db4`, `K3>1`, expoente `e>0` e `dlog B1/db4>0`.

Foram gerados **novos** `DOC1_repo_interval.json` e `DOC2_repo_interval.json` a partir desses endpoints, com status exato `REPO_INTERVAL_CERTIFIED_CONDITIONALLY`. Esses arquivos não importam `float`, constantes ajustadas nem envelopes relativos não demonstrados.

SHA-256 dos relatórios novos:
- `DOC1_repo_interval.json`: `d7902eb62aabbddeb35f71fe1c84c71efce82964b1aa0e4a39cdbfa75a8ff7bc`
- `DOC2_repo_interval.json`: `09c3d218a4def9403e4305f58212ce3cd3b666761ed685e9986af5f4ec179b87`.

**Partes não recertificadas:** raízes numéricas `b4*` do DOC1, a afirmação geral `G'(ε)<0` e os 28 gates ζ atribuídos ao DOC2. A conclusão sobre intervalos é rigorosa **dentro das hipóteses externas identificadas**, não uma validação integral de cada frase dos PDFs.

## 2. Tabela S2 robusta adotada de modo versionado

A tabela verbatim `NUMERIA_P76_COEFFICIENTS.table.TABLE` **permanece intacta**. Adotamos uma nova versão explícita no mesmo módulo: `TABLE_S2_PUBLIC_ROBUST`, e o consumidor `NUMERIA_P75_CS2.cs_v2.TABLE_V2_S2_PUBLIC_ROBUST` selecionável. B1 por linha:
- 10: 17.49
- 11: 19.67
- 12: 22.25
- 13: 25.48
- 14: 29.89
- 15: 36.59
- 16: 51.65.

Mínima margem rigorosa sobre a cota Arb S2 em 256 bits (linha 16): aproximadamente **0.0019044440953141922**, superior a 0.001.

O replay downstream **foi efetivamente executado** com a tabela S2 robusta (não a P79 de seis casas), em **256 bits, 500 zeros, 512 bandas**: 
- `status=conditional_finite_positive`;
- **`18 <= t <= 79`**, **4.284 células aceitas**, nenhuma célula inconclusiva;
- menor margem inferior de Suzuki: diádico exato `92196961544760134523615486393004446412848310998322068400591052010486925957521 * 2^(-274)` ≈ **3.037370682403981×10⁻⁶**.

O replay permaneceu condicionado às estimativas e teoremas citados; nenhuma ampliação da faixa para além de 79 é reivindicada.

## 3. Varrimento downstream real

`scan_downstream.py` examinou **717 arquivos** no checkout real e identificou **36 arquivos com referências** a B1/TABLE sob as expressões pesquisadas. Classificação:
- **2** tabelas/consumidores versionados;
- **16** fontes ou certificados históricos, mantidos sem alteração;
- **5** módulos da própria P84;
- **9** saídas geradas durante a CI;
- **1** workflow CI;
- **1** protótipo arquivado;
- **2** referências documentais;
- **0** possíveis consumidores de execução ainda sem classificação.

A contagem 0 vale **para as regras heurísticas de busca implementadas**, e não demonstra ausência de consumidores dinâmicos ocultos. O JSON recebido, com apenas 17 arquivos de `/tmp/kit_p78`, foi preservado como evidência antiga, mas não usado para alegar cobertura do repositório.

SHA-256 do novo varrimento: `827a55e722adfe9881497fb8e0bda8a2c91aede3bcb7a3998e2403e370d0fcec`.

## 4. Hipóteses científicas e PDFs

Os originais anexados DOC1/DOC2 também foram encontrados como PDFs e seus bytes verificados localmente, independentes dos JSONs:
- `p78_s2_proposicao_auxiliar_b1.pdf` — SHA-256 `acd532d7f48c828577f3d1080f19f6741ab23af6594412b34ac5dbe7ed246841`;
- `p78_s2_passagem_completa (3).pdf` — SHA-256 `a52ae9cac7777130c5383b4c70267663bbf0594c05798d6b3b001e32c06d90ab`.

Os PDFs não são arquivos presentes no checkout da CI e seus hashes **não foram verificados pelo runner GitHub**: sua conferência foi feita localmente, e devem ser fornecidos a auditores que queiram reproduzir essa parte da cadeia.

A P78 tem `certificate.json` para todas as 16 linhas e o cálculo P84 reconstruiu os respectivos B1, incluindo as linhas 1–9. Isso **não** equivale à demonstração independente das estimativas Ramaré–Zúñiga em todas as alturas.

**Ainda abertas:** `3×10¹⁰ <= X < 3×10¹²` sem novo limite público suficientemente forte, e as provas externas de Ramaré–Zúñiga para `S_ε(X)` e Montgomery–Vaughan/Preissmann para o valor médio. É possível obter uma certificação computacional condicional sem fechar essas premissas analíticas; não confundir ambas.

## Conclusão

**P84 — computação e integração versionada aprovadas.** Full REPO_CERTIFIED matemático dos dois documentos, incluindo todas as reivindicações analíticas, **não** está encerrado. RH global **não está demonstrada**. PR continua **rascunho**, sem merge e sem alterar os certificados P75–P83.
