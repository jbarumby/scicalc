# P79 — Relatório da reprodução Arb e da cobertura de Suzuki

**Execução verificada (GitHub Actions):** [45025753920](https://github.com/jbarumby/scicalc/actions/runs/45025753920) — job `audit` concluído com **success** em 22/11/2026.

**Artefato:** [numeria-p79-public-b1-and-suzuki](https://github.com/jbarumby/scicalc/actions/runs/45025753920/artifacts/14755746422), ID `14755746422`. Conteúdo: `p79-b1-192.json`, `p79-b1-256.json`, `p79-b1-256-512.json`, `p79-suzuki-cover.json`, `p79-suzuki-cover-512.json`, `p79-sha256sums.txt`.

## Coeficientes B1 substitutos

16 linhas comparadas com seus majorantes Arb em 192 e 256 bits; ambas as execuções deram **exatamente a mesma tabela decimal com arredondamento dirigido para cima**, SHA-256:

`d7ba517b4461e4791ad1ca92d69575d4f973b0d161b5a2a50b928ea2137af615`.

As linhas 1–9 retêm o B1 publicado e passam a desigualdade estrita. Os novos B1 das linhas 10–16 são:

| linha | b4 da P78 | B1 novo (limite superior decimal) |
|---|---:|---:|
| 10 | 557/1000 | 17.484560 |
| 11 | 557/1000 | 19.659253 |
| 12 | 279/500 | 22.239617 |
| 13 | 559/1000 | 25.472651 |
| 14 | 14/25 | 29.881218 |
| 15 | 561/1000 | 36.586170 |
| 16 | 303/500 | 51.648097 |

## Comparação das duas coberturas

| Faixas de densidade | Zeros | Bits | Resultado computacional | Teto da cobertura | Células aceitas |
|---:|---:|---:|---|---:|---:|
| 256 | 500 | 256 | `inconclusive` | `10308589/131072 ≈ 78.64829254` | 4248 |
| **512** | **500** | **256** | **`conditional_finite_positive`** | **79** | **4284** |

A primeira célula inconclusiva com 256 faixas foi
`[10308589/131072, 20617179/262144]`, largura `1/262144`.
Com 512 faixas **não restou célula inconclusiva** até 79.

A menor margem inferior estritamente positiva da cobertura com 512 faixas foi o diádico exato:

`46098626345208630096337756115931819003479026553077848153632446686159745792549 × 2^(-273)`
(aproximadamente `3.037380274×10^-6`), na célula `[2643/64, 661/16]`.

**Checksums úteis:** código `audit_public.py`:
`face1276ea7051964dc1ad96a72cf6e79fea35b7b15a09aea607b7c3c7e5b069`;
`p79-b1-256.json`:
`364b5cf0ce059ed733ac3eec8a7a46dcec8b770b05ccab8b75a671d2699f8e6a`;
`p79-suzuki-cover-512.json`:
`49596dde084b931ae2346d5a37a34c7d312c208557396ef5004ff88fb894d44b`.
O manifesto integral está no artefato.

## Leitura científica — não confundir os gates

**O cálculo intervalar passou**, e mostra que uma tabela B1 mais conservadora, baseada nos b4 públicos auditados na P78, **ainda permite** cobrir `18≤t≤79` ao refinar as faixas de densidade de 256 para 512.

**Não** constitui prova independente integral do Lema 6 de Chourasiya–Simonič: a P78 introduz hipóteses e desigualdades de Ramaré–Zúñiga Alterman **v1**, e toda a passagem até o Teorema 1 necessita auditoria matemática. O preprint Ramaré–Zúñiga possui versão v3 publicada em julho/2026, ainda não comparada formalmente linha a linha. Chourasiya–Simonič v2 permanece uma fonte independente com suas próprias hipóteses. A cobertura finita não prova RH global.

**Próxima prioridade:** cotejar v3 Ramaré–Zúñiga com v1 e justificar formalmente o transporte das desigualdades ao Lema 6; depois reproduzir todos os fatores da fórmula de Suzuki e as fontes externas.
