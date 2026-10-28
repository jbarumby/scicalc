# NUMERIA E7 — registro de incorporação e auditoria

Data de registro: 28/10/2026  
Branch: `integration/e8k-v14-package-2026-10-28`

## Artefato E7 publicado

O arquivo `numeria_e8k/v15/E7_original.bin` foi adicionado à branch pelo commit `a58281f8fedcc014a4f5a93ef7ca2a08c7e6a33e`.

- Tamanho local: 1.641.486 bytes.
- SHA-256 calculado no arquivo local: `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7`.
- Git blob SHA-1 informado pelo GitHub para o arquivo remoto: `5b6abb59aa09aa53b088d672a3c4e2c240e5ac4d`.
- O Git blob SHA-1 calculado localmente a partir dos mesmos bytes coincide com o blob SHA-1 remoto, confirmando a correspondência byte a byte do artefato enviado.
- O arquivo é um ZIP interno E7; o teste de integridade ZIP local passou e o manifesto interno confirmou 32/32 hashes.

## Artefatos recebidos e origem

| Artefato | SHA-256 | Estado |
|---|---|---|
| `E7_DENTRO_DO_ZIP.zip` (contêiner recebido) | `135c46704a009e1f5af90ee970fda880f518af5e8cc2a2206770dee13fb61c19` | ZIP íntegro |
| `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip` (ZIP interno; publicado como `E7_original.bin`) | `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7` | ZIP íntegro; manifesto: 32/32 hashes conferidos |

## Conteúdo e resultado reportado

O pacote inclui parecer E7, manifesto SHA-256, revisores/verificadores, resultado E7, resultados corrigidos B3 G11/G22, arquivos de referência da Proposição 6 para N=6, 8 e 10, registros B5 e código auxiliar para auditorias e pesos racionais.

O parecer declara aprovação condicionada às identidades matemáticas explicitadas e intervalos compostos sob lemas documentados. A auditoria inclui registros B5 e testes em raios de mesma paridade. O pacote ressalva que os testes de raios não verificam a posição do centro de cada integral.

## Limites científicos

- A integridade criptográfica comprova correspondência dos arquivos, não a verdade das alegações matemáticas.
- As verificações computacionais conferem consequências aritméticas das identidades adotadas; não formalizam por si sós todos os teoremas analíticos citados.
- A admissibilidade da função-teste e a interpretação/normalização da integral arquimediana na fórmula explícita global de Weil exigem demonstração independente completa.
- Este pacote não é uma prova formal completa da fórmula global de Weil e não prova a Hipótese de Riemann.

## Estado da incorporação

O binário E7 está versionado na branch de integração e o GitHub confirma um blob cujo identificador SHA-1 coincide com o calculado localmente para os mesmos bytes. O SHA-256 foi calculado localmente; o conector não disponibilizou o conteúdo binário remoto para recalcular diretamente o SHA-256 a partir do download remoto. O manifesto interno 32/32 foi validado localmente antes do envio.

A branch ainda precisa de revisão integrada antes de qualquer merge. Nenhum PR V15 foi aberto nesta etapa e a `main` não foi modificada por esta incorporação.
