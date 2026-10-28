# Status de incorporação — NUMERIA E8-K V15

Data: 28/10/2026

## Pacotes recebidos e validados localmente

- `E8K_V15_RESULTADOS.zip`: SHA-256 `bddb7cbc5c4d2216396731032a46090643057ee269a98754d8d7d04916bda374`; 4 entradas; teste de integridade ZIP passou.
- `E8K_V15_GITHUB_COMPLETO.zip`: SHA-256 `d22073ac6d0649517bf214dc5312b297781b97239a0a885c29e9f2e99b4197af`; 35 entradas; teste de integridade ZIP passou.
- `E7_DENTRO_DO_ZIP.zip` (arquivo externo enviado): SHA-256 `135c46704a009e1f5af90ee970fda880f518af5e8cc2a2206770dee13fb61c19`; integridade ZIP aprovada; contém `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip`.
- ZIP E7 interno: SHA-256 `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7`; integridade ZIP aprovada; verificação do manifesto interno: 32/32 arquivos conferidos.

O verificador V15 foi executado no pacote extraído e terminou com código 0: 11/11 gates PASS, 7/7 controles negativos detectados. Resultado declarado: `CONDITIONAL_2X2_ARCHIMEDEAN_NEGATIVE_DEFINITE`. Isso não prova a fórmula global de Weil nem a Hipótese de Riemann; Lean não foi compilado.

## Publicado nesta branch

Branch: `integration/e8k-v14-package-2026-10-28`

- `numeria_e8k/v15/README_V15.md`
- `numeria_e8k/v15/PARECER_TECNICO_V15.md`
- `numeria_e8k/v15/SHA256SUMS_V15.txt`
- `numeria_e8k/v15/workflows_staged/auditoria-e8k-v15.yml` (apenas staged; não é workflow ativo)
- Este registro de status.

## E7 — auditoria de integridade e limites

O pacote E7 enviado foi verificado localmente: ambos os níveis ZIP passaram no teste de integridade e o verificador do manifesto confirmou 32/32 hashes. O parecer E7 reporta verificações aritmético-geométricas, auditoria de 15.754 registros B5 e intervalos compostos sob lemas explicitados. Isso permanece condicionado às identidades matemáticas e às hipóteses documentadas. O próprio parecer ressalva que a normalização/interpretação da integral arquimediana na fórmula explícita global de Weil, a admissibilidade da função-teste e as normalizações precisam de prova separada. O teste de raios não verifica a posição do centro de cada integral. Portanto, E7 não constitui prova formal completa nem prova da Hipótese de Riemann.

**Importante:** o ZIP E7 e seus arquivos binários foram recebidos e validados nesta sessão, mas o binário ainda não foi carregado para o GitHub. Os hashes acima documentam o artefato recebido; não significam que o arquivo já esteja versionado no repositório.

## Ainda não incorporado ao repositório

Os artefatos binários/grandes e os dados operacionais do pacote não foram todos enviados ao repositório. Faltam, em particular, o script executável `v15/VERIFICAR_E8K_V15_BLOCO.py`, `v15/RESULTADO_E8K_V15.json`, os dois JSONL B5, o ZIP E7, a proveniência B3, a amostra do oráculo e o código Lean. O manifesto V15 foi preservado como referência, mas seus hashes só podem ser verificados contra os arquivos correspondentes depois de carregá-los.

O workflow está intencionalmente em `workflows_staged/`; não foi ativado porque a estrutura e as dependências que ele exige ainda não estão completas na raiz da branch. Nenhuma execução do GitHub Actions foi confirmada e nenhum PR foi aberto.

## Próxima etapa

Carregar o pacote completo preservando os caminhos esperados, conferir os manifestos SHA-256 no clone e executar V15 no GitHub Actions. Compilar Lean em etapa separada e verificar ausência de axiomas/`sorry` não autorizados antes de atribuir status de prova formal.
