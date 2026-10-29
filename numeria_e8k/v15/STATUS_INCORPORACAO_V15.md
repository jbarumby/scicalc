# Status de incorporação — NUMERIA E8-K V15

Data: 28/10/2026

## Pacotes recebidos e validados localmente

- `E8K_V15_RESULTADOS.zip`: SHA-256 `bddb7cbc5c4d2216396731032a46090643057ee269a98754d8d7d04916bda374`; 4 entradas; teste de integridade ZIP passou.
- `E8K_V15_GITHUB_COMPLETO.zip`: SHA-256 `d22073ac6d0649517bf214dc5312b297781b97239a0a885c29e9f2e99b4197af`; 35 entradas; teste de integridade ZIP passou.
- `E7_DENTRO_DO_ZIP.zip` (arquivo externo recebido): SHA-256 `135c46704a009e1f5af90ee970fda880f518af5e8cc2a2206770dee13fb61c19`; integridade ZIP aprovada; contém o ZIP E7 interno.
- ZIP E7 interno, publicado como `numeria_e8k/v15/E7_original.bin`: SHA-256 local `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7`; integridade ZIP aprovada; manifesto interno: 32/32 arquivos conferidos.
- O arquivo remoto está presente na branch. O Git blob SHA-1 remoto `5b6abb59aa09aa53b088d672a3c4e2c240e5ac4d` coincide com o Git blob SHA-1 calculado localmente para os mesmos bytes. O conector não permite recalcular SHA-256 diretamente a partir do conteúdo binário remoto.

## Validação V15 previamente executada localmente

O verificador V15 terminou com código 0: 11/11 gates PASS, 7/7 controles negativos detectados. Resultado declarado: `CONDITIONAL_2X2_ARCHIMEDEAN_NEGATIVE_DEFINITE`. Isso não prova a fórmula global de Weil nem a Hipótese de Riemann; Lean não foi compilado.

## Publicado na branch

Branch: `integration/e8k-v14-package-2026-10-28`

- `numeria_e8k/v15/README_V15.md`
- `numeria_e8k/v15/PARECER_TECNICO_V15.md`
- `numeria_e8k/v15/SHA256SUMS_V15.txt`
- `numeria_e8k/v15/STATUS_INCORPORACAO_V15.md`
- `numeria_e8k/v15/E7_INCORPORACAO_E_AUDITORIA.md`
- `numeria_e8k/v15/E7_original.bin` (ZIP E7 interno)
- `numeria_e8k/v15/workflows_staged/auditoria-e8k-v15.yml` (apenas staged; não é workflow ativo)

## Auditoria E7 e limitações

O parecer E7 reporta verificações aritmético-geométricas, auditoria de 15.754 registros B5 e intervalos compostos sob lemas explicitados. Isso permanece condicionado às identidades matemáticas e às hipóteses documentadas. A normalização/interpretação da integral arquimediana na fórmula explícita global de Weil, a admissibilidade da função-teste e as normalizações precisam de prova separada. O teste de raios não verifica a posição do centro de cada integral. Portanto, E7 não constitui prova formal completa nem prova da Hipótese de Riemann.

## Itens ainda pendentes

Os artefatos binários/grandes e os dados operacionais do pacote V15 não foram todos enviados ao repositório. Faltam, em particular, o script executável `v15/VERIFICAR_E8K_V15_BLOCO.py`, `v15/RESULTADO_E8K_V15.json`, os dois JSONL B5, a proveniência B3, a amostra do oráculo e o código Lean. O manifesto V15 foi preservado como referência, mas seus hashes só podem ser verificados contra os arquivos correspondentes depois de carregá-los.

O workflow está intencionalmente em `workflows_staged/`; não foi ativado porque a estrutura e as dependências que ele exige ainda não estão completas na raiz da branch. A consulta à API de Actions não retornou execuções para esta branch. Lean não foi compilado.

## Próxima etapa

Fazer revisão integrada do diff antes de abrir PR. A branch está 35 commits à frente de `main` e não está atrás; por isso, a revisão do escopo é importante antes de propor merge. Não mesclar automaticamente.
