# Status da incorporação E8-K V12

Data: 2026-10-28.

Esta branch é uma **área de staging**, não uma integração pronta para merge.

## Arquivos enviados nesta etapa
- Verificador de replay V12.
- Conferidor de evidências V2.
- Relatório V12 (cópia textual de trabalho).
- Documentação V12 e das melhorias.
- Manifestos de hashes do pacote original, renomeados/retidos como referências.
- Registro de origem dos arquivos B3.
- Workflows YAML preservados em `workflows/` dentro da área de staging; não são workflows ativos do GitHub Actions.

## Pendências bloqueantes
1. O ZIP binário `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip` (SHA-256 esperado `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7`) **não foi transferido para o repositório** por esta integração.
2. Os relatórios B3 precisam ser conferidos byte a byte contra os SHA-256 originais antes de serem tratados como evidência executável.
3. Os arquivos devem ser comparados com os hashes dos pacotes de origem; os manifestos originais não devem ser usados para certificar cópias de trabalho modificadas.
4. Só depois de incorporar o ZIP binário e verificar os hashes será seguro mover os workflows para `.github/workflows/` na raiz e executar CI.
5. Não abrir/mesclar PR de produção nem declarar CI verde antes dessa conferência.

## Escopo científico
A reprodução de 15.754 pontos alegada pelos relatórios é uma checagem computacional, não valida independentemente as cotas analíticas de quadratura ou cauda, a positividade global de Weil ou a Hipótese de Riemann.
