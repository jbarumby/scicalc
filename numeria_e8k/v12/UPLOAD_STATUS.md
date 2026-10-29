# Status da incorporação E8-K V12

Data: 2026-10-28.

Esta branch é uma **área de staging**, não uma integração pronta para merge.

## Arquivos disponíveis nesta branch

- `VERIFICAR_E8K_V12.py`: verificador de replay V12.
- `CONFERIR_E8K_V12_V2.py`: conferidor independente das evidências exportadas.
- `RESULTADO_E8K_V12.json`: relatório textual de resultado.
- `README_E8K_V12.md` e `README_MELHORIAS.md`: documentação de escopo e limitações.
- `proveniencia_b3/ORIGEM_DOS_ARQUIVOS.json`: metadados de proveniência.
- `proveniencia_b3/RESULTADO_B3_G11_WINDOWS.json`: relatório B3 G11 em cópia de trabalho.
- `workflows/auditoria-e8k-v12.yml` e `workflows/auditoria-e8k-v12-evidencias.yml`: YAMLs preservados em staging; não são workflows ativos do GitHub Actions.
- `VALIDACAO_LOCAL.md`: resultados da execução local do pacote de origem.
- `SHA256SUMS_ORIGINAL_RELEASE_REFERENCE.txt` e `SHA256_MELHORIAS_ORIGINAL_RELEASE_REFERENCE.txt`: referências aos hashes declarados pelo pacote original, não atestados de integridade da árvore staged.

## Pendências bloqueantes

1. O ZIP binário `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip` (SHA-256 esperado `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7`) **não foi transferido para o repositório**.
2. As cópias textuais staged precisam ser comparadas byte a byte com o ZIP original; os hashes de referência não certificam automaticamente esses arquivos.
3. O segundo JSON de proveniência B3 G22 e o relatório detalhado de conferência V2 ainda não foram adicionados.
4. Os arquivos precisam ser postos nos caminhos de raiz esperados pelo workflow. Só depois de incorporar o ZIP E7, confirmar todos os hashes e ajustar os caminhos será seguro mover os YAMLs para `.github/workflows/`.
5. Não abrir nem mesclar PR de produção e não declarar GitHub Actions verde antes de concluir essas verificações.

## Validação local do pacote original

O ZIP completo fornecido foi testado localmente: ZIP íntegro; `py_compile` passou; o replay V12 estrito retornou código 0 para 15.754 pontos e zero lacunas de proveniência; o conferidor V2 retornou 8/8 gates PASS e 6/6 mutações isoladas detectadas. Esses resultados referem-se ao pacote local completo, não à árvore staged nem ao GitHub Actions.

## Escopo científico

A reprodução de pontos e a consistência dos certificados exportados não demonstram independentemente as cotas analíticas de quadratura ou cauda, a positividade global de Weil ou a Hipótese de Riemann.
