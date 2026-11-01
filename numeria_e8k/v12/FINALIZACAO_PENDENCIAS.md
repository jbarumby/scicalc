# Finalização da integração E8-K V12

Data: 2026-10-28
Branch: `integration/e8k-v12-replay-2026-10-28`

## Estado confirmado

- A branch de integração existe; a branch `main` não foi alterada.
- Os scripts de verificação, relatórios de status e workflows em preparação estão na pasta `numeria_e8k/v12/`.
- A validação local do pacote V12 completo passou: compilação Python, reprodução de 15.754 pontos com hashes de entrada conferidos e 0 lacunas de proveniência; o autoteste do conferidor V2 passou em 8/8 gates e detectou 6/6 mutações de teste.
- Esses resultados são da cópia local do pacote, não de uma execução do GitHub Actions.
- O arquivo binário original `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip` ainda não está no repositório.
- Os arquivos de resultado transferidos manualmente não foram confirmados byte a byte contra os originais do pacote. Os manifestos SHA-256 estão preservados apenas como referência do release original.
- Os workflows permanecem em `numeria_e8k/v12/workflows/`; ainda não foram ativados em `.github/workflows/`.
- Nenhum PR foi aberto e nada foi mesclado.

## Gate para concluir

Não declarar a integração concluída, não ativar o workflow e não abrir/mesclar PR até que o ZIP binário original e os artefatos de evidência sejam adicionados sem alteração e os hashes SHA-256 sejam verificados contra o manifesto original.

## Passos finais no computador local

Execute em um clone local de `https://github.com/jbarumby/scicalc.git`:

```powershell
git fetch origin
git switch integration/e8k-v12-replay-2026-10-28
```

Copie o arquivo original `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip` e os artefatos originais do pacote para os caminhos correspondentes em `numeria_e8k/v12/`. Não recrie nem edite os JSONs de evidência. Em seguida:

```powershell
Get-FileHash .\numeria_e8k\v12\NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip -Algorithm SHA256
git status --short
git add numeria_e8k/v12
git commit -m "data(e8k): add original V12 package and evidence"
git push origin integration/e8k-v12-replay-2026-10-28
```

Compare o hash obtido com o manifesto original antes do commit. Depois de enviar os arquivos originais, execute a validação a partir da árvore clonada, ajuste os caminhos do workflow e só então considere ativá-lo em `.github/workflows/`.

## Integridade científica

A aprovação desses testes confirma apenas a reprodução computacional e os controles descritos. Não constitui prova da Hipótese de Riemann nem valida, por si só, uma nova demonstração matemática.
