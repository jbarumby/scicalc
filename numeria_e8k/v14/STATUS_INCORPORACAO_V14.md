# Status de incorporação — E8-K V14

Data: 2026-10-28  
Branch: `integration/e8k-v14-package-2026-10-28`

## O que foi feito

- O ZIP recebido `E8K_V14_GITHUB_COMPLETO.zip` foi aberto e sua estrutura foi inspecionada localmente.
- A verificação de integridade do ZIP terminou sem erro.
- Foi criada esta branch de integração, derivada da branch de trabalho E8K V12.
- Este README e este status foram publicados para registrar o escopo, os comandos e as limitações.

## O que ainda não foi incorporado como arquivos originais

O pacote contém 29 entradas, incluindo scripts Python, relatórios JSON, arquivos Lean, workflows, dois JSONL de pontos B5 (aprox. 6,2 MB combinados) e o ZIP E7 original (aprox. 1,64 MB). Esses arquivos binários/de evidência **não foram todos transferidos para esta branch por esta operação**. Portanto, a branch não deve ser tratada como cópia integral do ZIP V14.

O workflow V14 não foi ativado em `.github/workflows/` nesta etapa. Ativá-lo antes de enviar todos os caminhos referenciados causaria falha por arquivos ausentes.

## Como completar a incorporação sem alterar os bytes

1. Baixe o ZIP original anexado nesta conversa.
2. No clone local do repositório, mude para `integration/e8k-v14-package-2026-10-28`.
3. Extraia o ZIP na raiz de um diretório temporário e copie seu conteúdo para a raiz do clone, preservando os caminhos e sem recriar os JSONs.
4. Compare os hashes com `SHA256SUMS.txt` e investigue qualquer diferença antes de fazer commit.
5. Execute os cinco comandos listados em `README_E8K_V14.md`.
6. Só depois de todos os arquivos estarem presentes e os testes concluírem, mantenha ou mova os workflows para `.github/workflows/`, envie a branch e abra um PR para revisão.

Exemplo de comandos PowerShell, a partir da raiz do clone:

```powershell
git fetch origin
git switch integration/e8k-v14-package-2026-10-28
# Extraia E8K_V14_GITHUB_COMPLETO.zip em uma pasta temporária e copie os arquivos
# para a raiz do clone, preservando os caminhos relativos do ZIP.
Get-FileHash .\NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip -Algorithm SHA256
git status --short
# Após comparar os hashes e executar os testes:
git add .github NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip b5 proveniencia_b3 NUMERIA_E8K_FORMAL
git add CONFERIR_E8K_V14.py GERAR_E8K_V14.py GERAR_SHA256_V14.py MANIFESTO_SHA256.json ORACULO_RACIONAL_V14.py PARECER_E8K_V14.md README_E8K_V14.md REPLAY_E8K_V12_INTEGRAL.json RESULTADO_E8K_V13.json RESULTADO_E8K_V14.json SHA256SUMS.txt VERIFICAR_E8K_V12.py VERIFICAR_E8K_V13.py VERIFICAR_INVARIANTES_LEAN_V13.py referencias oraculo
git commit -m "feat(e8k): incorporate V14 replay and rational witnesses"
git push origin integration/e8k-v14-package-2026-10-28
```

## Não afirmar antes de comprovar

- Não afirmar que todos os arquivos V14 já estão no GitHub.
- Não afirmar que o GitHub Actions passou: nenhum run V14 foi executado nesta etapa.
- Não afirmar que Lean compilou: o parecer incluído no pacote diz que isso permanece por verificar.
- Não tratar estes testes como prova da Hipótese de Riemann.
