# Status de incorporação — E8-K V14

Data: 2026-10-28
Branch: `integration/e8k-v14-package-2026-10-28`

## Pacote completo recebido e testado localmente

O arquivo `NUMERIA_E8K_V14_TODOS_ARQUIVOS.zip` foi aberto e inspecionado. O ZIP passou na verificação de integridade. SHA-256 do ZIP:

`794ac407aa2278639fd2182a65f1ea36f02fa37044676ba90b8f71cdc5c29d12`

Contém 29 entradas, aproximadamente 10 MB descompactados, incluindo os workflows V13/V14, scripts Python, relatórios JSON, E7 ZIP original, código Lean, dois JSONL B5, proveniência B3, amostra de 128 oráculos e PDF de referência.

## Validação local concluída

- Compilação sintática dos sete scripts Python: passou.
- Replay V12: código de saída 0.
- Verificação V13: código de saída 0.
- Geração V14: código de saída 0.
- Conferidor V14 com autotestes: **9/9 gates PASS; 0 FAIL**.
- Verificação SHA-256: todos os 27 arquivos do manifesto passaram.

Um aviso de aquecimento interno de biblioteca de planilhas apareceu no stderr, mas não alterou os códigos de saída nem os gates aprovados. A compilação Lean não foi executada e continua não verificada.

## Estado real do GitHub

**A incorporação integral dos 29 arquivos ainda não está concluída.** Esta integração só publicou a documentação de acompanhamento e o manifesto de hashes; as limitações do conector atual não permitem transferir aqui os bytes do ZIP e de todos os arquivos de evidência como uma operação de upload binário. Não declarar que os JSONL, o ZIP E7 ou o PDF já estão no repositório.

- Workflow V14 ainda está em `workflows_staged/`, não ativo em `.github/workflows/`.
- Nenhuma execução de GitHub Actions V14 foi realizada.
- Nenhum PR foi aberto nem mesclado.
- `main` permanece sem alterações desta integração.

## Passo para completar a incorporação

No clone local de `https://github.com/jbarumby/scicalc.git`:

```powershell
git fetch origin
git switch integration/e8k-v14-package-2026-10-28
```

Extraia o ZIP completo na raiz do clone (preserve os caminhos de `.github/`, `b5/`, `proveniencia_b3/`, `oraculo/`, `referencias/`, `NUMERIA_E8K_FORMAL/` e os arquivos na raiz). Compare os hashes com `SHA256SUMS.txt`; execute os comandos do README; depois:

```powershell
git status --short
git add .github b5 proveniencia_b3 oraculo referencias NUMERIA_E8K_FORMAL
git add CONFERENCIA_E8K_V14.json CONFERIR_E8K_V14.py GERAR_E8K_V14.py GERAR_SHA256_V14.py MANIFESTO_SHA256.json NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip ORACULO_RACIONAL_V14.py PARECER_E8K_V14.md README_E8K_V14.md REPLAY_E8K_V12_INTEGRAL.json RESULTADO_E8K_V13.json RESULTADO_E8K_V14.json SHA256SUMS.txt VERIFICAR_E8K_V12.py VERIFICAR_E8K_V13.py VERIFICAR_INVARIANTES_LEAN_V13.py
git commit -m "feat(e8k): incorporate V14 replay and rational witnesses"
git push origin integration/e8k-v14-package-2026-10-28
```

Só após o push, confira os caminhos no GitHub, ative os workflows e abra um PR para revisão.

## Limites científicos

Os testes não constituem prova da Hipótese de Riemann. A compilação Lean e as premissas analíticas continuam precisando de verificação independente.
