# Validação local do pacote E8K V12

Data: 2026-10-28

O pacote de origem `E8K_V12_GITHUB_MELHORADO_COMPLETO.zip` foi extraído e executado localmente antes do staging no GitHub.

- ZIP íntegro: `ZipFile.testzip() == None`.
- `py_compile`: passou para `VERIFICAR_E8K_V12.py` e `CONFERIR_E8K_V12_V2.py`.
- Replay V12 estrito: código 0; `ALL_POINTS_REPRODUCED_B5_INPUT_HASHES_MATCHED`; 15.754 pontos; zero lacunas de proveniência.
- Conferidor V2 no modo `--root`: 8/8 gates PASS, incluindo 6/6 mutações isoladas detectadas.
- O aviso de warmup de planilha emitido pelo ambiente Python não impediu a compilação nem os dois comandos acima de concluírem com código 0.

## Interpretação e limite

Esses resultados foram executados contra o ZIP completo local, que contém o ZIP E7 binário e os artefatos originais. **Não são resultados de GitHub Actions** e não validam automaticamente os arquivos parciais atualmente staged na branch. A cópia binária do E7 continua pendente no GitHub; os workflows foram guardados sob `numeria_e8k/v12/workflows/` para não ativarem uma execução incompleta.

A execução comprova reprodução e consistência de evidências computacionais; não demonstra independentemente as cotas analíticas da quadratura ou da cauda e não prova a Hipótese de Riemann.
