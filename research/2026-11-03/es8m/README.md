# NUMERIA–RIEMANN ES8-M — Hankel, Stieltjes e momentos gaussianos

Data: 03/11/2026. **Status científico: verificação finita e parcial; não prova a Hipótese de Riemann.**

A ES8-M incorpora, após auditoria, dois manuscritos recebidos (DeepSeek e Qwen). Registra um novo resultado computacional: intervalos **racionais exatos** para os momentos completos `mu0` e `mu1` da fórmula explícita gaussiana em `a=1/50`, combinando componentes previamente cercadas em ES6/ES7 e uma nova certificação para a derivada do termo primo. Eles atestam somente a positividade escalar `H_0` e `H_0^{[1]}`; **Hankel/Stieltjes de dimensão ≥2 para a zeta NÃO estão certificados**.

Há ainda matrizes-modelo artificiais de dimensão 1 a 4 (LDL exato com frações) e uma comparação **não intervalar** de matrizes formadas com os primeiros 32 zeros críticos conhecidos, a 45 dígitos.

## Arquivos

- `es8m_moment_audit.py` — computação, certificados racionais, testes adversariais, diagnóstico numérico de momentos altos.
- `es6_weil_gaussian_baseline.py` — cópia da rotina ES6, **com correção no ramo par/ímpar da cota de arctangente**, sem modificar o original já versionado em ES6.
- `ES7_RESULTADOS_BASELINE.json` — intervalos ES7 usados como insumo, com exatidão decimal preservada nos inteiros escalonados.
- `test_es8m_moment_audit.py` — 14 testes, incluindo erros de arredondamento, fontes independentes e uso indevido da positividade.
- `RESULTADOS_ES8M.json` — intervalos racionais, controle de modelos, diagnósticos.
- `RELATORIO_ES8M.md` — hipóteses, provas de cauda, fontes e resultados.
- `sources/` — dois textos originais, **não validados** como teoremas.
- `MANIFESTO_SHA256.json` — integridade e versões.

## Reproduzir

```bash
python -m pip install mpmath==1.3.0 pytest==8.4.2
python -m pytest -q test_es8m_moment_audit.py
python es8m_moment_audit.py --zeros 32 --degree 5 --dps 45 --out RESULTADOS_ES8M.json
```

A computação usa o teorema clássico de Guinand–Weil como identidade **incondicional**: todos os zeros (inclusive eventuais zeros não críticos) pertencem à soma global. Os valores de `zetazero(n)` não certificam que o conjunto utilizado contenha todos os zeros. A estrutura de Stieltjes permanece uma **condição necessária sob RH** neste experimento; a recíproca sugerida nos textos não foi demonstrada.

**Prioridade seguinte:** certificar `mu2,mu3` da expressão analítica completa, depois fazer LDL intervalar de `H1` e `H1_shift`, incluindo análise de condicionamento e comparação normalizada com E7/V24.
