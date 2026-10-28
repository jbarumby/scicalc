# NUMERIA E8-K — CI revisado (sem superdeclarações)

## Instalação
Copie estes arquivos para a raiz do repositório, preservando `.github/workflows/`.
Acrescente **separadamente** o ZIP exato `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip`
(não está incluído neste pacote). Se usar Git LFS, o checkout já habilita `lfs: true`.
O pipeline verifica SHA-256 do ZIP (`4acf99a2...02a7`) e da fonte B5 (`a1753d30...3df9`).

## Por que o workflow original falharia
- V11 sem ZIP: retorna exit 2 com **1 PASS e 5 SKIP**, mesmo no cotejo independente aprovado.
- V11 com ZIP: o adaptador genérico não entende o ponto fixo B5;
  os gates `edge_cases`, `rational_tests` e `b5_pi_bounds_overlap` são SKIP.
  Por isso NÃO usamos `--strict` como gate padrão até consertar o próprio V11.
- O arquivo `fonte_b5.py` não é necessário: extraímos **somente** a fonte do ZIP pinado.
- O teste Hypothesis incluído está adaptado a extremos inteiros em escala 10^42.

## O que falha/bloqueia o CI
- ZIP ou B5 com hash incorreto;
- alteração inesperada no padrão PASS/SKIP do V11;
- falha na auditoria racional/estrutural/constantes do V10;
- falha property-based do B5.
Os SKIP restantes da V11 ficam registrados e não recebem o rótulo de certificação.

## Lean
O arquivo `lean_draft/E8K_Intervalos.lean` é o texto enviado pelo usuário,
**não compilado** e com `sorry`. O job Python registra isso sem apresentar o
rascunho como prova. Integração real com Lean+Mathlib exige `lean-toolchain`,
Lake, compilação por `leanprover/lean-action@v1` e remoção dos `sorry`.

## Segurança e alcance matemático
A fonte B5 é executada após checagem SHA-256; isso **não é sandbox**.
Não execute fontes não confiáveis. Os gates não certificam a quadratura contínua,
a cauda global de Weil nem a Hipótese de Riemann.

## Execução local
`python ci/verify_e7.py --zip NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip --extract-to audit/extracted/AUDITOR_PONTOS_B5.py --report audit/origem_e7.json`

`python ci/v11_gate.py independent --script VERIFICADOR_E8K_V11.py --out audit/resultado_v11_independente.json`

`python VERIFICAR_E8K_V10_UNIFICADO_CORRIGIDO.py --zip NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip --out audit/resultado_core_v10.json`

`E8K_SOURCE=audit/extracted/AUDITOR_PONTOS_B5.py pytest -q -x testes_property_e8k.py`
