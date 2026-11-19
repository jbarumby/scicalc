# ES8-P246 — auditoria independente dos novos certificados

## Resultado da auditoria

O ZIP recebido passou a auditoria do manifesto: todos os 14 arquivos listados em `SHA256SUMS_P246.txt` conferem, e os caminhos do arquivo são seguros. O verificador `p246_verify_uploaded.py` foi executado em modo autenticado. Seu resultado local coincidiu byte a byte com `P246_AUDIT_CERT.json`; os seis testes do pacote passaram em Python normal e otimizado.

O conteúdo do ZIP, isoladamente, não incluía os seis blobs-fonte P237 de aproximadamente 60 MB no total. Porém, os blobs estavam disponíveis no diretório local de trabalho `/home/ubuntu/NUMERIA_ES8P237_LOCAL/chunks`. Verifiquei os seis SHA-256 embutidos no replay P245D, bem como o hash dos brackets GL60. Todos coincidiram com os valores esperados.

Reexecutei o código P245D recebido, usando esses blobs, os brackets versionados em `numeria/p245/gate_c/P245B_ROOT_BRACKETS.json` e os módulos P245B/P237 do repositório. Cada um dos três intervalos `[25,50)`, `[500,525)` e `[975,1000)` terminou com 26 somas selecionadas independentemente envolvidas. As três saídas recém-geradas são **idênticas byte a byte** aos JSON P245D anexados. O total é 78/78 somas e 117.000 contribuições individuais de GL60 reexecutadas conforme a seleção do script.

Também regenerei os certificados analíticos P237. `P237_ANALYTIC_GATES.json` e `P237_GL60_ANALYTIC_CERT.json` foram estruturalmente idênticos aos arquivos recebidos. O verificador reconstruiu por frações os 754 limites de base; as somas e máximos por setor coincidiram com os valores arquivados. O limite complexo do símbolo melhora de `50` para `274/7`, reduzindo o orçamento de linha GL60 pelo fator exato `(274/7)/50 = 137/175`. Os novos limites ficam abaixo de `4×10⁻⁴⁷` em ambos os setores.

## O que estes resultados significam

A auditoria agora é mais forte que a checagem de inclusões armazenadas: os seis blocos-fonte foram autenticados por hash e o integrando GL60 selecionado foi efetivamente reexecutado. Ainda assim, a reprodução cobre apenas 13 pares matriciais por setor em três grupos de painéis. Os 18 casos ligados ao índice 376 usam o majorante conservador do modo alto descrito pelo pacote. Nenhuma das etapas refaz a matriz inteira de 377×377 por setor ou todas as 1.008 integrais contínuas com uma implementação independente.

O aperfeiçoamento `274/7` aperta o orçamento GL60 sob os lemas analíticos e a identificação do integrando declarados pelos certificados. Não altera por si só as caixas das entradas P237 nem prova a identidade total entre a matriz calculada e a forma de Weil original. A positividade global de Weil e a Hipótese de Riemann permanecem não provadas.

## Reprodução

Da raiz do repositório, com os seis arquivos P237 disponíveis em `NUMERIA_ES8P237_LOCAL/chunks`:

```bash
PKG=numeria/p246/uploaded/package
DATA=/home/ubuntu/NUMERIA_ES8P237_LOCAL/chunks
ROOTS=numeria/p245/gate_c/P245B_ROOT_BRACKETS.json
export PYTHONPATH="numeria/p245/gate_c:numeria/p237"
for START in 25 500 975; do
  STOP=$((START+25))
  python3 "$PKG/sources/p245d_additional_chunk_replay.py" \
    --references "$DATA" --roots "$ROOTS" --start "$START" \
    --out "/tmp/P246_REPLAY_${START}_${STOP}.json"
done
```

As saídas byte-idênticas e o resumo de hashes estão em `numeria/p246/uploaded/reproduced/`. O pacote recebido, sem alterações, está preservado em `numeria/p246/uploaded/package/`. Os seis blobs de aproximadamente 60 MB não foram adicionados ao repositório; a reprodução acima precisa dos arquivos locais descritos e o relatório não os apresenta como conteúdo do ZIP.

Elaborado por Manus AI para o projeto NUMERIA.
