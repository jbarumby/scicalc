# Continuação Codex — P75

Base preservada: [PR #74](https://github.com/jbarumby/scicalc/pull/74).

1. Faça o replay do código P74 com 192/256 bits e confirme zeros/SHAs.
2. Compare `cs_v2.py` com a versão **v2** do Corolário 1 e Tabela 1 de https://arxiv.org/html/2507.15184v2. A fonte é preprint; confirme também hipóteses e dependências de Platt–Trudgian/Bellotti–Wong.
3. Rode testes e o workflow P75. Corrija erros de código e formalização antes de confiar em `status`.
4. Execute o certificador completo com `--zeros 500 --precision 256 --segments 128 --end 79`. Exija todas as células positivas e cobertura desde 18.
5. Verifique fator de simetria 2/4, multiplicidades, uso `D(1/2+)`, índices nas bordas fechadas, gama incompleta, curvatura local e autenticação dos 500 zeros.
6. Se inconclusivo, preserve os resultados e isole a primeira célula problemática; teste malha menor ou novas cotas SEM inventar um certificado.
7. Compare com perfil de periódicos auditado sem reaproveitar um limite cuja proveniência esteja questionada.
8. Faça commits apenas nessa linha de pesquisa ou em nova branch descendente; não sobrescreva PR #74 nem mescle automaticamente.