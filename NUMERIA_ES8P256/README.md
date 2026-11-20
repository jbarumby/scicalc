# NUMERIA ES8-P256 — mapa da dependência `EveryCharacter → RH`

**Objetivo:** avançar a auditoria P255 do nível dos cinco enunciados finais para um mapa fonte-pinado dos elos imediatamente abaixo da prova Lean de RH. Isso identifica onde uma investigação matemática independente precisaria começar. Não produz nem reivindica uma prova independente da Hipótese de Riemann.

## Resultado matemático da leitura

A definição Lean `UpperTheory.EveryCharacter` afirma que, para todo módulo e todo caráter de Dirichlet, a função L não se anula quando `Re(s)>1/2`, exceto no polo trivial. `allChi` fornece essa proposição; `rhOfAll` aplica `Formalism.VW.rh_of_every`. Esse último passo usa o caráter trivial módulo 1 para obter a não anulação de ζ à direita da linha crítica; `rh_of_right` usa a equação funcional para refletir zeros e concluir RH.

O gargalo relevante não está nessa redução final. O teorema `UpperTheory.FinalTheorems.everyCharacter` deriva a não anulação universal a partir de `galoisLocalGlobalStrong_recParameter`. O passo intermediário `vwin_of_row` converte uma hipótese `Row0b` em `Top.VWin`; em seguida, `every_character_of_literature_of_vWin` usa um teorema de correspondência de integrais theta e matching local para obter `EveryCharacter`. Esse módulo também combina várias hipóteses e teoremas auxiliares de geometria aritmética, representações e teoria automórfica. O arquivo de quase 100 KB identifica o próximo alvo de leitura, mas o presente verificador não inspeciona nem reproduz seus corpos importados.

Isso coloca a ponte com Numeria em termos concretos: os certificados P249–P254 tratam de operadores, formas de Weil e janelas finitas; não demonstram `Row0b`, `VWin` nem a não anulação para todos os caracteres, módulos e pontos do semiplano. A próxima meta de pesquisa que pode contribuir à rota própria do Numeria é obter uma desigualdade analítica uniforme — para uma classe explícita de funções teste e todos os módulos/caracteres — que implique a mesma forma de não anulação. Uma identidade entre o determinante de Stepanov/`W(1_S)` e as matrizes de Weil do Numeria continua por provar; não se deve assumir essa ponte.

## Verificação externa do estado formal

Na última atualização visível da página de auditoria de Andrew Sutherland (MIT), datada de 20/11/2026 19:45 UTC, a página relata: compilação dos 175.187 módulos; auditoria de axiomas dos cinco teoremas finais reduzida aos três axiomas padrão; replay de todos os módulos pelo `leanchecker`; e aceitação da exportação completa pelo kernel independente `nanoda`. Naquele registro, o replay da exportação pelos kernels padrão/paranoico do Lean e uma execução ponta a ponta do Comparator ainda estavam em curso. [Fonte MIT](https://math.mit.edu/~drew/grh/)

Esse estado torna a formalização muito mais corroborada do que a mera alegação inicial. Ainda assim, o presente pacote só audita texto e proveniência de uma cadeia pequena: **não executa Lean, não reproduz o kernel do MIT, não confere os corpos de todos os teoremas importados e não constitui prova própria do Numeria**. O repositório de referência e os sete módulos copiados estão fixados ao commit `e50d09817ad81469f737e8f637561e6ccf86c670`; cada byte usado pelo verificador tem SHA-256 registrado no script.

## Reprodução local

Na raiz do repositório:

```bash
python -m unittest NUMERIA_ES8P256.test_dependency_chain -v
python -O -m unittest NUMERIA_ES8P256.test_dependency_chain -v
python -m NUMERIA_ES8P256.audit_dependency_chain --out /tmp/P256_DEPENDENCY_AUDIT.json
```

Os testes verificam hashes, elos específicos da redução e falham de forma fechada se houver alteração das fontes. Eles são uma auditoria estática de dependências, não uma checagem formal do Lean. Veja `P256_DEPENDENCY_AUDIT.json` e `lean_sources/` para o resultado e os trechos fonte completos.
