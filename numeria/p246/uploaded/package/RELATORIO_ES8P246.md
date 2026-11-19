# NUMERIA–RIEMANN — ES8-P246: revisão dos nove arquivos recebidos

**Escopo:** auditoria independente por inteiros/frações das inclusões *arquivadas* nos JSONs P245D e rederivação do orçamento analítico P237. **Não** é uma reexecução das 117.000 contribuições GL60 nem uma reprodução do operador de Weil–Zhu.

## 1. Materiais aproveitados e procedência

Foram vinculados por SHA-256 os nove arquivos enviados: três certificados de blocos `P245D_SELECTED`, três certificados analíticos P237, o código de replay P245D e dois testes. `p246_verify_uploaded.py` rejeita qualquer alteração nesses nove arquivos ao operar em modo autenticado.

Os três novos blocos são `[25,50)`, `[500,525)` e `[975,1000)`, com 25 painéis, 60 nós GL60, dois setores e 13 pares matriciais selecionados. Cada certificado registra **39.000 contribuições de nó × paridade × par de índices** e **26 inclusões de somas intervalares**. Nos três certificados são 117.000 contribuições *registradas* e 78 somas *arquivadas*.

Recalculei com `fractions.Fraction`, e não usando os booleans `source_Arb_interval_contained`, os dois extremos de cada bola Arb diádica, conferi sua representação como fração e verifiquei sua inclusão na caixa do replay. Resultado **78/78 inclusões válidas como operações racionais nos dados fornecidos**.

**18 dessas 78 comparações** envolvem o índice 376 (graus 752/753). As caixas registradas para esses casos são largas: aproximadamente `[-1500,+1500] · 2^-256`, ao contrário dos intervalos ajustados nos graus baixos/intermediários. Elas são limites conservadores, mas não representam recálculo preciso de Bessel para o grau alto. Restam **60 comparações com intervalos mais informativos**.

## 2. Melhoramento do orçamento GL60

O certificado `P237_SYMBOL_E3_RATIONAL_GATE.json` dá uma cadeia de desigualdades racionais, sujeita às identidades analíticas clássicas declaradas, que permite usar

`|Psi_L(t)-beta_504| < 274/7 < 48 < 50`, em todas as elipses Bernstein E3 de 1008 painéis.

A geometria `|Re t|≤3025/6`, `|Im t|≤1/3` e `|Lt|<406` foi recalculada em frações. Reconstituí, grau a grau, os **754 majorantes**

`a_n = ceil(sqrt(2 L (2n+1))) * min(2, 2*406^n/(2n+1)!!)`

com `L=161/200`, sem importar do certificado a soma nem o máximo dos `a_n`. Os 377 modos de cada setor produziram **somas racionais exatas idênticas** às de `P237_GL60_ANALYTIC_CERT.json`. A igualdade inteira do antigo orçamento também foi checada:

`R_old = (1008*50/3^120) * max(a_n) * sum(a_n)`.

Como o mesmo argumento admite `274/7` no lugar de `50`, um novo orçamento estritamente superior para o erro verdadeiro é obtido multiplicando o anterior por `137/175`:

| Setor | Antigo R (aprox.) | Novo R (aprox.) |
|:--|--:|--:|
| Par | 3,661198905758 × 10^-47 | **2,866195714793 × 10^-47** |
| Ímpar | 3,657549585347 × 10^-47 | **2,863338818243 × 10^-47** |

Todos os sinais/comparações foram decididos por frações exatas. Os decimais apenas apresentam os valores. Esse aperfeiçoamento é **condicional à validade dos lemas clássicos citados e à identificação da função quadraturada**; ele não altera por si só os arquivos Arb e hashes P237 preexistentes.

## 3. Limites de reprodução e problemas detectados

1. **Não estão presentes os seis blocos-fonte** `P237_CHUNK_{even,odd}_{0025_0050,0500_0525,0975_1000}.json.gz`. Os JSONs declaram os hashes desses arquivos, mas não é possível verificar os bytes originais nem reexecutar `evaluate_chunk` a partir deles apenas com os nove anexos. Uma consulta à árvore do branch `numeria-p228-arb08-rational-audit` do GitHub não localizou esses seis arquivos. Ausência naquela árvore não exclui outras branches.
2. O teste recebido `test_p245d_additional_chunk_replay.py` foi executado com as dependências importáveis e passou **4/4**. Contudo, testa somente os **registros JSON arquivados**; não chama `evaluate_chunk()` nem exige abrir os seis blocos P237. Logo seus quatro testes não constituem uma reexecução independente da quadratura.
3. O teste `test_p237_analytic_certificates.py` supõe uma versão de `p237_analytic_gates.py` que disponibiliza `need`. O arquivo recuperado de um pacote histórico `NUMERIA_ES8P237_GATE_C_ANALITICO.zip` não contém essa função e falha em `ImportError`; o arquivo atual da branch GitHub **contém** `need`. Trata-se de incompatibilidade de versões, não de refutação dos lemas. Fixar commit/manifesto é necessário para testar corretamente.
4. `p245d_additional_chunk_replay.py` resolve `P237_DIR` por posição relativa. Em uma instalação independente, o caminho pode não coincidir com a localização dos módulos. Recomenda-se parametrizar `--p237-dir` ou usar imports de pacote estáveis.
5. Os dois certificados `P237_SYMBOL_E3_RATIONAL_GATE.json` e `P237_GL60_ANALYTIC_CERT.json` foram ainda **regenerados com igualdade estrutural exata** usando os scripts históricos fornecidos em pacote anterior. Para `P237_ANALYTIC_GATES.json`, o pacote histórico tem schema `.1`, enquanto o arquivo atual é schema `.2`; o script exato da versão `.2` precisa ser vinculado para replay fiel.

## 4. Gate seguinte recomendado

**Prioridade 1:** obter os seis arquivos P237_CHUNK nomeados, verificar os seis SHA-256, e executar `p245d_additional_chunk_replay.py` para comparar *saída produzida* e certificados P245D existentes. Rodar também os testes com versões de código fixadas por commit.

**Prioridade 2:** substituir os graus 752/753 atualmente incluídos em caixas de largura de cerca de `2.6×10^-74` por uma recorrência/serie rigorosa de Bessel que informe a magnitude, não apenas a inclusão trivial.

**Prioridade 3:** ampliar a seleção para entradas adicionais e, por fim, reconstruir todas as 142.506 entradas triangulares, com orçamento de erro analítico e hashes de artefatos. Só então reavaliar de ponta a ponta a identificação com Weil–Zhu e a matriz de positividade finita.

### Alcance do resultado

**Validado neste pacote:** 78 inclusões exatas dos certificados arquivados, 754 majorantes independentes, fator rigoroso `137/175` no orçamento GL60 sob hipóteses analíticas clássicas, 6/6 testes P246 normais e otimizados. **Não validado:** cálculo original de todos os nós GL60, integridade dos seis chunks ausentes, matriz integral 377×377 completa, identidade de Weil–Zhu e RH. A Hipótese de Riemann permanece em aberto.
