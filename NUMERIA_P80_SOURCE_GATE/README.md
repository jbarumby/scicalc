# NUMERIA–RIEMANN P80 — versões v1/v3 e protocolo S3 seguro

**Base:** P79, branch `numeria-p79-b1-public-replay`. PR em rascunho; sem merge ou alteração retroativa dos certificados históricos.

## Evidência das fontes primárias

- [Ramaré–Zúñiga Alterman v1](https://arxiv.org/pdf/2603.25961v1), Teorema 1.2, p. 2, **e** limite adicional `S_ε(X) <= 0.526` para `10^12<=X<10^16`, `ε<=1/100`, explicitado nas contas do **§6, p. 21**. O `0.526` não está na tabela do Teorema 1.2!
- [Ramaré–Zúñiga Alterman v3](https://arxiv.org/pdf/2603.25961v3), Teorema 1.2, p. 2: para `10^12<=X<10^16`, `ε<=1/50`, passa de **0.5503 (v1)** para **0.5501 (v3)**. Em `10^16<=X<10^33`, `ε<=1/50`, passa de 0.5104 a **0.5103**. Mas para `ε<=1/25` na primeira zona **piora**, de 0.593 (v1) para 0.5936 (v3); por isso não adotar v3 globalmente sem comparação por zona!
- [Chourasiya–Simonič v2](https://arxiv.org/html/2507.15184v2), Lema 6, equações C3/C4. A hipótese do artigo é `X>=3*10^10`; nosso gate de substituição pública é explicitamente restrito a `X>=3*10^12`. O preprint CS cita a comunicação privada de Ramaré para a constante `b4=0.54`, ainda não demonstrada publicamente na forma e na totalidade do domínio do Lema 6.

## Resultado P80 reproduzível

O script `v3_gate.py` reconstrói **exatamente em frações racionais** o certificado P78 v1 para as 16 linhas, rejeitando qualquer divergência. Compara as melhores cotas publicadas em v1 e v3 dentro de cada subintervalo de `X` e restrição sobre `ε=2δ/log(X)`. Para a zeta, usa `ε ζ(1+ε)e^{-γε} <= 1+ε`.

Os limites inferiores racionais de log nas três faixas são `28,36,75`, demonstrados sem log flutuante a partir de `e<11/4` e potências inteiras.

**Ganho:** linha **10** admite `b4=0.556`, em lugar de `0.557`; nas demais 15 linhas, o maiorante escolhido da P78 é preservado por política (não requer reconstrução de B1). Em particular, o cálculo da linha 16 mantém a cota de v1 `0.593`, em vez da mais fraca v3 `0.5936`.

Como todas as cotas novas `b4_new <= b4_P78` e a P79 demonstrou monotonicidade de B1 em b4 na faixa pertinente, **a tabela de B1 da P79 continua sendo um majorante válido sob os mesmos insumos analíticos externos**. Não revendicar o ganho de b4 como nova extensão de t sem recalcular e certificar essa extensão.

## Arquivo S3 recebido

Original integral preservado em `source_input/p78_s2_full_protocol.py`. Ele é um **protótipo**, não o certificado da P80.

Fragilidades da versão recebida:
1. o `gate()` da bisseção de `b4*` é chamado **após encerrar** o `ctx.workprec`; a conversão de midpoints também ocorre fora do contexto. Os limites `b4*` são exportados em `float`, e o sucesso de um `gate()` não fornece sozinho uma cota superior para a raiz;
2. o `patch_all` mistura `action=keep_current` (sem campos `b1_new_*`) e `action=reprint`. O recertificador S2 antigo rejeita a linha `keep_current` quando passado diretamente ao parâmetro `--proposed`;
3. campos `dyadic` contêm valores `float`, não endpoints diádicos exatos; precisões distintas não são comparadas como gate obrigatório;
4. a opção `--rows` ignora silenciosamente índices fora da faixa quando combinados a válidos;
5. o texto original sugere editar a `TABLE`, o que contradiz a política do NUMERIA de preservar os certificados históricos.

A variante segura `s3_safe.py` reutiliza **os certificados Arb P79 em 192/256 bits** e classifica todas as 16 linhas sem bisseção de raiz, extraindo expoentes `e=2i/(48-i)` como frações exatas. Emite **overlay** separado dos valores `B1` e conserva a tabela histórica. Todas as 16 linhas podem ser verificadas com a mesma API.

## Reproduzir

```bash
python -m pip install python-flint==0.9.0 mpmath==1.3.0
python -m unittest -v NUMERIA_P80_SOURCE_GATE.test_p80
python -O -m unittest -v NUMERIA_P80_SOURCE_GATE.test_p80
python -m NUMERIA_P80_SOURCE_GATE.v3_gate --out p80-v3-source.json
python -m NUMERIA_P80_SOURCE_GATE.s3_safe --precision 192 --out p80-s3-192.json
python -m NUMERIA_P80_SOURCE_GATE.s3_safe --precision 256 --out p80-s3-256.json
```

A CI compara e preserva hashes das entradas e saídas e reproduz também a **P79 512 faixas**, sem trocar a tabela da cobertura.

**Limitações:** estimativas de ambos os artigos continuam importadas; reproduzir frações e intervalos não prova integralmente as demonstrações dos artigos. A análise pública continua restrita ao domínio parametrizado, em vez da abrangência original do Lema 6 CS. O resultado `18<=t<=79` permanece um certificado **finito, condicional**, e não constitui demonstração da Hipótese de Riemann global.
