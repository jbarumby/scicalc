# NUMERIA-Q — Marco de reprodução independente de divsteps

**Branch:** `research/divsteps-reproduction-2026-10`  
**Estado:** protocolo e testes de referência inicial; reprodução do circuito publicado ainda não executada.

## Objetivo

Criar uma trilha verificável para comparar o protótipo NUMERIA-Q 05.21 com o circuito reversível de inversão modular publicado por Gu, Kemeny, Yelin e Lukin (2026). Não reivindicamos descoberta, equivalência com o circuito publicado, nem vantagem quântica nesta etapa.

## Referência primária

- Artigo: https://arxiv.org/abs/2610.22367
- Código, circuitos e avaliador dos autores: https://github.com/andigu/uncompute
- Tarefa relevante: `tasks/modinv-point-add`
- Circuito de referência: `tasks/modinv-point-add/circuits/divsteps/`
- Prova de recuperação de decisões: `lean/divsteps`

O repositório dos autores informa para a construção divsteps padrão aproximadamente 4.006.673 portas Toffoli, 935 qubits de pico e custo Toffoli × qubits de 3.746.239.255. Esses valores são referências publicadas, não resultados reproduzidos por este branch.

## O que foi implementado neste marco

1. Testes independentes da referência clássica de inversão modular: resíduos de primos pequenos, entrada sem inverso, caso grande com 59 bits e exemplo pequeno com coeficientes intermediários negativos.
2. Uma especificação de reprodução, com os comandos do avaliador oficial.
3. Critérios explícitos para impedir que os testes clássicos sejam confundidos com a reprodução do circuito quântico reversível.

## Ambiente recomendado

O projeto upstream requer Python e `uv`; consulte o README upstream para a versão e dependências exatas. Execute em um clone separado, sem copiar arquivos para `main` nem alterar as versões congeladas 05.19.1 e 05.20.

```bash
git clone https://github.com/andigu/uncompute.git
cd uncompute
uv sync
uv run python tools/eval.py modinv-point-add --circuit tasks/modinv-point-add/circuits/divsteps/circuit.py
```

O avaliador pode levar minutos ou horas para circuitos grandes. Guarde o log integral, o commit upstream, a versão do Python/uv e o resultado certificado. Não declare reprodução concluída apenas porque o comando terminou: confira o veredito, as contagens certificadas e o hash do fluxo de operações.

## Critérios de aceitação

- [ ] Testes da referência clássica passam.
- [ ] O commit upstream e o hash dos arquivos usados estão registrados.
- [ ] O avaliador upstream aceita o circuito divsteps e emite as contagens certificadas.
- [ ] A prova Lean `lean/divsteps` é executada e seu resultado é registrado.
- [ ] O protótipo NUMERIA-Q e a referência são comparados com a mesma especificação e convenções de custo.
- [ ] Qualquer diferença de resultados é reproduzida e explicada antes de conclusões científicas.

## Limites e transparência

Os testes locais deste branch verificam aritmética modular clássica, incluindo casos de borda; não verificam o pior caso de iterações nem inspecionam o sinal do parâmetro interno `delta` de divsteps. Isso exige testar a implementação divsteps real e sua instrumentação, ainda não definida neste arquivo. A reprodução do circuito e a execução da prova Lean são tarefas pendentes. Não se afirma novidade científica, ganho de recursos, vantagem quântica ou validação independente do artigo.

