# NUMERIA–RIEMANN P259 — Gate B: setor ímpar completo

Base: PR #65, certificado de uma testemunha ímpar da P258 recuperada por piso escalonado P259. **Este Gate B não presume positividade.**

O código calcula, com Arb, os prefixos principais da matriz ímpar de 377 modos:

\[
R_{\rm P259}=R_{\beta}+(1/4-\beta)(I-G_{504})+(1/4)(I-G_{850}).
\]

A quadratura ponderada original vem do arquivo P253 com SHA fixado no código P258. As matrizes de concentração \(G_{504}\) e \(G_{850}\) são construídas pela **identidade analítica de extremos Bessel**, não apenas por soma numérica; seus arredondamentos são intervalares em Arb. A quadratura ponderada conserva sua cota analítica GL60, aplicada uma vez.

**Verificações obrigatórias:** cada tentativa de LDL retorna POSITIVO, NEGATIVO_CERTIFICADO ou INCONCLUSIVO. Um pivô negativo aproximado não constitui certificado; a testemunha precisa passar por um segundo cálculo de Rayleigh usando frações exatas e caixas intervalares. A busca é interrompida ao encontrar uma direção estritamente negativa ou um prefixo inconclusivo; não confundir falta de decisão com positividade.

**Escopo:** não reutilizar a constante de dilatação da P258 ao alterar a cauda para piso 1/2 após 850. Não reivindicar conclusão para todas as larguras acima do primo 7 a partir de um resultado numa largura fixa. Mesmo um bloco ímpar finito positivo necessita de teste Schur modal e de ambos os setores antes de justificar a transferência para a forma completa de Weil.

## Execução

No repositório, com python-flint conforme `NUMERIA_ES8P249/requirements.txt`, e com `P258_GL60.json` baixado do artefato GitHub Actions `44858303395`:

```bash
python -m unittest NUMERIA_ES8P259.test_gate_b -v
python -O -m unittest NUMERIA_ES8P259.test_gate_b -v
python -m NUMERIA_ES8P259.gate_b --p258 replay/P258_GL60.json \
   --max-dimension 96 --bits 1536 --out replay/P259_GATE_B_PREFIX.json
```

Se **todos** os prefixos até 96 forem positivos, continuar para `--max-dimension 377 --bits 2048`. Os resultados são JSON com proveniência e estados claros. O GitHub Actions executa um gate leve primeiro; o gate completo é uma etapa separada, para não confundir sucesso de testes com prova da positividade.

## Próxima decisão

- Direção negativa certificada: preservar vetor diádico e Rayleigh racional, investigando seus suportes e faixas.
- Positividade de 377 modos: usar limite de Schur e auditoria independente.
- Inconclusivo: aumentar precisão ou construir um precondicionador intervalar; não usar autovalores float para decidir sinais.
