# NUMERIA–RIEMANN E8-K V12 — replay integral de nós do B5

## Resultado reproduzido em 28/10/2026

A auditoria percorreu todos os **15.754 registros JSONL** (G11 9.524; G22 6.230) do ZIP E7 fixado pelo SHA-256 `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7`. A fonte B5 fixada é `a1753d3010e34f29371a9706cd6bb05b3b0c8491ca3cafb44bf98b7cd53d3df9`.

- Cada ponto médio foi confrontado com a malha B3, a série de registros B5, o valor reproduzido por `eval_point` do **mesmo** B5 e o produto intervalar com a largura da célula.
- Todas as somas das 11 faixas e ambas as somas finitas coincidiram exatamente com B5, inclusive por soma racional independente dos extremos exportados.
- A soma dos erros de quadratura e cauda **declarados** reproduziu os intervalos compostos de E7, em aritmética `Fraction` exata.
- **Limitação epistemológica**: regenerar `eval_point` com a mesma implementação B5 não constitui oráculo independente e não demonstra que as cotas analíticas originais sejam verdadeiras.

## Proveniência dos B3 originais

Os hashes dos B3 originais referidos pelo resultado B5 são:
- G11: `65a467c925b3824a23e814ab49d5ce456e528ce165e72a8f2b5dc8ab62368dc2`
- G22: `8f53deef567ab58134436fff3bab38d8049fedda32f88f8fc7f587723b4dd54f`

As versões B3 originais e corrigidas têm malhas iguais, mas bytes diferentes. A auditoria de proveniência identifica os arquivos originais e não confunde seus hashes com os das versões corrigidas.

## Próxima pendência (V13)

1. Validar independentemente as cotas globais da segunda derivada usadas nos erros da regra do ponto médio.
2. Auditar a hipótese, a normalização e a prova da cauda integral r>160; recertificar a equivalência da forma de Weil.
3. A formalização Lean permanece separada: código com `sorry` e sem compilação não é prova mecanizada.

**Não é prova da Hipótese de Riemann.**