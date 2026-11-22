# NUMERIA P80 — Relatório da execução reproduzida

**Workflow principal:** [GitHub Actions run #45030351631](https://github.com/jbarumby/scicalc/actions/runs/45030351631). **Job `certify`: completed / success**, sem falhas. Executado em 22/11/2026.

**Artefato:** [numeria-p80-v3-s3-source-gates](https://github.com/jbarumby/scicalc/actions/runs/45030351631/artifacts/14756893521), ID `14756893521`: `p80-v3-source.json`, `p80-s3-192.json`, `p80-s3-256.json`, `p80-p79-original-b1.json`, `p80-p79-512-cover.json`, `p80-sha256sums.txt`.

## Verificações aprovadas

- **Fontes matemáticas:** o limite de Ramaré–Zúñiga v1 `S_ε(X)<=0.526` para `10^12<=X<10^16`, `ε<=1/100` está no §6 (página impressa 21), apesar de não estar listado no Teorema 1.2. A comparação com v3 usa **a melhor cota por região**, não substitui globalmente v1 por v3.
- **Racionais exatos:** 16/16 linhas da P78 reproduzidas exatamente, com log lower bounds `28,36,75`, e comparação `max` entre as três regiões de X. O código produz **apenas uma redução conservadora**: linha 10, `b4=557/1000 -> 139/250 (=0.556)`.
- **Protocolo S3 seguro:** 16 linhas confrontadas com o auditor de B1 P79 a **192 e 256 bits**, sempre **9 preservadas / 7 substitutas**. As duas execuções deram a **mesma tabela P79** SHA-256 `d7ba517b4461e4791ad1ca92d69575d4f973b0d161b5a2a50b928ea2137af615`.
- **Reprodução Suzuki:** usando a tabela P79 **sem mudanças**, 500 zeros, 256 bits e 512 bandas de densidade, a cobertura condicional finita `18<=t<=79` foi integralmente reproduzida em **4.284 células**. Esta é **reprodução** da P79, não uma nova extensão certificada.
- **Integridade do arquivo S3 anexado:** SHA-256 `1faf315454a100a6a9850e843c795eca2bf734053ec9d32f202286709565b686` corresponde ao original fornecido, com a linha final preservada. A CI testa esse SHA e a sintaxe Python. Os arquivos históricos permanecem intocados.
- **CI normal e `python -O`:** testes de seleção de faixas, domínio de ε, comparação das versões e proveniência aprovados. O verificador falha explicitamente em situações de domínio inválido.

## Limitações importantes

O artigo Ramaré–Zúñiga Alterman (v1/v3) e Chourasiya–Simonič v2 permanecem fontes matemáticas **importadas**. Os scripts conferem domínio, implicações racionais e certificados numéricos, mas **não demonstram de forma independente os teoremas externos**.

A prova pública P80 de substituição foi construída **somente** para `X>=3×10^12` e `η=2δ/log X` com `0<δ<=0.303`. O Lema 6 original cobre `X>=3×10^10`; usar a P80 fora da janela parametrizada requer uma nova ponte formal.

A P79 oferece uma verificação **finita, condicional** de positividade. Isso **não** é prova da Hipótese de Riemann global.

### Próxima prioridade científica (P81, ainda não iniciada)

Auditar a transferência exata, termo a termo, das estimativas de `S_η(X)` ao Lema 6 CS, inclusive:
- identificação de `τ=1+η`, `η=2δ/log X`, os fatores `ζ(τ)e^{-γη}` e `C3/C4`;
- faixa `X>=3×10^12` que é suficiente para os parâmetros específicos da P79, versus o domínio original `X>=3×10^10`;
- a substituição simultânea das instâncias do Lema 6 na cadeia Teorema 1 → Corolário 1 → tabela P79;
- verificação das hipóteses de cada interpolação; revisão de fontes e de arredondamentos.

Todas as contas devem permanecer reprodutíveis e condicionais a insumos externos ainda não auditados independentemente.
