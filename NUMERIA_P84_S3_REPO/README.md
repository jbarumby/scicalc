# NUMERIA P84 — S3: adoção versionada de B1 e recertificação REPO por Arb

A P84 incorpora o patch S2 robusto no módulo `TABLE` da P76 e no consumidor `cs_v2` da P75 como uma **versão explicitamente selecionável**, sem sobrescrever a transcrição literal Chourasiya–Simonič v2.

- `NUMERIA_P76_COEFFICIENTS.table.TABLE`: **permanece** com os valores impressos de CS v2, necessários à reprodução histórica P76/P77/P79.
- `NUMERIA_P76_COEFFICIENTS.table.TABLE_S2_PUBLIC_ROBUST`: nova TABLE, com apenas sete valores `b1` alterados.
- `NUMERIA_P75_CS2.cs_v2.TABLE_V2_S2_PUBLIC_ROBUST`: a mesma tabela S2, pronta para seleção explícita pelos algoritmos de cobertura. `TABLE_V2` fica intacta.
- O programa `audit_s3.py` usa os dados **diretamente** da tabela S2, verifica versão/proveniência contra a P78 e executa o auditor P79 `row_audit` para **todas as 16 linhas** a 192, 256 e 384 bits.
- `build_docs.py` emite **novos** `DOC1_repo_interval.json` e `DOC2_repo_interval.json` baseados em extremos diádicos exatos de Arb. **Não importa resultados float nem envelopes arbitrários** dos documentos anexados.
- `scan_downstream.py` examina o repositório completo, ao contrário do JSON original, que examinou 17 arquivos em `/tmp/kit_p78`.
- A CI reexecuta a cobertura Suzuki com a **tabela S2 robusta de duas casas**, 500 zeros e 512 faixas. Não confundir com a cobertura P79 até `t=79`, que foi certificada com a tabela de seis casas.

## Valores robustos S2 (7 linhas)

| linha | CS v2 impresso | TABLE S2 robusta |
|---|---:|---:|
| 10 | 17.21 | 17.49 |
| 11 | 19.31 | 19.67 |
| 12 | 21.77 | 22.25 |
| 13 | 24.83 | 25.48 |
| 14 | 29.01 | 29.89 |
| 15 | 35.35 | 36.59 |
| 16 | 46.06 | 51.65 |

## Origem dos seis arquivos recebidos

**Originais preservados byte a byte e SHA-256 conferidos em CI:**
- `input_unverified/p78_s3_doc1_certificado.json`: `093dd70b4fc6ddf0385094cd112699034f2fad2c36a3c51758c7f7c044652137`.
- `input_unverified/p78_s3_doc2_certificado.json`: `ac356d2997e93945e8ce0751dd9e2f00054b45b9284a800d3ed690b57d8c6846`.
- `input_unverified/p78_s3_downstream_scan.json`: preservado para comparação; **apenas kit temporário**, não repositório.

O script de adoção e seu JSON, assim como o log recebidos, foram examinados, mas seus estados `REPO_CERTIFIED` não foram adotados automaticamente. O relatório original DOC1 contém `monotonia_deriv_lower_min` igual a uma margem de B1 na linha 10, apesar de as derivadas por precisão terem ordem de grandeza ≈16,5: portanto, a variável de resumo não é utilizável como cota de derivada. O DOC2 usa constantes ajustadas e envelopes `1e-11` e `1e-12`, sem demonstrar suas propriedades de inclusão em Arb.

## O que realmente é certificado na P84

**Certificação computacional condicional**:
- K1 refeito via quadratura `acb.integral` e cauda analítica; K2 por fórmula analítica; C3/C4/K3/B1 pela fórmula original do P79.
- Extremos diádicos exatos para cada quantidade relevante, margens racionais de `B1_S2-B1_upper>0`; nas linhas 10–16 exige margem `≥1/1000`.
- Estrutura monotônica `K3(b4)=αb4+β`, `α,β>0`, `K3>1`, `e=(2σ₂−1)/(2−σ₂)>0`, logo `B1(b4)=D(αb4+β)^e` e derivada positiva, no regime implementado.
- Validação por comparação entre três precisões, SHA-256 de entradas/código/saídas e CI fail-closed.
- Reutilização dos 96 testes racionais de transmissão `S_ε→C3/C4` da P81 e prova elementar de `b3=0.605` da P83.

**Não certificado por essa operação:**
- Limiar `b4*` exato do DOC1; seu número em float permanece diagnóstico.
- A alegada demonstração de monotonicidade de `G(ε)` e os 28 testes `arb.zeta` do DOC2; sua transcrição em float **não** fornece enclausuramento verificável.
- Teoremas externos Ramaré–Zúñiga sobre `S_ε(X)` e Preissmann/Montgomery–Vaughan sobre valor médio (a P83 já provou o fator `b3`).
- Cobertura `3×10^10 ≤ X < 3×10^12`. Testar os parâmetros nas fronteiras **não substitui** estimativas públicas para `S_ε(X)` no interior.
- Prova global da Hipótese de Riemann.

**Sobre as linhas 1–9:** a P78 contém `certificate.json` para todas as 16 linhas e a P79 já as testou em Arb. Assim, a P84 **pode** certificar o fator B1 das primeiras nove faixas sem o certificado S1 separado solicitado pelos anexos. Isso NÃO amplia a validade do Lema 6 a regiões não cobertas pelas fontes externas.

## Reproduzir

A partir da raiz do checkout GitHub:

```bash
python -m pip install python-flint==0.9.0 mpmath==1.3.0
python -m unittest -v NUMERIA_P84_S3_REPO.test_s3
python -m NUMERIA_P84_S3_REPO.audit_s3 --bits 192 --output p84-s3-interval-192.json
python -m NUMERIA_P84_S3_REPO.audit_s3 --bits 256 --output p84-s3-interval-256.json --replay-output p84-s2-robust-suzuki-512.json
python -m NUMERIA_P84_S3_REPO.audit_s3 --bits 384 --output p84-s3-interval-384.json
mkdir -p p84-grade
cp p84-s3-interval-*.json p84-grade/
# renomeie as cópias para s3-interval-{192,256,384}.json para o build_docs
python -m NUMERIA_P84_S3_REPO.build_docs --input-dir p84-grade --output-dir .
python -m NUMERIA_P84_S3_REPO.scan_downstream --repo-root . --out p84-real-downstream-scan.json
```

O workflow `.github/workflows/numeria-p84-s3-repo.yml` automatiza as tarefas e publica artefatos de saída. Nenhum resultado é apresentado como aprovado até que as execuções GitHub Actions pertinentes concluam com `success`.

**Situação científica:** a P84 certifica as novas contas intervalares, condicionalmente aos teoremas externos, sem aceitar conclusões de RH global.
