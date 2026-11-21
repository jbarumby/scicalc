# NUMERIA–RIEMANN — passagem ao Codex (P71 → P72)

> **Atualização 2026-11-21:** o gate Arb 20 zeros foi executado com sucesso no GitHub Actions [run 44939756190](https://github.com/jbarumby/scicalc/actions/runs/44939756190), e os dois certificados JSON foram publicados. Ver [CI_ARB_RESULT.md](CI_ARB_RESULT.md). Os passos de geração abaixo permanecem como *replay independente* e revisão analítica, não como ausência de certificação executada.

**Branch:** `numeria-p72-zero-side-codex`  
**Base:** PR #71 (`numeria-suzuki-prefix-balance`), que cobre a faixa aritmética até `t=log(10^8)`.  
**Estado P72:** prévias com 3, 10 e 20 zeros. A extensão até `t=44.35` é **condicional** aos intervalos dos 20 zeros: `preview_zeros.json` não é certificado por Arb. Não declarar esse resultado como teorema incondicional antes do gate.

## Sequência prioritária

1. Reproduzir `NUMERIA_P71_ZERO_SIDE/numeria_pr71_zero_side_audit.py` (Python e `-O`) e conferir SHA-256 do resultado da malha. Revisar a dependência das três alturas documentadas.
2. Reproduzir os testes da P72: `cd NUMERIA_P72_ZERO_SIDE && python -m unittest test_p72 -v && python -O -m unittest test_p72 -v`. Eles validam o mecanismo racional, **não** a proveniência dos 20 zeros.
3. Instalar `python-flint==0.9.0 mpmath==1.3.0` e executar `python NUMERIA_P72_ZERO_SIDE/export_zeros_flint.py --n 20 --precision 192 --out NUMERIA_P72_ZERO_SIDE/arb_zeros_20.json`. Conferir no código-fonte e na documentação a semântica de `acb.zeta_zeros`, numeração, multiplicidades, isolamento e conversão dos extremos dyádicos.
4. Executar `python NUMERIA_P72_ZERO_SIDE/verify_spectral_rational.py --zeros NUMERIA_P72_ZERO_SIDE/arb_zeros_20.json --n 20 --to 44.35 --out NUMERIA_P72_ZERO_SIDE/arb_20_t44.35.json`. Repetir com `-O`, conferir SHA-256 da malha e exigir `uniform_lower_rational > 0`.
5. **Endurecer a cadeia de confiança**: o campo `zero_certification: arb_enclosures` é metadado declarativo, não comprovação criptográfica. Evitar que o verificador marque `certified_positive: true` apenas porque uma entrada JSON contém tal string. Exigir origem autenticada/manifesto reprodutível ou uma etapa combinada que gere e valide os intervalos antes de reconhecer certificação.
6. Conferir em detalhes os lemas analíticos importados: identidade incondicional de Suzuki (2023, eq. 1.3), Platt–Trudgian até `H=3e12`, erro explícito `N(T)` de Bellotti–Wong, termo de fronteira na integração de Stieltjes e interpolação de curvatura `N h²/4`.
7. Somente após os gates anteriores, avançar para 100–3.000 zeros e investigar `t=47.48`. A prévia com 20 zeros para esse extremo é inconclusiva; não corresponde a um contraexemplo. Manter separada a cauda infinita (RH global continua aberta).

## Testes mínimos

```bash
python NUMERIA_P71_ZERO_SIDE/numeria_pr71_zero_side_audit.py
python -O NUMERIA_P71_ZERO_SIDE/numeria_pr71_zero_side_audit.py
(cd NUMERIA_P71_ZERO_SIDE && sha256sum -c NUMERIA_P71_ZERO_SIDE_SHA256SUMS.txt)
(cd NUMERIA_P72_ZERO_SIDE && sha256sum -c SHA256SUMS_P72.txt)
(cd NUMERIA_P72_ZERO_SIDE && python -m unittest test_p72 -v)
(cd NUMERIA_P72_ZERO_SIDE && python -O -m unittest test_p72 -v)
```

`NUMERIA_P72_ZERO_SIDE/SOURCES.tar.zlib.b64` é uma cópia de preservação de todos os 17 arquivos originais; os arquivos individuais estão publicados no mesmo diretório. As prévias em JSON **não são certificados de zeros**. O workflow `.github/workflows/numeria-p72-zero-side.yml` deve produzir os primeiros artefatos realmente baseados em Arb, se passar.

**Não mesclar por interpretação matemática automática.** Solicitar revisão das cotas, intervalos de zeros e integridade antes de anunciar um novo alcance certificado. 
