# P72 — Resultado do gate Arb no GitHub Actions

**Data (UTC):** 2026-11-21. **Run:** https://github.com/jbarumby/scicalc/actions/runs/44939756190  
**Resultado:** workflow concluído com sucesso. `python-flint==0.9.0`, `mpmath==1.3.0`, Arb `192` bits.

## Evidências persistidas

- `arb_zeros_20.json`: 20 ordenadas dos zeros da zeta isoladas por `acb.zeta_zeros(1,20)`, com extremos diádicos.
- `arb_20_t44.35.json`: reprodução racional (malha `h=0.005`) em `18<=t<=44.35`, `certified_positive: true`.
- Margem inferior publicada: `0.000176565011706`.
- SHA-256 do JSON de zeros: `c5993b75a6db7803251655955859cd935f4381b2c0bc20fd77b9f0258bdeddab`.
- SHA-256 do JSON da malha: `5fb6c2b838f7a6a7af88905869e8c8cd308f74875ddcd7ec70e122a066179fd3`.
- SHA-256 interno das linhas da malha: `4eda34dd6c306e15ace01ce62ac5c2a4862b1c7636ba22613572ab6216d88ae2`.

A execução local do verificador racional com o arquivo de zeros proveniente do CI produziu **JSON idêntico** ao do GitHub Actions. A pasta contém os dois arquivos no versionamento para revisão pelo Codex.

## Escopo e limites

É uma **certificação computacional finita baseada nos teoremas externos invocados**, não uma prova autossuficiente do domínio inteiro nem de RH. Precisa revisão independente da aplicação da fórmula de Suzuki, simetrias, Bellotti–Wong, Platt–Trudgian e da cadeia de confiança para `zero_certification`. A verificação de `t=47.48` ainda é inconclusiva com 20 zeros.

O workflow foi reforçado para exigir que novas execuções reproduzam exatamente os dois JSON versionados, além dos manifestos SHA-256 das fontes P71/P72. Se essa regressão falhar, o gate deve ser considerado aberto novamente.
