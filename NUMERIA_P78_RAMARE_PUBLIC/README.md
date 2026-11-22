# NUMERIA–RIEMANN P78 — Limite público, **restrito**, para o Lema 6

**Base:** PR #77, `numeria-p77-cs-v2-b1-quadrature`; **não** altera nem mescla P71–P77.
**Resultado:** um majorante alternativo para a soma dupla de Möbius **somente na faixa de parâmetros da Tabela 1 de Chourasiya–Simonič v2**, importando o preprint de Ramaré–Zúñiga Alterman **v1 (26 mar. 2026)**. Não afirma a demonstração completa do Lema 6 nem RH.

## Descoberta documental

O Lema 6 de [Chourasiya–Simonič v2](https://arxiv.org/html/2507.15184v2), entre as fórmulas de C3 e C4, utiliza `b4=0.54` e atribui uma desigualdade a comunicação privada de O. Ramaré. Posteriormente, Ramaré e Zúñiga Alterman publicaram [*On a Möbius double sum*, arXiv:2603.25961v1](https://arxiv.org/abs/2603.25961v1), examinando precisamente

\[S_\eta(X)=\sum_{d,e\le X}\frac{\mu(d)\mu(e)}{\operatorname{lcm}(d,e)^{1+\eta}}.\]

O artigo fornece explicitamente (Teoremas 1.1–1.2 e cálculos da seção 6) as seguintes cotas usadas aqui:

| Intervalo de X | 0 ≤ η ≤ | Cota superior de Sη |
|---|---:|---:|
| [10¹², 10¹⁶) | 0,01 | 0,526 (seção 6) |
| [10¹², 10¹⁶) | 0,02 | 0,5503 |
| [10¹², 10¹⁶) | 0,04 | 0,593 |
| [10¹⁶, 10³³) | 0,02 | 0,5104 |
| [10³³, ∞) | 0,01 | 0,4794 |

**Atenção à versão:** esses valores são de **v1**, fixada para tornar a auditoria reprodutível; o HTML experimental de v3 contém marcadores numéricos não resolvidos. Não substituir silenciosamente por v3. O novo trabalho dos autores também é um *preprint*, não uma fonte necessariamente revisada por pares.

## Lemma substituto com escopo limitado

Para `X≥3×10¹²`, `0<δ≤0,303` e `0<η≤2δ/log X`, definindo a função `λ_X(n)` como em Chourasiya–Simonič, temos a identidade

\[\sum_{n\ge1}\frac{\lambda_X(n)^2}{n^{1+\eta}}=\zeta(1+\eta)S_\eta(X).\]

A comparação integral para a zeta dá `ζ(1+η) ≤ 1+1/η`. Como a constante de Euler é positiva, `exp(γ η) ≥ 1`. Portanto, uma cota suficientemente grande `b4` segue de

\[\eta e^{-\gamma\eta}\zeta(1+\eta)S_\eta(X)\le(1+\eta)S_\eta(X).\]

Sem usar logaritmos numéricos, mostramos `e<11/4` por Taylor com resto geométrico e, com potências de racionais exatos,
`(11/4)^28 < 3×10¹²`, `(11/4)^36 < 10¹⁶`, `(11/4)^75 < 10³³`.
Logo, respectivamente, `log X>28`, `>36`, `>75` nas três faixas. Assim, em cada linha `η ≤ 2δ/28`, `2δ/36`, `2δ/75`, permitindo aplicar as cinco cotas publicadas acima.

O script `restricted_lemma6.py` testa as 16 escolhas de δ da tabela com **frações inteiras exatas**, em três faixas de X e sem ponto flutuante. Conclusão:

- **Linhas 1–9** (δ≤0,132): o valor original `b4=0,54` continua majorando a desigualdade **na faixa de parâmetros restrita**;
- **Linhas 10–16**: obtemos constantes substitutas dirigidas para cima, de `0,557` a `0,606`, sem reutilizar indevidamente `0,54`.

**Este resultado não demonstra `b4=0,54` para o enunciado completo do Lema 6** (`X≥3×10¹⁰` e `τ∈(1,1.1]`). O limite usa X≥3×10¹² e η no intervalo específico da P75/P77, inclusive as ocorrências `η=δ/log X` e `η=2δ/log X`.

**Proibição importante:** não aplicar Hölder diretamente a `S_η(X)`! A identidade fornece log-convexidade à série `ζ(1+η)S_η(X)`, e não a `S_η` isolada. Uma interpolação anterior desse tipo foi descartada.

## Limites e próximo gate

Este certificado **não** reconstrói os 16 `B1` com os novos valores `b4`, não fecha novamente `18≤t≤79` e não declara RH global. A fonte nova substitui uma comunicação privada por uma desigualdade **publicamente escrita, porém mais fraca** e ainda dependente do teorema importado.

1. Refatorar o auditor `NUMERIA_P77_B1/audit_b1.py` para aceitar b4 por linha, sem alterar seu arquivo histórico.
2. Calcular os novos `C3,C4,K3,L1,L2,A1,A2,B1` com Arb 192/256 bits; se ultrapassarem B1 impresso, criar majorantes dirigidos para cima.
3. Atualizar um **perfil alternativo** da P75 com os novos coeficientes, sem substituir os dados originais; repetir toda a certificação de `t=18` a `t=79` e a junção com P71, sem esconder células inconclusivas.
4. Auditar as demais dependências analíticas, publicar certificados e considerar alternativas publicadas em periódicos.

## Reprodução

```bash
python -m NUMERIA_P78_RAMARE_PUBLIC.restricted_lemma6 --verify NUMERIA_P78_RAMARE_PUBLIC/certificate.json
python -m unittest NUMERIA_P78_RAMARE_PUBLIC.test_p78 -v
python -O -m unittest NUMERIA_P78_RAMARE_PUBLIC.test_p78 -v
sha256sum NUMERIA_P78_RAMARE_PUBLIC/*.py NUMERIA_P78_RAMARE_PUBLIC/certificate.json
```

**Ponto de prudência:** o certificado prova uma **implicação exata a partir das cotas publicadas**, não rederiva as cotas originais de Ramaré e Zúñiga Alterman. A avaliação matemática completa de sua demonstração permanece independente desta verificação.