> **Replay 2026-11-21:** o workflow [P75/256 faixas](https://github.com/jbarumby/scicalc/actions/runs/44947512003) passou e cobriu \(18\\le t\\le79\) com 500 zeros, 256 bits e 4.284 células, com limite inferior \(3,0358387546\\times10^{-6}\). Ver [CI_RESULT.md](CI_RESULT.md). A dependência de Chourasiya–Simonič v2 continua sendo um **preprint** e requer revisão externa.

# NUMERIA–RIEMANN · P75: Chourasiya–Simonič v2 em Arb + curvatura local

**Base de continuidade:** [PR #74](https://github.com/jbarumby/scicalc/pull/74) (`numeria-p74-zero-side-certified`, `13e67b6...`). **Esta etapa preserva integralmente a P74.**

## Escopo

Implementar em `python-flint==0.9.0` um majorante da cauda dos zeros potencialmente fora da reta crítica **sem utilizar os coeficientes de Kadiri–Lumley–Ng empregados na P74**. O novo perfil utiliza o **preprint Chourasiya–Simonič v2**; sua revisão documental e matemática deve ser considerada um gate independente, mesmo que o código e a CI passem.

- Fonte: [arXiv:2507.15184v2](https://arxiv.org/html/2507.15184v2), **Corolário 1 e Tabela 1**, 30 set. 2025.
- Verificação finita de RH até `H=3×10^12`: Platt–Trudgian (2021), https://doi.org/10.1112/blms.12460.
- Cota de contagem total `W(H)`: Bellotti–Wong (2025), https://doi.org/10.1090/mcom/4133.
- Fórmula de Suzuki: https://doi.org/10.1112/jlms.12785.

## Cota de densidade

Para `sigma ∈ (1/2,1)`, `N(sigma,H)=0` e, com `L=log H`,

`alpha=3(1−sigma)/(2−sigma)`, `p=(7−5sigma)/(2−sigma)`, `lambda=2−alpha`,

`D_sigma = Σ_{gamma>H, beta>=sigma} 1/gamma² ≤ 2B1*Gamma(p+1,lambda*L)/lambda^(p+1) + [B2*(L²+L+1/2)+B3*(L+1/2)]/H²`.

**Importante:** os decimais `B1, B2, B3` são os 16 valores **da versão v2**; em uma borda entre duas faixas ambas são válidas, e toma-se o menor limite superior. As operações são realizadas em Arb com arredondamento para fora. A avaliação via `mpmath` na auditoria anterior era *apenas diagnóstica*.

Para `1/2 = sigma_0 < sigma_1 < ... < sigma_K = 1`, utilize `D(1/2+) ≤ W(H)/2` e a monotonicidade de `D`: `U(t) ≤ 2W(H) + 4Σ D_upper(sigma_j+) [cosh((sigma_{j+1}-1/2)t)-cosh((sigma_j-1/2)t)]`. Essa forma é conservadora e evita supor RH acima de `H`.

## Curvatura local

Para `S_N(t)=2Σ_j (1−cos(gamma_j t))/gamma_j²`, temos `S_N''(t)=2Σcos(gamma_j t)` e `|S_N'''(t)|≤2Σgamma_j`. Para célula `[a,b]` de comprimento `h`, `m=(a+b)/2`:

`M_ab = min(2N, upper(S_N''(m))+h*Σupper(gamma_j))`.

A cota inferior em cada célula é a parábola do secante `secant − M_ab*h²*u(1−u)/2`. O algoritmo minimiza esta parábola a partir de intervalos dos valores nos extremos, usando apenas arredondamento para fora. Uma falha gera célula *inconclusiva*, não um sinal negativo de Suzuki.

## Executar

```sh
python -m pip install python-flint==0.9.0 mpmath==1.3.0
python -m unittest NUMERIA_P75_CS2.test_math NUMERIA_P75_CS2.test_flint -v
python -O -m unittest NUMERIA_P75_CS2.test_math NUMERIA_P75_CS2.test_flint -v
python -m NUMERIA_P75_CS2.certify --start 18 --end 79 --zeros 500 \
  --precision 256 --segments 128 --mesh 1/64 --out /tmp/p75-cs2-79.json
```

O comando acima pode **retornar inconclusivo**. Só se poderá relatar positividade até `79` após saída `status: conditional_finite_positive`, margem diádica estritamente positiva em todas as células, SHA-256 da cobertura, replay independente e revisão da fonte. Nenhuma conclusão global de RH decorre de uma janela finita.

## Trabalho futuro e preservação

1. Reproduzir o certificado P74 e autenticar 500 zeros sem depender de uma etiqueta JSON.
2. Verificar as constantes e condições **do preprint v2** diretamente contra a fonte original.
3. Executar a cobertura desde `t=18`, verificando a sobreposição com P71 (aritmética até `10^8`).
4. Comparar 16, 32, 64, 128 faixas, com e sem curvatura local, incluindo custo computacional e margem.
5. Separar perfil **somente literatura em periódicos** do perfil que usa o preprint; não mascarar um pelo outro.
6. Não mesclar os PRs automaticamente. Registrar limites positivos **e inconclusivos** sem inflar o resultado.