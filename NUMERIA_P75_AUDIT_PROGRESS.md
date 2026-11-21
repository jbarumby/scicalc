# Auditoria independente da P75 — progresso até 21/11/2026

## Resultado da comparação até t=79

O replay independente da P75 produziu o mesmo arquivo JSON da CI do PR #75, byte por byte: SHA-256 `292cd84348575ef0f1328d58d82a0334f82e526ab9f7c38e52d05ff6616ea1cb`. Com os coeficientes de Chourasiya–Simonič v2, a execução certifica positividade finita em `18 ≤ t ≤ 79`, com 500 zeros isolados, 256 bits, 4.284 células aceitas, nenhuma célula inconclusiva e margem inferior mínima positiva. Isso continua sendo um resultado condicional às cotas da cauda e à fórmula importada, não uma prova global da Hipótese de Riemann.

Como controle, substituí a cota de densidade do preprint pelo perfil de Kadiri–Lumley–Ng (KLN), mantendo a mesma fórmula de Suzuki, os mesmos 500 zeros, a mesma precisão, a malha e a curvatura local da P75. A execução baseada somente nos resultados publicados usados pelo projeto — KLN para densidade, Bellotti–Wong para contagem total, Platt–Trudgian para a verificação finita da RH e Suzuki para a representação — certifica até `t = 4.344.693/65.536 ≈ 66,294754`. A primeira célula não resolvida é `[4.344.693/65.536, 17.378.773/262.144]`. “Não resolvida” significa que este majorante não provou positividade nessa célula; não significa que a função seja negativa ali.

A convenção de contagem da KLN merece cuidado: o artigo define `N(σ,T)` contando zeros com `β>σ`, ao passo que a cauda usa `β≥σ`. Para evitar depender de ausência de zeros exatamente sobre um limiar racional, o baseline conservador usa, para cada limiar-alvo, a linha publicada imediatamente abaixo dele. Pela monotonicidade, os zeros com `β≥σ` ficam então incluídos na contagem estrita `β>σ₀`. A nova saída JSON, incluindo a cobertura ordenada e suas margens, acompanha esta nota.

## Hipóteses do preprint Chourasiya–Simonič

A página do arXiv registra a v2 em 30/09/2025 e a classifica no repositório como manuscrito de 31 páginas, com o comentário “Comments are most welcome”; não apresenta referência editorial de periódico. Até esta consulta, não confirmei publicação em periódico, portanto trato o resultado como preprint, sem inferir que ele esteja incorreto. [1]

A aplicação do Teorema 1 no Corolário 1 usa `T₀=H₀=H=H_RH=3×10¹²`, `h=1` e `d∈{0,28; 0,29}`. Os demais parâmetros tabulados são positivos (`κ∈[0,07;0,34]` e `δ≤0,303`). Para esses valores, as condições gerais do Teorema 1 são satisfeitas: `T≥H≥H₀≥1002`; `T≥T₀≥max(e^{2d},e^e)`; `h≥3×10¹⁰/T₀`; e `0<δ≤(1/20)log(hT₀)`. A hipótese finita de RH até `3×10¹²` é atribuída a Platt–Trudgian. [1] [2]

Há uma sutileza no primeiro intervalo da Tabela 1. O Teorema 1 exige `σ₁>1/2+d/log T`, enquanto o Corolário estende a estimativa à parte muito próxima da reta crítica por continuidade/monotonicidade e cobre o pequeno intervalo restante com a cota de contagem total. Assim, a linha da tabela com `σ₁=1/2` não deve ser lida como aplicação direta do Teorema 1 nessa fronteira; a justificativa é a prova do Corolário 1 no próprio artigo. [1]

O preprint define `N(σ,T)` com a convenção `β≥σ`, apropriada para a camada cumulativa usada no código. As linhas da Tabela 1 cobrem intervalos fechados, e a implementação escolhe o menor entre as duas cotas disponíveis nas fronteiras comuns. A cota de cauda é obtida por integração de Stieltjes/partes: se `Aσ(u)` conta os zeros com `H<γ≤u` e `β≥σ`, então

`Dσ = Σ_{γ>H, β≥σ} γ⁻² = 2∫_H^∞ Aσ(u)u⁻³du ≤ 2∫_H^∞ N(σ,u)u⁻³du`.

Substituindo `N(σ,u)≤B₁u^α(log u)^p+B₂(log u)²+B₃log u`, com `α=3(1−σ)/(2−σ)`, `p=(7−5σ)/(2−σ)` e `λ=2−α>0`, a integral fica exatamente

`Dσ ≤ 2B₁ Γ(p+1,λ log H)/λ^(p+1) + [B₂((log H)²+log H+1/2)+B₃(log H+1/2)]/H²`.

A forma de gama incompleta, os dois termos logarítmicos e a monotonicidade usada na soma por camadas estão coerentes com a identidade analítica; o código avalia as quantidades com arredondamento para fora em Arb. O teste independente de identidades e integrais em mpmath é apenas diagnóstico, enquanto a saída certificadora usa Arb. [1]

A passagem da massa total para a parte fora da reta usa as simetrias conjugada e funcional `β↔1−β`: para `σ>1/2`, no máximo metade da massa positiva total `W(H)=Σ_{γ>H}γ⁻²` fica no semiplano direito. A camada em `θ=β−1/2` é majorada por somas de incrementos positivos `cosh(θ_{j+1}t)−cosh(θ_jt)`, o que justifica a monotonicidade e a segurança do uso das linhas fechadas `β≥σ_j`. A alternativa grosseira `2W cosh(t/2)` permanece como limitante superior. Os fatores são os mesmos implementados em `cs_v2.py` e na cauda legada P74; a fonte para a representação de zeros é o artigo publicado de Suzuki. [3]

Os 500 zeros retidos são intervalos Arb disjuntos, com sinais opostos de Hardy Z nas pontas e uma contagem `zeta_nzeros` numa altura separadora. O artigo de Platt–Trudgian cobre o intervalo até `3×10¹²`; esses passos verificam a extração finita utilizada pela execução e não substituem a hipótese de densidade acima dessa altura. A fórmula de Suzuki usada também é publicada e sua positividade nesta faixa continua sendo apenas o objeto finito testado. [2] [3]

## Pendência principal e decisão

A verificação dos parâmetros do Teorema 1 não basta para reproduzir independentemente todos os coeficientes `B₁,B₂,B₃` da Tabela 1. Os valores publicados são decimais arredondados e, no texto consultado, não encontrei um certificado de arredondamento externo ou uma tabela de intervalos dirigidos que demonstre, por si só, que cada decimal impresso é uma majoração para todas as linhas. O preprint afirma as desigualdades e o replay testa o efeito computacional da tabela; isso não equivale a recalcular as constantes a partir das fórmulas dos `A_i` do Teorema 1. A próxima prioridade é reproduzir esses cálculos com Arb a partir dos parâmetros tabulados, ou inflar os coeficientes por uma margem explicitamente provada e refazer o replay.

O baseline conservador baseado em periódicos não cobre o trecho de `66,294754` a `79`, ao passo que o perfil do preprint cobre até `79`. Portanto, a contribuição específica do preprint é substancial para esta execução. A decisão recomendada é **não tentar ainda `t>79`**: primeiro reproduzir independentemente as constantes tabuladas, estabelecer limites dirigidos para seu arredondamento e preservar uma variante jornal-only reproduzível. O PR #75 permanece aberto e em rascunho; esta nota está em branch independente baseada na P74, sem alteração ou merge dos PRs #74/#75.

## Referências

[1]: https://arxiv.org/abs/2507.15184v2 "An explicit form of Ingham's zero density estimate, versão 2"
[2]: https://doi.org/10.1112/blms.12460 "The Riemann hypothesis is true up to 3·10^12"
[3]: https://doi.org/10.1112/jlms.12785 "Aspects of the screw function corresponding to the Riemann zeta-function"
[4]: https://doi.org/10.1016/j.jmaa.2018.04.071 "Explicit zero density for the Riemann zeta function"
[5]: https://doi.org/10.1090/mcom/4133 "Improved estimates for the argument and zero-counting function of the Riemann zeta-function"


### Verificação complementar de sensibilidade

Como teste de robustez computacional — não como prova de arredondamento — aumentei exatamente todos os coeficientes `B₁,B₂,B₃` do preprint por um fator `1,001`. A cobertura ainda chega a `t=79`, com as mesmas 4.284 células aceitas e nenhuma célula inconclusiva. Isso sugere que o resultado não depende de uma variação relativa minúscula na tabela, mas não substitui a derivação intervalar das constantes a partir do Teorema 1.


### Artefatos reproduzíveis

- [Baseline jornal-only com limiar KLN estrito corrigido](P75_JOURNAL_ONLY_BASELINE.json)
- [Sensibilidade com inflação uniforme exata de 0,1%](P75_CS_V2_PLUS_0P1PCT_SENSITIVITY.json)


O baseline tem um script autocontido, executável no branch de auditoria baseado na P74: [P75_JOURNAL_ONLY_BASELINE.py](P75_JOURNAL_ONLY_BASELINE.py). Reexecutei-o com `python3 P75_JOURNAL_ONLY_BASELINE.py --out /tmp/p75-journal.json`; os campos de status, célula final, contagens, margem mínima, primeira célula inconclusiva e SHA-256 da cobertura coincidiram exatamente com `P75_JOURNAL_ONLY_BASELINE.json`.
