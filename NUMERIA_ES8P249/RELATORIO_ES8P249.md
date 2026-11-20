# NUMERIA–RIEMANN · ES8-P249

A soma oscilatória dos primos permite antecipar de **504 para 309** o início de uma cauda com símbolo residual positivo na janela `L=161/200`, com o mesmo `β` da P237. Arb certificou a parte finita da cauda em 24.960 subintervalos. A fórmula clássica de Binet e a cota triangular dão a extensão a frequências arbitrariamente altas. Um novo cálculo com esse corte teria **618 painéis GL60**, 390 a menos que na P237, uma redução de 38,7%. O certificado do bloco de baixa frequência para o corte 309 ainda precisa ser construído.

A P249 foi construída sobre o último commit P247 recuperado no histórico público. A revisão da P240 e da P241 usou o relato da P244. Os arquivos originais dessas duas etapas e os da P243 e da P248 seguem pendentes de localização. A positividade global da forma de Weil e a Hipótese de Riemann permanecem em aberto. [Base P247](https://github.com/jbarumby/scicalc/commit/6adabcca7cc09f60df740934114d49966893d3d7)

## Cancelamento diádico

Na janela fixa, `5 < exp(2L) < 6` e as potências primas relevantes são `q = 2, 3, 4, 5`, com bases `p = 2, 3, 2, 5`. Escrevemos

```text
a_q = 2 log(p) / sqrt(q)
ω_q = log(q)
S(t) = Σ_q a_q cos(ω_q t)
A = Σ_q a_q
M₂ = Σ_q a_q ω_q²
```

Em cada célula de frequência de largura `2^-j` e centro `c`, seja `S_j(t) = S(c) + (t-c) S'(c)`. Taylor dá `|S(t)-S_j(t)| ≤ M₂ 2^(-2j)/8`. A decomposição até o nível `J` é

```text
S = S₀ + Σ_{j=1}^J (S_j − S_{j−1}) + (S − S_J).
```

Os coeficientes `S(c)` e `S'(c)` são somas com sinal, de modo que os cancelamentos precedem a aplicação de valores absolutos. Para cada célula de raio `r`, o certificado usa `|S(c)| + r|S'(c)| + r² M₂/2` como majorante de `|S|`. A tabela apresenta cotas superiores para `sup|S|/A`, arredondadas para cima em relação aos extremos diádicos do JSON.

| Faixa de frequência | Cota superior |
|---|---:|
| [0, 1] | 1 |
| [1, 2] | 0,604166 |
| [2, 4] | 0,679301 |
| [4, 8] | 0,261330 |
| [8, 16] | 0,642019 |
| [16, 32] | 0,783348 |
| [32, 64] | 0,868085 |
| [64, 128] | 0,931277 |
| [128, 256] | 0,955843 |
| [256, 504] | 0,989148 |

A igualdade `S(0)=A` impede reduzir a cota em toda a faixa desde zero. Também vale `limsup_{t→∞} S(t)=A` para qualquer conjunto fixo e finito desses termos. A aproximação simultânea de Dirichlet aplicada a `ω_q/(2π)` produz inteiros positivos com todas as fases arbitrariamente próximas de múltiplos de `2π`. Se os inteiros assim obtidos permanecerem limitados, um deles é um período comum exato e seus múltiplos produzem retornos arbitrariamente tardios. Portanto a melhora útil para a estratégia P246 combina cancelamento em faixas finitas e crescimento arquimediano. [Diagnóstico cofinal P246](../numeria/README_ES8P246_COFINAL_SCALING.md)

## Sinal da cauda com o β da P237

Fixe `B(t)=log(t/(2π))−1/t` e `β=B(504)−A`. Em cada célula, a monotonicidade de `B` e a cota de Taylor para `S` dão

```text
B(t) − β − S(t) ≥ B(c−r) − β − S(c) − r|S'(c)| − r² M₂/2.
```

A busca verifica células consecutivas a partir de 504 e informa o primeiro corte inteiro com cobertura pelo refinamento escolhido. A largura `1/32` e a largura `1/64` certificaram `[481,504]`. Com largura `1/128`, as 24.960 células de `[309,504]` deram extremo inferior comum `250055138131793279450921429863262946377 × 2^-140 > 1,794 × 10^-4`. A primeira célula inconclusiva na busca foi `[39521/128,39522/128]`. O certificado dá um corte suficiente com esse procedimento. A partir de 504, `S≤A` e `B(t)≥B(504)`.

A ponte analítica usa a fórmula clássica de Binet para `z=1/4+it/2`. Para `t>0`, temos `|x²+z²|≥t/4`, `Re(1/(2z))≤1/(2t)` e `log|z|≥log(t/2)`. Como `∫₀∞ x/(exp(2πx)−1) dx = 1/24`, o resto integral de Binet tem módulo no máximo `1/(3t)`. Segue que `Re ψ(z)−log π ≥ log(t/(2π))−5/(6t) ≥ B(t)`. Assim, na normalização usada pela P237/P239, o resultado de sinal aplica-se a `Ψ_L(t)−β` para `t≥309`. A identidade de Binet e a avaliação clássica da integral são insumos analíticos. [Revisão da redução de Zhu](../numeria/README_ES8P239_ZHU_EXTENSION.md) · [NIST DLMF §5.9](https://dlmf.nist.gov/5.9)

Com `β` fixo, o bloco contínuo de baixa frequência com corte 309 é menor ou igual, na ordem de formas quadráticas, ao bloco com corte 504. A diferença é a contribuição positiva de `[309,504]`. Para usar o corte menor, a próxima etapa precisa remontar o bloco dos 618 painéis, recalcular o orçamento GL60, incorporar `β` e os termos de polo, testar a positividade do bloco e reaplicar as cotas de cauda modal e acoplamento de P239. O resultado atual é para `L=161/200`.

## Operadores matriz-livre

Sejam `N` células uniformes `I_i` de largura `h=2L/N`, com funções ortonormais `e_i=1_{I_i}/√h` estendidas por zero. A P249 implementa:

- **Termo primo em frequência completa.** O multiplicador `S(t)` corresponde à soma `Σ_q (a_q/2)(T_{ω_q}+T_{−ω_q})` de translações simétricas. O lag `d` é `Σ_q (a_q/2)[H(dh−ω_q)+H(dh+ω_q)]`, com `H(u)=max(1−|u|/h,0)`. A fórmula é a interseção geométrica exata das células na aritmética real. Cada translação toca até dois lags por direção.
- **Termo primo de quadratura finita.** Para frequências `t_l` e pesos `w_l` que já incluem o multiplicador e a normalização, o lag é `h Σ_l w_l sinc²(t_l h/2) cos(t_l h d)`. Em um painel P236, `w_l=w_GL,l S(t_l)/(4π)`; a contribuição entra em `C` com sinal negativo. A escala grossa, os detalhes e o resíduo atuam por linearidade da quadratura. Para um painel inteiro, omitir o resíduo do nível `J` custa no máximo `L M₂ 2^(-2J)/(8π)` em norma de operador, pois a norma da contribuição unitária de cada nó é limitada por `2L` e a soma dos pesos GL positivos após normalização é `1/(2π)`.

Ambos são Toeplitz na base de células. Uma imersão circulante com preenchimento por zeros e comprimento de pelo menos `2N−1` permite aplicar o operador em `O(N log N)` operações e `O(N)` memória, sem formar a matriz. A montagem dos lags de quadratura é `O(NM)` para `M` nós e a FFT do núcleo é reutilizada. A soma direta de poucos deslocamentos custa `O(N·#q)` e pode ser mais barata para as quatro potências desta janela. A FFT intervalar usa bolas complexas ACB e raízes da unidade calculadas com arredondamento externo. Os testes exercitam preenchimento, dimensões pequenas, operadores complexos e oráculos diretos independentes.

A matriz P236/P237 de coeficientes de Legendre mantém a sua própria estrutura. Para compará-la com a matriz de células, a P249 calcula os coeficientes da projeção ortogonal pela primitiva exata de cada polinômio de Legendre. Em uma célula, a identidade de covariância em integral dupla dá um erro de Fourier de no máximo `h³ |t| sup|φ'_n|/12`. Somar sobre as células dá `2L h² |t| sup|φ'_n|/12`. A cota de Markov dá `sup|φ'_n|≤sqrt((2n+1)/(2L)) n²/L`. Cauchy–Schwarz limita a amplitude antes e depois da projeção por `sqrt(2L)`, de modo que o erro de um produto de duas amplitudes é limitado pela soma dos erros individuais vezes `sqrt(2L)`. A ponderação pelos módulos dos pesos primos dá a cota de discretização por entrada calculada no replay.

## Erros e verificações executadas

| Verificação | Resultado |
|---|---|
| FFT binária de translações, 65 células e vetor diádico fixo | Erro máximo por coordenada `<8,749 × 10^-15`, incluindo a montagem dos lags, contra produto direto Arb |
| FFT binária de quadratura, 128 células, painéis 0, 503, 1007 e vetor diádico fixo | Erro máximo `<6,34 × 10^-17` por coordenada, incluindo nós, pesos, lags e FFT, contra a ação ACB |
| Decomposição de quadratura até o nível 5 | Diferença observada entre somar os componentes e aplicar o operador inteiro `<1,39 × 10^-16`. A cota de norma ao descartar o resíduo é `<2,210 × 10^-4` por painel |
| Projeção do termo primo, 36 entradas de baixo grau | Maior limite Arb para o erro absoluto de grade `<3,109 × 10^-5`, entre três painéis e duas paridades com 128 células. A cota a priori pode ser mais larga |
| Replay da quadratura completa `C` dos pilotos P236 | `3.204/3.204` intervalos arquivados contidos nos replays em Arb a 900 bits. São 534 entradas de cada um de seis arquivos P245, incluindo graus altos. A maior largura do replay é `<2^-320` |

Os erros de FFT binária se referem aos vetores executados. A FFT ACB fornece um intervalo para a ação sobre cada entrada fornecida. O replay da quadratura `C` inclui a digama e a subtração de `β` no integrando, como nos pilotos P236. Os pilotos representam somas GL60 de painel. O controle de erro da integral contínua e o certificado LDL pertencem às etapas históricas abaixo. A verificação do erro de grade da P249 cobre as entradas de baixo grau explicitadas no JSON. [Dados e extremos diádicos](P249_TRIAL.json)

## Relação com os certificados anteriores

| Etapa | Evidência revisada e relação com P249 |
|---|---|
| P237 | 1.008 painéis por setor, 377 modos por paridade. O LDL direto foi inconclusivo nos índices 60 e 66. LDL do centro e transferência de Weyl certificaram as duas caixas finitas. A P249 preserva `L` e `β` e estabelece um corte suficiente para o sinal da cauda. [Fonte](../numeria/README_ES8P237.md) |
| P238 | Margem finita comum de `10^-18` e condição de Schur para o complemento da base. A margem é do bloco avaliado com o corte antigo. [Fonte](../numeria/README_ES8P238_P207_BRIDGE.md) |
| P239 | Margem de `5 × 10^-19` na janela fixa, no domínio da redução de Zhu, com controle de perda diagonal e acoplamento modal. A P249 usa a mesma convenção de primos, suporte e paridades. [Fonte](../numeria/README_ES8P239_ZHU_EXTENSION.md) |
| P240–P241 | A revisão P244 registra o controle do sinal arquimediano com Binet/Zhu para `t≥504`. A P249 verifica uma faixa adicional com a soma assinada. [Revisão](../numeria/p244/package/RELATORIO_ES8P244.md) |
| P242 | Replay independente C/MPFR a 1.024 bits, com 754 pivôs de centro e transferência de Weyl sobre as caixas P237. O teste de pivôs é distinto da verificação da quadratura isolada. [Fonte](../numeria/README_ES8P242_INTEGRATION.md) |
| P243 | Arquivos e resultados pendentes de localização. |
| P244 | Autenticação dos replays MPFR e formulação de um critério de extensão cofinal com corte suficiente `max(504,4π exp(A_L))`. O resultado P249 tem uma janela fixa. [Fonte](../numeria/README_ES8P244_GLOBAL_EXTENSION.md) |
| P245B/C/D | Replays selecionados em ponto fixo e orçamento analítico GL60. A P249 amplia o replay Arb dos pilotos P236 para todas as amostras arquivadas nos seis arquivos usados, com o mesmo integrando e uma implementação escalar independente. [Gate B](../numeria/README_ES8P245_GATEB_AUDIT.md) · [Gate C](../numeria/README_ES8P245_GATEC_AUDIT.md) · [Gate D](../numeria/README_ES8P245_GATE_D_ANALYTIC_FOUNDATIONS.md) |
| P246 | A auditoria reproduziu 78 somas selecionadas de blocos e baixou a cota do símbolo complexo de 50 para `274/7`. O diagnóstico cofinal dá `A_L∼4e^L` e explica o custo do corte triangular. As cotas de P249 são para frequência real. [Replay](../numeria/README_ES8P246_REPLAY_AUDIT.md) · [Escala cofinal](../numeria/README_ES8P246_COFINAL_SCALING.md) |
| P247 | 26 inclusões selecionadas de dois setores em `[25,50)`, em ponto fixo inteiro a 512 bits. Mantém o majorante analítico de Bessel e os brackets GL de 224 bits. [Fonte](../numeria/p247/README_P247_512_BIT.md) |
| P248 | Arquivos e resultados pendentes de localização. |

Os manifestos checados para P237, P242, P244, P245B e P246 conferiram no commit base. O README P238 e o código P245C conferiram nos commits de origem de seus manifestos e têm revisões posteriores no P247. [Origem do README P238](https://github.com/jbarumby/scicalc/commit/0ed9b58864411cb4e909a9e452be9e43b875aeab) · [Origem do pacote P245C](https://github.com/jbarumby/scicalc/commit/433015ce37da8941086e66aab093791bfc4c311e)

O replay autenticado P242/P244 foi executado novamente em C/MPFR a 1.024 bits. Os 754 pivôs e ambos os traços completos de pivôs e inversas coincidiram com os arquivos de referência. Passaram os 13 testes P249 tanto na execução normal quanto em `python -O`, além de 81 testes históricos selecionados de P151 e P237–P247.

## Reprodução

Na raiz do repositório, com Python 3.11 ou superior:

```bash
python -m pip install -r NUMERIA_ES8P249/requirements.txt
python -m unittest NUMERIA_ES8P249.test_p249 -v
python -O -m unittest NUMERIA_ES8P249.test_p249 -v
python -m NUMERIA_ES8P249.run_trial --out /tmp/P249_TRIAL.json
```

Cada par `[m,e]` do JSON representa o diádico exato `m × 2^e`. A parte intervalar registra extremos dirigidos e identifica os arquivos de referência por SHA-256. Os diagnósticos binários são da execução registrada e podem variar nos últimos bits com a plataforma. As cotas matemáticas e os testes independentes fornecem verificações por execução.
