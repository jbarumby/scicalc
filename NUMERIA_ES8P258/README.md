# NUMERIA–RIEMANN · ES8-P258 · Auxiliar reforçado no limiar de 7

A P258 parte do commit [bb6ace41](https://github.com/jbarumby/scicalc/commit/bb6ace41f4e9739cde488d3952803b78e5a92118) do PR #63. O cálculo testa o modelo auxiliar truncado de quatro potências primas, com corte 504, piso espectral 1/4 e meia-largura racional 0,972955074527656652552676371721590. Essa largura fica cerca de 1,35×10⁻³⁴ acima de log(7)/2. Há 377 modos por setor de paridade.

## Resultado executado

**O piso 1/4 certifica o setor par, mas o setor ímpar ainda tem uma direção negativa.** Os dois cálculos de G concordaram nas 71.253 entradas triangulares de cada setor. O C/MPFR repetiu os 377 pivôs do setor par. No ímpar, a testemunha usa 29 coordenadas diádicas, preenchidas com zeros até a dimensão 377. O seu sinal foi conferido em Arb e por frações exatas, com as duas cotas GL60 incluídas.

| Modelo na largura testada | Setor par | Setor ímpar |
| --- | --- | --- |
| P253, piso antigo, corte 504 | Rayleigh certificado < −1,12×10⁻²⁹ | Rayleigh certificado < −5,30×10⁻²³ |
| P258, piso 1/4, corte 504 | Menor autovalor > 4,53×10⁻²⁸ | Rayleigh certificado < −1,7300×10⁻²⁴ |

As testemunhas usadas em etapas distintas podem ser diferentes. Os valores negativos da tabela são cotas de Rayleigh. O diagnóstico de condicionamento provou κ₂>6,69×10²⁷ no setor par e κ₂>8,69×10²³ no ímpar. A negatividade sobrevive à dilatação para a largura racional inferior a log(7)/2 registrada no JSON. A extensão por zero então impede a positividade do mesmo auxiliar, com corte 504 e piso 1/4, em toda janela maior ou igual ao limiar. Uma congruência invertível preserva a direção negativa. [Resumo durável e testemunha](P258_RESULTADO.json) · [Certificados GL60 e de extremos](https://github.com/jbarumby/scicalc/actions/runs/44858303395/artifacts/14688468811) · [Replay independente MPFR e frações](https://github.com/jbarumby/scicalc/actions/runs/44859003923/artifacts/14688504136)

### Obstrução para qualquer piso constante no corte 504

Em t=58933/100, a avaliação da digama e dos quatro cossenos em Arb provou Ψ₄(t)<7/25. Um piso constante válido para toda a cauda precisa ser menor que esse teto. O teste em 7/25 serve como **majorante de forma pelo ordenamento em β**, pois, para β₁≤β₂, R_{β₂}−R_{β₁}=(β₂−β₁)(I−G₅₀₄)≥0 no espaço da janela. O replay encontrou, para a forma majorante com piso 7/25, um limite superior de Rayleigh ímpar menor que −4,04×10⁻²⁵ já na largura racional abaixo do limiar. Logo, nenhum piso constante admissível em |t|≥504 torna positivo esse auxiliar de quatro termos acima do limiar. O teto 7/25 é usado para comparação e a certificação do piso 1/4 na cauda continua válida. [Dados, digests e largura inferior](P258_RESULTADO.json)

A janela de Weil certificada pelo NUMERIA permanece em 0,95. O obstáculo se refere ao auxiliar truncado, que descarta a parte do símbolo acima do piso. Uma próxima tentativa é reter pisos por faixas de frequência. Um candidato já tem validade analítica. Use 1/4 entre 504 e 850 e 1/2 a partir de 850. De fato, Binet e a soma de amplitudes dão Ψ₄(t)≥β+log(850/504)+1/504−1/850>1/2 para t≥850, pois log(850/504)>2(850−504)/(850+504)>1/2. Em relação à P258, isso acrescentaria (1/4)(I−G₈₅₀). A margem dessa nova matriz e a reserva para 7 ainda precisam ser certificadas.

## Certificados que entram no cálculo

A P249 separou a parte oscilatória de primos e controlou o complemento modal. A P250 certificou a cauda antiga a partir de 309. A P251 verificou o bloco finito com o β antigo e sua passagem para Weil na janela de meia-largura 0,805. A P252 ampliou essa janela para 0,95, com margem superior a 3,41×10⁻²⁷ na forma de Weil. A P253 certificou direções negativas no auxiliar de quatro termos com corte 504 e piso antigo na largura usada agora. Seus limites superiores de Rayleigh eram −1,12×10⁻²⁹ no setor par e −5,30×10⁻²³ no ímpar. Esses limites incluíam o erro analítico GL60 e não se transportam automaticamente ao piso elevado. [P251 Gate C](../NUMERIA_ES8P251/gate_c/README.md) · [P252](../NUMERIA_ES8P252/README.md) · [P253](../NUMERIA_ES8P253/README.md)

A P254 derivou a identidade de atualização e documentou um teste de Schur. Ao carregar a soma GL60 da P253, a P258 exige o SHA-256 fixado no código, assim como os dados de corte, largura e paridade. A verificação de sinal da P253 para Ψ₄(t)−1/4 em |t|≥504 é calculada de novo em Arb a 192 bits e comparada byte a byte ao certificado arquivado. O piso maior é aplicado só no complemento de frequência após 504. [Atualização e Schur P254](../NUMERIA_ES8P254/README.md) · [Cauda com piso maior](../NUMERIA_ES8P253/data/FOUR_PRIME_RAISED_FLOOR.json)

## Construção e erro da matriz

Na convenção angular, a transformada do polinômio de Legendre ortonormal de grau n tem amplitude √(2L(2n+1)) jₙ(Lt) e fase iⁿ. Um fator comum de fase cancela dentro de cada setor. Aos modos da paridade oposta corresponde integral nula no intervalo simétrico. O bloco G₅₀₄ integra a multiplicação das amplitudes reais no semieixo de 0 a 504, com fator 1/π. Pela identidade de Plancherel, G₅₀₄ é a compressão da projeção sobre a banda e satisfaz 0≤G₅₀₄≤I.

A mudança de piso acrescenta (1/4−β)(I−G₅₀₄) ao bloco antigo. A quadratura da parcela arquimediana e de primos ponderada da P253 e os polos originais entram uma vez. O termo I−G₅₀₄ mede a massa de Fourier fora da banda.

Cada setor é calculado por 1.008 painéis GL60 de meia unidade, 60.480 nós, com produtos UᵀU em Arb. A cota de erro da concentração usa diretamente a elipse de Bernstein de parâmetro 3. Nela, a normalização uₙ(t) é limitada por aₙ=ceil(√(2L(2n+1)))·2 min(1,Yⁿ/(2n+1)!!), com Y≥|Lt|. Na parte analítica, a exatidão até grau 119 e a série de Chebyshev dão erro por entrada menor que

    2T · 106/(111·3¹²⁰) · aₙ aₘ.

Aqui 106/111 é um majorante de 3/π, e 3 vem da cota 12 para quadratura na elipse combinada com o jacobiano 1/4. A soma máxima por linha limita a norma 2 da matriz de erro. Ao somar a atualização, o auditor multiplica essa cota pelo limite superior intervalar de 1/4−β e soma a cota analítica da parte ponderada. O erro de arredondamento continua incluído nas caixas Arb.

A segunda construção integra jₙjₘ exatamente por identidades de extremos. O Wronskiano x²(jₘjₙ′−jₙjₘ′) tem derivada [n(n+1)−m(m+1)]jₙjₘ. Para a diagonal, a derivada de x(jₙ²+jₙ₊₁²) é (2n+1)jₙ²−(2n+3)jₙ₊₁² e a semente da integral diagonal é Si(2x)−sin²(x)/x. Isso permite conferir todas as entradas GL60, acrescidas de seus erros, mediante avaliações Arb independentes nos extremos.

## Decisões de sinal e condicionamento

O teste LDL usa pivôs em Arb e a cota da inversa do centro com a perturbação total. Um pivô inconclusivo não determina o sinal. Para uma direção de tentativa, as coordenadas do vetor são diádicos exatos. A avaliação de Rayleigh inclui os erros analíticos e o mesmo sinal precisa aparecer no replay de inteiros e frações sobre os intervalos da matriz. O diagnóstico de condicionamento usa um vetor de resíduo pequeno para limitar σ_min por cima e uma entrada diagonal para limitar ‖A‖₂ por baixo.

Se o limite superior de Rayleigh for inferior à perda por dilatação, a negatividade também vale numa largura racional logo abaixo de log(7)/2. O modelo auxiliar conserva símbolos, polos e corte em todas essas larguras e a extensão por zero transporta essa direção a janelas maiores. Para a dilatação unitária de g em [−1,1], F_L(t)=√L Ĝ(Lt). Em |t|≤504, Plancherel dá ‖∂_L F_L‖≤1/(2L)+504<505 para 0,95≤L≤1, com ‖g‖=1. O residual na banda tem módulo menor que 21 pela cota P239 e pela diferença entre os pisos. Isso limita a derivada da contribuição integral por 21210. As amplitudes dos polos e suas derivadas têm módulo menor que 3. A regra do produto dá mais 36. Usamos 21250 como constante total. A extensão por zero refere-se ao operador limitado auxiliar.

Para o teste de Schur no caso positivo, o orçamento modal de P249 continua válido. A P239 mostra |Ψ₄−β|<1384/75 na banda real até 504. Como |1/4−β|<1/4, |Ψ₄−1/4|<20, o mesmo majorante usado no complemento. Com a≤λ_min do bloco finito, d=1/4−‖D−(1/4)I‖ e b²≥‖B‖², a desigualdade (a−m)(d−m)>b² verifica a margem m no auxiliar. Essa forma de separar a parte limitada corresponde à fatoração de Schur de Dijksma e Langer. [Estimativa P239](../numeria/README_ES8P239_ZHU_EXTENSION.md) · [Dijksma e Langer](https://doi.org/10.1007/s43036-020-00050-0)

## Limiar do primo 7

A perda 2225δ da P257 usa o piso antigo. Ao elevar para 1/4, há nova reserva em |t|≥2520. Binet dá Ψ₄(t)−1/4 > log(5)+1/630−1/4 > 51/40+1/12. A P258 também separa as duas faixas de largura δ=2L−log(7). Para P a projeção de Fourier em |t|≤Ω, o núcleo de 1_E P 1_E é sin(Ω(x−y))/(π(x−y)). Cada bloco diagonal tem norma no máximo δΩ/π. As faixas têm distância superior a 1 quando 2L<log(8), de modo que cada bloco cruzado tem norma no máximo δ/π. O teste de Schur ou a soma de normas dos dois blocos dá ‖1_E P‖²≤δ(Ω+1)/π.

Usando Ω=2520, κ=log(7)/√7<37/50, π>157/50 e τ=51/40, completa-se o quadrado na contribuição da frequência alta. A perda fica abaixo de (Ω+1)κτδ/[π(τ−κ)]<1416δ e resta 1/12 na cauda alta. Assim,

    Q_Weil(f) ≥ R_{1/4}(f) − 1416δ‖f‖² + (1/12)‖P_{|t|>2520} f‖².

Uma margem positiva m independente para o auxiliar permitiria uma extensão local escolhendo δ≤m/2832 dentro da janela auxiliar já certificada. O auditor emite esse resultado só quando os dois setores e o complemento modal passam. [P257](../NUMERIA_ES8P257/README.md) · [Donoho e Stark, projeções de concentração](https://doi.org/10.1137/0149053)

## Reprodução

O comando de extremos testa o sinal da atualização usando G contínuo. O caminho GL60 também confronta todos os intervalos com esse primeiro método.

~~~bash
python -m pip install -r NUMERIA_ES8P249/requirements.txt
python -m unittest NUMERIA_ES8P258.test_p258 -v
python -O -m unittest NUMERIA_ES8P258.test_p258 -v
mkdir -p replay/chunks
python -m NUMERIA_ES8P258.audit --out replay/P258_ENDPOINT.json
seq 0 84 924 | xargs -P 2 -I@ bash -c '
  python -m NUMERIA_ES8P258.concentration \
    --start "$1" --stop "$(($1+84))" --out "replay/chunks/$1.json.gz"
' _ @
python -m NUMERIA_ES8P258.concentration \
  --chunks replay/chunks/*.json.gz --out replay/P258_G504.json.gz
sha=$(sha256sum replay/P258_G504.json.gz | cut -d ' ' -f 1)
python -m NUMERIA_ES8P258.audit --gram replay/P258_G504.json.gz \
  --gram-sha "$sha" --out replay/P258_GL60.json
~~~

O workflow restrito à branch P258 gera os dois relatórios JSON, o agregado GL60 comprimido e o manifesto SHA-256. A investigação textual de Alpöge e da soma ponderada de von Mangoldt nas P255 e P256 continua como linha separada. A identidade de Guinand–Weil é usada na passagem condicional das faixas. Nenhum passo usa a Hipótese de Riemann ou um teorema da formalização de Alpöge como premissa de positividade. [P255](https://github.com/jbarumby/scicalc/pull/61) · [P256](../NUMERIA_ES8P256/README.md)

A execução principal passou 58 testes em Python normal, nas etapas P253, P254, P256, P257 e P258, e passou os cinco testes novos também com -O. A segunda execução conferiu novamente os três SHA-256 da execução principal, recompôs o bloco par e executou o LDL em C/MPFR com arredondamento dirigido. Os arquivos da quadratura **ponderada** de base vieram da P253 autenticada. A P258 regenerou integralmente a quadratura **não ponderada** de G₅₀₄ e repetiu o sinal da cauda. [Execução principal](https://github.com/jbarumby/scicalc/actions/runs/44858303395) · [Execução independente](https://github.com/jbarumby/scicalc/actions/runs/44859003923)

Depois de gerar ou extrair o pacote principal em replay/, o replay independente pode ser repetido em Linux com cc, GMP e libmpfr.so.6.

~~~bash
python -m NUMERIA_ES8P258.postcheck --replay replay --out P258_POSTCHECK.json
sha256sum P258_POSTCHECK.json
~~~
