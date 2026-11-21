# NUMERIA–RIEMANN · ES8-P258 · Auxiliar reforçado no limiar de 7

A P258 parte do commit [bb6ace41](https://github.com/jbarumby/scicalc/commit/bb6ace41f4e9739cde488d3952803b78e5a92118) do PR #63. O cálculo testa o modelo auxiliar truncado de quatro potências primas, com corte 504, piso espectral 1/4 e meia-largura racional 0,972955074527656652552676371721590. Essa largura fica cerca de 1,35×10⁻³⁴ acima de log(7)/2. Há 377 modos por setor de paridade.

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

Para o teste de Schur no caso positivo, o orçamento modal de P249 continua válido. A P239 mostra |Ψ₄−β|<1384/75 na banda real até 504. Como |1/4−β|<1/4, |Ψ₄−1/4|<20, o mesmo majorante usado no complemento. Com a≥λ_min do bloco finito, d=1/4−‖D−(1/4)I‖ e b²≥‖B‖², a desigualdade (a−m)(d−m)>b² verifica a margem m no auxiliar. Essa forma de separar a parte limitada corresponde à fatoração de Schur de Dijksma e Langer. [Estimativa P239](../numeria/README_ES8P239_ZHU_EXTENSION.md) · [Dijksma e Langer](https://doi.org/10.1007/s43036-020-00050-0)

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
