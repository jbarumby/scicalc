# Gate D — fundamentos analíticos revisados do Gate C

**Autor: Manus AI**

Esta nota fecha uma lacuna de reprodutibilidade do NUMERIA: os scripts P237 que justificam a cota complexa do símbolo e o orçamento de erro GL60 existiam no diretório de trabalho, mas não na árvore versionada. Agora os certificadores racionais, os portões geométricos E3 e os testes correspondentes podem ser reproduzidos a partir de `numeria/p237/`.

O resultado é uma prova analítica **condicional às identidades clássicas citadas abaixo**, combinada com verificações por `Fraction` exatas. Não se trata de um novo critério para RH, nem de prova da positividade global de Weil. O ponto de valor é transformar os fatores numéricos do Gate C em hipóteses explícitas, rastreáveis e testáveis, em vez de deixá-los como constantes opacas.

## Geometria exata das elipses

Para `rho=3`, a elipse de Bernstein tem semieixos `5/3` e `4/3`. Sob o mapa de cada painel `t=(2p+1)/4+s/4`, `p=0,...,1007`, as meias larguras tornam-se `5/12` e `1/3`. Assim `|Re(t)| <= 3025/6` e `|Im(t)| <= 1/3`.

Para os dois argumentos holomorfos da digama, `z_±=1/4±it/2`, resulta `1/12 <= Re(z_±) <= 5/12` e `|z_±| <= |Re(z_±)|+|Im(z_±)| <= 505/2 < 253`. Em particular, os argumentos permanecem no semiplano direito e não encontram os polos da digama.

Para `L=161/200`, o retângulo que contém `Lt` satisfaz

`|Re(Lt)| <= 19481/48`, `|Im(Lt)| <= 161/600 < 1/2`,

`|Lt|^2 <= (19481/48)^2 + (161/600)^2 < 406^2`.

Essas desigualdades são recalculadas exatamente em `p237_analytic_gates.py`; as comparações falham com exceção explícita, inclusive quando Python é executado com `-O`.

## Cota do símbolo complexo

A identidade de frações parciais da digama é

`psi(z)=-gamma-1/z+sum_{n>=1} z/[n(n+z)]`.

Para `Re(z)>0`, `|z|<253` e `|z|>=1/12`, vale `|n+z|^2=|z|^2+2n Re(z)+n^2>|z|^2`; logo cada termo é estritamente menor que `1/n`. Nos termos `n>253`, também `|n+z|>=n+Re(z)>n`, então `|z/[n(n+z)]|<253/n^2` e `sum_{n>253} n^-2 < integral_253^infinity x^-2 dx=1/253`. Como `H_253<7` (verificado por soma racional exata) e `gamma<1`, resulta `|psi(z)|<1+12+7+1=21`. A identidade e seu domínio estão no DLMF, eq. 5.7.6 [1].

No prolongamento holomorfo do símbolo, a parte de digama é a média de `psi(z_+)` e `psi(z_-)`. Ambos obedecem ao mesmo limite, logo a média também é menor que 21 em módulo. Para os quatro termos oscilatórios indexados por `q=2,3,4,5` — com `q=4` correspondendo ao segundo harmônico de 2 —, a soma dos coeficientes é menor que 11. Além disso, `|cos(t log q)|<=cosh(2/3)<=9/7`; a última desigualdade vem da série de `cosh` e de `(2k)!>=2^k`. Com `|log(pi)|<2` e `0<beta_504<2`, certificado racionalmente por `p237_beta_tight.py`, obtém-se

`|Psi_L(t)-beta_504| < 21+2+11(9/7)+2 = 274/7 < 40 < 50`.

`p237_symbol_majorant.py` verifica as desigualdades finitas racionais e chama o certificado exato de beta. As identidades analíticas — incluindo frações parciais da digama, conjugação, `gamma<1` e séries de `cosh` — continuam explicitamente listadas como insumos clássicos; o script não afirma formalizá-las em Lean.

## Cotas complexas dos fatores de Bessel

A cota uniforme não precisa supor `|P_n(x)|<=1`. Pela representação plano-de-onda da Bessel esférica, DLMF 10.54.2 [2],

`j_n(z)=(-i)^n/2 integral_{-1}^1 exp(izv) P_n(v) dv`.

Pela desigualdade de Cauchy–Schwarz e pela norma de Legendre `integral_{-1}^1 P_n(v)^2 dv=2/(2n+1)` [3],

`|j_n(z)| <= exp(|Im(z)|)/sqrt(2n+1) <= exp(|Im(z)|)`.

Como `|Im(Lt)|<=161/600<1/2`, a série exponencial dá `exp(|Im(Lt)|)<=exp(1/2)<2`. A segunda representação DLMF 10.54.1 [4], junto com `|cos(zv)|<=exp(|Im(z)|)` para `v` real e a integral beta de `sin^(2n+1)`, dá

`|j_n(z)| <= exp(|Im(z)|) |z|^n/(2n+1)!! < 2*406^n/(2n+1)!!`.

Combinando os dois limites e o teto inteiro para a normalização, chega-se exatamente ao majorante de amplitude usado no Gate C:

`a_n = ceil(sqrt(2L(2n+1))) * min(2, 2*406^n/(2n+1)!!)`.

## Erro de quadratura e orçamento de linha

Para `g(s)=u_i(t(s))u_j(t(s))(Psi_L(t(s))-beta)/(4*pi)` em cada painel, o majorante do integrando é `M_ij=(50/12)a_i a_j`. O Teorema 8.1 de Trefethen dá `|c_k|<=2M_ij*3^-k` para os coeficientes de Chebyshev quando a função é holomorfa na elipse `E_3` [5]. A cauda a partir do grau 120 é no máximo `3M_ij/3^120`.

A regra Gauss–Legendre com 60 nós é exata até o grau 119, tem pesos positivos, e a soma dos pesos é 2; as propriedades gerais de Gauss quadrature constam em DLMF 3.5(v) [6]. Como tanto a integral quanto a quadratura têm norma de operador 2, o erro por painel é no máximo `4` vezes a cauda de Chebyshev. Isso resulta em

`|integral(g)-GL60(g)| <= 50*a_i*a_j/3^120`.

Somando os 1008 painéis, obtém-se a cota de entrada `(1008*50/3^120)*a_i*a_j`. O certificador P237 recalcula os 754 majorantes inteiros/racionais e confirma, para cada paridade, erro de soma de linha inferior a `10^-46`; também verifica o limite de erro diagonal do último modo. O Gate C amostral P245 continua sendo outra coisa: seus nós e entradas selecionadas não são uma avaliação completa de todos os elementos de todos os blocos.

## Replay independente de três chunks adicionais

O replay P245D recalculou as contribuições GL60 em três janelas de 25 painéis: `[25,50]`, `[500,525]` e `[975,1000]`. Em cada janela foram avaliados 13 pares da matriz, nas duas paridades, em todos os 60 nós por painel. Isso totaliza 117.000 avaliações de par-painel e 78 inclusões de intervalos de soma de chunk: em todos os casos, a bola Arb P237 de 900 bits ficou contida no intervalo racional de ponto fixo P245D de 256 bits.

Para a primeira janela, os graus de Bessel baixos foram calculados por uma série inteira com cauda geométrica racional; para as duas janelas de frequência alta, usou-se a recorrência esférica intervalar até grau 33. O modo índice 376 não foi calculado ponto a ponto: foi contido por uma majorante analítica. Cada resultado JSON inclui os intervalos P245D, a bola Arb de referência original, seus extremos racionais e o SHA-256 do chunk completo usado como fonte. Os seis hashes esperados estão fixados no runner.

Esse é um cruzamento independente dos **13 pares escolhidos**, não uma reconstrução dos 71.253 elementos triangulares por setor em cada chunk. Os seis chunks completos de aproximadamente 10 MB cada não foram duplicados no repositório; para reproduzir o cálculo, forneça-os ao argumento `--references`. Os certificados JSON guardam as referências dos pares selecionados e permitem auditar as inclusões sem abrir os arquivos comprimidos completos.

## Artefatos e reprodução

Os arquivos versionados nesta etapa são:

- `numeria/p237/p237_analytic_gates.py`: geometria E3 e domínios complexos em aritmética racional exata;
- `numeria/p237/p237_symbol_majorant.py`: majorante do símbolo e beta racional;
- `numeria/p237/p237_gl60_analytic_certificate.py`: pesos, erro GL60 e orçamento por paridade;
- `numeria/p237/test_p237_analytic_certificates.py`: testes dos três portões, incluindo execução fail-closed;
- `numeria/p245/gate_c/p245d_additional_chunk_replay.py`: replayer racional de pares selecionados para os três blocos adicionais;
- `numeria/p245/gate_c/test_p245d_additional_chunk_replay.py`: testes de indexação triangular e pins dos dados.

Os certificados gerados são `numeria/p237/P237_ANALYTIC_GATES.json`, `P237_SYMBOL_E3_RATIONAL_GATE.json` e `P237_GL60_ANALYTIC_CERT.json`, além de `P245D_SELECTED_25_50.json`, `P245D_SELECTED_500_525.json` e `P245D_SELECTED_975_1000.json` em `numeria/p245/gate_c/`.

Do diretório raiz do repositório:

```bash
python3 numeria/p237/p237_analytic_gates.py --out /tmp/P237_ANALYTIC_GATES.json
python3 numeria/p237/p237_symbol_majorant.py --out /tmp/P237_SYMBOL_E3.json
python3 numeria/p237/p237_gl60_analytic_certificate.py --out /tmp/P237_GL60.json
python3 -m unittest discover -s numeria/p237 -p 'test_p237_analytic_certificates.py' -v
python3 -O -m unittest discover -s numeria/p237 -p 'test_p237_analytic_certificates.py' -v
python3 numeria/p245/gate_c/p245d_additional_chunk_replay.py --references /path/to/P237/chunks --roots numeria/p245/gate_c/P245B_ROOT_BRACKETS.json --start 25 --out /tmp/P245D_25_50.json
python3 numeria/p245/gate_c/p245d_additional_chunk_replay.py --references /path/to/P237/chunks --roots numeria/p245/gate_c/P245B_ROOT_BRACKETS.json --start 500 --out /tmp/P245D_500_525.json
python3 numeria/p245/gate_c/p245d_additional_chunk_replay.py --references /path/to/P237/chunks --roots numeria/p245/gate_c/P245B_ROOT_BRACKETS.json --start 975 --out /tmp/P245D_975_1000.json
```

## O que permanece por provar

A derivação analítica do orçamento de quadratura não autentica, por si só, que a implementação P237 coincide com a redução de Weil–Zhu em todos os termos e convenções. Também não transforma um certificado de bloco finito em positividade global em todo o espaço de funções-teste. Essa passagem exige uma identidade funcional completa e uma desigualdade quantitativa uniforme para a cauda e seu acoplamento, sem perda da margem espectral. Logo, o resultado continua sendo um certificado de janela compacta, não uma prova da Hipótese de Riemann.

Este documento não reivindica novidade bibliográfica para as identidades clássicas de análise numérica. A contribuição técnica aqui é a derivação específica no parâmetro `L=0,805`, o fechamento exato das constantes do Gate C e uma trilha reproduzível que liga essas hipóteses às margens finitas já auditadas.

## Referências

[1]: https://dlmf.nist.gov/5.7.E6 "NIST DLMF, Digamma partial-fraction expansion, equation 5.7.6"
[2]: https://dlmf.nist.gov/10.54.E2 "NIST DLMF, Plane-wave integral for spherical Bessel functions, equation 10.54.2"
[3]: https://dlmf.nist.gov/18.3 "NIST DLMF, Legendre orthogonality and squared norm, Table 18.3.1"
[4]: https://dlmf.nist.gov/10.54.E1 "NIST DLMF, spherical-Bessel integral representation, equation 10.54.1"
[5]: https://people.maths.ox.ac.uk/trefethen/trefethen_sample.pdf "Nick Trefethen, Approximation Theory and Approximation Practice, Theorem 8.1"
[6]: https://dlmf.nist.gov/3.5 "NIST DLMF, Gaussian quadrature, exactness and positive weights, Section 3.5(v)"
