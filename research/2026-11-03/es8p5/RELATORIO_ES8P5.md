# NUMERIA–RIEMANN — ES8-P5: reconciliação dos arquivos históricos E7/V16

**Data:** 3 de novembro de 2026. **Escopo:** validação de proveniência e compatibilidade matemática *finita* das matrizes projetadas. **Não é demonstração da positividade universal de Weil nem da Hipótese de Riemann.**

## 1. Evidência primária recuperada

O usuário disponibilizou cinco arquivos:

- `E7_original (1).bin`, `E7(1).zip` e `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3(1).zip`: exatamente 1.641.486 bytes cada, SHA-256 idêntico `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7`.
- `E7_DENTRO_DO_ZIP(1).zip`: contêiner externo com o mesmo E7, de SHA-256 `135c46704a009e1f5af90ee970fda880f518af5e8cc2a2206770dee13fb61c19`.
- `E8K_V16_GITHUB_COMPLETO(1).zip`: SHA-256 `aee5bfedf977f7a0a77d66949753faaab673bb9fa18b01b4425160e464936826`. Contém os ZIPs E7 e Fase4, com SHA-256 `4acf99...` e `b0d7d335534239a6863e655295a658282551d334b0808ee801710a4b5699d8fe`, respectivamente.

Executou-se `ZipFile.testzip()` nos arquivos, a verificação de **7/7** hashes internos do `SHA256SUMS.txt` da V16, de **32/32** hashes do manifesto E7 e de **30/30** hashes individuais da Fase4. Identificadores Git de objetos e SHA-256 de conteúdo **não são a mesma coisa**.

A execução **sem modificação** do código `v16/AUDITAR_E8K_V16_FORMULA.py` forneceu 12/12 gates, com pivôs positivos **condicionais** relatados em N=6,8,10 e um JSON final exatamente igual, byte a byte, ao `v16/RESULTADO_E8K_V16.json` que veio no ZIP.

## 2. Teste técnico: 30 integrais arquimedianas

O E7/V16 constrói para modos centrados

\[
F_k(r)=\int_{-L}^{L}\phi_k(x)e^{irx}\,dx,\quad
\Omega(r)=\operatorname{Re}\psi(1/4+ir/2)-\log\pi,\quad L=9/20,
\]

\[
G_{ij}=\frac1\pi\int_0^\infty
\operatorname{Re}(F_i(r)\overline{F_j(r)})\Omega(r)\,dr.
\]

A Fase4 certifica **condicionalmente a seus majorantes analíticos** as integrais até R=160; o E7 fornece uma cota simétrica do módulo da cauda em [160,∞). A ES8-P3 reconstrói o termo completo por uma série exata de núcleos exponenciais com caudas Euler–Maclaurin intervalares. Comparamos, para cada par de mesma paridade, os intervalos **já acrescidos da cauda** E7 com os intervalos novos ES8-P3. Resultado: **30/30 interseções**.

A necessidade de comparar as mesmas integrais é mensurável: sem agregar a cauda analítica, apenas **2/30** janelas Fase4 em [0,160] interceptam os intervalos ES8-P3 para [0,∞). Isso não é uma falha da Fase4: as integrais finita e infinita são quantidades diferentes.

## 3. Primo 2, cancelamento de polos e matrizes

O suporte \(2L=0,9\), com \(\log2<0,9<\log3\), deixa apenas a potência prima 2. A convenção comparada é

\[
P_{ij}=-\frac{2\log2}{\sqrt2}\,C_{ij}(\log2),
\]

com a convolução real e simétrica `C`. Verificaram-se **30/30** interseções dos intervalos da V16 com os da ES8-P2.

Cada coluna de projeção indexada por \(k\ge3\) tem pivô \(p=1\) para k ímpar e \(p=2\) para k par, e a entrada de cancelamento

\[
B_{p,k-2}=-\frac{k}{p}\frac{p^2\pi^2+L^2}{k^2\pi^2+L^2}.
\]

As **8/8** entradas históricas desses coeficientes intersectam as reconstruídas, preservando a identidade racional de anulação dos dois momentos de polos.

Por fim, compararam-se as matrizes reais simétricas

\[
A_N=B^{\mathsf T}(G+P)B,\quad N\in\{6,8,10\}.
\]

| N | Dimensão | Entradas da parte triangular superior comparadas | Interseções |
|---:|---:|---:|---:|
| 6 | 4×4 | 10 | 10 |
| 8 | 6×6 | 21 | 21 |
| 10 | 8×8 | 36 | 36 |
| **Total** | | **67** | **67** |

A ES8-P3 reproduz a reconstrução infinita usando outra estratégia de integração. A Fase4 usa 1.600 células em [0,160]; sua projeção é refeita a partir dos dados originais e não inferida do JSON de pivôs. A comparação mantém os sinais, as projeções e a convenção de paridade. O código exige **todos** os encontros de intervalos, e qualquer falha é reportada.

## 4. Interpretação rigorosa

O resultado desta fase é **compatibilidade intervalar**, com proveniência histórica verificável. Interseção de intervalos não garante identidade das integrais reais sem provar que ambos os intervalos realmente as contêm. A autenticidade de arquivo não autentica cada derivação matemática.

Permanecem hipóteses ou obrigações não encerradas aqui:

1. Demonstrar independentemente cada grandeza e majorante complexo da Fase4 (interpolação, Cauchy e arredondamento), a cauda analítica E7 e as cotas Euler–Maclaurin da ES8-P3.
2. Identificar plenamente a forma \(B^\top(G+P)B\) com o critério de Weil para as funções de teste admissíveis e os termos de polos, inclusive todos os fatores e convenções de Fourier. A fórmula explícita de Guinand–Weil é um teorema clássico externo; sua especialização ao código deve ser demonstrada.
3. Não confundir **pivôs LDL** com menores autovalores na norma original; a base B não é ortonormal.
4. Não estender certificados em N=6,8,10 a N arbitrário, a todo L ou a todas as funções-teste. **RH permanece em aberto.**

## 5. Reprodutibilidade e GitHub Actions

Há **dois níveis de reprodução**:

- **Reabertura integral**: `--v16-upload` calcula hashes de todos os dados dentro dos arquivos ZIP, carrega o código V16 original e exige um replay byte a byte de seu JSON. Foi realizado localmente nesta execução, juntamente com os quatro arquivos E7 redundantes/aninhados.
- **Reprodução por fixture**: `--fixture` usa os campos racionais extraídos do V16 após a verificação integral, recompõe com a fonte da ES8-P3 os 30 intervalos arquimedianos, os 30 primos, as 8 colunas e as 67 entradas. **Sem ZIP binário, o GitHub Actions executa somente este segundo nível**, não deve alegar reabertura do arquivo original.

A fixture, os hashes e o resultado completo ficam sob `research/2026-11-03/es8p5`; os arquivos binários originais permanecem no ZIP entregue para reprodução integral, não são publicados novamente no GitHub. Os testes incluem mutações de sinal, projeção e intervalos que deverão falhar.

**Conclusão:** a divergência potencial entre a reconstrução ES8-P3 e as fontes históricas E7/V16 foi transformada em **97 verificações de entradas** (30 arquimedianas + 67 projetadas), mais 30 verificações da parcela de primos e 8 coeficientes de projeção, todas compatíveis nos intervalos analisados. Isso encerra a lacuna *documental e de comparação finita* relativa ao material agora recebido; não encerra a prova matemática global.
