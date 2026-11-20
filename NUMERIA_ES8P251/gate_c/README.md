# NUMERIA–RIEMANN ES8-P251 · Gate C

A fórmula de Guinand–Weil transfere a margem do Gate B à forma clássica de Weil para funções complexas suaves de suporte compacto na janela. Para `L₀=161/200`, a margem é superior a `1,6186×10^-18`. Uma estimativa unilateral por dilatação dá ainda

```text
Q_W(f) ≥ 5,486×10^-19 ||f||₂²
f ∈ C_c^∞(ℝ;ℂ), supp(f) ⊆ [−L,L],  L₀ ≤ L ≤ L₀ + 5×10^-22.
```

A extensão reutiliza o mesmo corte 309, β da P237 e quatro potências primas. Trata-se de uma extensão local, de tamanho muito pequeno. A extensão cofinal continua em aberto. [Certificado Gate C](data/P251_GATE_C_CERT.json) · [Gate B](../gate_b/README.md) · [P250](../../NUMERIA_ES8P250/README.md)

## Identificação com a forma de Weil

Para uma função da classe indicada, use as convenções

```text
F(z) = ∫ f(x)e^(izx) dx
g(u) = ∫ f(x+u) conjugado(f(x)) dx
H(z) = F(z) conjugado(F(conjugado(z))) = ∫ g(u)e^(izu) du
Q_W(f) = Σ_ρ H((ρ−1/2)/i)
```

A soma percorre os zeros não triviais, com multiplicidades. A função `H` é inteira e decai rapidamente em faixas horizontais limitadas, por integração por partes em `F`. Isso a torna admissível para a fórmula explícita. O valor espectral é `H(z)` nos argumentos complexos dos zeros. Para `t` real, `H(t)=|F(t)|²`. Na variável multiplicativa, `φ(u)=u^(-1/2) f(log u)` satisfaz `Mφ(s)=F((s−1/2)/i)` com `Mφ(s)=∫₀∞ φ(u)u^(s−1)du`. A involução `u^(-1) conjugado(φ(1/u))` produz a conjugação em `1−conjugado(s)`, como na forma hermitiana de Weil. [Guinand–Weil, Lema 2.2](https://par.nsf.gov/servlets/purl/10637681)

Na convenção do lema citado, a transformada de `H` com núcleo `e^(−2πitξ)` vale `2πg(2πξ)`. A parcela de primos é portanto `−Σ_q Λ(q)/√q · [g(log q)+g(−log q)]`. O suporte de `g` está em `[−2L,2L]`. Pelo conjunto de potências primas na janela, a inversão de Fourier reescreve a parcela como o multiplicador `−S(t)` abaixo. Os dois termos arquimedianos da fórmula dão `Re ψ(1/4+it/2)−log π`. Assim,

```text
S(t) = Σ_{q∈{2,3,4,5}} [2Λ(q)/√q] cos(t log q)
Ψ(t) = Re ψ(1/4+it/2) − log π − S(t)
ℓ_±(f) = ∫ f(x)e^(±x/2) dx
Q_W(f) = 2 Re(ℓ_+(f) conjugado(ℓ_−(f))) + (1/(2π))∫ℝ Ψ(t)|F(t)|² dt.
```

A parcela dos polos é `2(|∫ f cosh(x/2)|²−|∫ f sinh(x/2)|²)`. O setor par tem o primeiro sinal e o ímpar, o segundo. No termo integral, a interferência entre as transformadas par e ímpar é ímpar em `t` e integra a zero contra `Ψ`. Em cada setor puro, o módulo ao quadrado da transformada é par, e a integral bilateral vira a integral no semieixo com fator `1/π`. Esses argumentos valem para coeficientes complexos.

Na base ortonormal de Legendre de `L²([-L,L])`, a transformada do grau `n` é `i^n√(2L(2n+1))j_n(Lt)`. A fase comum do setor ímpar cancela no produto hermitiano, dando as amplitudes reais usadas pelo Gate B. Avaliar a transformada em `−i/2` dá o vetor de polo com entradas `√(2L(2n+1))i_n(L/2)`, onde `i_n` é Bessel esférica modificada. [DLMF §10.54](https://dlmf.nist.gov/10.54.E2) · [Zhu, fórmula geométrica e paridades](https://arxiv.org/html/2608.24827v2)

Ponha `w(t)=Ψ(t)−β` e `T=309`. O operador-modelo usa os mesmos polos, a diagonal `βI` e o integral de `w` em `|t|≤T`. Por Plancherel, sua diferença para a forma de Weil é exatamente

```text
Q_W(f) − R_L(f) = (1/(2π))∫_{|t|>T} w(t)|F(t)|² dt ≥ 0.
```

O sinal usa a P250 e Binet com o β fixado em 504. O operador-modelo é limitado em `L²([-L,L])`, pois é composto de uma diagonal constante, uma compressão de multiplicador limitada a uma banda finita e os polos de posto finito. Por isso o complemento da base de Legendre pode ser tratado em `L²` pelo Gate B. A passagem à forma de Weil acima usa diretamente sua classe suave admissível. [Replay do sinal P250](../../NUMERIA_ES8P250/README.md)

## Extensão unilateral da janela

A dilatação unitária leva `h∈L²([-1,1])` a `f(x)=L^(-1/2)h(x/L)`. Com a mudança `ξ=Lt`, o termo de frequência finita vira

```text
C_L(h) = (1/(2π))∫_{|ξ|≤309L} w(ξ/L)|ĥ(ξ)|² dξ.
```

Na faixa até o novo extremo, `5<exp(2L)<7`. As potências primas são as mesmas, pois 6 tem von Mangoldt nulo. Logo `w` e β são comuns às duas janelas. As desigualdades racionais do auditor dão `Σ_q a_q log(q) < 11/2`. Para o termo de digama, a série `ψ'(z)=Σ_{k≥0}(z+k)^−2` dá `|(Re ψ(1/4+it/2))'|<9` por comparação com `Σ_{k≥1}k^-2<2`. Para `t≥2`, com `y=t/2`, a soma multiplicada por `y` é limitada pelo primeiro termo e a integral monótona, dando `t |(Re ψ)'| ≤ 1/y+π/2<3`. Segue `|w'|<29/2` em toda a banda e `|t w'(t)|<1702,5` em `[0,309]`. [DLMF §5.15](https://dlmf.nist.gov/5.15.E1)

No trecho espectral comum, derivar `w(ξ/L)` em relação a `L` e usar Plancherel limita a perda por `(1702,5/L₀)(L−L₀)||h||₂²`. Na faixa acrescentada, o argumento de `w` dista de 309 por no máximo `309(L−L₀)/L₀`. A margem P250 em 309 excede `1794/10^7`. A cota de `|w'|` preserva o sinal nessa faixa em todo o incremento proposto.

Os vetores de polo sob a dilatação são `a_±(x)=√L e^(±Lx/2)`, em `[-1,1]`. Cada derivada tem módulo pontual limitado por `(1+L)/(2L)` vezes o módulo de seu vetor. A regra do produto no operador de posto dois dá limite estrito `4(1+L)e^L<24` para a derivada em norma. O β fica fixo. Juntando os termos, a perda estrita para `L>L₀` é menor que `2140(L−L₀)||h||₂²`. O auditor confirma em frações

```text
1,6186×10^-18 − 2140·(5×10^-22) = 5,486×10^-19.
```

## PDF e circuitos quânticos

O capítulo sobre spin relativístico e a equação de Dirac ajuda a interpretar a decomposição por paridade e a simetria de operadores autoadjuntos. A decomposição usada aqui decorre da reflexão `x↦−x` e do símbolo par `Ψ`. A identidade de Weil permanece a da fórmula explícita, com seus fatores arquimedianos e primos. [Livro de Gordon Baym](https://www.routledge.com/link/link/p/book/9780367091873)

O GitHub já contém notebooks de um qubit para IBM e Google, uma comparação de circuitos reversíveis de inversão modular, a reprodução automatizada do divsteps e a auditoria Spin–Zeta. O simulador do NUMERIA-Q 05.21 passou em oito testes e 3.200 verificações exaustivas da base computacional na reexecução. O divsteps automatiza verificações de aritmética modular. [Google](https://github.com/jbarumby/scicalc/pull/3) · [IBM](https://github.com/jbarumby/scicalc/pull/4) · [NUMERIA-Q 05.21](https://github.com/jbarumby/scicalc/pull/2) · [Divsteps](https://github.com/jbarumby/scicalc/pull/5) · [Spin–Zeta](https://github.com/jbarumby/scicalc/blob/83c58fd41df3748f3ba0d27be21266b0813fda55/research/2026-11-03/es3/README.md)

A avaliação para o Gate B inclui novos limites variacionais em Arb. Os vetores de tentativa são diádicos exatos. O quociente de Rayleigh é encerrado com o erro de montagem e a cota GL60. O menor autovalor par é inferior a `3,2374×10^-18` e o maior supera `4,65`. O diagnóstico certifica condicionamento `κ₂>1,43×10^18` no setor par. Isso mede a escala real do problema. Nove qubits bastam para etiquetar os 259 modos. A estimação de fase direta ainda exige implementar a ação do operador, preparar sobreposição com seu menor modo e resolver sua escala em relação à norma. A construção padrão com potências controladas tem custo de resolução proporcional ao inverso dessa separação normalizada. Uma eventual vantagem com outra codificação ou com precondicionamento exigiria uma análise própria. Para o replay atual, foi implementada a automação clássica com limites dirigidos. [Dados espectrais](data/P251_GATE_C_SPECTRAL.json) · [Estimação por fase e preparação de estado](https://doi.org/10.1103/PRXQuantum.3.010318)

## Reprodução automatizada

Na raiz do repositório, com as versões da P249 instaladas:

```bash
python -m pip install -r NUMERIA_ES8P249/requirements.txt
python -m NUMERIA_ES8P251.gate_c.bridge_audit --out /tmp/P251_GATE_C_CERT.json
python -m NUMERIA_ES8P251.gate_c.spectral_diagnostic --out /tmp/P251_GATE_C_SPECTRAL.json
python -m unittest NUMERIA_ES8P251.gate_c.test_bridge_audit -v
python -O -m unittest NUMERIA_ES8P251.gate_c.test_bridge_audit -v
sha256sum -c NUMERIA_ES8P251/gate_c/SHA256SUMS
cmp /tmp/P251_GATE_C_CERT.json NUMERIA_ES8P251/gate_c/data/P251_GATE_C_CERT.json
cmp /tmp/P251_GATE_C_SPECTRAL.json NUMERIA_ES8P251/gate_c/data/P251_GATE_C_SPECTRAL.json
```

O primeiro comando de auditoria autentica o agregado e o certificado do Gate B, reexecuta a certificação Arb e o complemento modal, recompõe a P250 com aritmética dirigida e verifica as novas desigualdades racionais. A derivação analítica acima usa teoremas clássicos. Passaram os seis testes do Gate C e os 41 das etapas anteriores, tanto em Python normal quanto com `-O`. A próxima questão quantitativa é obter margens úteis e cortes controlados para janelas macroscópicas e, depois, para uma sequência cofinal. A Hipótese de Riemann permanece em aberto.
