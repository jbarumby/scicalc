# NUMERIA P90 — Resultado confirmado e escopo científico

**Execução CI:** [GitHub Actions #45051710194](https://github.com/jbarumby/scicalc/actions/runs/45051710194), **completed/success**, 22/11/2026.

**Artefato:** [numeria-p90-sharp-neighbor-and-one-interior](https://github.com/jbarumby/scicalc/actions/runs/45051710194/artifacts/14766240834), ID 14766240834. Inclui `p90-nearby-sharp.json` e `p90-SHA256SUMS.txt`.

SHA-256 exato do JSON produzido pela CI:
`be73fd1f1d24f99579ebfc36f5ba8175411a21daf43bb3a0d051aed72d3af9c4`.

## 1. Teorema analítico independente novo

Demonstramos a **desigualdade afiada de Preissmann Lema 6 com a constante original** para qualquer configuração de pontos reais finitos e distintos e pares `a<b` com:
- `b=a+1`, sem pontos intermediários;
- `b=a+2`, exatamente um ponto intermediário.

**Os pontos exteriores são completamente arbitrários** em número e espaçamento.

A demonstração, no [README](README.md), é uma prova matemática geral por convexidade e comparação com as caudas exteriores da P89 — **não é generalização dos 96 casos finitos**.

Com `h=x_b-x_a`, `s_a=delta_a/h`, `s_b=delta_b/h`, `Z=zeta(2)=pi²/6`,
`S(s)=sum_{j>=1} 1/[j²(1+js)²]`, `D(s)=(2Z-S(s))/s−3s`:
- `D(s)` é decrescente (a derivada do núcleo é o negativo do polinômio `1+3z+6z²+2z³` sobre denominador positivo);
- `S(0)=Z`, `S(1/2)=2Z−11/4`, `S(1)=2Z−3`;
- pares adjacentes: `D(s)>=(1−s)(Z/s+3)>=0`;
- pares com um nó interno em `0<t<=1/2`: a estimativa via cordas convexas reduz o orçamento interior ao fator `(1−2t)(1+t)/2>=0`, usando apenas `Z>3/2`.

O termo afiado `−3(delta_a+delta_b)/h⁴` foi **preservado sem enfraquecimento**.

## 2. Testes e certificados exatos

- **8 testes unitários em Python normal e 8 em `-O`**, todos aprovados.
- **96 verificações racionais de pares próximos**, por `fractions.Fraction`, usando a cota **INFERIOR** de π derivada rigorosamente por Machin. As margens foram todas positivas.
- **40 identidades de malha uniforme** em `m=1,...,40`, que confirmam `M_uniform(m)=2D(1/m)` em frações exatas. A prova analítica geral usa somas harmônicas e não depende dos 40 exemplos.
- **8 verificações racionais de orçamento interior** em malhas irregulares com `m>=3` lacunas; verificações **apenas finitas**, não prova universal.
- Derivada piecewise exata para o caso de um ponto móvel no intervalo e checagem dos sinais ao redor de `1/3`, `1/2` e `2/3`.
- A CI compara os relatórios, exibe SHA-256 dos códigos/insumos e bloqueia `full_Preissmann_Lemma6_independently_reproved=true`.

## 3. Dados do usuário P87, agora completos

Foram examinados:
- `p87_teorema_preissmann (1).tex` SHA `d16771dab750ef4f00e00bf841ddcf5f179e412ef6392914ea819a0e7b9d706d`;
- `p87_verificacao (1).py` SHA `5089f7016827dffb5a6202cd783ca70ab099ba5a87db37e954e2192a1b8ff14a`;
- `p87_ledger.json` SHA `f593dc80743a5510c064aaa6ad3d4b2149f9950389c2294e16e425121b728da7`;
- `p87_execucao.log` SHA `1b594dda439681ab78a286dca2aa810c3dab4b2149f9950389c2294a16`.

O JSON confirma `gates_total=24`, `gates_failed=0` e **`imported.proved_here=false`**; o log corresponde ao mesmo SHA de certificado. Essas são evidências de execução dos gates programados, **não 24 provas analíticas independentes**.

A errata P89 continua válida: representação de Laplace na nota original sem `1/t`, desigualdade exponencial com sinal invertido e uso de **limite superior de π** no bloco B5. O uso de `numpy/mpmath` em B3/B4/B6/B7 constitui diagnóstico, não prova universal.

## 4. Região ainda aberta (P91)

Com `m>=3` lacunas `g_1+...+g_m=1`, `g_j>0`, e `t_k=sum_{j<=k}g_j`, definimos o orçamento
\[
\mathcal M(g):=\sum_{k=1}^{m-1}
\frac{\min(g_k,g_{k+1})}{t_k^2(1-t_k)^2}.
\]

**O seguinte é apenas um CRITÉRIO SUFICIENTE CONJECTURADO nesta fase:**
\[
\mathcal M(g)\stackrel{?}{\le}D(g_1)+D(g_m)\qquad
\text{para todo }m\ge3\text{ e todo }g>0,\ \sum g=1.
\]

A desigualdade é **exata com igualdade** em malhas internas uniformes (`g_j=1/m`), mas essa igualdade não prova o caso geral. A comparação variacional de Preissmann Lema 5 distingue regimes de crescimento/decrescimento das folgas; `variational.py` mostra que um único ponto interno já apresenta mudanças de sinal na derivada de seu termo, portanto um argumento ingênuo de suavização falha.

**Próxima prioridade P91:** provar ou refutar o critério suficiente universal na região `m>=3`, dividindo o espaço das folgas nos cones definidos por `min(g_k,g_{k+1})` e reconstruindo as mudanças de espaçamento do Lema 5. Uma eventual refutação do **critério suficiente** não refuta Preissmann Lema 6: significaria apenas que o envelope de caudas exteriores da P89 é demasiado frouxo.

## 5. Limites do projeto NUMERIA–RIEMANN

A desigualdade ponderada geral de Preissmann continua teorema publicado usado como entrada na cota de `F_X`, pois **falta o caso com dois ou mais pontos entre a e b**. Ramaré–Zúñiga `S_epsilon(X)` e `3e10<=X<3e12` seguem externos. A janela Suzuki até `t=79` permanece finita/condicional. **Nenhuma prova global da RH**.

PR segue **draft**, sem merge automático.
