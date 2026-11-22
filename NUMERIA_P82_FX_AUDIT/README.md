# NUMERIA–RIEMANN P82 — Auditoria independente da estrutura de F_X

**Objeto principal:** desigualdade para `F_X(1+δ/log X,0,T,1/2)` no Lema 6 de Chourasiya–Simonič v2, a partir da definição de `f_X` e das referências originais de Kadiri–Lumley–Ng (KLN) e Ramaré–Zúñiga.

**Escopo das conclusões:** prova matemática direta de identidades aritméticas e transformações algébricas, mais **derivação condicional** da desigualdade de F_X com três teoremas externos claramente apontados. Não afirma auditoria independente completa desses três teoremas.

## Arquivos recebidos, incorporados como fontes para revisão

- `p78_s2_passagem_completa (3).pdf`: SHA-256 `a52ae9cac7777130c5383b4c70267663bbf0594c05798d6b3b001e32c06d90ab`;
- `p78_s2_proposicao_auxiliar_b1.pdf`: SHA-256 `acd532d7f48c828577f3d1080f19f6741ab23af6594412b34ac5dbe7ed246841`;
- `p78_s2_proposicao_auxiliar_capa.html`: SHA-256 `a81879bc525fb4dd6ea4ed2d930dfbc2d2b3d30e011a07e0f6a9d2ce9a096994`.

A segunda nota corrige a classificação da "Proposição auxiliar" e a condição **estrita** `b4 < b4_star`. Isso é proveitoso e deve ser preservado. A primeira nota ainda contém a prova não fechada de `G'(\varepsilon)<0` por termo `O(\varepsilon²)` não majorado e considera como "exatamente" a igualdade `L*=logT0` de um ajuste. **Não usar esses argumentos na prova P82**. O ramo ζ foi refeito no PR #81 usando `ζ(1+ε)<=1+1/ε` e reservas racionais de b4.

## 1. Coeficientes de Dirichlet e um termo constante essencial

Para `X≥1` inteiro, `M_X(s)=Σ_{d≤X} μ(d)d^(−s)` e `f_X(s)=ζ(s)M_X(s)−1`.
Para `Re(s)>1` temos

```
a_X(n) = Σ_{d|n, d≤X} μ(d)
λ_X(n) = a_X(n) − 1_(n=1).
λ_X(n)=0 para 1≤n≤X.
```

Como `Σ_{d|n} μ(d)=0` para `n>1` e `μ(1)=1`, a igualdade acima é **exata**.
Para `τ>1`, a expansão do quadrado e a troca de somas absolutamente convergentes dão

```
ζ(τ) S_(τ−1)(X)
 = Σ_{d,e≤X} μ(d) μ(e) ζ(τ)/lcm(d,e)^τ
 = Σ_{n≥1} a_X(n)^2/n^τ
 = 1 + Σ_{n≥1} λ_X(n)^2/n^τ.
```

**Consequência:** `Σ λ_X(n)^2/n^τ = ζ(τ)S_(τ−1)(X) − 1`. O artigo CS usa apenas a desigualdade segura `Σλ²/n^τ≤ζ(τ)S`; a identificação anterior da P81 sem `−1` era uma **abreviação para um majorante**, e não igualdade literal para a série de λ. Os scripts P81 são conservadores, por isso **não precisam de correção numérica**.

O código `fx_audit.py` reproduz a versão truncada da convolução em números racionais exatos, inclusive o termo `1`, em dez casos separados.

## 2. Desigualdade de F_X: separar identidade do teorema importado

Para `σ>1`, os coeficientes reais `u_n=λ_X(n)n^(-σ)` são admissíveis no teorema de valor médio (KLN 2018, Lema 3.6, equação 3.14). Este teorema externo diz:

```
∫_0^T |Σ u_n n^(it)|² dt
    ≤ Σ u_n² [T + πm0(n+1)],
m0 = sqrt(1+(2/3)*sqrt(6/5)).
```

O grupo `n+1` deve ser preservado. Reorganizando seus termos, **sem usar aproximações**:

```
F_X(σ,0,T,1/2)
 ≤ πm0 Σ λ_X(n)²/n^(2σ−1)
   +(T+πm0) Σ λ_X(n)²/n^(2σ).
```

A passagem de `T+πm0(n+1)` à segunda linha é puramente algébrica, testada em seis exemplos racionais exatos. **Não** demonstramos de modo independente a desigualdade de valor médio KLN; ela é importada com hipótese de coeficientes reais e somabilidade.

## 3. Decomposição em X e 5X

Escreva `a=δ/log X>0`, `σ=1+a` e `q=2+2a`. Como `λ_X(n)=0` para `n≤X`:

```
Σ λ_X(n)²/n^(2+2a)
 = Σ_(X<n≤5X) λ_X(n)²/n^(2+2a)
   + Σ_(n>5X) λ_X(n)²/n^(2+2a).

X<n≤5X: n^(−2−2a)≤X^(−2a)n^(−2).
n>5X:   n^(−2−2a)≤(5X)^(−1−a)n^(−1−a).
```

São desigualdades de potências exatas, verificadas adicionalmente em seis exemplos racionais (para expoentes inteiros). Agora são importados, **distintamente**:

- **KLN 2018, Lema 3.4 (3.8)**: `Σ_(X<n≤5X) λ²/n²≤b3/X`, `b3=0.605`. A fonte KLN atribui essa cota a Ramaré, [24, Lema 5.6].
- **Ramaré–Zúñiga v1/v3, P81**: `Σ λ²/n^(1+a)≤b4/a` sob `X≥3×10^12`, importando as cotas públicas de `S_a(X)` e `ζ(1+a)`. Este limite inclui a margem `1+a` usada pela P80.

Com `a=δ/log X`, obtém-se

```
Σ λ²/n^(2+2a)
 ≤ 1/X * [ b3 exp(−2δ)
           + (b4*logX)/(5δ) exp(−δ) 5^(−a)].
```

**Importante:** o fator `5^(−a)≤1` controla o segundo termo. Não simplesmente remover `e^(γa)>1` de um limite superior, como na versão anterior do PDF.

## 4. Resultado para F_X e compatibilidade com CS v2

Usando o limite de `Σλ²/n^(1+2a)≤b4/(2a)` da P81 e a estimativa de cauda precedente:

```
F_X(1+a,0,T,1/2)
 ≤ πm0*b4*(logX)/(2δ)
  + (T+πm0)/X * [
       b3 exp(−2δ) + b4*(logX)*exp(−δ)/(5δ)].
```

Se `X≥X0≥3×10^12`, `log X≥log X0>0`, e `γ>0\), o majorante acima é no máximo

```
(log X) * [ C3(X0,δ) + (T+πm0)/X * C4(X0,δ) ],

C3 = π*m0*b4*exp(2γδ/log X0)/(2δ),
C4 = b4/(5δ expδ) + b3*exp(−2δ)/log X0.
```

Essa é precisamente a **forma estrutural** de Chourasiya–Simonič v2, Lema 6, mas com b4 recalculado por P78–P81 e **somente para a região pública X≥3×10^12**. É mais fraca ou igual ao majorante intermediário (não troca a tabela certificada P79).

Na prova do Corolário 1 de CS, `X=hT`, `h=1`, `T≥3×10^12`; assim a região usada nas faixas da tabela P79 está dentro da restrição válida. Entretanto, o Lema 6 **enunciado originalmente** cobre `X≥3×10^10`, que não foi reproduzido integralmente aqui.

## 5. O que permanece em aberto

**Três insumos externos explicitamente separados**:

1. KLN Lema 3.6, equação (3.14) — a desigualdade de valor médio para séries com coeficientes reais (derivada de Ramaré/Montgomery–Vaughan).
2. KLN Lema 3.4, equação (3.8) — a cota da parte curta `X<n≤5X` com `b3=0.605` (via Ramaré).
3. Ramaré–Zúñiga Alterman v1/v3 — o majorante da soma dupla de Möbius `S_ε(X)` dentro de seus domínios, importado na P81.

**Status:** a identidade e a passagem algébrica para F_X estão fechadas, sob os três teoremas externos. Não há prova independente desses três insumos nesta etapa, nem demonstração do Lema 6 em todo seu domínio original, nem prova da Hipótese de Riemann global.

## Reproduzir

```bash
python -m unittest -v NUMERIA_P82_FX_AUDIT.test_fx
python -O -m unittest -v NUMERIA_P82_FX_AUDIT.test_fx
python -m NUMERIA_P82_FX_AUDIT.fx_audit --out p82-fx-audit.json
```

Fontes primárias: [CS v2 Lemma 6](https://arxiv.org/html/2507.15184v2), [KLN (2018), Lemmas 3.4/3.6](https://www.cs.uleth.ca/~kadiri/articles/Explicit-zero-density-for-the-Riemann-zeta-function-JMAA-May2018.pdf), [RZ v3](https://arxiv.org/html/2603.25961v3).
