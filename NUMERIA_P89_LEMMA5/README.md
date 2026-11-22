# NUMERIA–RIEMANN P89 — Preissmann Lema 5: controle uniforme das caudas exteriores

**Base:** PR #88, que reduziu o Lema 6 à positividade do compensador misto `E_ab`. A P89 não reexecuta S3, DOC1/DOC2, a tabela S2 nem a janela finita de Suzuki. A nota externa P87 fornecida pelo usuário foi auditada separadamente em [P87_UPLOAD_ERRATA.md](P87_UPLOAD_ERRATA.md).

## 1. Novo teorema para TODOS os espaçamentos reais: caudas exteriores

Considere uma configuração FINITA arbitrária de pontos reais distintos, ordenada `x_1<...<x_N`, com \(\delta_i=\min_{j\ne i}|x_i-x_j|\). Fixe índices `a<b` e `h=x_b-x_a>0`. Defina

\[
L_{a,b}=\sum_{k<a}\frac{\delta_k}{(x_a-x_k)^2(x_b-x_k)^2},\quad
R_{a,b}=\sum_{k>b}\frac{\delta_k}{(x_k-x_a)^2(x_k-x_b)^2}.
\]

**Teorema P89 (caudas exteriores):** para TODO `N`, TODO espaçamento real distinto e qualquer par `a<b`,

\[
\boxed{L_{a,b}\le \delta_a\sum_{j=1}^\infty
\frac1{(j\delta_a)^2(h+j\delta_a)^2},\quad
R_{a,b}\le \delta_b\sum_{j=1}^\infty
\frac1{(j\delta_b)^2(h+j\delta_b)^2}.}
\tag{P89-T}
\]

É um limite **independente do Lema 5 de Preissmann** (mas usa o **teorema geral de Shan**, já reconstruído na P87). O kernel pode ser escrito \(\varphi_h(t)=1/[t^2(t+h)^2]\), `t>0`. Como \(\varphi_h(t)=t^{-2}(t+h)^{-2}\), é positivo, estritamente decrescente e estritamente convexo. Mais formalmente, as derivadas primeiras de \(\log\varphi\) são negativas e a segunda derivada de \(\varphi\) é estritamente positiva, pela regra do produto de duas funções positivas com primeiras derivadas negativas e segundas derivadas positivas.

### Demonstração analítica completa de (P89-T)

Olhe primeiro para os pontos à **esquerda** de `x_a`. Escreva suas distâncias à âncora em ordem crescente `0<y_1<...<y_M`, com `y_n=x_a-x_{a-n}`. Defina os intervalos `g_1=y_1`, `g_n=y_n-y_{n-1}` para `n>1`. Temos \(\delta_a\le g_1\) e, para cada ponto à esquerda, sua folga mínima global é exatamente \(\min(g_n,g_{n+1})\), com a convenção `g_(M+1)=g_M` no extremo mais à esquerda.

**Aplique o teorema unidirecional Shan/P87** à sequência positiva `a_n=g_n/g_1`, com `a_1=1`, e à função
\[
f(t)=g_1^2\varphi_h(g_1t).
\]
Ela satisfaz todas as hipóteses de Shan: `f>=0`, decrescente, convexa e somável nos inteiros. O teorema diz
\[
\sum_{n=1}^{M}\min(a_n,a_{n+1})f\left(\sum_{j=1}^n a_j\right)
\le\sum_{j=1}^\infty f(j).
\]
Usando \(\min(a_n,a_{n+1})=\delta_{a-n}/g_1\) e \(y_n=g_1\sum_{j\le n}a_j\), isso é
\[
g_1 L_{a,b}\le g_1^2\sum_{j=1}^\infty\varphi_h(jg_1),
\quad\text{logo}\quad
L_{a,b}\le g_1\sum_{j=1}^\infty\varphi_h(jg_1).
\]
Para todo `j>=1`, a função
\[
g\mapsto g\varphi_h(jg)
=\frac1{g\,j^2(jg+h)^2}
\]
é estritamente decrescente. Como \(\delta_a\le g_1\), resulta
\[
L_{a,b}\le \delta_a\sum_{j=1}^\infty\varphi_h(j\delta_a).
\]
A mesma demonstração, por reflexão ao redor de `x_b`, fornece a desigualdade para `R_{a,b}`.

**A prova funciona para todo conjunto finito de reais distintos.** Não é generalização a partir de exemplos.

### Cota racional finita de cada cauda

Para `K>=1` inteiro, \(\varphi_h(t)\) é decrescente, então
\[
\delta\sum_{j>K}\varphi_h(j\delta)
\le\int_{K\delta}^\infty \varphi_h(t)\,dt
\le\int_{K\delta}^\infty t^{-4}\,dt
=\frac1{3(K\delta)^3}.
\]
Assim, cada limite da P89-T possui **certificado de cota superior racional**:
\[
\boxed{\quad
S_{\rm ext}\le
\delta\sum_{j=1}^K\frac1{(j\delta)^2(h+j\delta)^2}
+\frac1{3(K\delta)^3}\quad}
\]
para \(\delta,h\in\mathbb Q_{>0}\), calculado por `fractions.Fraction` em `exterior.py`. Nos testes tomamos `K=48`, mas o teorema vale para qualquer `K`.

## 2. Separação rigorosa: exterior positivo, interior negativo

Com `a<b` e `h>0`, o termo cruzado que apareceu na P88 é
\[
C_{a,b}=\sum_{k\ne a,b}\frac{\delta_k}
{(x_k-x_a)(x_k-x_b)}.
\]
Separe em `C_left + C_mid + C_right`. Termo a termo,
\[
C_\mathrm{left}\ge0,\qquad
C_\mathrm{right}\ge0,\qquad
C_\mathrm{mid}\le0,
\]
pois os dois fatores do denominador têm sinais iguais **fora** de `(x_a,x_b)` e sinais opostos **dentro**.

A P89 registra os três termos exatamente para todos os pares diagnosticados. Isso esclarece o obstáculo de `E_ab>=0` da P88: controlar as **duas caudas exteriores** não resolve a possível contribuição negativa do trecho **entre a e b**.

A soma de L6 separa-se de modo exato:
\[
S_{a,b}=
\underbrace{\sum_{k<a} \cdots}_{L_{a,b}}
+\underbrace{\sum_{a<k<b}\cdots}_{M_{a,b}}
+\underbrace{\sum_{k>b}\cdots}_{R_{a,b}}.
\]
As duas parcelas exteriores agora possuem majorações universais auditadas sem o Lema 5. O desafio é a **compensação afiada do termo central**, junto aos déficits diagonais \(\Delta_a,\Delta_b\).

## 3. Relação com Preissmann, Lema 5

A fonte original [Preissmann, *Sur une inégalité de Montgomery–Vaughan*, Enseign. Math. 30 (1984), pp. 104–110](https://doi.org/10.5169/seals-53823) usa o **Lema 5** para uniformizar parte dos espaçamentos exteriores e comparar somas de recíprocos na demonstração do Lema 6. A digitalização textual apresenta índices e coeficientes degradados: não tomamos uma conjectura tipográfica como se fosse um enunciado rederivado.

O resultado (P89-T) resolve **uma parte exterior uniforme** pela rota já fundamentada em Shan/P87. Não afirma reproduzir as duas alternativas variacionais, a dependência dos extremos ou os termos de fronteira exatos do Lema 5 na formulação de 1984. A hipótese `full_Lemma5_variational_proof_reproduced=false` é intencional e CI a verifica.

## 4. Documentos recebidos: auditoria e correção prioritária

O PDF/TeX P87 do usuário apresenta uma rota alternativa por Bellman. Ela contém um **erro efetivo no Lema 5.3** (falta `1/t` na mudança de variável da trigama e há uma desigualdade exponencial invertida), apesar de sua monotonicidade final ser verdadeira. A correção completa está em [P87_UPLOAD_ERRATA.md](P87_UPLOAD_ERRATA.md).

O arquivo de CI e `p87_verificacao.py` fazem comparação com `pi_upper²/3` no bloco B5: para provar `Q_i<pi²/3` por certificados racionais é necessário `pi_lower²/3`. As margens devem ser persistidas em texto racional, não reduzidas a `float`.

Os blocos de cálculos espectrais `numpy` são diagnósticos; a nota reconhece que a hipótese cruzada `(H)` permanece importada. A P89 não reexecuta indevidamente todos esses testes nem os classifica como provas globais.

## 5. Procedimento reproduzível

```bash
python -m unittest -v NUMERIA_P89_LEMMA5.test_p89
python -O -m unittest -v NUMERIA_P89_LEMMA5.test_p89
python -m NUMERIA_P89_LEMMA5.exterior --out p89-exterior.json
```

A CI publica JSON exato e SHA-256 das fontes utilizadas. Os exemplos com pontos racionais são **regressão**, não a prova para todos os reais, que é o argumento analítico acima.

## 6. Próximo alvo específico: P90

Desenvolver uma desigualdade **mista/interior** para `a<k<b` capaz de compensar as parcelas negativas de `C_mid`, preservando a **constante afiada** e o termo \(-3(\delta_a+\delta_b)/h^4\) do Lema 6. Reconstituir cuidadosamente do Lema 5 original a passagem com duas alternativas de espaçamentos, evitando dados obtidos apenas por OCR. Até lá, o Lema 6 geral e a constante espectral \(\pi m_0\) permanecem **teoremas externos publicados**, não resultados independentes completos.

Outros limites abertos na cadeia NUMERIA–RIEMANN: estimativas Ramaré–Zúñiga para `S_epsilon(X)` e o intervalo `3e10 <= X < 3e12`. A janela Suzuki até `t=79` é **finita e condicional**; não há prova global de RH.
