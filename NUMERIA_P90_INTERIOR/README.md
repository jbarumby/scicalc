# NUMERIA–RIEMANN P90 — Lema 6 afiado para pares com 0 ou 1 ponto intermediário

**PR base:** P89, que demonstrou as caudas exteriores do Lema 6 para **todos** os espaçamentos reais. Esta fase não refaz os certificados S2/Arb nem amplia artificialmente Suzuki.

**Resultado novo desta fase:** demonstração **analítica uniforme**, independente do Lema 5 ou do Lema 6 original de Preissmann, da desigualdade afiada
\[
\sum_{k\notin\{a,b\}}\frac{\delta_k}{(x_k-x_a)^2(x_k-x_b)^2}
\le\frac{\pi^2}{3h^2}
     \left(\frac1{\delta_a}+\frac1{\delta_b}\right)
 -\frac{3(\delta_a+\delta_b)}{h^4},
\quad h=x_b-x_a>0
\tag{L6}
\]
**para toda configuração finita de reais distintos, sempre que \(b-a\le2\)**. Os pontos fora do par podem ser arbitrariamente numerosos e irregulares. O caso \(b-a=1\) significa **nenhum ponto interior**; \(b-a=2\), **um ponto interior**.

**A P90 NÃO demonstra L6 quando \(b-a\ge3\)**, nem reconstitui integralmente o Lema 5 de 1984. A desigualdade geral continua **teorema publicado externo**, necessária para preservar \(m_0\) no argumento de \(F_X\).

## 1. Fontes matemáticas e os cinco uploads

- [E. Preissmann, *Sur une inégalité de Montgomery–Vaughan* (1984), Lemas 5–6, pp. 104–110](https://www.e-periodica.ch/cntmng?pid=ens-001:1984:30::29). A extração textual identifica **dois regimes variacionais** no Lema 5 (sequências de lacunas decrescentes ou com variações). As fórmulas escaneadas estão parcialmente ilegíveis e seu desenvolvimento completo NÃO é atribuído a esta fase.
- [W. Yangjit, *On the Montgomery–Vaughan weighted generalization of Hilbert's inequality*, arXiv:2203.14950v1, Lema 3.2, equação (3.6), p. 9](https://arxiv.org/pdf/2203.14950): enunciado exato de (L6), citado de Preissmann.
- P87 [PR #87](https://github.com/jbarumby/scicalc/pull/87): reconstrução independente da cota diagonal \(Q_i\le\pi^2/3\) via método de Shan.
- P89 [PR #89](https://github.com/jbarumby/scicalc/pull/89): cota uniforme das duas caudas exteriores com núcleo \(\varphi_h(t)=1/[t^2(t+h)^2]\).

**Novos arquivos P87 nesta rodada:** `p87_teorema_preissmann (1).tex` e `p87_verificacao (1).py` têm os **mesmos SHA-256** já auditados em P89: `d16771dab750ef4f00e00bf841ddcf5f179e412ef6392914ea819a0e7b9d706d`, `5089f7016827dffb5a6202cd783ca70ab099ba5a87db37e954e2192a1b8ff14a`. O `p87_execucao.log` (SHA `1b594dda439681ab78a286dca2aa810c3dab4b2149f9950389c2294a16`) e `p87_ledger.json` (SHA `f593dc80743a5510c064aaa6ad3d4b2149f9950389c2294e16e425121b728da7`) foram agora efetivamente examinados. Confirmam **24 gates e zero falhas reportadas**, com a dependência cruzada explicitamente importada. Contudo, o B5 desses arquivos usa **π superior** e margens finais `float`, e os blocos B3/B4/B6/B7 são diagnósticos; portanto o relatório antigo NÃO certifica por si só o teorema espectral geral. Ver [P89_UPLOAD_ERRATA](../NUMERIA_P89_LEMMA5/P87_UPLOAD_ERRATA.md).

## 2. Preparação: envelope externo e orçamento interior exato

Fixe \(a<b\) e \(\delta_i\) a distância entre o ponto \(x_i\) e seu vizinho mais próximo. Escale \(h=x_b-x_a\) para 1. Escreva \(s_a=\delta_a/h\) e \(s_b=\delta_b/h\), com \(0<s_a,s_b\le1\), pois \(\delta_a,\delta_b\le h\).

Defina
\[
Z=\zeta(2)=\pi^2/6,\quad
S(s)=\sum_{j=1}^\infty\frac1{j^2(1+js)^2},\qquad
D(s)=\frac{2Z-S(s)}s-3s.
\]
A P89 dá a cada cauda exterior o limite
\[
L_{\rm ext}\le h^{-3}\frac{S(s_a)}{s_a},
\qquad
R_{\rm ext}\le h^{-3}\frac{S(s_b)}{s_b}.
\]
O lado direito de (L6), multiplicado por \(h^3\), é
\[
2Z\left(\frac1{s_a}+\frac1{s_b}\right)-3(s_a+s_b).
\]
Escrevendo
\[
M=h^3\sum_{a<k<b}
 \frac{\delta_k}{(x_k-x_a)^2(x_k-x_b)^2}\ge0,
\]
temos um **critério suficiente universal**:
\[
\boxed{M\le D(s_a)+D(s_b)\quad\Longrightarrow\quad\text{(L6)}.}
\tag{P90-BUDGET}
\]
É uma redução rigorosa, sem supor o Lema 6. **Ainda não demonstramos P90-BUDGET para qualquer quantidade de pontos internos.** A P90 demonstra os casos de zero e um ponto interno.

### 2a. Valores exatos de S nos nós de interpolação

A função \(S\) é decrescente e **convexa** para \(s\ge0\) (cada parcela \((1+js)^{-2}/j^2\) é convexa, e a série converge uniformemente em compactos).
Por frações parciais e somas harmônicas, para qualquer inteiro \(m\ge1\),
\[
\boxed{S(1/m)=2Z-H_m^{(2)}-\frac2mH_m.}
\]
Em particular,
\[
S(0)=Z,\quad S(1/2)=2Z-\frac{11}{4},\quad S(1)=2Z-3.
\tag{P90-KNOTS}
\]
Todos os coeficientes de \(Z\) e os termos racionais são identidades exatas, não aproximações decimais.

### 2b. D(s) estritamente decrescente — comparação variacional uniforme

Escreva
\[
D(s)=\sum_{j\ge1}\frac1j\,q(js)-3s,
\qquad q(z)=\frac{2-(1+z)^{-2}}{z}.
\]
A derivada racional é
\[
q'(z)=-\frac{1+3z+6z^2+2z^3}{z^2(1+z)^3}<0\quad(z>0).
\]
Portanto **cada termo da soma é estritamente decrescente** com \(s\), e \(D(s)\) também. A comparação termo a termo evita qualquer troca injustificada de derivadas e somas infinitas.

## 3. Teorema P90-A — L6 afiado quando a e b são vizinhos

Se \(b=a+1\), não há termos interiores: \(M=0\). Basta demonstrar \(D(s)\ge0\) para \(0<s\le1\).

Por convexidade de \(S\) e pela corda entre os pontos 0 e 1:
\[
S(s)\le(1-s)S(0)+sS(1)=Z+s(Z-3).
\]
Consequentemente,
\[
\boxed{
D(s)\ge(1-s)\left(\frac Zs+3\right)\ge0.
}
\]
Assim \(0=M\le D(s_a)+D(s_b)\), e (L6) está demonstrada para **todos os N e todos os espaçamentos reais** quando \(b=a+1\). A igualdade-limite para redes uniformes infinitas é compatível com \(D(1)=0\).

## 4. Teorema P90-B — L6 afiado com exatamente um nó intermediário

Agora seja \(b=a+2\). Na escala \(h=1\), o único ponto interior está em \(t\in(0,1)\). Por reflexão podemos supor \(0<t\le1/2\). Suas distâncias aos extremos são \(t\) e \(1-t\). A distância ao vizinho mais próximo é exatamente \(\min(t,1-t)=t\) e, portanto,
\[
M=\frac{t}{t^2(1-t)^2}=\frac1{t(1-t)^2}.
\]
Como \(\delta_a/h\le t\), \(\delta_b/h\le1-t\) e \(D\) é decrescente, é suficiente mostrar
\[
D(t)+D(1-t)\ge\frac1{t(1-t)^2}.
\tag{P90-ONE}
\]

Pela convexidade e os três valores exatos (P90-KNOTS):
\[
\begin{aligned}
S(t)&\le(1-2t)Z+2t(2Z-11/4)
      =Z+2tZ-\frac{11t}{2},\\
S(1-t)&\le2t(2Z-11/4)+(1-2t)(2Z-3)
      =2Z-3+\frac t2.
\end{aligned}
\]
Isso implica
\[
D(t)+D(1-t)
\ge\frac{Z(1-2t)}t+\frac{11/2-3t}{1-t}.
\]
Multiplicando a diferença entre esse majorante inferior e \(M\) por \(t(1-t)^2>0\), obtemos
\[
Z(1-2t)(1-t)^2
 +t(1-t)\left(\frac{11}{2}-3t\right)-1.
\]
Como \(Z=\pi^2/6>3/2\), a expressão é no mínimo
\[
\frac32(1-2t)(1-t)^2
+t(1-t)\left(\frac{11}{2}-3t\right)-1
=\boxed{\frac{(1-2t)(1+t)}2\ge0}.
\]
Logo (P90-ONE) e (L6) valem, **independentemente da quantidade e distribuição dos pontos exteriores**, para \(\boxed{b=a+2}\). No caso \(t=1/2\) os limites tornam-se iguais ao orçamento da malha uniforme infinita, mostrando que a constante \(-3\) foi preservada exatamente.

**Este resultado é uma demonstração analítica para todos os números reais admissíveis**, não extrapolação de matrizes finitas ou amostras racionais.

## 5. Igualdade exata para qualquer malha uniforme, incluindo múltiplos pontos interiores

Na escala \(h=1\), tome \(m\) lacunas internas iguais a \(1/m\), com \(m\ge1\). A contribuição do trecho interior é
\[
M_{\mathrm{unif}}(m)=m^3\sum_{k=1}^{m-1}\frac1{k^2(m-k)^2}
=2mH_{m-1}^{(2)}+4H_{m-1}.
\]
Pela identidade (P90-KNOTS),
\[
D(1/m)=mH_m^{(2)}+2H_m-\frac3m
\]
e, portanto, **exatamente para todo \(m\)**,
\[
\boxed{M_{\mathrm{unif}}(m)=2D(1/m).}
\]
Isso recupera a rigidez da igualdade uniforme da P88 e evidencia que **não há folga na constante disponível** para uma prova de todos os espaçamentos.

O programa `nearby_sharp.py` confere a identidade por `fractions.Fraction` em \(m=1,\dots,40\). Essa regressão não substitui a identidade analítica geral.

## 6. Problema variacional residual e uma advertência concreta

Para \(m\ge3\) lacunas interiores, escreva
\[
g_1+\cdots+g_m=1,\quad g_j>0,\quad
t_k=\sum_{j=1}^{k}g_j.
\]
O núcleo interior normalizado é
\[
\mathcal M(g)=\sum_{k=1}^{m-1}
\frac{\min(g_k,g_{k+1})}{t_k^2(1-t_k)^2}.
\]
O **alvo suficiente, ainda aberto**, é
\[
\boxed{\mathcal M(g)\le D(g_1)+D(g_m)
\quad\text{para TODOS os }m\text{ e }g>0,\ \sum g=1.}
\tag{P91-GATE}
\]
Se fosse demonstrado, então como \(s_a\le g_1\), \(s_b\le g_m\) e \(D\) decresce, fecharia o Lema 6 afiado para **toda separação**, eliminando o elo externo da reconstrução. Não supomos (P91-GATE) como fato.

O cálculo variacional exige os diferentes regimes da função \(\min(g_k,g_{k+1})\). Mesmo com **um só ponto móvel** \(t\) e extremos fixos em 0 e 1, o núcleo interior é
\[
M(t)=
\begin{cases}
1/[t(1-t)^2],&t\le1/2,\\
1/[t^2(1-t)],&t\ge1/2.
\end{cases}
\]
Sua derivada é
\[
M'(t)=
\begin{cases}
(3t-1)/[t^2(1-t)^3],&t<1/2,\\
(3t-2)/[t^3(1-t)^2],&t>1/2.
\end{cases}
\]
Ela muda de sinal em \(t=1/3\) e \(t=2/3\), e apresenta cúspide em \(t=1/2\). Assim, a mera afirmação de que *regularizar espaçamentos sempre melhora cada termo* é **falsa**. O Lema 5 original distingue regimes variacionais por um motivo substantivo. `variational.py` guarda exemplos exatos dessas trocas.

**Próxima fase P91:** investigar (P91-GATE) no interior com \(m\ge3\), por divisões em cones \(g_k\le g_{k+1}\)/\(g_k\ge g_{k+1}\), transformações de adjacência e comparação convexa entre pontos de suporte \(s=1/m\). Isso pode exigir refazer parte das duas alternativas do Lema 5 de Preissmann; a prova geral continua aberta.

## 7. Computação, escopo e reprodutibilidade

```bash
python -m unittest -v NUMERIA_P90_INTERIOR.test_p90
python -O -m unittest -v NUMERIA_P90_INTERIOR.test_p90
python -m NUMERIA_P90_INTERIOR.nearby_sharp --out p90-nearby-sharp.json
```

A auditoria executa:
- comparações diretas do Lema 6 por frações exatas e **limite INFERIOR** rigoroso para π, nos pares com 0 ou 1 ponto interior;
- comparações de orçamento **estritamente positivas** em oito sequências interiores irregulares com \(m\ge3\), usando um **limite racional inferior** de \(D\) via caudas explícitas da série; isso é diagnóstico, NÃO demonstra P91-GATE;
- igualdades de malhas uniformes para todo \(m\ge1\), deduzidas de identidades de somas harmônicas;
- verificação de fatorações algébricas e sinais dos regimes de uma variável, sem `float`.

**Não comprovado na P90:** Lema 5 variacional completo de Preissmann, Lema 6 para qualquer distância de índices, desigualdade espectral de Preissmann \(\pi m_0\) sem o teorema publicado, extensões de Ramaré–Zúñiga para \(\,3\times10^{10}\le X<3\times10^{12}\), nem RH global.

Os PRs P84–P89 e a cobertura Suzuki \(18\le t\le79\) continuam intactos, finitos e condicionais. Não fazer merge automático.
