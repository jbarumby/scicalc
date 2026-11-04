# NUMERIA–RIEMANN — ES8-P6: auditoria independente das cotas analíticas

**Data:** 03/11/2026. **Escopo:** os 30 pares de mesma paridade da Fase 4, `L=9/20`, `R=160`, 1.600 células de largura `h=1/10`, 16 nós por célula, raio de Cauchy `ρ=3/10`. Não é prova da Hipótese de Riemann.

## 1. Fontes e níveis de evidência

- Fonte primária: `NUMERIA_FASE4_30_ENTRADAS_AUDITORIA_PROJECAO_2026-10-24.zip` extraída de `E8K_V16_GITHUB_COMPLETO(1).zip`, obtido do usuário. O programa P6 lê o arquivo original, confere todos os quatro hashes dos scripts e todos os 30 hashes dos certificados e reproduz a expressão fechada do erro. O ZIP de origem **não é** incorporado novamente ao repositório.
- E7 original: `NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3(1).zip`, SHA-256 `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7`. O código compara, sob aritmética racional, as cotas mais fortes para `G11` e `G22` com `RESULTADO_INTEGRAIS_E7.json` dentro do pacote.
- `source_snapshot/` preserva, por hashes SHA-256, os dois implementadores históricos dos valores nodais: uma via inteiros fixos e outra via `fractions.Fraction`. Eles são usados em **nove** controles de interseção de intervalos; não representam 768 mil reavaliações integrais (30 × 1.600 × 16).

## 2. Majorante de Cauchy: demonstração analítica independente

Tomemos `d=ρ+h/2=7/20`, `α=1/4-d/2=3/40` e `1-2Ld=137/200>0`. Para o modo real compactamente suportado

\[
\phi_k(x)=\frac{1}{\sqrt L}\sin\frac{k\pi(x+L)}{2L}\;\mathbf1_{[-L,L]}(x), \qquad
F_k(z)=\int_{-L}^{L}\phi_k(x)e^{izx}\,dx,
\]

`F_k(z)` é inteira. Para a função de entrada real, a extensão *holomorfa* usada na fórmula de Cauchy é

\[
H_{ij}(z)=\frac{F_i(z)F_j(-z)+F_j(z)F_i(-z)}2,
\quad
\Omega_\mathbb C(z)=\frac{\psi(1/4+iz/2)+\psi(1/4-iz/2)}2-\log\pi,
\quad q_{ij}(z)=H_{ij}(z)\Omega_\mathbb C(z)/\pi.
\]

A igualdade de `H` com `Re(F_i \overline{F_j})` vale **na reta real**; a operação `Re` não deve ser aplicada a funções de `z` complexo ao empregar Cauchy. Os polos mais próximos da extensão `Ω` estão em `z=±i/2` e o disco completo `|z-c|≤d=0,35` está livre deles. Para `|Im z|≤d`,

\[
|F_k(z)|\le2\sqrt L\,e^{Ld}\le\frac{2\sqrt L}{1-Ld},\qquad
|H_{ij}(z)|\le\frac{4L}{(1-Ld)^2}\le\frac{4L}{1-2Ld}.
\]

Da série clássica de frações parciais de `ψ`, para `Re w≥α>0`,

\[
\psi(w)=-\gamma-\frac1w+\sum_{n\ge1}\frac{w}{n(n+w)},\qquad
|\psi(w)|\le1+\frac1\alpha+2|w|,
\]

pois `0<γ<1`, `|n+w|≥n`, `∑_{n≥1} n^{-2}<2`. Com `w=1/4±iz/2` e `|z|≤|c|+d`, mais `log π<2`,

\[
|\Omega_\mathbb C(z)|\le\frac72+\frac1\alpha+|c|+d
\le 3+\left(1+\frac1\alpha\right)\left(1+\frac{|c|+d}{2}\right).
\]

A última desigualdade é exata pois a diferença entre os lados é `1/2+(1/α−1)(|c|+d)/2≥0`. Como `π>3`, conclui-se a cota da Fase 4 sem invocar seu gerador original:

\[
\boxed{M(c)=\frac{4L}{3(1-2Ld)}\left[3+\left(1+\frac1\alpha\right)\left(1+\frac{|c|+d}{2}\right)\right].}
\]

Para os centros não negativos da integração presente, ela se reduz à **função racional afim**

\[
M(c)=\frac{2381}{137}+\frac{860}{137}\,c.
\]

As identidades da digama e os lemas de Cauchy usados aqui são resultados matemáticos clássicos, não teoremas formalizados em Lean neste pacote. Fonte primária: [DLMF §5.7.6](https://dlmf.nist.gov/5.7.E6).

## 3. Erro de interpolação: expressão fechada

Para os nós interpolatórios equiespaçados `u_j=j/15`, os pesos são as integrais exatas dos polinômios de Lagrange, calculadas independentemente por produtos de polinômios. O novo teste verifica **16/16** condições de exatidão dos momentos, até o grau 15.

Para um ponto real da célula, o resto da interpolação (teorema de Rolle) e a desigualdade de Cauchy implicam

\[
|f(x)-p_{15}(x)|\le M(c)\left(\frac h\rho\right)^{16}
\prod_{j=0}^{15}\max\{j/15,1-j/15\},
\]

com o produto exato

\[
P_{16}=\frac{16416784384}{1601806640625}.
\]

Somando as integrais por célula (`c_k=(k+1/2)h`), utiliza-se **algebricamente**, sem loop original,

\[
\sum_{k=0}^{1599} c_k =\frac{h\,1600^2}{2}=128000.
\]

Logo a estimativa total é

\[
\begin{aligned}
E&=hP_{16}(h/\rho)^{16}\sum_{k=0}^{1599} M(c_k)\\
&=\boxed{\frac{4154891126177792}{209922127267236328125}}
\approx 1,979253535711606\times10^{-5}.
\end{aligned}
\]

O código demonstra por frações que essa fórmula fechada é igual à soma direta das 1.600 células e exatamente igual aos **30 erros históricos**. Também verifica os limites externos registrados para as 30 integrais **sem recalcular todos os valores nodais**.

## 4. Cauda arquimediana — Proposições 5 e 6

Para `r≥1`, a fórmula de duplicação e reflexão da digama fornece

\[
\Omega(r)=-(\gamma+3\log2+\log\pi)+S(r)-\frac{\pi}{2\cosh(\pi r)},
\qquad S(r)=\sum_{n≥0}\frac{r^2}{(n+1/2)((n+1/2)^2+r^2)}.
\]

Por comparação integral de funções positivas decrescentes, `0≤S(r)−1/2 log(1+4r²)≤2`. A cota simples `0<γ+3ln2+lnπ<6` e `π/(2cosh πr)<2` prova conservadoramente `|Ω(r)|≤8+ln(3r)`. Esse é o **majorante da Proposição 5**, suficiente para a auditoria antiga.

Reforço independente: com a expansão de Euler–Maclaurin de `γ`, as séries racionais de `atanh`, limites para `π` pela fórmula de Machin e `cosh(π)>∑_{k=0}^8 (157/50)^{2k}/(2k)!`, verificamos *apenas por desigualdades racionais*:

\[
3,8013<\gamma+3\log2+\log\pi<3,8014,\quad
0,8047<\tfrac12\log5<0,8048,\quad
\frac{\pi}{2\cosh\pi}<0,14.
\]

Deduções: `Ω(r)>−3,14` e `Ω(r)<ln r−0,99` para `r≥1`; portanto `|Ω(r)|≤ln r+3,14`. Pela fórmula fechada dos modos e `R=160≥2w_{max}`:

\[
|A_i(r)|\le\frac{2w_i}{\sqrt L\,r^2(1-(w_i/R)^2)},\quad
\left|\frac1\pi\int_R^\infty A_i(r)A_j(r)\Omega(r)\,dr\right|
\le \frac{4w_iw_j}{\pi L(1-\varepsilon_i)(1-\varepsilon_j)}
\frac{\log R+3,14+1/3}{3R^3},
\quad\varepsilon_i=(w_i/R)^2.
\]

Essa é a **Proposição 6 mais forte**, agora rederivada e avaliada com extremos racionais externos. A cota estrita é menor que a conservadora nos 30 pares. A comparação com dois valores efetivos do E7 confirmou:

| Entrada | Cota E7 histórica | Cota nova independente |
|:--|--:|--:|
| G11 | 0,0000240068359093 | 0,0000240068359093 |
| G22 | 0,0000963022970370 | 0,0000963022970370 |

A diferença racional absoluta entre cada par é inferior a `10^-14`. As identidades de reflexão/duplicação de `ψ` usadas nesta seção constam da [DLMF §5.5](https://dlmf.nist.gov/5.5).

## 5. Auditoria nodal delimitada

Duas implementações históricas cujos hashes foram autenticados foram executadas com argumentos de teste e comparadas por `Fraction` exato:

- cinco avaliações de `Ω` em `r∈{0,1,4,10,160}`;
- quatro avaliações de amplitudes em pares `(k,r)={(1,0),(1,4),(2,4),(3,10)}`.

Resultado: **9/9 interseções**. A checagem não pressupõe confiança em `mpmath` e não arredonda a decisão para `float`. Porém as duas implementações partilham a origem de algumas identidades analíticas; essa concordância **não certifica exaustivamente cada nó de cada uma das 30 integrais**. O código original utiliza rotinas de seno, sinc, logaritmo, constante de Euler e digama cujas provas de resto deverão ser verificadas linha a linha ou formalizadas antes de um certificado inteiramente independente.

## 6. Resultado e limitações

| Gate | Status |
|:--|:--|
| Majorante holomorfo usado no erro da Fase 4 | **REDERIVADO** a partir da série clássica da digama |
| Pesos de interpolação e soma exata do erro | **VERIFICADOS** em aritmética racional independente |
| Cotas de cauda conservadora e forte | **REDERIVADAS** com desigualdades racionais e identidades clássicas |
| Autenticidade e coerência dos 30 certificados históricos | **VERIFICADAS** por hashes e metadados |
| Reavaliação de todos os pontos nodais | **NÃO EXECUTADA** (9 controles por intervalos) |
| Reprovação/refutação de RH ou positividade universal de Weil | **NÃO DEMONSTRADA** |

Esta etapa fortalece diretamente duas premissas da Fase 4 e as caudas E7, mas não é um atestado completo sobre os valores nodais. Um ciclo no GitHub Actions reproduz os 22 testes e os dados previamente autenticados; não recebe novamente os ZIPs binários históricos e, portanto, não deve alegar que reabriu os arquivos originais em execução remota.

**Próxima obrigação ES8-P7:** auditoria exaustiva da inclusão intervalar nodal com divisão independente das 30×1.600×16 avaliações, além da verificação global das premissas da fórmula explícita para a classe suavizada.
