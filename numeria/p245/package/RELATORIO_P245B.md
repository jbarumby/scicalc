# NUMERIA–RIEMANN — ES8-P245 Gate B: reprodução intervalar independente de amostras GL60

**Data:** 18/11/2026 · **Configuração:** `L = 161/200`, `T♯ = 504`, GL60, 377 modos por paridade.  **Estado:** certificado intervalar selecionado; não é uma certificação da matriz integral completa, do critério global de Weil ou da RH.

## Conclusão reproduzível

- **60/60 raízes** de `P_60` isoladas por intervalos racionais disjuntos de largura `2⁻²²⁴`, com verificação exata das mudanças de sinal. Uma raiz por intervalo e 60 intervalos disjuntos para um polinômio de grau 60 significam que são todas as raízes; a fórmula clássica fornece os pesos GL60, com intervalos dirigidos. A soma dos pesos contém exatamente 2.
- **48/48 comparações de amostras:** 8 pares de índices `(0,0),(0,1),(1,1),(0,16),(16,16),(0,376),(1,376),(376,376)` × 3 painéis `{0,503,1007}` × 2 setores `{par,ímpar}`. Para cada entrada, o **intervalo Arb original do piloto P236 está inteiramente contido** no intervalo independente calculado, e o arquivo piloto possui SHA-256 fixo.
- Das 48, **30** usam Bessel calculada por recorrência ou série em aritmética de intervalos; **18** que envolvem o grau 752/753 usam um majorante de Bessel centrado em zero, mais largo em precisão absoluta. Todas incluem digama, logaritmo, seno/cosseno, raízes e pesos GL60 controlados independentemente. Para os 18 casos extremos foram ainda emitidos limites de magnitude mais fortes, por frações exatas.
- **12/12 testes** passaram com Python convencional e novamente **12/12** com `python -O`, incluindo rejeição de uma raiz adulterada e de arquivo P236 cujo hash foi alterado. Uma segunda execução do auditor gerou o mesmo JSON byte a byte.

Os hashes e a lista exata das caixas aparecem em `P245B_CERT.json`. O auditor foi implementado sem `python-flint`, Arb, MPFR nem SciPy. Na reprodução final, não chama sequer mpmath: lê as 60 propostas de intervalos racionais de `P245B_ROOT_BRACKETS.json` e **prova novamente** que cada uma isola uma raiz real.

## 1. Erros e arredondamento com inteiros exatos

O sistema de intervalos representa `[a,b]/2²⁵⁶`, com inteiros Python arbitrariamente grandes. Soma, produto, divisão e raízes quadradas sempre arredondam para fora. Um teste de inclusão `a ≤ x ≤ b` é decidido por produtos inteiros e `Fraction`; nenhum `float` determina um sinal ou uma inclusão. `mpmath` foi utilizado antes apenas para sugerir os 60 intervalos iniciais, posteriormente certificados por sinais exatos do polinômio de Legendre.

Para uma raiz `x` de `P_60`, o peso é calculado a partir da identidade exata

\[
 w(x)=\frac{2(1-x^2)}{60^2 P_{59}(x)^2},
\]

avaliando `P_59` por recorrência intervalar. A regra usada é a mesma do P237,

\[
C_{ij}^{(p)}=\sum_{k=1}^{60}\frac{w_k}{4\pi}\bigl(\Psi_L(t_{pk})-\beta_{504}\bigr)U_i(t_{pk})U_j(t_{pk}),\qquad
 t_{pk}=\frac{2p+1+x_k}{4}.
\]

A constante `π` usa o intervalo de Machin `16 arctan(1/5)−4 arctan(1/239)` com restos alternados; os logaritmos usam a série racional de `atanh` após redução binária, e seno/cosseno, Taylor com argumento reduzido por um múltiplo **certificadamente incluído** de `2π`. Os restos algébricos de todas essas séries são comparados com um ulp fixo `2⁻²⁵⁶` e incluídos no resultado.

## 2. Digama complexa por outra representação

Não é utilizada a função `acb.digamma`. Parte-se da identidade clássica de Binet (DLMF 5.9.15), válida para `Re z > 0`:

\[
 \psi(z)=\log z-\frac1{2z}-2\int_0^\infty\frac{u\,du}{(u^2+z^2)(e^{2\pi u}-1)}.
\]

Primeiro desloca-se o argumento `1/4+it/2` para `z=512+1/4+it/2` pela recorrência da digama; os **512 termos** são computados como frações de expressões intervalares. O termo de Binet é expandido em 30 termos `-B_{2k} z^{-2k}/(2k)` com números de Bernoulli calculados como frações exatas.

Para `0≤t≤504`, `a=512+1/4` e `|Im z|≤252`, o resto depois de `K=30` termos satisfaz o limite absoluto

\[
 |R_K|\le
 \frac{4(2K+1)!}{6^{2K+2}\,512^{2K}\,((512+1/4)^2-252^2)}
 <2^{-256}.
\]

Essa cota segue diretamente da expansão finita do denominador na integral de Binet, de `π>3` e de `ζ(s)<2` para `s≥2`. **Cuidado corrigido durante a auditoria:** coeficientes de Bernoulli numericamente iguais podem aparecer em graus distintos; é obrigatório enumerá-los pelo índice, não buscá-los pelo valor. A correção eliminou uma diferença prévia de aproximadamente `10⁻²⁵` no teste-piloto.

## 3. Bessel e graus extremos

As funções `j_0`, `j_1`, `j_2` e `j_3` são avaliadas a partir de seno/cosseno. Os graus 32 e 33 usam recorrência de três termos nos painéis 503 e 1007 e série de potências alternada com resto explícito no painel 0. Para os graus 752 e 753, a desigualdade clássica

\[
 |j_n(x)|\le \frac{x^n}{(2n+1)!!}\quad (x\ge0),
\]

implica, no painel `p`, com `x ≤ 161(p+1)/400`, um limite racional `|U_{376}| < 50x_{\max}^n/(2n+1)!!`. Essa cota é válida para os testes de magnitude, mas **não resolve o sinal de cada Bessel de alto grau**.

| Painel | Setor | `|C_(0,376)|` e `|C_(1,376)|` estritamente menores que | `|C_(376,376)|` estritamente menor que |
|---:|---|---:|---:|
| 0 | Par | `10^-2360` | `10^-4722` |
| 0 | Ímpar | `10^-2364` | `10^-4730` |
| 503 | Par | `10^-328` | `10^-658` |
| 503 | Ímpar | `10^-329` | `10^-660` |
| 1007 | Par | `10^-102` | `10^-206` |
| 1007 | Ímpar | `10^-103` | `10^-207` |

São **cotas para as somas GL60 nos três painéis selecionados**, não para a matriz inteira nem para o erro de discretização da integral. Elas usam avaliações intervalares de `Ψ_L−β` nos 60 nós e os majorantes positivos de Bessel, sem reusar os valores Arb dos índices altos para calculá-las.

## 4. Escopo e próximos trabalhos

**Verificado neste gate:** cálculo intervalar independente das 48 entradas de amostra do **funcional quadraturado GL60**, e inclusão dos intervalos Arb publicados dentro das novas caixas; raiz/peso/digama/trigonometria/Bessel com limites explícitos; 18 cotas adicionais de magnitude por frações exatas.

**Não verificado:** a fórmula geral de Weil–Zhu como identidade de operadores; o limite analítico GL60 da diferença *integral verdadeira − soma quadraturada*; as outras entradas de cada painel; a auditoria dos outros 1.005 painéis por setor; a matriz P237 inteira; todo o espectro e a Hipótese de Riemann. As identidades clássicas de Binet, Legendre, Bessel, Machin, os momentos de Bernoulli e a regra GL60 são entradas matemáticas declaradas, não teoremas formalizados em Lean.

**Próximo Gate C:** reutilizar essa biblioteca intervalar para um subconjunto estratificado dos painéis 1–1006, aumentar a precisão adaptativa nos graus extremos e derivar/validar separadamente o resto analítico global GL60 antes de falar em prova integral independente.

## Reprodução

```sh
cd NUMERIA_ES8P245_GATE_B_INTERVALAR
python3 p245b_certified_gl60.py --pilot-dir pilot_sources --roots P245B_ROOT_BRACKETS.json --out /tmp/P245B_REPLAY.json
python3 -m unittest -v test_p245b_certified_gl60
python3 -O -m unittest -q test_p245b_certified_gl60
sha256sum /tmp/P245B_REPLAY.json P245B_CERT.json
```

**Fontes matemáticas:** [DLMF §5.9.15, Binet](https://dlmf.nist.gov/5.9.E15); [DLMF §10.54, representação integral da Bessel esférica](https://dlmf.nist.gov/10.54); [DLMF §10.53, série de potência da Bessel esférica](https://dlmf.nist.gov/10.53).
