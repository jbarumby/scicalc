# NUMERIA–RIEMANN — Fase NL1/ES1 (2026-11-03)

**Resultado:** um certificado finito EXATO de sistema não linear **apenas para um modelo espectral 2×2**; três zeros conhecidos da zeta encontrados **exploratoriamente** com Newton 2D; confronto numérico com cotas registradas no **modelo não projetado V24**; análise de 44 referências deduplicadas.

**Não há demonstração da Hipótese de Riemann neste pacote.**

## 1. Origem e rastreabilidade

- Repositório consultado (via conexão GitHub): `jbarumby/scicalc`.
- Origem da fase V24: `audit/e8k-v24-global-arch-prime-pole-taper-2026-11-02` (commit `3ffe3799b238a02ad8ec5c0573a16a294affaa20`).
- Arquivos examinados: `v24/README_V24.md`, `v24/RESULTADO_V24_CONTINUIDADE.json` e `v23/README_V23.md`.
- `EXTRATO_V24_COTAS.json` contém **somente quatro linhas de cotas** transcritas da fonte acima; NÃO é o JSON V24 original, nem uma reprodução independente de suas demonstrações.
- Bibliografia: os CSVs fornecidos `References_LeapSpace (5).csv` (15 registros) e `(6).csv` (31 registros); dois registros em comum (IDs 314, 320); 44 referências únicas. A classificação decorre dos títulos e resumos disponibilizados, não de uma leitura integral dos 44 artigos.

## 2. Exemplo de sistema não linear e certificação exata

Escolheu-se a matriz de brinquedo **autoadjunta**

\[
H(a,b)=\begin{pmatrix}a&1\\1&b\end{pmatrix}.
\]

O sistema é

\[
F(a,b)=\begin{pmatrix}\det H(a,b)\\\operatorname{tr}H(a,b)-5/2\end{pmatrix}
=\begin{pmatrix}ab-1\\a+b-5/2\end{pmatrix}=0.
\]

Newton 2D com as sementes (1,8;0,7) e (0,7;1,8) convergiu para (2;1/2) e (1/2;2). Em torno de cada solução, usou-se a caixa com semilargura 1/1000 em cada eixo.

O operador de Krawczyk é

\[
K(X)=c-CF(c)+(I-CJ(X))(X-c),\qquad C=J(c)^{-1}.
\]

Todas as operações de intervalos no `numeria_nl1.py` usam **frações racionais exatas** da biblioteca padrão (`fractions.Fraction`), e não arredondamento usual de ponto flutuante. Em ambas as caixas, verificou-se `K(X)` estritamente contido no interior de `X`, e `||I-CJ(X)||_∞ ≤ 1/750 < 1`. Portanto há **exatamente uma solução por caixa**. A asserção é matemática rigorosa **nesse modelo algébrico**.

Isso **não** certifica nenhuma propriedade de uma matriz de Weil ou de um operador Hilbert–Pólya. Uma matriz de ordem 2 que sabemos construir não resolve o problema de identificar o operador global associado à zeta.

## 3. Zeros da zeta — descoberta exploratória

Foi resolvido o sistema real de duas componentes

\[
\big(\Re\zeta(\sigma+it),\Im\zeta(\sigma+it)\big)=(0,0),
\]

com derivada complexa da zeta via mpmath (`zeta(s, derivative=1)`), equivalendo ao Jacobiano real

\[
J(\sigma,t)=\begin{pmatrix}\Re\zeta'(s)&-\Im\zeta'(s)\\\Im\zeta'(s)&\Re\zeta'(s)\end{pmatrix}.
\]

As condições iniciais têm partes reais 0,62; 0,56; 0,44 (a reta crítica **não** foi imposta como restrição), respectivamente.

| semente (σ,t) | raiz obtida σ | raiz obtida t | |ζ(s)| aproximado | iterações |
|---|---:|---:|---:|---:|
| (0,62;14,10) | 0,5 | 14,13472514173469379046 | 3,08×10⁻⁸⁶ | 7 |
| (0,56;21,06) | 0,5 | 21,02203963877155499263 | 8,35×10⁻⁸⁶ | 7 |
| (0,44;25,05) | 0,5 | 25,01085758014568876321 | 2,11×10⁻⁸⁶ | 7 |

Os resíduos são **aproximações numéricas**, NÃO limites superiores certificados. Esses três exemplos são zeros conhecidos. Não foi implementado um operador intervalar complexo para a zeta, nem foi excluída a presença de zeros fora da reta crítica.

## 4. V24 — gate de perturbação condicional

Da V24, apenas para o **modelo não projetado de 10 modos**, foram recuperados os seguintes limites superiores de norma por soma de linha:

| δ | cota V24 sobre ||E||∞ |
|---:|---:|
| 10⁻³ | 95,495337151444... |
| 10⁻⁴ | 0,954953371514... |
| 10⁻⁵ | 0,009549533715... |
| 10⁻⁶ | 0,000095495337... |

Os quatro dados são compatíveis com um fator constante vezes δ²; a razão foi verificada por `Decimal` de alta precisão. Essa é uma **conferência de consistência de valores transcritos**, não a validação independente da fórmula V24.

Se a matriz exata **projetada** P tiver limite inferior de autovalor `μ > 0`, `||E||₂ ≤ ε` e `||B||₂ ≤ β`, então, depois de demonstrar que a projeção e os sinais são os mesmos,

\[
\lambda_{\min}(P+B^*EB)\;\geq\;\mu-\beta^2\varepsilon.
\]

O utilitário `perturbation_gate_exact` verifica a desigualdade racional `μ > β² ε`, mas **não fornece μ, ε ou β a partir da matriz de Weil real**. A margem de `0,001` utilizada em um exemplo no JSON foi escolhida artificialmente para demonstrar a operação; é expressamente marcada como NÃO verificada em V24.

**Lacuna decisiva:** cotejar as entradas E7/V16/V24, sinais dos termos de polos/primos, normalização de Fourier e a projeção `Bᵀ(G+P)B`. Sem isso não é lícito transferir os resultados da V24 aos pivôs LDL certificados de outra fase.

## 5. Referências incorporadas: decisões concretas

A planilha CSV `REFERENCIAS_CLASSIFICADAS.csv` contém todas as **44 referências únicas** e a origem de cada uma. Prioridades propostas:

1. **344 — Semenov & Semenova (2019), DOI `10.1007/s10559-019-00159-y`:** referência central sobre Krawczyk e zeros analíticos, com exemplo de função zeta. Para dar o próximo passo real, implementar extensão complexa com valores e derivadas *intervalares*, não só Newton mpmath.
2. **341 — G. A. Hiary (2011), DOI `10.4007/annals.2011.174.2.4`:** acelerar avaliações de zeta em alturas grandes. Sua complexidade não equivale, por si só, a um certificado.
3. **342 — Borwein, Bradley & Crandall (2000), DOI `10.1016/S0377-0427(00)00336-8`:** seleção e cotejo independente de algoritmos de avaliação da zeta.
4. **316 — Alain Connes (1999), DOI `10.1007/s000290050042`:** fórmula de traço e geometria não comutativa. Corrigir a grafia **A. Cormes** do CSV para **Alain Connes**; a hipótese espectral global ainda exige demonstração.
5. **314 — Gilles Lachaud (2003), DOI `10.1016/S0377-0427(03)00621-6`:** critérios espectrais/torais; há correspondências condicionais que não dispensam provar as hipóteses.
6. **330 — Xiao & Yin (2015):** Newton de alta ordem, útil depois que Newton clássico já produzir raízes estáveis; precisão e velocidade não são certificação.
7. **328 — Lazurchak, Makarov & Sytnyk (2008):** aproximações por ambos lados para equações de operadores não lineares; candidato quando se tiver operador definido e cotas apropriadas.
8. **345 — Dino Vlahek (2026):** localização multiescala Gauss–wavelet de zeros em funções oscilatórias; potencial descoberta, não um substituto para provar unicidade por intervalos.
9. **349 — Kendall Atkinson (1992):** métodos de Galerkin/colocação para equações integrais não lineares; possível ponte ao sistema finito, se a aproximação do operador for provada.
10. **332–340, 353–356:** métodos de evolução numérica de Dirac; não devem ser confundidos com a construção rigorosa de um operador autoadjunto cujo espectro corresponda exatamente aos zeros não triviais da zeta.

Referências públicas conferidas:
- https://doi.org/10.1007/s10559-019-00159-y
- https://annals.math.princeton.edu/2011/174-2/p04
- https://doi.org/10.1016/S0377-0427(00)00336-8
- https://doi.org/10.1007/s11075-026-02484-8
- https://python-flint.readthedocs.io/en/latest/acb.html

## 6. Matriz de evidências e próximos gates

| Obrigação | Resultado | Status |
|---|---|---|
| Newton 2D em sistema matricial de exemplo | duas soluções | EXECUTADO |
| Krawczyk exato (frações racionais) no exemplo | duas caixas únicas certificadas | CERTIFICADO NO MODELO |
| Newton 2D da função zeta, três zeros conhecidos | resíduos de ~10⁻⁸⁶ | NUMÉRICO, NÃO CERTIFICADO |
| Auditoria de quatro cotas V24 | relação quadrática consistente | VERIFICAÇÃO DE CONSISTÊNCIA |
| Testes unitários | 10/10 passaram | EXECUTADO |
| Intervalos rigorosos para ζ e ζ′ em retângulos complexos | ainda não implementado; depende de Arb/FLINT ou avaliador equivalente | PENDENTE |
| Identificação completa do operador de Connes/Lachaud/Sierra com matriz de Weil implementada | não demonstrada | PENDENTE |
| Transporte rigoroso V24 → Bᵀ(G+P)B e LDL | não demonstrado | PENDENTE |
| Positividade para todas as funções-teste admissíveis | não demonstrada | PENDENTE |

**Gate NL2 proposto**: (1) instalar e versionar `python-flint==0.9.0` em ambiente CI isolado; (2) especificar avaliação intervalar de ζ(s) e ζ′(s) numa caixa complexa admissível; (3) implementar um Krawczyk complexo com teste positivo e negativos; (4) cotejar independentemente zeros conhecidos com provas de contagem de zeros/argument principle quando aplicável; (5) só depois integrar com E7/V24 e a forma de Weil.

## 7. Reproduzir

Requisitos Python 3.10+ com `mpmath`, `pytest`:

```bash
python numeria_nl1.py --out RESULTADO_NL1_ES1.json --v24 EXTRATO_V24_COTAS.json
python -m pytest -q test_numeria_nl1.py
python classificar_referencias.py --files '/caminho/References_LeapSpace (5).csv' '/caminho/References_LeapSpace (6).csv' --out REFERENCIAS_CLASSIFICADAS.csv
```

A aritmética intervalar racional usa apenas a biblioteca padrão. Na sessão desta execução, `python-flint` **não estava instalado**; por isso **não** foi emitido certificado intervalar para a função zeta.

**Interpretação final:** os resultados demonstram que a arquitetura Newton + Krawczyk pode operar e ser auditada; o projeto precisa agora de uma ponte analítica, rigorosamente justificada, das matrizes-testes para os objetos genuínos da fórmula de Weil. Sem essa ponte, seria incorreto alegar avanço demonstrativo sobre RH.
