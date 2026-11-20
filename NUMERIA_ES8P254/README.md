# NUMERIA–RIEMANN — ES8-P254 Gate A: piso espectral e Schur por blocos

**Base revisada:** `jbarumby/scicalc`, P251 B/C, P252 espectral, P253 limiar de 7. Esta etapa **não** altera certificados históricos e **não** certifica uma nova janela.

## Achado e correção da estratégia

A P253 documenta testemunhas negativas para o **modelo auxiliar truncado de quatro potências primas** perto de `L = log(7)/2`, em cortes 309, 350 e 504; essas direções não são testemunhas negativas da forma completa de Weil. A P253 também registra, com auditoria Arb, que no **modelo de quatro potências** o residual satisfaz `Psi_4(t) - 1/4 > 0` para `|t| >= 504`. Não é lícito usar esse piso em `|t| >= 309` com base no certificado existente.

Na convenção `F(t)=∫f(x)e^(itx)dx`, Plancherel implica

```
G_T[i,j] = (1/(2π)) ∫_{|t|<=T} F_i(t) overline(F_j(t)) dt,
0 <= G_T <= I,
R_{β'}(f)-R_β(f) = (β'-β) <f,(I-G_T)f>.
```

Essas fórmulas usam a **mesma base ortonormal**, a mesma normalização angular, o mesmo cutoff T, polos idênticos e os mesmos quatro termos primos nos dois modelos. Para `T=504`, `β'=1/4`, `β` do P250, e o certificado **separado** de sinal para o espectro complementar, o novo modelo é limitado inferiormente pela forma de Weil e o acréscimo ao modelo velho é semidefinido positivo. O novo bloco finito se monta por

```
A_new = A_old + (1/4-beta)*(I-G_504).
```

**A positividade do acréscimo não demonstra, isoladamente, que o novo bloco seja positivo.** A P253 continua exigindo um novo certificado de margem, incluindo a quadratura para G, seus erros analíticos e uma fatoração LDL/Schur rigorosa. A presença de autovetores de Rayleigh negativos no modelo antigo não decide o novo.

## Critério matricial estrito

Se um operador hermitiano tem blocos `[[A,B],[B*,D]]`, dados limites **rigorosos** `A >= a I`, `D >= d I` e `||B||_2^2 <= b²`, então

```
A_total >= m I  se  a>m, d>m e (a-m)*(d-m)>b².
```

A prova vem de `2|<u,Bv>| <= 2b||u||||v||`, mais o critério de Sylvester da matriz escalar `[[a-m,-b],[-b,d-m]]`. Para `B` com intervalos de entradas, a norma de Frobenius é um majorante seguro: `b² <= Σ_ij max(|B_ij^-|,|B_ij^+|)^2`. O código usa **frações exatas**, nunca autovalores binário64 para tomar decisões. Este critério pode ser mais útil do que uma cota global de Weyl se as bases forem separadas em bloco de baixa frequência, modos de fronteira e cauda, mas requer margens `a,d` verificadas.

## Pacote e execução

`floor_schur.py` propaga intervalos exatos da correção de piso, verifica entradas simétricas, oferece um teste Gershgorin conservador e o teste escalar de Schur; falha se os dados não fecham. `test_floor_schur.py` contém testes de casos adversariais, intervalos e uma matriz de exemplo em que Weyl não fecha, mas Schur fecha. A matriz de exemplo **não** é um resultado NUMERIA.

```
python -m unittest discover -s NUMERIA_ES8P254_GATE_A -p 'test_*.py' -v
python -O -m unittest discover -s NUMERIA_ES8P254_GATE_A -p 'test_*.py' -v
```

Sem `python-flint` disponível neste ambiente, **não reexecutei Arb/MPFR dos blocos da P253** e não inferi nova margem real. O próximo gate precisa:

1. Regerar intervalarmente `G_504` na base de Legendre com a normalização Fourier correta e um resto GL60 próprio, não reaproveitar sem prova o resto de `Ψ u_i u_j`.
2. Montar o novo bloco `A_new` com intervalos dirigidos e calcular por LDL ou Schur as margens dos blocos.
3. Aplicar ao operador contínuo a reserva de alta frequência e o termo de 7 com a estimativa de faixa estreita da P253.
4. Se um gate falhar, fornecer o índice, a desigualdade e o intervalo inconclusivo, sem inferir indefinição da forma de Weil.

## Springer — triagem dos três periódicos indicados (novembro de 2026)

- **Javed e Maji**, *Toeplitz algebra and symbol map via Berezin transform on H²(Dⁿ)*, *Advances in Operator Theory*, 18/11/2026, DOI [10.1007/s43036-026-00546-1](https://doi.org/10.1007/s43036-026-00546-1). Estrutura de operadores Toeplitz e mapa do símbolo em Hardy; **não** fornece automaticamente uma desigualdade para a compressão por células NUMERIA, que é de outra classe.
- **Wang**, *Self-adjoint extensions with compact resolvent*, mesmo periódico, 03/11/2026, DOI [10.1007/s43036-026-00545-2](https://doi.org/10.1007/s43036-026-00545-2). Relevante se futuramente construir um operador autoadjunto explícito; sem aplicação direta aos certificados atuais.
- **Alrimawi, Moradi e Sababheh**, *Optimized numerical radius bounds*, *Annali dell'Università di Ferrara*, 16/10/2026, DOI [10.1007/s11565-026-00771-9](https://doi.org/10.1007/s11565-026-00771-9). O resumo trabalha matrizes operatórias 2×2 com termos mistos e otimização paramétrica. **Inspira** o uso de limites relativos por blocos; não copiamos fórmulas específicas do artigo porque o conteúdo integral não foi consultado.
- **Boggiatto, Boiti e Oliaro**, *Umbrella theorems for time-frequency representations*, mesmo periódico, 16/11/2026, DOI [10.1007/s11565-026-00781-7](https://doi.org/10.1007/s11565-026-00781-7). Restrições de concentração de bases ortonormais e suas transformadas Fourier/Mellin são relevantes conceitualmente, sem fornecer uma estimativa numérica do projeto.
- **Kchaou e Krichen**, *On generalized essential pseudospectra ... operator matrices*, *Applied Mathematics—A Journal of Chinese Universities*, 19/10/2026, DOI [10.1007/s11766-026-4908-y](https://doi.org/10.1007/s11766-026-4908-y). Pseudoespectro de matriz 2×2 e perturbações, mas em Banach não reflexivo; exige comparação de hipóteses para ser usado no espaço Hilbert NUMERIA.
- **Shariati, Yaghouti e Alipanah**, *An operational matrix method based on scaling function of Daubechies wavelet to solve fractional differential equations*, mesmo periódico, 19/10/2026, DOI [10.1007/s11766-026-5005-y](https://doi.org/10.1007/s11766-026-5005-y). A base de wavelets de suporte compacto pode inspirar uma discretização local alternativa; o artigo resolve EDOs fracionárias e não identifica a forma de Weil.

Estas análises baseiam-se nas páginas e **resumos disponíveis**, não em provas completas dos textos pagos. A origem das identidades usadas no algoritmo P254 é Plancherel e o critério elementar de Schur, **não** os teoremas específicos destes trabalhos.

## Limites

P254 Gate A confirma apenas uma estratégia de reconstrução matematicamente correta e uma infraestrutura de validação fail-closed. O modelo auxiliar, o operador contínuo, uma ponte Weil–Zhu global, a exaustão cofinal e a hipótese de Riemann **não** recebem novos certificados de positividade nesta etapa.