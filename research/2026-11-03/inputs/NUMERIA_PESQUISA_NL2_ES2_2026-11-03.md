# NUMERIA–RIEMANN — Radar bibliográfico NL2/ES2

**Data da pesquisa:** 03/11/2026  
**Escopo:** localizar literatura adicional e identificar uma rota de reprodução e auditoria do projeto NUMERIA. Complemento ao pacote `NUMERIA_NL1_ES1_2026-11-03` e às 44 referências classificadas nos dois CSVs anteriores.  
**Status:** pesquisa bibliográfica e projeto metodológico. **Nenhuma nova prova da Hipótese de Riemann, certificação de zeros, ou auditoria de matrizes foi executada nesta pesquisa.**

## 1. Achado principal — ponte direta para matrizes finitas de Weil

**Masatoshi Suzuki, “Spectral convergence associated with Weil's quadratic form” (preprint depositado em 23/10/2026),** https://arxiv.org/abs/2610.32704 . No Teorema 1.1, **assumindo RH**, define operadores de diferenciação autoadjuntos em espaços baseados na forma de Weil localizada e prova aproximação espectral, com valores próprios positivos `t_(a,n)` aproximando as ordenadas `γ_n` dos zeros distintos à medida que a expansão de suporte `a→∞` prossegue. No Teorema 1.2 demonstra convergência de funções inteiras espectrais normalizadas — **também sob RH**.

A Seção 6 fornece um experimento reproduzível a partir dos dados aritméticos da fórmula de Weil, sem inserir zeros da zeta como insumos:

- `χ_(n,a)(x) = exp(iπnx/a) 1_{[-a,a]}(x) / sqrt(2a)`, para `n=-N,...,N`.
- `G_(a,N) = (Q_W(χ_(m,a), χ_(n,a)))_(m,n=-N,...,N)`, matriz complexa Hermitiana de ordem `2N+1`.
- A fórmula explícita de Weil exige apenas potências de primos até `p^r ≤ exp(2a)` por causa do suporte da convolução. A construção das entradas da matriz **não exige RH**; as demonstrações de convergência espectral apresentadas no artigo exigem.
- Para `N=18` (matriz 37×37) o autor apresenta experimentos com `a=(1/2) log X`, `X=3,5,7,11`, 140/180 dígitos, e encontra aproximações para os seis primeiros zeros. Observa números de condição elevados: até `2,33×10^34` no último caso. **A concordância numérica não representa uma cota de erro rigorosa.**

**Implicação para o NUMERIA:** não se deve confundir essa base trigonométrica complexa de 37 modos com a base de senos/cossenos e matrizes projetadas E7/V24. Primeiro deve-se implementar as duas definições separadamente e cotejar normalização de Fourier, sinais, polos, termos primos, projeções e índices. Após isso, produzir testes de regressão e apenas então tentar provas de equivalência.

## 2. Trabalhos ligados diretamente à construção espectral

1. **Masatoshi Suzuki (2026), “Weil's quadratic form via the screw function”**, https://arxiv.org/abs/2606.09096 (versão 3 de 23/09/2026). Reescreve a forma de Weil por funções contínuas e métodos de espaços de de Branges; contém uma **conjectura** espectral envolvendo limites de operadores de diferenciação de intervalos finitos. A construção do quadro é apresentada sem pressupor RH; a conjectura permanece em aberto.
2. **Masatoshi Suzuki (2025), “On the Hilbert space derived from the Weil distribution”**, DOI https://doi.org/10.4153/S0008414X25101739, publicação on-line 03/11/2025. Relaciona espaço de Weil a espaços de de Branges e estabelece condições equivalentes a RH por igualdades de normas. O isomorfismo de espaços da introdução é descrito **sob RH**. Não confundir equivalência com prova incondicional.
3. **Alain Connes, Caterina Consani, Henri Moscovici (2026), “Zeta spectral triples”**, https://doi.org/10.4171/ELM/37/3 . Produz operadores autoadjuntos com espectros numericamente próximos aos zeros de zeta a partir de produtos de Euler finitos. O próprio resumo diz que **uma prova rigorosa de convergência** para o espectro global estabeleceria RH; ela ainda não foi demonstrada no artigo.
4. **Alain Connes e Caterina Consani (2023), “Spectral triples and ζ-cycles”**, https://doi.org/10.4171/LEM/1049 . Estudos de pequenos autovalores da forma de Weil localizada, funções prolatas e evidência numérica de aproximação espectral. Referência-chave para interpretar autovalores muito pequenos e riscos de condicionamento.
5. **H. Yoshida (1992), “On Hermitian forms attached to zeta functions”**, citado por Suzuki 2026, seção 6, fórmulas (5.15)–(5.16). Fonte a obter integralmente antes de transcrever qualquer fórmula de matriz. **Ainda não cotejado integralmente nesta pesquisa.**
6. **Kohji Matsumoto e Masatoshi Suzuki (2026), “M-functions and screw functions: Applications to Goldbach's problem and zeros of the Riemann zeta-function”**, https://doi.org/10.1016/j.jnt.2025.09.013 , *Journal of Number Theory* 280, 918–946. O resumo anuncia critérios necessários e suficientes para RH em termos de funções screw; a aplicação operacional à base atual ainda não foi demonstrada.

## 3. Certificação de zeros e solução de sistemas não lineares

7. **D. Platt e T. Trudgian (2021), “The Riemann hypothesis is true up to 3·10^12”**, https://doi.org/10.1112/blms.12460 . Certificação rigorosa de todos os zeros até altura `3×10^12` com aritmética intervalar; referência para validação externa de métodos, não uma demonstração global.
8. **Biblioteca FLINT — zeros da zeta**, https://flintlib.org/doc/acb_dirichlet/zeta-zeros.html e https://fredrikj.net/fdoc/html/acb_dirichlet/zeta-zeros.html . A API nativa `acb_dirichlet_hardy_z_zeros` realiza isolamento e refinamento com métodos de Riemann–Siegel ou Platt, combinando contagem de zeros pelo método de Turing. As ligações Python documentam `acb.zeta_zero(n)`, `acb.zeta_zeros(n,num)` e `acb.zeta()` em https://python-flint.readthedocs.io/en/latest/acb.html . **A instalação atual não contém `python-flint`; nenhuma execução dessas funções foi feita neste levantamento.**
9. **V. Semenov (2011), “Method for the Calculation of All Non-Multiple Zeros of an Analytic Function”**, https://doi.org/10.2478/cmam-2011-0003 . Propõe subdivisão de retângulos com exclusão de zonas sem raízes ou localização de uma raiz; seus resultados tratam zeros não múltiplos. Complementa o Krawczyk para evitar a perda silenciosa de zeros em uma caixa finita.
10. **V. Semenov e E. Semenova (2019), “Method for Localizing the Zeros of Analytic Functions Based on the Krawczyk Operator”**, https://doi.org/10.1007/s10559-019-00159-y . Referência já classificada na NL1; relata aplicação exemplificativa à zeta. Para certificação de zeta faltam avaliar rigorosamente ζ e ζ′ em retângulos complexos e garantir coberturas de exclusão/contagem quando se afirma completude.
11. **Jean-François Burnol (2000), “Sur les formules explicites I: analyse invariante”**, https://doi.org/10.1016/S0764-4442(00)01687-6 . Verifica critério de positividade sob **condições de suporte**. Exige determinar explicitamente as condições antes de afirmar que `L=0,45` está protegido por tal resultado.
12. **J. Brian Conrey e Xian-Jin Li (1998), “A note on some positivity conditions related to zeta- and L-functions”**, https://arxiv.org/abs/math/9812166 . Analisa dificuldades de propostas de de Branges por condições de positividade. Útil como teste adversarial lógico contra inferência de positividade geral a partir de subconjuntos finitos.

## 4. Rota técnica sugerida

### Trilha NL2 — certificação de zeros, como controle de qualidade

1. Preparar CI com versão de FLINT e python-flint fixadas; confirmar API disponível e versão na execução.
2. Comparar os três zeros aproximados pela NL1 com `acb.zeta_zero(1..3)`; registrar intervalos e identidades de refinamento. Essa reprodução testa a engenharia, mas **não agrega nova prova matemática** a fatos já certificados.
3. Construir `K(X)=c-Rζ(c)+(1-Rζ′(X))(X-c)` em aritmética complexa intervalar, numa caixa livre de polos. Com prova de inclusão estrita e de contração, certificar existência/unicidade **local**. Zeros múltiplos e completude global exigem métodos adicionais.
4. Se o objetivo é excluir *todos* os zeros em uma região, cobrir toda a região por caixas com testes de exclusão e inclusões validadas, ou usar princípio do argumento/Turing para certificar contagem. **Somente Newton convergente e resíduos pequenos não bastam.**

### Trilha ES2 — reprodução da forma de Weil em matrizes de Suzuki

1. Usar primeiro a definição do autor para `G_(a,N)`: entradas aritméticas independentes das ordenadas dos zeros, com `a=log(X)/2` e `X ∈ {3,5,7,11}`.
2. Testar `N=3,6,10,18` (matrizes 7×7, 13×13, 21×21 e 37×37). Registrar a condição numérica `κ`, uso de precisão intervalar, sinais e marcos de projeção.
3. Implementar resolução **verificada** dos sistemas lineares `G_(a,N) c = φ_(a,N)(z)` e obter funções espectrais aproximadas `W_N`, seguindo as fórmulas (6.3)–(6.5) do artigo; rastrear erros de inversão, truncamento e sensibilidade.
4. Comparar as raízes finitas detectadas com valores de referência certificados; guardar **apenas como observação** até haver prova de isolamento das raízes de `W_N` e cotas de erro de Galerkin.
5. Em paralelo, cotejar E7/V24: Fourier, digama, fatores π, polos, primos, Fourier conjugada, base trigonométrica, projeções B*GB e caudas. Sem equivalência demonstrada, relatórios de ambos os modelos devem permanecer separados.
6. Não usar a convergência condicionada a RH para concluir RH: isso seria raciocínio circular.

## 5. Decisões de prioridade

- **P1 (imediato):** Suzuki, arXiv:2610.32704, seção 6; Connes–Consani–Moscovici 2026; Suzuki 2606.09096. A combinação apresenta uma rota concreta de matriz aritmética, sistemas lineares, funções espectrais e análise de condicionamento.
- **P2:** FLINT + Platt–Trudgian + Semenov para certificar zeros conhecidos e testar ferramentas de validação.
- **P3:** Suzuki 2025 e Matsumoto–Suzuki 2026 para critérios equivalentes, sem assumir que eles removem as hipóteses exigidas pelas provas.
- **P4:** Yoshida 1992, suporte de Burnol e análise de Conrey–Li para uma auditoria rigorosa das premissas e da extrapolação finita→infinita.

## 6. Pendências honestas

- Nenhuma reprodução dos experimentos de Suzuki foi executada nesta pesquisa.
- `python-flint` não está instalado no ambiente desta sessão.
- O artigo 2011 de Semenov foi identificado por título/resumo/metadados; nenhuma prova do corpo completo foi inspecionada aqui.
- Nem a conexão aritmética E7/V24⇔G de Suzuki nem a cota uniforme na dimensão N foram demonstradas.
- Uma certificação de qualquer número finito de zeros **não prova RH**; o passo infinito continua sendo a dificuldade central.
