# NUMERIA–E8K V16 — elo global e plano isolado da lacuna analítica

**Base de trabalho:** commit V15 `fe6bdb7463172ea7107cd050410fe0f94196de29`.  
**Estado:** análise de dependências e especificação do próximo lema; este documento não afirma que o lema analítico já foi provado.

## 1. Qual seria o teorema que liga o bloco arquimediano ao critério de Weil?

Não é um teorema de sinal do bloco 2×2 isoladamente. O elo necessário tem, pelo menos, duas partes distintas:

1. **Identidade global de Weil (com normalização e domínio de teste fixados):** reescrever a forma quadrática do critério de Weil como a soma/contribuição arquimediana, termos sobre primos e potências de primos, contribuições de polos e quaisquer termos de fronteira/projeção previstos pela formulação adotada.
2. **Teorema de transferência de positividade:** demonstrar que a forma quadrática global, com todos esses termos, satisfaz a desigualdade exigida para toda função de teste admissível (ou demonstrar uma equivalência com a positividade do critério de Weil).

A V15 demonstra, condicionalmente, uma propriedade da **contribuição arquimediana 2×2**. Ela não demonstra a identidade global, não limita os outros termos e não prova a positividade da forma quadrática global. De fato, uma contribuição arquimediana negativa não implica por si só que a forma global seja negativa: as outras contribuições ainda precisam ser calculadas e combinadas. Também não se pode deduzir RH de uma verificação de dimensão finita.

## 2. Separação do estado lógico

### A. Demonstrado em Lean (arquivo-alvo `NUMERIA_E8K_FORMAL/NumeriaE8K.lean`)

O arquivo define intervalos fechados racionais `RatI` e demonstra propriedades de contenção para operações elementares:
- singleton contém o ponto;
- contenção para soma, negação, subtração e produto por quatro cantos;
- contenção de quadrado;
- contenção do recíproco quando o intervalo está estritamente positivo ou estritamente negativo;
- identidades/cotas de largura para soma, subtração, negação e produto;
- arredondamento dirigido abstrato e propriedades de `floor/ceil` para o tick racional (10^{-42}).

O workflow V13 reporta que `lake build` terminou com sucesso no commit fixado. Isto verifica a compilação dessa formalização Lean, não a validade das afirmações analíticas fora dela.

### B. Verificado por aritmética exata / programas

- V13: as cotas racionais publicadas dominam as cotas recalculadas a partir das premissas analíticas codificadas; 10 gates numéricos.
- V14: integridade/proveniência declarada, replay de 15.754 pontos, agregação das 11 bandas e 128 oráculos racionais, além de controles negativos.
- V15: SHA-256 e manifesto, agregação exata das bandas, recomposição dos intervalos compostos, sinal dos extremos diagonais, teste racional de robustez para (|G_{12}|le 1/10), 7 controles negativos e reprodutibilidade do JSON.

Estes programas provam que as verificações aritméticas codificadas passam sobre os dados fixados. Não provam automaticamente que as fontes numéricas representam as integrais pretendidas nem que as premissas analíticas são verdadeiras.

### C. Premissas/condições ainda abertas

- Identidade de digama e representação de Fourier usadas na construção de (F_k).
- Cotas analíticas da derivada segunda que justificam o erro de ponto médio.
- Validade da majorante da cauda infinita atribuída à Proposição 6.
- Justificação dos valores de intervalo nodais (Omega(0)) importados da proveniência B3.
- Identidade de Fourier/paridade usada para concluir (G_{12}=0).
- Fórmula explícita global de Weil, normalização, termos primos, polos e condições do espaço de funções de teste.
- Passagem da fórmula global para a positividade exigida pelo critério de Weil em todas as funções de teste.

## 3. Próximo lema escolhido: validade analítica da majorante da cauda (Proposição 6)

**Por que este lema:** V13 verifica numericamente, com frações exatas, a desigualdade entre a cota de cauda publicada e a expressão racional calculada. Porém, a justificativa analítica de que a expressão realmente majoriza a cauda integral infinita permanece externa ao verificador. É um ponto de dependência explícito que afeta diretamente os intervalos compostos de (G_{11}) e (G_{22}).

### Enunciado a formalizar — sem ocultar hipóteses

A formulação final precisa ser transcrita da fonte matemática E7/Proposição 6, com todas as definições de (F_k), domínio, normalização e hipóteses. Em formato de obrigação de prova:

> Para cada índice (k) admissível e parâmetros (L,R) que satisfaçam as hipóteses da Proposição 6, a integral de cauda real correspondente, em (|r|ge R), é absolutamente convergente e seu valor absoluto é limitado pela expressão de cauda declarada no E7/V13.

A expressão computada pelo script contém os fatores
[
w=rac{kpi_{m hi}}{2L},qquad
D=1-(w/R)^2,
]
e exige, no mínimo, (Rge1), (2wle R) e (D>0). O código V13 usa ainda uma cota superior racional para (pi) e para o logaritmo. Esses testes algébricos são condições necessárias da avaliação, **não substituem** a demonstração de convergência nem a desigualdade integral.

### Plano isolado de formalização e testes

1. Extrair e citar literalmente o enunciado completo da Proposição 6 e suas definições a partir do pacote E7 fixado pelo SHA-256 `4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7`.
2. Dividir a prova em três lemas verificáveis separadamente: (i) positividade dos denominadores e domínio; (ii) majorante pontual do integrando; (iii) integração da majorante e limite da cauda.
3. Formalizar primeiro (i) em Lean sobre (mathbb Q), sem integrais nem constantes transcendentais; provar (ii) no ambiente analítico adequado somente depois de tornar explícitas as definições e hipóteses; provar (iii) usando um teorema de comparação de integrais e convergência absoluta.
4. Manter um script independente de teste racional para os pré-requisitos (Rge1), (2wle R), (D>0), positividade da expressão majorante e estabilidade dos limites. Incluir casos de fronteira e mutações negativas que devem falhar.
5. Só promover o lema a “demonstrado” quando o enunciado corresponder à fonte original, todas as dependências estiverem provadas no Lean ou citadas como hipóteses explícitas, e `lake build` + testes negativos passarem no mesmo commit.

### Critério de conclusão V16

- [ ] Enunciado matemático original e definições recuperados e registrados.
- [ ] Hipóteses e constantes nomeadas; nenhuma hipótese analítica implícita.
- [ ] Lema de domínio/denominadores provado em Lean.
- [ ] Majorante pontual e comparação integral provadas, ou cada dependência restante rotulada como premissa.
- [ ] Testes racionais reproduzíveis e mutações negativas.
- [ ] Workflow executado no commit que contém a formalização.
- [ ] Nenhuma afirmação de prova global de Weil/RH sem o elo global completo.

## 4. Conclusão atual

O próximo passo matemático não é ampliar a dimensão da matriz nem inferir positividade global a partir do bloco arquimediano. É fechar — ou declarar formalmente condicionada — a validade analítica da cauda que sustenta os intervalos, e em seguida enfrentar a identidade global de Weil com todos os termos. Este documento é um roteiro de auditoria/formalização, não o certificado de uma prova analítica já concluída.
