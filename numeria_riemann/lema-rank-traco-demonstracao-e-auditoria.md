# Lema rank–traço (Lema 3.2): demonstração e auditoria

**Enunciado.** Sejam $P, Q$ matrizes hermíticas $d \times d$ com $P \succeq 0$, $\operatorname{rank} P \leq r$ e $n_+(Q) \leq b$, onde $n_+(Q)$ é o número de autovalores estritamente positivos de $Q$. Então

$$r \;\geq\; 2\operatorname{tr} P + 4\operatorname{tr} Q - 4b - \|P + Q\|_{\mathrm{HS}}^2, \tag{3.1}$$

onde $\|M\|_{\mathrm{HS}}^2 = \operatorname{tr}(M^2)$. A igualdade ocorre exatamente quando $\operatorname{rank}P=r$, $n_+(Q)=b$ e, após conjugação unitária, $P=\Pi_1$ e $Q=2\Pi_2$, para projeções ortogonais $\Pi_1 \perp \Pi_2$ de postos $r$ e $b$, respectivamente, ambas nulas no complemento ortogonal comum.

O lema é a forma matricial das desigualdades escalares $(x-1)^2 \geq 0$ e $(x-2)^2 \geq 0$; com $Q = 0$ e otimização do coeficiente $2$, recupera-se $\operatorname{rank} P \geq (\operatorname{tr} P)^2 / \|P\|_{\mathrm{HS}}^2$.

**Status das seções deste documento.**

| Seção | Conteúdo | Status |
|---|---|---|
| Prova 1 | Rota de von Neumann (a prova da fonte) | Provada aqui em detalhe; preservada com revisão matemática |
| Prova 2 | Rota variacional (expansão da Observação 3.3 da fonte) | Provada aqui com os passos omissos preenchidos; auditoria dos pontos delicados incluída |
| Aplicação à $\zeta$ | A cadeia $N_0^s \geq (2 - R(\psi))N$ | **Não provada aqui.** Reproduzida com as hipóteses e estimativas externas enumeradas; protocolo de auditoria |
| Evidência computacional | Monte Carlo | **Não é prova nem validação independente.** Apenas teste de falseamento |

A prova pura é a Parte I. Nada na Parte I depende da função zeta: é álgebra linear sobre $\mathbb{C}$ (vale verbatim sobre $\mathbb{R}$ com matrizes simétricas reais e conjugação ortogonal).

---

## Parte I — O lema

### Prova 1 (rota de von Neumann) — preservada, com revisão

**Linha 1 — decomposição canônica de $Q$.** Seja $Q = Q_+ - Q_-$ a decomposição espectral canônica: $Q_+$ retém os autovalores estritamente positivos, $Q_-$ os estritamente negativos, e o núcleo de $Q$ não entra em nenhum dos dois. Então $Q_\pm \succeq 0$, $Q_+ Q_- = Q_- Q_+ = 0$, $\operatorname{rank} Q_+ = n_+(Q) \leq b$, e

$$\operatorname{tr} Q = \operatorname{tr} Q_+ - \operatorname{tr} Q_-. $$

**Linha 2 — expansão do quadrado de Hilbert–Schmidt.** Como $\|M\|_{\mathrm{HS}}^2 = \operatorname{tr}(M^2)$ e $\operatorname{tr}((P+Q)^2) = \operatorname{tr}(P^2) + 2\operatorname{tr}(PQ) + \operatorname{tr}(Q^2)$, com $\operatorname{tr}(PQ) = \operatorname{tr}(PQ_+) - \operatorname{tr}(PQ_-)$ e $Q^2 = Q_+^2 + Q_-^2$ (os termos cruzados são nulos):

$$\|P + Q\|_{\mathrm{HS}}^2 = \|P\|_{\mathrm{HS}}^2 + \|Q_+\|_{\mathrm{HS}}^2 + \|Q_-\|_{\mathrm{HS}}^2 + 2\operatorname{tr}(PQ_+) - 2\operatorname{tr}(PQ_-).$$

**Linha 3 — descarte dos termos não negativos.** *Justificativa de $\operatorname{tr}(PQ_+) \geq 0$* (adicionada nesta revisão): como $P \succeq 0$, existe $S \succeq 0$ com $P = S^2$; então $\operatorname{tr}(PQ_+) = \operatorname{tr}(S Q_+ S) \geq 0$, pois $S Q_+ S \succeq 0$ e o traço de uma matriz semidefinida positiva é não negativo. Logo

$$\|P + Q\|_{\mathrm{HS}}^2 \;\geq\; \underbrace{\|P\|_{\mathrm{HS}}^2 - 2\operatorname{tr}(PQ_-) + \|Q_-\|_{\mathrm{HS}}^2}_{E} \;+\; \|Q_+\|_{\mathrm{HS}}^2.$$

**Linha 4 — desigualdade de traço de von Neumann.** Sejam $p_1 \geq \cdots \geq p_d \geq 0$ os autovalores de $P$ e $n_1 \geq \cdots \geq n_d \geq 0$ os de $Q_-$. A desigualdade de von Neumann (hermíticas, ambas ordenadas de forma decrescente) dá $\operatorname{tr}(P Q_-) \leq \sum_i p_i n_i$. Portanto

$$E \;\geq\; \sum_i p_i^2 - 2\sum_i p_i n_i + \sum_i n_i^2 \;=\; \sum_i (p_i - n_i)^2.$$

**Linha 5 — desigualdades escalares, coordenada a coordenada.** Como $\operatorname{rank} P \leq r$ e os $p_i$ estão em ordem decrescente, $p_i = 0$ para todo $i > r$. Então, para $i \leq r$:

$$\sum_{i \leq r} (p_i - n_i)^2 \;\overset{(x-1)^2 \geq 0}{\geq}\; \sum_{i \leq r} \bigl(2(p_i - n_i) - 1\bigr) = 2\operatorname{tr} P - 2\sum_{i \leq r} n_i - r \;\geq\; 2\operatorname{tr} P - 4\operatorname{tr} Q_- - r,$$

onde o último passo usa $n_i \geq 0$, logo $-2n_i \geq -4n_i$; e, para $i > r$:

$$\sum_{i > r} (p_i - n_i)^2 = \sum_{i > r} n_i^2 \;\overset{x^2 \geq 0}{\geq}\; 0 \;\geq\; -4\sum_{i > r} n_i.$$

No bloco positivo, seja $q_1, \ldots, q_{b'}$ ($b' = \operatorname{rank} Q_+ \leq b$) a lista de autovalores não nulos de $Q_+$:

$$\|Q_+\|_{\mathrm{HS}}^2 = \sum_{j \leq b'} q_j^2 \;\overset{(x-2)^2 \geq 0}{\geq}\; \sum_{j \leq b'} (4q_j - 4) = 4\operatorname{tr} Q_+ - 4b' \;\geq\; 4\operatorname{tr} Q_+ - 4b.$$

**Linha 6 — soma e rearranjo.**

$$\|P + Q\|_{\mathrm{HS}}^2 \;\geq\; \bigl(2\operatorname{tr} P - 4\operatorname{tr} Q_- - r\bigr) + \bigl(4\operatorname{tr} Q_+ - 4b\bigr) = 2\operatorname{tr} P + 4\operatorname{tr} Q - r - 4b,$$

que é (3.1) rearranjada. $\blacksquare$

**Revisão matemática da Prova 1 (nota do editor).** Pontos conferidos nesta revisão: (i) a expansão da Linha 2, incluindo a anulação dos termos cruzados $Q_+Q_-$, que exige a decomposição *canônica* da Linha 1 (uma repartição arbitrária do núcleo também funciona, mas então $Q_+Q_-$ é apenas semidefinida e o termo extra $-2\operatorname{tr}(Q_+Q_-) \leq 0$ se descarta igualmente); (ii) a aplicação de von Neumann exige o emparelhamento decrescente–decrescente — o mesmo emparelhamento usado em $\sum(p_i - n_i)^2$; (iii) cada desigualdade escalar é usada no único lugar onde é necessária: $(x-1)^2 \geq 0$ nas coordenadas $i \leq r$ (carga 1 por autovalor de $P$), $(x-2)^2 \geq 0$ no bloco positivo de $Q$ e, de forma oculta, no rebaixamento $-2n_i \geq -4n_i$ e em $n_i^2 \geq 4n_i - 4$ na cauda; (iv) a prova prova o enunciado com $r := \operatorname{rank} P$; como o lado direito não contém $r$, vale para qualquer $r \geq \operatorname{rank} P$; (v) **nitidez**: a igualdade em (3.1) força igualdade em von Neumann (autovetores comuns), $p_i = 1$ nos $i \leq r$ e $q_j = 2$ nos $j \leq b$ — consistente com o caso de igualdade classificado na Prova 2.

### Prova 2 (variacional) — corrigida e validada

Esta é a expansão completa da Observação 3.3 da fonte. A versão anterior deste documento comprimia dois passos delicados; eles estão agora explícitos e com a direção dos argumentos verificada.

Defina, no grupo unitário $U(d)$ (compacto),

$$g(U) \;=\; \bigl\|P + U Q U^*\bigr\|_{\mathrm{HS}}^2.$$

**Passo 0 — redução ao minimizador (a direção do argumento, explicitada).** Os invariantes: $\operatorname{tr}(UQU^*) = \operatorname{tr} Q$, o espectro (logo $n_+$ e o posto), e $\operatorname{rank} P$ são os mesmos ao longo da órbita. Como $g$ é contínua em conjunto compacto, atinge um mínimo em algum $U_0$; ponha $Q_0 := U_0 Q U_0^*$. **A redução funciona em uma única direção:** se provarmos

$$r \;\geq\; 2\operatorname{tr} P + 4\operatorname{tr} Q_0 - 4b - \|P + Q_0\|_{\mathrm{HS}}^2, \tag{$*$}$$

então, como $g(I) = \|P + Q\|_{\mathrm{HS}}^2 \geq g(U_0) = \|P + Q_0\|_{\mathrm{HS}}^2$ (definição de mínimo) e $2\operatorname{tr} P + 4\operatorname{tr} Q - 4b = 2\operatorname{tr} P + 4\operatorname{tr} Q_0 - 4b$ (invariância), a desigualdade $(*)$ implica (3.1) para $(P, Q)$. O mínimo *não* precisa ser classificado — basta a existência e a anulação da derivada nele.

**Passo 1 — condição de primeira ordem.** Seja $U(t) = e^{itK}U_0$ com $K$ hermítica arbitrária. Então $Q(t) := U(t)Q_0U(t)^* = e^{itK} Q_0 e^{-itK}$ satisfaz $\dot Q(0) = i[K, Q_0]$, que é hermítica. Como $\frac{d}{dt}\operatorname{tr}(M(t)^2) = 2\operatorname{tr}(M \dot M)$ para caminhos hermíticos,

$$\frac{d}{dt}\bigg|_0 g(U(t)) = 2i\,\operatorname{tr}\bigl((P + Q_0)\,[K, Q_0]\bigr).$$

*Os dois traços, separadamente:* por ciclicidade,

$$\operatorname{tr}\bigl(Q_0[K, Q_0]\bigr) = \operatorname{tr}(Q_0 K Q_0) - \operatorname{tr}(Q_0^2 K) = \operatorname{tr}(K Q_0^2) - \operatorname{tr}(Q_0^2 K) = 0,$$

$$\operatorname{tr}\bigl(P[K, Q_0]\bigr) = \operatorname{tr}(P K Q_0) - \operatorname{tr}(P Q_0 K) = \operatorname{tr}((Q_0 P - P Q_0)K) = -\operatorname{tr}\bigl([P, Q_0]\,K\bigr).$$

*Por que isso cobre todas as direções:* as curvas $e^{itK}U_0$ geram todo o espaço tangente em $U_0$ (a álgebra de Lie de $U(d)$ é $\{iK : K \text{ hermítica}\}$). No mínimo, a derivada se anula ao longo de todas elas:

$$\operatorname{tr}\bigl([P, Q_0]\,K\bigr) = 0 \quad \text{para toda } K \text{ hermítica}.$$

*Conclusão:* $[P, Q_0]$ é anti-hermítica; escreva $[P, Q_0] = iA$ com $A$ hermítica. A condição vira $\operatorname{tr}(AK) = 0$ para toda $K$ hermítica; escolhendo $K = A$ obtemos $0 = \operatorname{tr}(A^2) = \|A\|_{\mathrm{HS}}^2$, logo $A = 0$ e

$$[P, Q_0] = 0.$$

**Passo 2 — diagonalização simultânea.** Hermíticas que comutam diagonalizam juntas: existe base ortonormal comum com autovalores $(p_i, b_i)$, $i = 1, \ldots, d$. Nessa base, $\operatorname{rank} P = \#\{i : p_i > 0\}$, $n_+(Q_0) = \#\{i : b_i > 0\}$, e

$$\|P + Q_0\|_{\mathrm{HS}}^2 = \sum_i (p_i + b_i)^2.$$

Basta então provar $\sum_i \Delta_i \geq 0$ com o excedente por coordenada

$$\Delta_i = \mathbf{1}_{p_i > 0} + 4\,\mathbf{1}_{b_i > 0} + (p_i + b_i)^2 - 2p_i - 4b_i.$$

**Os três casos de sinal.**

*Caso 1: $p_i > 0$, $b_i > 0$.* Com $s = p_i + b_i$, usando $-2p_i - 4b_i = 2p_i - 4s$:

$$\Delta_i = s^2 - 2p_i - 4b_i + 5 = (s - 2)^2 + 2p_i + 1 \;\geq\; 1 \;>\; 0.$$

Nunca há igualdade neste caso.

*Caso 2: $p_i > 0$, $b_i \leq 0$.* Com $n = -b_i \geq 0$ e $y = p_i - n$:

$$\Delta_i = 1 + (p_i - n)^2 - 2p_i + 4n = (y - 1)^2 + 2n \;\geq\; 0,$$

com igualdade se, e somente se, $n = 0$ e $y = 1$, isto é, $(p_i, b_i) = (1, 0)$.

*Caso 3: $p_i = 0$.* Se $b_i > 0$: $\Delta_i = (b_i - 2)^2 \geq 0$, com igualdade se, e somente se, $b_i = 2$. Se $b_i \leq 0$: $\Delta_i = b_i^2 - 4b_i \geq 0$, com igualdade se, e somente se, $b_i = 0$.

Somando: $\sum_i \Delta_i \geq 0$, o que prova $(*)$ e, pelo Passo 0, (3.1). $\blacksquare$

**Caso de igualdade (completo).** A soma é zero se, e somente se, toda coordenada está em um dos três tipos de igualdade: $(1, 0)$, $(0, 2)$ ou $(0, 0)$. Isto é: $P = \Pi_1$ e $Q = 2\Pi_2$ com $\Pi_1 \perp \Pi_2$ projeções ortogonais de postos $r$ e $b$, ambas nulas no complemento ortogonal comum. Em particular, o caso 1 (mistura de sinais na mesma coordenada) *nunca* é gratuito.

**Auditoria da Prova 2 (nota do editor).** Os quatro pontos delicados, e como foram resolvidos: (1) *a direção da redução ao minimizador* — o Passo 0 agora mostra explicitamente que $(*) \Rightarrow (3.1)$ usando apenas $g(U_0) \leq g(I)$; não há necessidade de classificar minimizadores; (2) *cobertura do espaço tangente* — justificada via álgebra de Lie de $U(d)$; (3) *a condição de primeira ordem* — a separação dos dois traços e a identificação $[P,Q_0] = iA$ com a escolha $K = A$ estão agora completas (na versão anterior, a passagem era comprimida); (4) *o caso de igualdade* — a coordenada $(0,0)$ estava omitida; incluída. **Limitação honesta:** esta auditoria é uma revisão interna feita nesta edição, com verificação numérica de apoio (Parte III); não substitui um referee independente. A fonte relata (Apêndice A) que a desigualdade de traço de von Neumann e as duas direções da lei de inércia de Sylvester foram formalizadas em Lean no repositório que acompanha o artigo.

### Leitura no vocabulário da função zeta (interpretativa)

No teorema da fonte, $\widetilde{G} = P + Q$ comprime a forma de Weil numa família de Gabor com $d \sim N(T, 2T)$ elementos. Os três casos fotografam as espécies de zeros:

| Coordenada | Espécie | Tarifa |
|---|---|---|
| $(1, 0)$ — igualdade no caso 2 | zero simples na reta crítica | 1 (gratuíto) |
| $(0, 2)$ — igualdade no caso 3 | par fora da reta $\{\rho, 1-\bar\rho\}$, ou ponto múltiplo no reagrupamento $(P_1, Q')$ | 4 |
| $p_i > 0$, $b_i > 0$ — caso 1 | mistura de "na reta" e "fora da reta" na mesma direção | sempre sobra $\geq 1$ |

**Dois níveis de integrabilidade** (Observação da fonte, §7.2(c)): com a partição $(P, Q)$ o lema recupera $m^2 \geq 2m - 1$; com o reagrupamento $(P_1, Q')$ (múltiplos migram para o lado cobrado a 4) recupera também $m^2 \geq 3m - 2$, i.e. $(m-1)(m-2) \geq 0$: zeros simples e duplos são gratuitos, multiplicidade $m \geq 3$ paga a mais.

---

## Parte II — A cadeia até $\tfrac23$: hipóteses, estimativas e protocolo de auditoria

**O que esta parte é e o que não é.** O teorema $\tfrac23$ (zero simples e na reta crítica, $N_0^s \geq (\tfrac23 - o(1))N$) é **resultado da fonte**, que relata formalização em Lean (Apêndice A; `#print axioms` devolvendo apenas os três axiomas padrão, conforme o relato da fonte). Esta parte **não o reprova**: enumera a cadeia completa, com cada hipótese e cada termo de erro com sua localização na fonte, para que um auditor independente possa verificá-los um a um. Os únicos itens verificados diretamente nesta revisão estão marcados com ✓.

### A cadeia

Com $P_1$ = parte de $\widetilde{G}$ vinda dos zeros simples na reta e $Q' = \widetilde{G} - P_1$:

$$N_0^s + o(N) \;\overset{\text{Prop. 4.1}}{\geq}\; \operatorname{rank} P_1 \;\overset{\text{Lema 3.2}}{\geq}\; 4\operatorname{tr} \widetilde{G} - 2N - \|\widetilde{G}\|_{\mathrm{HS}}^2 \;\overset{\text{Prop. 4.2 + T. 5.7}}{=}\; \bigl(2 - R(\psi) - o(1)\bigr)N.$$

Os três insumos analíticos, com suas hipóteses e erros:

| Insumo | Enunciado | Hipóteses que usa | Erro |
|---|---|---|---|
| Prop. 4.1 (lado dos zeros) | $\operatorname{rank} P_1 \leq s_1$, $\operatorname{tr} P_1 \leq s_1$, $n_+(Q') \leq s_2 + p$, $\operatorname{tr} P_1 + 2n_+(Q') \leq N$ | equação funcional ($\rho \leftrightarrow 1 - \bar\rho$, daí $v_{1-\bar\rho} = \overline{v_\rho}$); Lema 2.1 em $z = z'$; Lema 3.1 (inércia sob pull-back) | exato (identidades) |
| Prop. 4.2 (traço) | $\operatorname{tr} \widetilde{G} = N(I') + O_\chi(T^{1/2}L^2)$ | Lema 2.1 (Poisson–Gabor, densidade crítica $h = 2\pi/L$); decaimento (2.8) de $\widehat\phi$; densidade de zeros $\ll L$; suavização $C^2$ | $O_\chi(T^{1/2}L^2) = o(N)$ |
| T. 5.7 (norma HS) | $\|\widetilde{G}\|_{\mathrm{HS}}^2 = (R(\psi) + O_\chi(L^{-1}))N$ | fórmula explícita de Weil (2.1)–(2.2); Chebyshev–Mertens (Lema 5.1); desigualdade de Montgomery–Vaughan (Lema 2.2) com $\delta_n^{-1} \leq 2n$, válida para $X \leq T$; Stirling para $\Gamma'/\Gamma$; cauda Prop. 4.3 com $N(t,t+1) \ll \log t$ | $O_\chi(L^{-1})N$ |

Taxa final (Remark 6.1 da fonte): $N_0^s \geq \bigl(2 - R(\psi) - c_\psi \frac{\log\log T}{\log T}\bigr)N$ para $T \geq T_0(\psi)$.

### Constantes das janelas

- $R(\psi_0) = \tfrac{4}{3}$ — **verificado simbolicamente nesta revisão** ✓ ($\int\psi_0^2 = 1$, $\iint|u - v|\,\psi_0\psi_0 = \tfrac13$, $\int\psi_0 = 1$), dando $\tfrac23$.
- $R(\psi_{\mathrm{MT}}) = \tfrac12 + \tfrac{1}{\sqrt2}\cot\tfrac{1}{\sqrt2} = 1{,}3274993\ldots$ — **verificado numericamente** ✓, dando $0{,}67250$ (simples na reta) e $0{,}83625$ (distintos).

### Protocolo de auditoria (o que um auditor independente deve verificar na fonte)

1. **Lema 2.1 (Poisson–Gabor):** o retículo dual de $h = 2\pi/L$ é $L\mathbb{Z}$, só o termo $m = 0$ contribui, e $\operatorname{supp}(\phi * \phi) \subset [-L, L]$ garante ausência de aliasing. É a única identidade que faz $\operatorname{tr}\widetilde{G}$ "contar" zeros.
2. **Prop. 4.2 e 4.3:** a soma de cauda $\sum_k \vartheta(\cdot)^2$ com decaimento $|r|^{-4}$; a cobertura dos zeros com $D_\rho$ pequeno; a necessidade da suavização $C^2$ (Remark 4.4 da fonte: com corte abrupto a cauda **não** é $o(N)$).
3. **Prop. 5.4 (diagonal dos primos):** a avaliação do termo $n = m$; o controle de $\mathcal{O}_1$ (fora-diagonal) pela desigualdade de Montgomery–Vaughan com $\delta_n^{-1} \leq 2n$ — o único ponto onde a restrição $X \leq T$ entra.
4. **T. 5.7:** a redução à integral dupla (Prop. 5.2, erro $O(N\log L / L^2)$ — origem do $\log\log T/\log T$ da taxa); os termos cruzados (Prop. 5.5, $\ll_\chi L^2\sqrt{X}$).
5. **Constantes de janela:** $R(\psi_0) = 4/3$ ✓ e $R(\psi_{\mathrm{MT}})$ ✓ (verificados nesta revisão, simbólico/numérico).
6. **A regulação $(P_1, Q')$** e o segundo nível de integrabilidade $m^2 \geq 3m - 2$.
7. **Livro-contábil do $o(N)$:** a versão acumulada por janelas diádicas e a dedução de $N_d \geq \tfrac12(3 - R(\psi))N$.
8. **Formalização:** conferir no repositório da fonte (tag v1.0) os enunciados de topo (`two_thirds_on_critical_line` e correlatos) e o resultado de `#print axioms` relatado.

### Limites inerentes (para dimensionar a auditoria)

- A desigualdade (3.1) é **sharp** nas cargas declaradas: a configuração $\tfrac23 N$ zeros simples ortogonais na reta $+ \tfrac16 N$ duplos realiza $\operatorname{tr} = N$, $\|\cdot\|_{\mathrm{HS}}^2 = \tfrac43 N$, $s_1 = \tfrac23 N$. Nada na cadeia distingue "dois terços" de "todos" com informação de ordem 1.
- A classe de certificados de banda 1 tem teto $\leq 0{,}6819$ (teorema Lean da fonte, §7.2); ambos os insumos valem para as funções de Davenport–Heilbronn e Epstein, nas quais o análogo da RH é **falso**. Logo esta cadeia, por construção, não pode ultrapassar afirmações de densidade — a auditoria acima valida a cadeia, não a RH.

---

## Parte III — Evidência computacional: status e limites

**O que os testes são.** Um experimento de *falseamento*: procurar contraexemplos numéricos ao enunciado (3.1). Parâmetros: matrizes simétricas reais $7 \times 7$; conjugações ortogonais aleatórias via Gram–Schmidt; gerador congruente-linear determinístico (semente 123456789); Box–Muller para gaussianas.

**Resultados.** 4.000 pares aleatórios $(P, Q)$ com posto e índice controlados: **0 violações** (menor excedente $+5{,}89$). Configuração de igualdade $P = \Pi_1$, $Q = 2\Pi_2$ (50 instâncias): erro $0$ em ponto flutuante. Formas fatoradas dos três casos: coincidem com a definição bruta com erro máximo $1{,}8 \times 10^{-15}$ em 2.000 pontos.

**O que os testes não são.** (i) **Não são prova**: as demonstrações da Parte I são a única validação matemática aqui; (ii) **não são validação independente** dos resultados da fonte sobre $\zeta$: nenhum teste abaixo envolve a função zeta, os primos ou qualquer estimativa analítica da Parte II; (iii) usam aritmética de ponto flutuante com tolerâncias ($10^{-8}$ para violações, $10^{-15}$ para identidades) — erros de arredondamento não são controlados como em aritmética de intervalos certificada; (iv) cobrem apenas matrizes simétricas **reais**, um subconjunto próprio das hermíticas; (v) são amostrais: ausência de contraexemplo em 4.000 casos não implica ausência em geral. Em suma: plausibilidade e checagem de consistência, no espírito de um teste de unidade; o estatuto de prova reside exclusivamente nas Partes I (aqui) e II (na fonte, com a formalização Lean lá relatada).

---

## Referências

- **Fonte primária:** *More than two thirds of the zeros of the Riemann zeta function are simple and on the critical line* (paper-v5), §2–§7: Lema 3.2 e Observação 3.3 (as duas provas), Proposições 4.1–4.3 e Teorema 5.7 (a cadeia), §7.2 (nitidez e teto), Apêndice A (formalização Lean).
- J. von Neumann, desigualdade de traço para matrizes hermíticas.
- E. Beltrami / J. Sylvester, lei de inércia de formas hermíticas.
- H. L. Montgomery (1973); F. Aryan (2022); S. A. C. Baluyot, D. A. Goldston, A. I. Suriajaya, C. L. Turnage-Butterbaugh (2024) — o lado dos primos.
- J. B. Conrey, A. Ghosh, S. M. Gonek (1998) — o dispositivo de cargas $m^2 \geq 3m - 2$ sob RH, tornado incondicional pela cadeia da fonte.