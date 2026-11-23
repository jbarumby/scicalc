# NUMERIA–RIEMANN — Auditoria do PDF “Reconstrução Detalhada: Passagem do Lema 5 para o Lema 6”

**Fonte:** PDF de cinco páginas enviado pelo usuário. **Veredito:** relevante e recuperável como **passagem condicional do Lema 5 ao Lema 6**. **Não é prova independente do Lema 5**, pois assume sua desigualdade com termo uniforme `O(z1/zM)` na equação (1). Três pontos do texto exigem correção antes de qualquer atribuição de rigor à demonstração.

## 1. Página 2, equação (6) — sinal invertido numa igualdade

Para r,s,t distintos, escrevendo h=s-t, a identidade correta é

\[
\frac1{(r-s)(r-t)}=\frac1h\left(\frac1{r-s}-\frac1{r-t}\right)
=-\frac1h\left(\frac1{s-r}+\frac1{r-t}\right).
\]

No documento falta o sinal menos antes da última expressão. Contraexemplo racional: r=3,s=2,t=0, os três valores relevantes são, respectivamente, `1/3`, `1/3`, `-1/3` na fórmula impressa. **Importante:** ao elevar ao quadrado, esse sinal desaparece; a equação (7), a decomposição correta (8) e o restante da prova algébrica podem ser preservados. A decomposição (8) é

\[
\frac1{(r-s)^2(r-t)^2}=\frac1{h^2}\left(\frac1{(r-s)^2}+\frac1{(r-t)^2}\right)
+\frac2{h^3}\left(\frac1{r-t}-\frac1{r-s}\right).
\]

## 2. Página 1, “propriedades imediatas de y” — desigualdade falsa

A afirmação `δ_r^(y)>=δ_r^(x) para todo r` é **falsa** depois de substituir os pontos exteriores por progressões aritméticas. Exemplo: pontos originais `x=(0,1,101,201)`, índices centrais t=0,s=1. Então `δ_s=1` e o primeiro ponto exterior à direita possui `δ_2^(x)=100`, ao passo que a substituição regular `y=(0,1,2,3)` possui `δ_2^(y)=1`.

O que é necessário para a prova não é monotonicidade das folgas ponto a ponto; é a **comparação das somas inteiras das caudas exteriores**, que pode ser obtida pelo teorema P89-T do repositório: para o núcleo decrescente e convexo `φ_h(z)=1/[z²(z+h)²]`, cada cauda original é majorada por sua cauda uniforme de passo `δ_s` ou `δ_t`. Esse reparo depende da demonstração independente de P89; não deve ser apresentado como consequência da falsa desigualdade de folgas.

## 3. Página 4, equação (21) — falta parte positiva

Ponhamos `A=N+t-s>0` e `B=(y_(s+L)-y_t)/(y_(t+1)-y_t)`. A desigualdade impressa compara uma soma `sum_(A<m<=B) 1/m` com `(B-A)/A`. Quando B<A, a soma é **vazia (zero)** e o lado direito é negativo; portanto a desigualdade como escrita é falsa.

Contraexemplo explícito compatível com a construção: `t=0,s=2`, pontos centrais `y_t=0,y_1=100,y_s=101`, ambos os prolongamentos externos com passo `δ_s=δ_t=1`, `N=500`. Daí `y_(s-N)=-498`, `L=floor((101-(-498))/1)=599`, `A=498`, `B=(101+599)/100=7`. A equação (21) diz falsamente `0<=-491/498`.

**Correção:**
\[
\sum_{A<m\le B}\frac1m\le\frac{(B-A)_+}{A},\qquad (z)_+=\max(0,z).
\]

Como a geometria do prolongamento fornece
\[
B-A\le 2h/\delta_t,\quad h=y_s-y_t,
\]
obtemos sempre a cota necessária à etapa seguinte:
\[
\frac{(B-A)_+}{A}\le\frac{2h}{A\delta_t}.
\]
Logo a conclusão assintótica de (23) pode ser preservada com o reparo.

## 4. Parte da demonstração que pode ser aproveitada

Assumindo a **cota quantitativa do Lema 5** (equação (1) do PDF), e após substituir os pontos exteriores por caudas uniformes com a justificação independente de P89, as comparações reparadas fornecem, no limite do truncamento:

\[
\sum_{r\ne s,t}\delta_r\left(\frac1{y_r-y_t}-\frac1{y_r-y_s}\right)
\le -\frac{\delta_s+\delta_t}{h},\qquad t<s,\quad h=y_s-y_t.
\tag{A}
\]

Trata-se de uma desigualdade de **compensação de primeira ordem**, que preserva o termo negativo da expansão por frações parciais. A identidade algébrica (8) então implica

\[
\sum_{r\ne s,t}\frac{\delta_r}{(y_r-y_s)^2(y_r-y_t)^2}
\le\frac1{h^2}\sum_{r\ne s,t}\delta_r\left[\frac1{(y_r-y_s)^2}+\frac1{(y_r-y_t)^2}\right]
-\frac{2(\delta_s+\delta_t)}{h^4}.
\]

As cotas de energia de um índice (P87/Lema 1), removendo o termo cruzado entre s,t, implicam:

\[
\sum_{r\ne s,t}\frac{\delta_r}{(y_r-y_s)^2}\le\frac{\pi^2}{3\delta_s}-\frac{\delta_t}{h^2}
\]

e uma fórmula simétrica para t. Assim se recupera o Lema 6 afiado **condicionalmente à etapa (A)**:

\[
\sum_{r\ne s,t}\frac{\delta_r}{(y_r-y_s)^2(y_r-y_t)^2}
\le\frac{\pi^2}{3h^2}\left(\frac1{\delta_s}+\frac1{\delta_t}\right)
-\frac{3(\delta_s+\delta_t)}{h^4}.
\]

A comparação de caudas de P89 transfere a conclusão da configuração estendida y para x. **Essa passagem não prova sozinha o Lema 5**, pois a inequação variacional quantitativa (1) é tomada como hipótese; portanto a reconstrução integral independente ainda está aberta.

## 5. Próxima tarefa prioritária

Extrair da fonte primária e reconstruir de modo independente o **Lema 5 com erro uniforme explícito**, incluindo as convenções dos pontos de fronteira e os índices que determinam o último espaçamento. Priorizar uma demonstração por variações locais das lacunas, `min(g_i,g_(i+1))`, somas harmônicas e cancelamentos. Depois, verificar se a prova do pdf, com os três reparos, fecha a ponte L5→L6 sem apelos não demonstrados.

**Certificação:** `verificar_reconstrucao.py` usa `Fraction` para os erros concretos em (6), monotonicidade das folgas e (21), e para conferir exemplos do reparo da cauda. A matemática universal remanescente é analítica e condicional; nenhuma amostragem valida os Lemas 5 ou 6 gerais.