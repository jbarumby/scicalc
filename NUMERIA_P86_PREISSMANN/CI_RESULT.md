# NUMERIA P86 — Resultado confirmado: Preissmann, matriz de Hilbert e LDL racional

**CI concluída:** [GitHub Actions #45040745790](https://github.com/jbarumby/scicalc/actions/runs/45040745790) — **completed/success**, 22/11/2026.

**Artefato:** [numeria-p86-preissmann-spectral-rational](https://github.com/jbarumby/scicalc/actions/runs/45040745790/artifacts/14761039785), com `p86-spectral-ledger.json` e `p86-SHA256SUMS.txt`.

**SHA-256 do JSON produzido pelo runner:**
`75c1c1e5ca163b3adf618a4c12fe2da57b0ea5f2187009d572adb1952dabb7ff`.

## Progresso matemático efetivo

**A. Prova analítica uniforme em dimensões pequenas (sem Preissmann):** para qualquer `2<=N<=5`, quaisquer `x_i` reais distintos, `||B||_2<=N−1<=4<pi*m0`, onde `B_ij=sqrt(delta_i delta_j)/(x_i−x_j)`. A prova geral usa apenas `delta_i<=|x_i−x_j|` e `||B||_2<=sqrt(||B||_1||B||_infty)`. **Não é mera extrapolação dos testes finitos.**

**B. Prova analítica uniforme do controle de energia diagonal com constante 4 (sem Preissmann):** para qualquer `N>=2`,

```
Q_i := sum_(j!=i) delta_i delta_j/(x_i-x_j)^2 <= 4.
```

A prova geral está no [README](README.md), §3a: em cada lado de `x_i`, o vizinho mais próximo contribui no máximo 1; as demais parcelas têm majorantes telescópicos com soma no máximo 1. **Esta constante é mais fraca que a cota pi²/3 usada por Preissmann**, e por isso não encerra seu Lema 1 na precisão necessária para m0.

**C. Certificados espectrais para configurações finitas explícitas:** identidade
`H(rho)=rho D^(-1)+K D K=D^(-1/2)(rho I+B²)D^(-1/2)`, de entradas racionais, e `||B||²<rho iff H(rho)>0`. O algoritmo `Fraction` produziu e reconstruiu por LDL todos os pivôs exatos de **oito conjuntos explícitos** de tamanhos `[2,3,5,9,12,12,6,7]`. Foi usado o limite mais exigente
`rho=((157/50)*(263/200))² < (pi*m0)²`.
Esses resultados são verdadeiros **para as oito configurações certificadas**, não para toda dimensão nem todos os espaçamentos.

## Testes e garantias

- **8 testes unitários aprovados** em Python normal e **8** sob `python -O`.
- Auditoria alternativa no workflow reconstituiu **todos os pivôs LDL**, as energias diagonais `Q_i` e a igualdade matricial em frações exatas.
- Entradas adversariais: pontos repetidos ou insuficientes geram erro; um `rho` intencionalmente pequeno produz falha de positividade.
- Nenhum `float`, heurística numérica ou tolerância é usado nas decisões dos certificados finitos.
- Os dois teoremas uniformes novos (casos N<=5, Q_i<=4 para todos N) estão **demonstrados por argumentos gerais**, e não por amostras.

## Reconstrução do texto primário e lacunas

Fonte: [Preissmann, L'Enseignement Mathématique 30 (1984), pp. 95–113](https://www.e-periodica.ch/cntmng?pid=ens-001:1984:30::29). A digitalização permite localizar, mas tem equações com OCR parcialmente degradado; evitar aprovar uma fórmula ilegível por inferência.

- Lema 1, pp. 96–99: maximização por pontos equidistantes e limite mais forte `Q_i<=pi²/3` **não provados independentemente**.
- Lemas 2–5, pp. 99–108: passos variacionais e maximização de padrões de espaçamento **não reconstruídos integralmente**.
- Lema 6, pp. 108–110: desigualdade de três denominadores e termo cruzado **não reconstruídos integralmente**.
- Prova final, pp. 111–113: a parte algébrica da matriz de Hilbert e seu quadrado foi auditada; o **limite uniforme para qualquer N** continua a depender dos lemas anteriores.

A fonte contemporânea [Carneiro–Littmann (2023)](https://ems.press/content/serial-article-files/44548), p. 3, confirma explicitamente o `m0=sqrt(1+(2/3)sqrt(6/5))` de Preissmann; seu limite alternativo `2*pi` é mais fraco e **não** substituiria `pi*m0` na P84.

## Conclusão e próximo gate

**P86 fechou o núcleo de álgebra espectral e avançou dois controles gerais independentes**, mas **não** provou a desigualdade ponderada geral de Preissmann. A cota `F_X` das etapas anteriores continua condicional a esse teorema externo publicado.

**P87 prioritária:** demonstrar a cota de energia mais forte `Q_i<=pi²/3` diretamente do Lema 1 para qualquer espaçamento, auditando todos os passos variacionais; depois investigar o Lema 6 e a passagem espectral final. A fonte de Ramaré–Zúñiga e a faixa `3e10<=X<3e12` permanecem abertas.

Sem nova extensão da cobertura de Suzuki e sem prova global da Hipótese de Riemann.
