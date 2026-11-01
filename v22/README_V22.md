# NUMERIA–RIEMANN: Fase V22 — cotejo analítico e quadratura

**Status científico:**
- V22-A: aritmética do erro de quadratura B3 conferida, com uma correção de arredondamento; **majorantes de segunda derivada permanecem condicionais**.
- V22-B: contagem exata de pontos e Frobenius de curva elíptica sobre F5 e F25.
- V22-C: hipóteses de suporte/suavidade de Burnol conferidas; **aplicabilidade à função-teste NUMERIA não demonstrada**.
- V22-D: exemplo racional finito de núcleo screw positivo; **não é o núcleo que caracteriza RH**.

## V22-A — principal lacuna numérica

O revisor independente lê diretamente o ZIP E7 histórico (SHA-256 fixo), carrega os resultados corrigidos de B3 para G11/G22 e recompõe, com Fraction:

    |erro_ponto_medio| <= M * (b-a)^3/(24*n^2),

desde que M seja um limite válido para sup |f''| em toda a banda.

**Achado:** os seis erros de banda de G11 envolvem as expressões exatas. Na banda **G22 r∈[96,160]**, o erro historicamente exportado fica **abaixo** da cota reavaliada por aproximadamente 2,39×10^-46. A nova cota total de G22 adiciona exatamente **1/10^45**, arredondando essa banda para cima. A correção de E2, que exige soma dos erros de banda para cima, é preservada.

O revisor valida 11 bandas contínuas, 9524+6230=15754 registros B5. **Não valida independentemente a desigualdade M≥sup|f''|**: trata os M históricos como hipóteses e expõe essa lacuna. O teste também rejeita E7 adulterado.

## V22-B — Frobenius como um MODELO provado

E:y²=x³−x sobre F5, com F25=F5[u]/(u²−2). Contagem exata:
 #E(F5)=8, #E(F25)=32.
Daí a_5=5+1−8=−2, e o polinômio X²+2X+5 possui raízes −1±2i, módulo √5. Esse exemplo ilustra a RH para curvas sobre corpos finitos tratada por Katz–Sarnak; **não é a zeta de Riemann clássica**.

## V22-C — Burnol e a hipótese de suavidade

J.-F. Burnol (2000), *Sur les formules explicites I: analyse invariante*, C. R. Acad. Sci. Paris Série I 331, 423–428, DOI: https://doi.org/10.1016/S0764-4442(00)01687-6. Arquivo original SHA-256:
 e110e9c541d7727647de91bc687725548e618e2211306e4f0397136350c96823.

Se g tem suporte logarítmico em [-L,L], com L=9/20, a autocorrelação multiplicativa tem suporte em [-2L,2L]. Verificação exata:
    L=0.45 < log2 < 2L=0.90 < log3.
Assim, **o termo primo 2 entra**; não estamos na região suficientemente curta livre de primos. O Teorema 3 de Burnol garante positividade para algum suporte pequeno, não certifica o valor L=0.45.

Há outra obrigação: a extensão por zero de modos seno deixa salto na derivada nas bordas do suporte, logo **não é C∞ de suporte compacto** como as funções do teorema de Burnol. É necessário demonstrar aproximação por suavização com cota de erro e compatibilidade de formas.

O artigo constrói o operador condutor local Hν=log|x|ν+log|y|ν (onde y é a variável de Fourier), com termos finitos e arquimedianos da fórmula explícita, não um operador de Hilbert–Pólya com espectro dos zeros.

## V22-D — novo critério de positividade

K. Matsumoto e M. Suzuki (2026), *M-functions and screw functions: Applications to Goldbach's problem and zeros of the Riemann zeta-function*, Journal of Number Theory 280, 918–946; DOI https://doi.org/10.1016/j.jnt.2025.09.013.
**Teorema 1.3:** a função g_H1(t)=H1(exp(t))−H1(1) é uma screw function em R se e somente se RH vale. Sua hipótese é uma positividade de núcleo em **todos** os números de pontos e **todos** os tempos reais, além da função específica H1.

O módulo de teste verifica somente um **núcleo auxiliar ilustrativo** g(t)=Σ aω(cos(ωt)−1), aω>0, em tempos múltiplos de π/2. Os 31 menores principais são não negativos por aritmética exata; uma mutação com peso negativo é rejeitada. Nenhum dado de H1(zeta) é testado. Uma quantidade finita de amostras PSD **não demonstra** RH.

## Outros artigos examinados do ZIP

| Fonte | Relevância e limite |
|---|---|
| Debrouwere (2021), representações sequenciais de espaços de Fréchet de funções inteiras | Média para fundamentar funções-teste e espaços de Xian-Jin Li (2015); não demonstra sua pertinência automática aos modos seno. |
| Artigo sobre nuclearidade e Grothendieck–Lidskii em espaços quaterniônicos (2024) | Contexto funcional; trata operadores quaterniônicos, não a identidade de traço de Weil. |
| Mehta, Saha, Viswanadham (2016), funções zeta múltiplas | Referência de continuação analítica; não é diretamente a zeta clássica. |
| Kirila (2021), combinações de funções L | Contraexemplos e advertências na região Re(s)>1; não substitui positividade Weil. |
| Hochfilzer e Oliver (2022), razões de funções L de Artin | Aspectos de cancelamento de zeros; ligação indireta. |

ZIP original recebido SHA-256:
 7f55b110e5953b3f8c7cad01c1d49756aff6a7bb66deadefa8ae3311e8c8ac56.

PDF de Matsumoto–Suzuki SHA-256 (dentro do ZIP):
 b8ca61e522f32deb5cfc615510a9a2dd122f038e617a1b9d1ad80a73de443342.

## Executar

    python v22/AUDITAR_V22_QUADRATURA.py --e7 NUMERIA_E7_REVISAO_FORMAL_2026-10-26_3.zip --out v22/RESULTADO_V22_QUADRATURA.json
    python v22/AUDITAR_V22_FROBENIUS.py --out v22/RESULTADO_V22_FROBENIUS.json
    python v22/AUDITAR_V22_SUPORTE.py --out v22/RESULTADO_V22_SUPORTE.json
    python v22/AUDITAR_V22_SCREW.py --out v22/RESULTADO_V22_SCREW.json
    python v22/TESTAR_V22.py

O workflow executa a auditoria, repete os JSONs byte a byte, confere os hashes fixados quando presentes e publica testemunhas. Não introduz arquivos PDF integrais potencialmente protegidos ao GitHub.

**Próximas obrigações:** provar os majorantes M de segunda derivada por intervalos independentes; construir uma função-teste suave ou demonstrar extensão adequada da fórmula explícita; cotejar as normalizações local/arq/primos de Burnol e Li (2015/2023). Nenhum resultado V22 prova a positividade global da forma de Weil ou RH.
