# P92-R — parecer complementar sobre auditoria tripla enviada (23/11/2026)

**Objeto:** auditoria dos três pontos A/B/C da reconstrução R1–R7 de Preissmann, materiais recebidos:
\`AUDITORIA_MATEMATICA_INDEPENDENTE_TRES_PONTOS.md\`,
\`p92r_auditoria_tres_pontos.py\`,
\`p92r_auditoria_resultados.json\`,
\`SHA256SUMS (1).txt\` e
\`NUMERIA_P92R_AUDIT_TRIPLE.zip\`.

**Estado:** Parecer favorável à incorporação COMO EVIDÊNCIA DE AUDITORIA E REPRODUÇÃO; não substituir a prova analítica por testes, não reivindicar formalização, publicação revisada por pares ou demonstração da RH. PR permanece **draft**.

## 1. Procedência e reprodução, em ambiente distinto

Os três hashes SHA256 do manifesto da auditoria coincidem byte a byte com os arquivos individuais entregues e também com os incluídos no ZIP. O ZIP passa na verificação de integridade. No pacote anterior \`NUMERIA_P92R_PREISSMANN5_bundle.zip\`, os **9 arquivos** referidos por seu próprio manifesto SHA256 conferiram byte a byte.

Os **10 testes unitários** do pacote P92-R foram reexecutados: 10/10 OK.

A implementação recebida foi copiada para ambiente isolado, com seu caminho absoluto \`BUNDLE = /home/z/my-project/p92r/NUMERIA_P92R_PREISSMANN5\` substituído pelo caminho local de desempacotamento; a parte D (mpmath) foi temporariamente desativada. Resultados reproduzidos com ZERO falhas:
- S: 13 verificações simbólicas;
- A: 4.812 verificações racionais;
- B: 3.960 verificações racionais;
- C: 1.852 verificações racionais;
- **Subtotal efetivamente reexecutado em ambiente isolado: 10.637 verificações.**

O JSON enviado relata ainda 273 verificações multiprecisão em D, com zero falhas. **Essas 273 não foram reexecutadas nesta revisão complementar.**

O relatório anuncia em um trecho 10.911 verificações; a soma aritmética das próprias contagens S + A + B + C + D é **10.910**, não 10.911.

## 2. Achados matematicamente substantivos

**A — Bellman e rampas convexas.** A documentação identifica corretamente a medida positiva de segunda derivada para \`f_X(u)=max(1/u-1/X,0)\`, em particular a massa pontual \`X^-2\` em \`T=X\`; omiti-la invalida a representação do truncamento, que é essencial para o termo explícito do R3. A prova de R1 usa fórmulas exatas em cada ramo e continuidade nos pontos \`h=L/n\`. As verificações simbólicas, e os casos exatos, sustentam a dedução; a universalidade decorre do argumento algébrico e de Tonelli, não da quantidade de exemplos.

**B — Caudas exteriores.** \`phi_h(u)=1/[u²(u+h)²]\` satisfaz
\[
\phi_h''(u)=\frac{2(10u^2+10uh+3h^2)}{u^4(u+h)^4}>0.
\]
Adicionalmente \`phi_h(u),phi_h'(u)->0\` quando \`u->+infty\`. Integração por partes fornece
\[
\phi_h(u)=\int_u^\infty(T-u)\phi_h''(T)\,dT,
\]
mostrando explicitamente a associação à classe R1 de misturas positivas de rampas. Isso justifica usar R2 nas caudas.

**C — Compensador misto.** O relatório confirma que R5 requer o termo \`+2h/(A L)\`, e fornece contraexemplo exato à versão SEM o termo:
\[
g=(1,1),\quad A=L=R=1,\quad h=2;\qquad
\Sigma_t^A-\Sigma_s^A=17/60>-7/15.
\]
Com o reparo, o limite é \`53/15\`. O limite \`A->infty\` deve ser interpretado como convergência da DIFERENÇA PAREADA assinada; as somas individuais \`Sigma_t\` e \`Sigma_s\` não convergem em geral.

## 3. Correções/qualificações adicionais na auditoria recebida

1. **Erro dimensional no microlema C-iii.** O relatório apresenta \`RHS*h² >= (pi²/3)(1/delta_s+1/delta_t)-6 >= 2pi²/3-6\`, que é falso para h arbitrário (por exemplo h=delta_s=delta_t=100). A forma correta, obtida multiplicando pelo fator dimensional \`h³\`, é
\[
h^3\,RHS=\frac{\pi^2}{3}\left(\frac h{\delta_s}+\frac h{\delta_t}\right)
-3\frac{\delta_s+\delta_t}{h}
\ge \frac{2\pi^2}{3}-6>0,
\]
pois \`0<delta_s,delta_t<=h\`. Esse erro não aparece no argumento principal R1–R7; afeta apenas o microlema extra da auditoria.

2. **Teste D não é aritmética intervalar dirigida.** A parte D usa \`mpmath.mp.mpf\` com \`mp.mp.prec=192/256/384\` e comparações numéricas ordinárias; alta precisão não é, por si, intervalo de inclusão rigorosa. Há majorantes analíticos de cauda, mas não certificados de arredondamento em cada operação. Classificar D como **ensaios multiprecisão auxiliares**, não como certificação Arb/MPFI.

3. **Portabilidade do auditor.** \`p92r_auditoria_tres_pontos.py\` usa um caminho local absoluto específico do ambiente original; parametrizar por \`--bundle\`/relativo ao script para que a auditoria possa ser reexecutada diretamente na CI.

4. **Terminologia de B2.** O \`guard_phi\` produzido pelo script é um majorante racional que inclui uma cota por integral; não é a igualdade com a integral exata (em geral ela contém termos logarítmicos). A majorização é justificável porque
\[
\frac1{\delta^2x^2(\delta x+h)^2}\le\frac1{\delta^2h^2x^2};
\]
logo o integral da direita a partir de M é \`1/(delta² h² M)\`, que é menor que o guard implementado. O uso conservativo do guard não invalida os testes.

5. **Limite de escopo.** A verificação independente reproduz identidades-chave e desigualdades racionais, mas não é formalização computacional da prova para todos os reais (Lean/Coq/Isabelle). A cadeia geral continua exigindo leitura da demonstração analítica completa, atenção aos domínios de Tonelli, sinal das somas e limite pareado.

## 4. Esclarecimentos obrigatórios no README antes do merge

- \`Sigma_t^A\` e \`Sigma_s^A\` são somas assinadas de \`delta_r/(y_r-y_i)\), **ambas excluindo os índices t,s**; \`P_t/P_s\` incluem o nó oposto ANTES da remoção. Nos truncamentos os deltas das bordas vêm da extensão bilateral completa, passos L/R.
- \`Sigma_t\` e \`Sigma_s\` isoladamente são geralmente divergentes; só a diferença pareada converge absolutamente (caudas com decaimento \`O(r^-2)\`).
- Incluir explicitamente \`phi_h=\int_u^\infty(T-u)phi_h''(T)dT\` sob as hipóteses \`phi_h,phi_h'->0\`.
- Manter classe de evidência distinta: provas analíticas; verificações racionais/simbólicas; ensaios multiprecisão não intervalares.
- **Não fazer merge automático; não inferir RH.**

**Conclusão:** a auditoria enviada é útil e confirma muitos passos relevantes da P92-R; há reparos documentais pequenos e um erro dimensional em um microlema adicional. Não foi identificado nesta revisão um contraexemplo às desigualdades R1–R7. Isso é evidência favorável, não equivalência a uma checagem formal da prova geral.
