# NUMERIA–RIEMANN — ES8-P250: redução adaptativa e corte inteiro ótimo

**Base:** branch `numeria-p249-prime-multiscale`, commit `707185541303393401798f0bb46e1371a3e6a291`.

## Resultado

O Codex P249 registrou 24.960 subintervalos Arb para comprovar o sinal residual na faixa [309,504], com a mesma constante beta da P237. Repetimos a demonstração com uma cadeia aritmética **independente de Arb**, apoiada na biblioteca P245B fixada por SHA-256: frações exatas, seno, cosseno, log e digama com arredondamento dirigido em uma grade de 256 bits.

O algoritmo adaptativo fecha todos os **195 intervalos unitários** entre 309 e 504 em **230 células terminais**, após **265 tentativas** e profundidade máxima **7**. A menor cota racional verificada é 20773768130915242101238847830408510328436019544188775636079480328315169623 / 2^256 ≈ **1.79405763103e-4**.

Além disso, para a frequência racional **t=6169/20=308,45**, o símbolo real exato `Re psi(1/4+it/2) − log(pi) − S(t) − beta` está em um intervalo dirigido estritamente negativo: aproximadamente **−0.3179321079788687**. Logo nenhum corte inteiro ≤308 pode satisfazer a condição ponto a ponto para todas as frequências futuras **com este beta fixo e este símbolo**. Junto à positividade para t≥309, concluímos: **309 é o menor corte inteiro possível sob essa escolha**. Isso não demonstra que 309 seja o melhor corte fracionário nem otimiza beta ou outra decomposição.

## Método

Escreva `S(t)=Σ_q a_q cos(ω_q t)` com `(q,p)=(2,2),(3,3),(4,2),(5,5)`, `a_q=2log(p)/√q`, `ω_q=log(q)`, `A=Σ a_q`, `M2=Σ a_q ω_q²`, `B(t)=log(t/(2π))−1/t` e `β=B(504)−A`.

Em cada `[a,b]`, centro `c`, raio `r`, certificamos

```
B(t)−β−S(t) ≥ B(a)−β−S(c)−r |S'(c)|−(r²/2) M2 > 0.
```

Todos os sinais são verificados por inteiros exatos e intervalos dirigidos. Binet fornece `Re ψ(1/4+it/2)−log π ≥ B(t)` para `t≥309`, sob as identidades matemáticas clássicas da P249. Para `t≥504`, `B(t)≥B(504)` e `S(t)≤A`; a desigualdade é estrita para `t>504`. A amostra negativa usa a avaliação intervalar independente da digama por recorrência e resto de Binet da P245B.

## Reprodução

A partir da raiz do repositório, onde já existe `numeria/p245/package/p245b_certified_gl60.py`:

```bash
python3 -m NUMERIA_ES8P250.adaptive_signed_tail --out /tmp/P250_CERT.json
python3 -m unittest NUMERIA_ES8P250.test_p250 -v
python3 -O -m unittest NUMERIA_ES8P250.test_p250 -v
sha256sum /tmp/P250_CERT.json
```

SHA-256 esperado do certificado gerado: `2bde93ccce25b8d2a1c6e59b49bcce43475fb2c7bc83366edaa1a19b7d4efe38`. O certificado completo com as 230 células integra o ZIP de reprodução da P250. O auditor usa a base P245B já auditada e suas identidades clássicas; não é uma formalização autônoma em Lean.

## Próxima etapa

Não gastar mais recursos procurando um **corte inteiro inferior** com o mesmo beta e a mesma condição de sinal. **Prioridade P251:** reconstruir/certificar o bloco de baixa frequência com 259 modos por paridade e 618 painéis GL60, com avaliação matriz-livre ou blocos Arb. É necessário verificar a margem finita, a identificação Weil–Zhu e a ponte de Schur. **Nenhuma afirmação de positividade global ou RH está demonstrada.**
