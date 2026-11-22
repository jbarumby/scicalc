# NUMERIA P85 — Resultado de GitHub Actions: Hilbert real e raízes b4*

**Execução confirmada:** [GitHub Actions #45039702283](https://github.com/jbarumby/scicalc/actions/runs/45039702283), job `source-and-roots` concluído **success** em 22/11/2026.

**Artefato:** [numeria-p85-hilbert-b4star-intervals](https://github.com/jbarumby/scicalc/actions/runs/45039702283/artifacts/14761515605), id `14761515605`; arquivos:
- `p85-hilbert-source.json`;
- `p85-b4star-192.json`, `p85-b4star-256.json`, `p85-b4star-384.json`;
- `p85-SHA256SUMS.txt`.

**SHA-256 dos certificados:**
- Hilbert: `0600bd7808cb6b7c286824df13eedce90df350a2638efad3681890faf5f32ab6`;
- raízes 192 bits: `bbce683516fcb6c093f0f3df7e8ef08f26a62bd4e483fccba8bb13359a1ac948`;
- raízes 256 bits: `349ab4a54a7d17769c5205674026e05fd0c724c383989dba37f28256cf63df1a`;
- raízes 384 bits: `94c4423093eb0e64f78b936003264de2e6ccd87d35aa6c52175dbfdcb3ebe986`.

## Porta A — transferência de Hilbert

O Lema 6.5 de Ramaré afirma a estimativa ponderada
`int_0^T |sum u_n n^(it)|² dt <= sum |u_n|²(T+pi*m0(n+1))` para coeficientes **reais**, onde `m0=sqrt(1+(2/3)*sqrt(6/5))`. A fonte mais profunda é a desigualdade ponderada de **Preissmann (1984)**, e esse teorema é **ainda um insumo externo**.

**Algebra independentemente reconstruída:** (i) integração de `exp(it d)` exige denominador `i d` — o `i` não está impresso no display do artigo de Ramaré; (ii) cancelamento exato dos termos de fronteira em `t=0` por simetria quando `u_n` são reais; (iii) `log(1+1/n)>=1/(n+1)` por comparação da integral; (iv) limite `N→∞` aplicado aos coeficientes de `f_X`, pois `sigma>1` garante convergência absoluta e somabilidade dos pesos. Também `263/200<m0<329/250` é verificado sem floats.

**Testes:** seis testes unitários aprovados em Python normal e `-O`; quatro exemplos independentes de antissimetria validados com `fractions.Fraction`. A transição de Preissmann para `F_X` está justificada; não se afirma uma nova demonstração independente da desigualdade ponderada geral.

## Porta B — limiares b4* com Arb sem ponto flutuante

**42 limites de raiz com sinais rigorosos**: 7 linhas × 2 impressões de B1 × 3 precisões `workprec` = **192, 256, 384 bits**. A quadratura `K1` da P79/P77 foi recomputada, e a fatoração `B1(b4)=D(αb4+β)^e` foi reconstruída algébrica e intervalarmente com `α,β,D,e>0` e `K3>1`.

Para cada limiar impresso `P`, o programa verificou com intervalos Arb:

```
B1(b4*_lower).upper < P < B1(b4*_upper).lower
```

e, em **todas** as sete linhas,

```
b4*_old.upper < b4_public < b4*_robust.lower.
```

As raízes do cálculo com **256 bits** são cercadas pelos seguintes pares **racionais exatos** de 10 casas decimais:

| Linha | b4 público | antigo: raiz em [inf, sup] | S2 robusto: raiz em [inf, sup] |
|---:|---:|---:|---:|
| 10 | 0.557 | [0.5404552002, 0.5404552005] | [0.5573302735, 0.5573302738] |
| 11 | 0.557 | [0.5404084366, 0.5404084369] | [0.5575138485, 0.5575138488] |
| 12 | 0.558 | [0.5403577756, 0.5403577759] | [0.5583922273, 0.5583922276] |
| 13 | 0.559 | [0.5400253328, 0.5400253331] | [0.5592179943, 0.5592179946] |
| 14 | 0.560 | [0.5401486588, 0.5401486591] | [0.5602007869, 0.5602007872] |
| 15 | 0.561 | [0.5400784332, 0.5400784335] | [0.5610649509, 0.5610649512] |
| 16 | 0.606 | [0.5400681104, 0.5400681107] | [0.6060224697, 0.6060224700] |

São **raízes do majorante explícito B1**, não limiares de validade da Hipótese de Riemann. A segurança do cálculo é a aritmética dirigida Arb: o texto decimal acima é apenas a apresentação de pares racionais cujos sinais foram comprovados; o certificado JSON inclui os diádicos originais e os sinais comparados.

O gate de 192/256/384 passou com **42/42** intervalos e intersecção consistente de cada trio de precisão. Nenhum dos resultados float reivindicados no DOC1/S3 foi usado como limite de prova.

## Conclusão e pendências

**Fechado nesta P85:** fator `pi*m0` para a passagem Hilbert→F_X *condicional à desigualdade de Preissmann*, verificação da fatoração, monotonicidade de `B1(b4)` e todos os limiares `b4*` de substituição das sete linhas do S2, agora com intervalos rigorosos.

**Ainda não fechado:**
- demonstrar do zero a desigualdade ponderada de Preissmann / Montgomery–Vaughan;
- demonstrar independentemente os limites Ramaré–Zúñiga de `S_ε(X)`, incluindo uma fonte pública que permita `3×10¹⁰≤X<3×10¹²`;
- o argumento global sobre `G'(ε)<0` e os alegados 28 gates adicionais zeta do DOC2.

A janela de Suzuki P84 até `t=79` permanece **finita e condicional**, com a TABLE robusta preservada; não constitui prova da RH global. Não alterar certificados históricos, não fazer merge automático.
