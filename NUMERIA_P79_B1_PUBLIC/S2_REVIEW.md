# Revisão P79-S2 — três arquivos recebidos (22/11/2026)

**Escopo:** avaliação de material complementar. Não altera os certificados originais P77–P79, não reivindica extensão nova e não substitui o gate de fonte externa ainda em aberto.

## Arquivos preservados literalmente

Em `prototype_input/`:
- `gemini-code-1795353291902.py` — SHA-256 `d6edcd2b02b278f17999a34102a33f2733b7231c08f0e3984676c9dc160d912f`;
- `p78_b1_resolution.py` — SHA-256 `10be6eb754ec5a834b97cb1654533f53e77000722e05e7a2d600ac9e1e6ca42d`;
- `p78_b1_resolution_results.json` — SHA-256 `c8f568359a72252aa3e18f52a470a4685606f9f0c861259a8adf5e4de6ef665d`.

A integridade dos três arquivos é verificada pelo teste `test_s2_exact.py`. O terceiro arquivo **declara `mode=STANDALONE_RECONSTRUCTION`**, portanto seus valores de ajuste e de limiar `b4_star` em float NÃO são certificados Arb.

## Ganho matemático independente: monotonicidade exata em b4

No auditor da P77, com os demais parâmetros fixos:

`C3 = (pi*m0/(2*delta))*exp(2*delta*gamma/L)*b4`;
`C4 = b4/(5*delta*exp(delta)) + 0.605*exp(-2*delta)/L`.

Assim, `K3 = alpha*b4 + beta`, com `alpha>0`, `beta>0`. Todos os parâmetros da tabela P78 satisfazem `delta<=303/1000` e `b4>=27/50`. Como `pi>3`, `m0>1`, `gamma>0` e `L>0`, obtemos **analiticamente**:

`K3 > C3 > 3*(27/50)/(2*(303/1000)) = 270/101 > 1`.

Portanto, os fatores `max(1,K3^{-r})` nos canais L1/L2/comum valem exatamente 1, pois `r>0`. Como `sigma1<=sigma2` e `e(sigma)=(2sigma-1)/(2-sigma)` é crescente (`e'(sigma)=3/(2-sigma)^2>0`), o canal positivo restante é `K3^e` para `e=e(sigma2)`.

Logo, **nesta estrutura algébrica da P77**,

`B1(b4) = D*(alpha*b4+beta)^e`, `D>0` independente de `b4`, e
`dB1/db4 = B1*e*alpha/(alpha*b4+beta)>0` nas faixas auditadas.

Para as linhas 10–16, `sigma2=(linha+16)/32`, `e=2*linha/(48-linha)`, exatos:

| linha | e exato |
|---|---|
| 10 | 10/19 |
| 11 | 22/37 |
| 12 | 2/3 |
| 13 | 26/35 |
| 14 | 14/17 |
| 15 | 10/11 |
| 16 | 1 |

O teste com `fractions.Fraction` verifica esses expoentes, a desigualdade racional `270/101>1` e os coeficientes correspondentes da tabela P78, sem aproximar `sigma` por float.

## Reprints de duas casas: válidos como majorantes alternativos

A tabela P79 já certificou em Arb 192/256 bits os seguintes **limites superiores estritos**:

| linha | B1 P79 em 6 casas | S2 2 casas | Diferença exata S2 − P79 |
|---|---:|---:|---:|
| 10 | 17.484560 | 17.49 | 0.005440 |
| 11 | 19.659253 | 19.66 | 0.000747 |
| 12 | 22.239617 | 22.24 | 0.000383 |
| 13 | 25.472651 | 25.48 | 0.007349 |
| 14 | 29.881218 | 29.89 | 0.008782 |
| 15 | 36.586170 | 36.59 | 0.003830 |
| 16 | 51.648097 | 51.65 | 0.001903 |

Logo, **por transitividade**, as propostas de duas casas também são majorantes rigorosos, *assumindo os gates e arredondamentos do P79*. Mas são **menos precisas**: não substituir silenciosamente a tabela P79 na execução que obteve `18<=t<=79` com 512 faixas. Essa cobertura não foi reexecutada com a tabela S2 de duas casas.

O título `reprint_minimal` do protótipo S2 indica mínimo sob ajuste float; a comparação acima prova que as propostas dominam os majorantes P79, **não** a minimalidade rigorosa. Prova de minimalidade exigiria também um limite inferior certificado de B1.

## Auditoria do arquivo Gemini

Contribuições aproveitáveis:
- explicita uma cauda `D = (1+c/(2L))^(4/3) * (1+c1*sqrt(c)+c2)^(2/3)`, que pode ser comprovada a partir da integranda de K1 para `u>=64` (com fatores positivos). Essa escolha tende a ser mais ajustada que a forma conservadora `(1+c/(2L))^2*(1+c1*sqrt(c)+c2)`.
- explicita as decisões `strict` em Arb e separa a exportação de decimais para diagnóstico.

Pontos que **não** devem entrar no pipeline como prova sem revisão:
- o cabeçalho menciona `ARF_RND_CEIL/ARF_RND_FLOOR`, mas a implementação usa `float(x.str(17))` e `math.nextafter`: isso não equivale a fornecer extremos diádicos exatos com direção de arredondamento formalmente demonstrada;
- `arb_bounds_str` usa `.str(digits, radius=False)` e rotula os extremos como rigorosos, mas a formatação decimal pode arredondar para qualquer lado; não reutilizar a string como intervalo certificador;
- `ceil_2dp_strict` inicializa via `int(float(...))` e termina com `f"{n/100:.2f}"`; preferir a aritmética de frações com `B1_upper_dyadic` já implementada na P79;
- a variante da integral K1 e as operações `arb.max`/`arb.union` não foram reproduzidas no CI oficial; não substituir módulos anteriores até testar;
- o script S2 captura qualquer exceção do modo `REPO_CERTIFIED` e retorna ao modo `STANDALONE_RECONSTRUCTION`: para gates de prova, **falhar explicitamente** em vez de fazer fallback silencioso.

## Conclusão e próxima prioridade

A S2 dá uma justificativa algébrica simples da monotonicidade de B1 e uma tabela alternativa mais legível, mas **não supera o resultado já obtido pela P79**, que cobre finitamente `18<=t<=79` em 512 faixas sob as fontes importadas.

O ponto de pesquisa que continua faltando é confrontar Ramaré–Zúñiga v1/v3 e **justificar integralmente** a transferência de suas desigualdades à cadeia do Lema 6 de Chourasiya–Simonič, antes de elevar a conclusão matemática. A verificação numérica finita não demonstra RH global.
