# P76 — decomposição diagnóstica de B1: caminho de K2 e K3

**Status:** análise numérica exploratória `mpmath` (70 dígitos), **não é certificado Arb** e **não demonstra B1**.

Fonte: [Chourasiya–Simonič arXiv:2507.15184v2, Teorema 1 e Eqs. (35)–(40)](https://arxiv.org/html/2507.15184v2).

Com `T0=H0=H=3×10^12`, `h=1`, foi avaliado o ramo de `A2`, que não depende da integral de `K1`:

`K3=C3+C4*(Gamma(3/2)/sqrt(2*kappa)+pi*m0/T0)` (pela integração exata de Eq. 36);

`K2=2^(1/3)*sqrt(3/(4*kappa))*Gamma(3/2)+2^(1/3)*(4*kappa/3)*(A0/T0)^3` (Eq. 38).

Em seguida aplicou-se Eq. (35) para `L0`, Eq. (40) para `L2` e a fórmula de `A2` no Teorema 1. A parcela normalizada em Eq. (10) é:

`B1_from_A2 = A2 * (log H)^[-6(1-sigma2)/(2-sigma2)]/(2*pi*d)`.

**Atenção:** para este diagnóstico apenas, considerou-se `A0=10^(log10 A0 impresso)`. A tabela arredonda `log10 A0` em duas casas; portanto o cálculo NÃO garante que o parâmetro real otimizado dos autores seja exatamente esse. A derivação rigorosa deverá usar A0 com proveniência e arredondamento justificáveis.

| Faixa sigma | B1 impresso | Parcela B1_from_A2 calculada | Percentual aproximado |
|---|---:|---:|---:|
| [0,50000;0,53125] | 5,069 | 0,0071716 | 0,14% |
| [0,62500;0,65625] | 9,387 | 0,0466128 | 0,50% |
| [0,75000;0,78125] | 15,35 | 0,339336 | 2,21% |
| [0,84375;0,87500] | 21,77 | 1,74536 | 8,02% |
| [0,87500;0,90625] | 24,83 | 3,10508 | 12,51% |
| [0,93750;0,96875] | 35,35 | 10,51444 | 29,74% |
| [0,96875;1,00000] | 46,06 | 19,96006 | 43,33% |

**Prioridade matemática:** Enclausurar rigorosamente todos os componentes `K1,K2,K3,L0,L1,L2,A1,A2` em Arb e verificar os 16 `B1`; a integral de `K1` precisa de decomposição em intervalo finito e cauda explícita (não apenas `mpmath.quad`). Separadamente, resolver a ambiguidade do parâmetro A0 tabulado. Não usar as percentagens acima para concluir positividade ou atualizar o certificado de t=79 sem nova execução.

A prova computacional validada nesta branch continua restrita a **B2/B3 em Arb 192/256 bits**. Consulte [CI_B23_RESULT.md](CI_B23_RESULT.md).
