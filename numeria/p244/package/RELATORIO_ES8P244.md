# NUMERIA–RIEMANN — ES8-P244
## Auditoria analítica independente e repetição do certificado MPFR P242

**Data:** 18/11/2026. **Janela:** L = 161/200 = 0,805, T♯ = 504, 377 modos em cada paridade.

## 1. Resultado principal da revisão dos arquivos enviados

Recompilei o código C fornecido `p242_mpfr_interval_ldl.c`, exportei novamente os 142.506 centros diádicos a partir do arquivo P237 autenticado e executei as duas fatorações LDL e inversões intervalares com MPFR em 1024 bits. O **replay produziu exatamente os mesmos bytes** dos dois arquivos enviados:

| Setor | Pivôs positivos | SHA-256 da saída repetida | ν = limite superior para ‖C⁻¹‖²_F (valor apenas ilustrativo) |
|---|---:|---|---:|
| par | 377/377 | `0080384607c00fd63d3b4c481c17555118791b36b19372b91ea47cb40984fd2b` | 4,3561885896 × 10³⁴ |
| ímpar | 377/377 | `d61f312cd24c3f2b518692798022ea90dfcc81b0f4d42d9fad321d16dd1ba5a5` | 3,6842757790 × 10²⁸ |

O revisor racional `p242_verify_exact.py` voltou a confirmar por frações exatas as somas das entradas dos dois intervalos matriciais e a desigualdade `(δ + 10^-18)^2 ν < 1` em ambas as paridades. **Conclusão finita:** sob a correção do algoritmo C/MPFR, toda matriz real simétrica nas duas caixas P237 montadas é definida positiva com margem estrita >10⁻¹⁸.

### Falha de *fail-closed* encontrada

`p242_verify_exact.py` apenas lê o número `ν` declarado na linha `# inverse_frobenius_squared_upper` da saída MPFR. Ele não calcula ν nem vincula esse número a uma execução de `p242_mpfr_interval_ldl.c`. Prova adversarial: ao substituir `ν` por `1` nos dois arquivos de rastreamento **sem alterar os 377 pivôs**, a rotina ainda anuncia certificado completo. Isso **não** invalida o resultado original, que foi agora repetido pelo programa C; invalida o uso do verificador Python **isolado** como autenticador de qualquer rastreamento arbitrário.

**Correção P244:** `p244_replay_fail_closed.py` confere o SHA dos códigos C e do verificador, autentica os dados P237, extrai os centros, compila e executa novamente MPFR, exige a igualdade byte a byte entre saídas regeneradas e originais, e só então chama o auditor racional. Uma tentativa de alteração de ν é rejeitada. A exigência de igualdade byte a byte é propositalmente conservadora: diferentes versões de MPFR podem produzir intervalos corretos mas não idênticos, caso em que o verificador falha sem declarar um falso sucesso.

## 2. P237 — majorante complexo em elipses e quadratura GL60

A rederivação parte da extensão holomorfa `D(t)=[ψ(1/4+it/2)+ψ(1/4-it/2)]/2`, *não* de `Re ψ` para argumento complexo. Nas elipses E₃ em todos os 1008 painéis, `|Im t|≤1/3`, `|Re t|≤3025/6`, `1/12≤Re z±≤5/12`, `|z±|<253`. A série de ψ usada no Gate C dá `|ψ(z±)|<21`; as estimativas conservadoras dos primos, de log π e de β dão `|Ψ(t)−β|<274/7<50`.

Da identidade Legendre–Bessel e de Rodrigues, `|j_n(Lt)|≤e^{L/3} min(1,|Lt|^n/(2n+1)!!)`; como `L/3<1/2`, o majorante `a_n` da P237 é válido sob esses lemas clássicos. No painel, `|g(s)|<50 a_i a_j/12`. O erro de aproximação Chebyshev para graus ≥120 é `3M/3^120`; o operador erro integral menos GL60 tem norma ≤4 na norma do supremo, de onde decorre por painel `50 a_i a_j/3^120`. Em 1008 painéis, o fator exato é `1008·50/3^120`.

O auditor P244 recalculou **os 754 pesos modais e as duas cotas racionais por soma de linha** sem reutilizar os certificados JSON Gate C. As estimativas resultantes permanecem menores que `4×10⁻⁴⁷` (par e ímpar). Os testes não formalizam os teoremas clássicos importados.

## 3. P239 — fechamento de uma justificativa analítica de `M<20`

A prova escrita em P239 usa Binet no argumento `z=5/4+iy`, com `0≤y≤252`. A estimativa `log|z|<6`, isoladamente, **não** implica `|log z|<6`, já que `z` é complexo. A cota é, no entanto, recuperável e foi explicitamente justificada:

- `|z|<253` e `|arg z|<π/2<11/7`.
- A soma finita de nove termos positivos da série de `exp(28/5)` excede `253`, logo `log|z|<28/5`.
- Assim `|log z|²=(log|z|)²+(arg z)²<(28/5)²+(11/7)²<36`.

Com `|u²+z²|≥(5/4)²` (verifica-se minimizando em `u²≥0`, em dois casos `y≤5/4` e `y≥5/4`), Binet fornece `|Re ψ(1/4+iy)|<6+2/5+4/75+4`. Os limites racionais independentes `A_L<5`, `|log π|<2` e `0<β<1` importado da P237 implicam o majorante

`|Ψ_L(t)−β|<1384/75<20` para `0≤t≤504`.

## 4. P240–P241 — sinal do resto de alta frequência

O Lema 3.1 do preprint de Xuefeng Zhu fornece, para `t≥15/4`, a estimativa `Re ψ(1/4+it/2)−log π ≥ log(t/(2π))−1/t`, sustentando a desigualdade one-stroke. A derivação mais forte do pacote P240 pode ser auditada separadamente, usando a fórmula clássica de Binet em `z=1/4+it/2`, para `t≥504`:

- Para `0≤u≤t/4`, `|u²+z²|≥t²/8` e a parcela da integral de Binet é ≤`2/(3t²)`.
- Para `u≥t/4`, `|u²+z²|≥t/4`, e a cauda exponencial é <`1/(3t²)`; o confronto final é racional e uniforme em `t≥504`.
- `Re log z ≥ log(t/2)` e `Re[1/(2z)]≤1/(2t²)`.

Portanto `Re ψ(1/4+it/2)>log(t/2)−3/(2t²)` e

`Ψ_L(t)−β₅₀₄ > 1/504−3/(2·504²) = 335/169344 > 0` para `t≥504`, assumindo a fórmula de Binet e os coeficientes aritméticos da fonte clássica.

O controle de alta frequência é **uma desigualdade de sinal**, não um resto pequeno em norma.

## 5. Verificação diagnóstica adicional da montagem

Realizei uma quadratura Gauss–Legendre 60×1008 independente em `float64` com NumPy/SciPy para 12 entradas da matriz C de baixa ordem (seis pares e seis ímpares). Todas diferiram dos respectivos centros de P237 por menos de `1,7×10⁻¹⁵`. **Esta comparação não é rigorosa**: não serve para substituir a quadratura certificada Arb, apenas para flagrar possíveis erros grosseiros de integrando/sinal/normalização.

## 6. Estado dos gates após P244

| Gate | Conclusão |
|---|---|
| P242 centros, pivôs e inversa | **Reproduzidos integralmente com C/MPFR**; duas saídas idênticas às enviadas |
| Verificador P242 isolado | **Vulnerabilidade demonstrada**; fluxo rígido corrigido pela P244 |
| P237 geometria, símbolo E3, GL60 | Dedução revisada; desigualdades racionais reexecutadas; teoremas clássicos importados |
| P239 símbolo `M<20` | Justificativa complexa de `|log z|<6` complementada |
| P240–P241 cauda de alta frequência | Sinal revisado contra Binet e o Lema 3.1 de Zhu |
| Montagem integral P237 por segunda biblioteca rigorosa | **Não** reproduzida integralmente; apenas 12 testes numéricos diagnósticos |
| Redução one-stroke / fórmula explícita | Concordância de estrutura com Zhu, mas o teorema externo é importado; não formalizado em Lean |
| Todos os suportes / RH | **Não demonstrados** |

**Fontes:** Zhu, *Weil positivity in compact windows* (arXiv:2608.24827v2), especialmente Teorema 1.1, Lema 3.1, §4, §5; DLMF §5.9 (Binet para digama). Os resultados sobre `L=0,805` dependem da montagem P237 e do limite finito independente, não do certificado do artigo em `L=0,8`.

**Próxima etapa útil (P245):** conferir *rigurosamente* e de forma independente uma amostra estratificada da própria integração GL60 (baixa e alta ordem, duas paridades e painéis de borda) usando uma segunda biblioteca de bolas complexas e, depois, ampliar para a reconstituição de todos os termos primos e polos no nível de entradas. Só chamar a positividade de Weil na janela fixa de certificada após validar a cadeia inteira de inclusões analíticas e computacionais.
