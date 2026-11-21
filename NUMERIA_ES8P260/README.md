# NUMERIA–RIEMANN — ES8-P260 | Ahlfors / digama racional

Base: PR #66 (P259 Gate B), sem modificar a branch principal.

## Resultado

O auditor P260 certifica com **frações exatas**, sem Arb nem ponto flutuante, a parte arquimediana na frequência 504 e os pisos uniformes do símbolo de quatro potências primas:

| Faixa | Piso |
| --- | ---: |
| |t| ≥ 504 | 1/4 — já certificado pela P253, não refeito aqui |
| |t| ≥ 850 | 1/2 — confirmação independente da P258 |
| |t| ≥ 1100 | **3/4** — novo |
| |t| ≥ 1500 | **1** — novo |
| |t| ≥ 2000 | **5/4** — novo |

A referência bibliográfica é Lars V. Ahlfors, *Complex Analysis*, terceira edição, capítulos 4–5 (resíduos, Gama, Stirling, zeta). A íntegra do livro não foi fornecida, e o desenvolvimento presente usa a série clássica de Euler–Weierstrass, com referências NIST/DLMF. [Prova matemática](PROVA_P260.md) · [Fontes](SOURCES.md).

## Reprodução (somente Python padrão)

Execute na raiz do repositório:

    mkdir -p replay
    python -m unittest NUMERIA_ES8P260.test_ahlfors -v
    python -O -m unittest NUMERIA_ES8P260.test_ahlfors -v
    python -m unittest NUMERIA_ES8P260.test_pole_parity -v
    python -O -m unittest NUMERIA_ES8P260.test_pole_parity -v
    python -m NUMERIA_ES8P260.ahlfors_digamma --out replay/P260_RESULTADO.json
    python -m NUMERIA_ES8P260.pole_parity --out replay/P260_POLE.json
    sha256sum replay/P260_RESULTADO.json

O workflow P260 compara os resultados com o JSON de referência e guarda os hashes. **26 testes passaram localmente tanto no Python normal quanto em -O**, incluindo 240 identidades racionais de normalização dos polos sob a hipótese explícita do funcional 2ℓ₊ℓ₋.

## Próximo Gate B

A forma auxiliar reforçada é

    R_P260 = R_P259 + (1/4)(I-G_1100) + (1/4)(I-G_1500) + (1/4)(I-G_2000).

A prioridade é calcular essas três contribuições no **vetor negativo de 57 coordenadas** da P259 com Arb e Rayleigh racional, sem refazer 377×377 prematuramente.

**A P260 não certifica positividade de toda a matriz, não transfere automaticamente o resultado para Weil e não prova RH.**
