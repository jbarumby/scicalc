# NUMERIA–RIEMANN: reprodução e comparação dos certificados (N = 6, 8, 10)

**Data:** 2026-10-21. Parâmetros: L = 0,45, R = 160, primos ativos {2}. Dados brutos: `results/` (192 bits, tolerância de integração 2⁻⁵⁰) e `results/reproduction-256bits/` (256 bits, tolerância 2⁻⁷⁰). Resumo legível por máquina: `results/verification-summary.json`.

## O que foi executado

1. **Verificador racional exato** (`verify_certificate.py`) em cada um dos seis relatórios: lê só os extremos decimais exportados (arredondados para fora), refaz o teste tipo Rump em `Fraction`, sem Arb nem ponto flutuante.
2. **Comparação 192 × 256 bits** (`--compare`): todos os pares de enclosures correspondentes devem se intersectar; medidas também a maior diferença entre pontos médios e as larguras.
3. **Checagem independente em ponto flutuante** (`check_enclosures_float.py`, numpy/scipy): recalcula a matriz projetada por outro caminho (Gauss–Legendre composta, `digamma` complexo do scipy, overlaps e restrições por quadratura direta, sem fórmulas fechadas nem Arb) e confere se cada valor cai no enclosure exportado.

## Resultados

| N | Bits | LDL | Rump | Verificador exato | δ | Menor pivô deslocado | Checagem float |
|---|---:|---|---|---|---:|---:|---|
| 6  | 192 | passou | passou | verificado | 0,012676 | 0,21892 | dentro |
| 6  | 256 | passou | passou | verificado | 0,012676 | 0,21892 | dentro |
| 8  | 192 | passou | passou | verificado | 0,027929 | 0,20366 | dentro |
| 8  | 256 | passou | passou | verificado | 0,027929 | 0,20366 | dentro |
| 10 | 192 | passou | passou | verificado | 0,051784 | 0,17981 | dentro |
| 10 | 256 | passou | passou | verificado | 0,051784 | 0,17981 | dentro |

Comparação entre 192 e 256 bits:

| N | Enclosures se intersectam (ambos os sentidos) | Maior \|Δ ponto médio\| | Maior largura de entrada (192 / 256) |
|---|---|---:|---|
| 6  | sim | 9,2 × 10⁻³⁰ | 1,448 × 10⁻² / 1,448 × 10⁻² |
| 8  | sim | 2,0 × 10⁻²⁹ | 2,354 × 10⁻² / 2,354 × 10⁻² |
| 10 | sim | 3,5 × 10⁻²⁹ | 3,523 × 10⁻² / 3,523 × 10⁻² |

A checagem em double concorda com os pontos médios do Arb até 4,4 × 10⁻¹⁵ (N = 6, 8) e 7,3 × 10⁻¹⁵ (N = 10). O menor autovalor do ponto médio em double (0,1326; 0,1053; 0,1030) coincide com o diagnóstico do relatório.

## Leitura dos números

- A largura das entradas é a mesma a 192 e a 256 bits: ela é dominada pela **cauda analítica** (raio ∝ R⁻³), não pela precisão aritmética nem pela tolerância do integrador. Para ganhar margem em N maior, aumente `--R` ou use o refinamento do Lema 4 (`TAIL_PROOF.md`, §7).
- A margem do teste é folgada: δ é de 6 % a 29 % do menor autovalor do ponto médio.

## O que esta reprodução não é

- É **o mesmo código** rodado em outra precisão e tolerância, na mesma máquina e com a mesma versão do `python-flint` (0.9.0). Não é reprodução em ambiente independente, nem com outra versão do Arb, nem com uma segunda implementação de aritmética intervalar.
- O verificador exato valida o **passo final** (enclosure exportado ⇒ positividade). Que os enclosures contenham a matriz verdadeira depende do Arb, de `TAIL_PROOF.md` e da derivação em `WEIL_FORM_DERIVATION_PT.md`.
- A checagem em double é evidência, não prova: ela detecta erros de álgebra, sinal e normalização, mas não certifica nenhum intervalo.
- Os relatórios antigos (2026-10-19 e 2026-10-20) guardam intervalos como texto com ~3 dígitos; não dá para compará-los rigorosamente com os atuais. Só se verificou que os limites inferiores antigos ficam abaixo dos novos.

## Para encerrar a auditoria faltam

1. Rodar o CI no GitHub e comparar com estes números (ambiente independente).
2. Uma segunda implementação intervalar (por exemplo, `mpmath.iv` ou outra versão do Arb) para os enclosures da matriz.
3. Revisão humana independente de `TAIL_PROOF.md`, de `WEIL_FORM_DERIVATION_PT.md` e da leitura das hipóteses do teorema de Guinand–Weil citado.

Nada aqui demonstra a hipótese de Riemann nem a positividade global de Weil; o certificado vale só para esta família finita.
