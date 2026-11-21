# NUMERIA–RIEMANN — ES8-P259 · Gate A (piso por faixas e testemunha ímpar)

**Continuação exata do PR #64 / P258.** O auxiliar P258 usa quatro potências primas, corte 504, piso 1/4 e 377 modos por paridade, na meia-largura racional `0.972955074527656652552676371721590`. O setor par está certificado positivo e o setor ímpar possui testemunha negativa com **29 coordenadas diádicas**.

## Resultado da P259

Mantendo o modelo até 850, eleve o piso apenas acima de 850:

```text
Psi4(t) >= 1/4  para 504 <= |t| < 850
Psi4(t) >= 1/2  para |t| >= 850
R_P259 = R_P258 + (1/4)(I-G_850)
```

A segunda desigualdade usa o mesmo limite Binet/P250 de P258, `beta>0`, e a desigualdade racional `log(850/504)>346/677>1/2`.

**Novo resultado rigoroso:** para **a testemunha ímpar particular da P258**, uma prova por 23 períodos completos de Fourier mostra

```text
massa normalizada de frequência >850   > 9.5780e-24
Rayleigh original da testemunha P258     ≈ -1.73007525684909e-24
Rayleigh no novo auxiliar P259            > 6.6443e-25
```

As decisões de sinal usam `fractions.Fraction` (não `float64`). O auditor usa os **extremos racionais exatos** do Rayleigh P258, sem substituir esses números por decimais. No resultado publicado, `float` aparece **somente para apresentar** valores após a decisão exata.

### Método

Escreva cada Bessel esférica `j_n(x)` como combinação racional de `sin(x)/x` e `cos(x)/x`, usando a recorrência com **coeficientes inteiros**. Para a função ímpar formada pelos 29 modos de Legendre, isto dá

```text
F(t) = sqrt(2L)/(Lt) * [A(1/(Lt)) sin(Lt) + B(1/(Lt)) cos(Lt)].
```

Os radicais `sqrt(2n+1)` são envolvidos por raízes quadradas inteiras a 320 bits. Os extremos de `pi` usam o somatório alternado exato de `pi=16 atan(1/5)-4 atan(1/239)`. Para cada período `[850+2πk/L,850+2π(k+1)/L]`, `k=0,…,22`, o auditor demonstra com frações exatas que a variação conjunta das amplitudes é menor que `1/32` da amplitude de referência. A desigualdade triangular em `L²` preserva pelo menos `(15/16)^2` da energia trigonométrica; o extremo do último período é menor que 1000.

Somar os 23 limites racionais dá a massa de Fourier estritamente positiva acima. O certificado prova que a correção `(1/4)(I-G_850)` **elimina a antiga testemunha negativa da P258**.

### Limitação decisiva

**Esta é uma decisão de sinal para um único vetor, NÃO um certificado de positividade do bloco ímpar 377×377.** Um novo vetor negativo pode existir. A negatividade do auxiliar de piso constante da P258 e seu diagnóstico espectral não são apagados retroativamente. Não se reivindica nova extensão da forma completa de Weil nem prova da Hipótese de Riemann.

## Reprodução

O auditor P259 usa só **Python 3.11+ padrão**. Baixe o artefato `numeria-p258-certificates` da [execução P258](https://github.com/jbarumby/scicalc/actions/runs/44858303395/artifacts/14688468811), extraia seu arquivo `P258_GL60.json` para `replay/`, e execute na raiz do repo:

```bash
mkdir -p replay
gh run download 44858303395 --repo jbarumby/scicalc -n numeria-p258-certificates -D replay
sha256sum replay/P258_GL60.json
python -m NUMERIA_ES8P259.p259_witness --p258 replay/P258_GL60.json --out replay/P259_RESULT.json
python -m unittest NUMERIA_ES8P259.test_p259 -v
python -O -m unittest NUMERIA_ES8P259.test_p259 -v
```

O programa **rejeita** fonte alterada se o SHA-256 diferir de `1aa765357cfe1fa46f864d5b57fc407157375e75b3fa8b23b77f0db5c2e4985c`. Um workflow dedicado executa esses comandos e publica o JSON como artefato.

Leia [a demonstração completa](PROVA_P259.md) e o [resumo certificado](P259_RESULTADO.json). Para execução imediata, há também um pacote local ZIP com a entrada P258 completa e 12 testes que passaram tanto normalmente quanto com `-O`.

**Próxima prioridade — Gate B:** certificar `G_850` e a fatoração intervalar dos **377×377** modos ímpares. Se o bloco ainda for indefinido, obter uma **nova** direção negativa em aritmética dirigida, não reutilizar a testemunha agora reparada. Não mesclar os PRs anteriores antes dessa auditoria.
