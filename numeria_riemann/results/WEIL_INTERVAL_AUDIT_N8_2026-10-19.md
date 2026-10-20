# NUMERIA–RIEMANN — auditoria intervalar ampliada N=8

**Data:** 2026-10-19  
**Branch:** `numeria-riemann/explicit-formula-prototype`  
**Execução CI:** https://github.com/jbarumby/scicalc/actions/runs/39180626906  
**Artefato N=8:** https://github.com/jbarumby/scicalc/actions/runs/39180626906/artifacts/12156246209

## Configuração executada

- (L=0{,}45), (N=8); dimensão projetada (N-2=6).
- (R=160), passo (0{,}02), 8.000 células de quadratura.
- 128 termos da série de digamma; precisão Arb de 192 bits.
- Aritmética intervalar Arb nas entradas, restrições, sobreposições e cauda analítica.

## Resultado observado

O job e os testes auxiliares passaram tecnicamente, mas o estado matemático reportado foi **`inconclusive_interval_ldl`**. Os limites inferiores dos pivôs LDL divulgados foram:

| Pivô | Limite inferior |
|---:|---:|
| 1 | 0,2149252225123455007633539770 |
| 2 | 0,6591348761638142069767114403 |
| 3 | 0,3643216324439690834660359699 |
| 4 | 0,9835229029954672328707609370 |
| 5 | -0,2238731886330354285625576508 |

Os quatro primeiros limites reproduzem os valores positivos da auditoria (N=6), enquanto o quinto deixa de fornecer a condição necessária para certificar positividade. Isso **não demonstra**, por si só, que a matriz exata seja indefinida: o teste intervalar só certifica positividade quando cada pivô tem limite inferior estritamente positivo. Também não é correto descartar o sinal negativo como mero erro numérico sem análise adicional.

Somente a potência prima (2) ficou ativa para (L=0{,}45) na rotina de sobreposição de primos.

## Próximo teste já adicionado ao CI

Foi adicionada uma segunda execução para (N=8), com passo (0{,}01), 256 termos e precisão de 256 bits. A comparação entre os dois relatórios ajudará a distinguir um intervalo largo por sobre-enclausuramento de um pivô que permanece incompatível com uma certificação positiva. Mesmo que o segundo cálculo permaneça inconclusivo, o resultado continua restrito a essa matriz finita e não prova nem refuta a Hipótese de Riemann.

## Escopo e cautelas

- O resultado é uma auditoria finita, intervalar e exploratória da matriz projetada.
- `inconclusive_interval_ldl` não é um certificado de positividade nem uma prova de indefinição.
- Nenhuma conclusão global sobre positividade de Weil ou a Hipótese de Riemann decorre deste cálculo.


Atualizacao: o workflow refinado 39180676723 terminou com sucesso; consulte os resultados detalhados no registro desta auditoria.
