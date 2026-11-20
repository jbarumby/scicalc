# NUMERIA–RIEMANN · ES8-P253 · Limiar do primo 7

A coleção *Spectral Theory and Operators in Mathematical Physics* motivou duas escolhas. A localização por suporte de Gabel et al. sugere isolar as faixas em que o novo deslocamento pode atuar. A decomposição fora da diagonal de Seelmann sugere tratar a ligação entre essas faixas separadamente. O limite abaixo é demonstrado diretamente pela fórmula explícita, Plancherel e uma conta de forma quadrática. Os resultados dos artigos tratam dos seus próprios operadores. [Coleção](https://link.springer.com/collections/hijacgibha) · [Gabel et al., Lema 3.1](https://doi.org/10.1007/s11785-023-01448-3) · [Seelmann, Teorema 1.4](https://doi.org/10.1007/s11785-022-01209-8)

## Testes do modelo auxiliar

O modelo auxiliar usa quatro potências primas, β da P250 e os polos originais. Foram montados os dois setores GL60 e verificadas direções negativas em Arb com vetores diádicos exatos. As cotas superiores para seus quocientes de Rayleigh já incluem o erro analítico GL60.

| Corte | Meia-largura auxiliar | Setor par | Setor ímpar |
| ---: | --- | ---: | ---: |
| 309 | `0,973` | `< −1,85×10^-26` | `< −1,79×10^-17` |
| 350 | `0,972955074527656652552676371721590` | `< −4,70×10^-29` | `< −8,68×10^-20` |
| 504 | `0,972955074527656652552676371721590` | `< −1,12×10^-29` | `< −5,30×10^-23` |

A segunda largura supera `log(7)/2` em aproximadamente `1,35×10^-34`. Na primeira tentativa de ampliação do corte, sementes de Bessel a 640 bits geraram caixas excessivamente largas nos painéis acrescentados. Recalcular com 960 bits estreitou as caixas. Para o corte 504, usamos 1.008 painéis, 377 modos por setor e sementes a 1.536 bits. O replay MPFR do centro par interrompeu no pivô de índice 42, em acordo com Arb. As direções negativas da tabela se referem ao modelo auxiliar truncado. A forma de Weil contém ainda a cauda e o termo de 7. [Vetor 309](data/T309_NEGATIVE_WITNESS.json) · [Vetor 350](data/T350_NEGATIVE_WITNESS.json) · [Vetor 504](data/T504_NEGATIVE_WITNESS.json) · [Auditoria do corte 504](data/P253_THRESHOLD_504.json)

## Critério para atravessar o limiar

Suponha que o modelo auxiliar de quatro termos tenha margem certificada `m>0` numa largura `L_aux>log(7)/2`, com corte entre 309 e 504. A extensão por zero transfere essa margem a qualquer janela menor para o mesmo modelo auxiliar. Escolha `0<δ≤m/5040` com `log(7)+δ<2L_aux<log(8)`, e ponha `L_*=(log(7)+δ)/2`. Para toda `f ∈ C_c^∞(ℝ;ℂ)` com suporte em `[-L_*,L_*]`, a forma de Weil satisfaz

```text
Q_W(f) ≥ (m−2520δ) ‖f‖₂² + (1/12) ‖P_{|t|>2016} f‖₂².
```

Em particular, a margem garantida seria de pelo menos `m/2`, numa janela local que inclui o novo deslocamento.

## Deslocamento nas bordas

Ponha `d=log(7)` e `κ=d/√7`. Para `d<2L<2d`, a translação comprimida `J f(x)=f(x−d)`, com extensão por zero, leva isometricamente a faixa `E_−=[−L,−L+δ_L]` à faixa `E_+=[L−δ_L,L]`, onde `δ_L=2L−d`. A translação simétrica `J+J*` troca as faixas e zera no meio. Sua norma é exatamente 1 para toda `δ_L>0`. Em `2L=d` as faixas têm medida zero. A fórmula explícita acrescenta exatamente `−κ⟨f,(J+J*)f⟩`, pois a amplitude `2κ` em frequência multiplica `cos(td)`. [Normalização de Guinand–Weil](../NUMERIA_ES8P251/gate_c/README.md)

O bloco integral de banda finita é compacto no intervalo e os polos têm posto finito. Logo, aplicar a correção de 7 diretamente ao antigo modelo limitado `βI+K` coloca `β−κ<0` no espectro essencial. De fato, uma sequência ortonormal de vetores da forma `(u_n+Ju_n)/√2`, com `u_n` na faixa esquerda, tem `(J+J*)f_n=f_n` e `Kf_n→0`. Se uma malha uniforme de pelo menos duas células tem largura de célula `h>δ_L`, a compressão do deslocamento às constantes por célula tem apenas o par de entradas opostas nos cantos, de módulo `δ_L/h`. A malha pode esconder a obstrução. Na forma completa, a cauda de alta frequência oferece a reserva usada a seguir.

## Absorção pela cauda de alta frequência

Escreva `E=E_−∪E_+` e `Ω=2016`. A desigualdade `2ab≤a²+b²` dá `|⟨f,(J+J*)f⟩|≤‖1_E f‖₂²`. Na convenção angular de Fourier da P251, a projeção `P_Ω` sobre `|t|≤Ω` tem `‖1_E P_Ω‖²≤|E|Ω/π`, pela norma de Hilbert–Schmidt de seu núcleo no par posição-frequência. Separe `f=P_Ω f+P_{>Ω}f` e aplique Young com pesos `5/2` e `5/3`.

```text
κ ‖1_E f‖₂² ≤ (5/2)κ(2δ_L Ω/π) ‖f‖₂² + (5/3)κ ‖P_{>Ω}f‖₂²
             ≤ 2520δ_L ‖f‖₂² + (5/4) ‖P_{>Ω}f‖₂².
```

O resíduo `w₄` de quatro potências primas é positivo em `|t|≥309` pela P250. Para `|t|≥2016`, Binet e a soma máxima das quatro amplitudes dão `w₄(t)≥B(2016)−B(504)>4/3`, onde `B(t)=log(t/(2π))−1/t`. Essa reserva absorve o termo `5/4` e deixa `1/12` no resíduo. O auditor verifica `κ<3/4`, `π>3` e `2L_*<2L_aux<log(8)`, de forma que 7 é a única nova potência prima no extremo. [Sinal P250](../NUMERIA_ES8P250/README.md)

## Sinal do símbolo com cinco termos

Em uma auditoria independente, o resíduo completo com cinco potências primas e o mesmo β foi certificado por Arb a 192 bits para `|t|≥1551`. Os 704 intervalos diádicos cobrem `[1551,2200]`, e Binet com a soma das amplitudes cobre `t≥2200`. Em `t=15501/10`, a avaliação direta intervalar da digama dá um resíduo em `(−0,088524106,−0,088524104)`. Assim, **1551 é o menor corte inteiro para o sinal uniforme do novo símbolo com esse β**. A estimativa de borda trata 7 em posição e trabalha com o modelo auxiliar. A auditoria de sinal abre outra rota para janelas maiores na faixa `7<exp(2L)<8`, que exige certificar o bloco retido com cinco termos até 1551. [Certificado de sinal](data/PRIME7_TAIL.json)

Para o modelo auxiliar de quatro termos, um replay adicional certificou o piso espectral `1/4` a partir de 504, usando 162 intervalos até 650 e a cota de massas depois de 650. Substituir o β antigo por esse piso adicionaria ao modelo `(1/4−β)P_{|t|>504}`, comprimido à janela. O piso da cauda já está verificado. A próxima etapa é certificar a matriz auxiliar com essa atualização e recalcular a reserva usada na absorção de 7. [Certificado do piso maior](data/FOUR_PRIME_RAISED_FLOOR.json)

## Reprodução

Na raiz do repositório, instale as dependências da P249. O replay MPFR requer compilador C, GMP e MPFR no Linux.

```bash
python -m pip install -r NUMERIA_ES8P249/requirements.txt
sha256sum -c NUMERIA_ES8P253/SHA256SUMS
work=$(mktemp -d)
matrix=NUMERIA_ES8P253/data/P253_CLOSE504.json.gz
hash=$(sha256sum "$matrix" | cut -d' ' -f1)
width=0.972955074527656652552676371721590
python -m NUMERIA_ES8P253.threshold_audit \
  --matrix "$matrix" --sha256 "$hash" --cutoff 504 --modes 377 \
  --half-width "$width" --out "$work/P253_THRESHOLD_504.json"
python -m NUMERIA_ES8P253.low_cutoff_witness \
  --matrix "$matrix" --sha256 "$hash" --cutoff 504 --modes 377 \
  --half-width "$width" --out "$work/T504_NEGATIVE_WITNESS.json"
python -m NUMERIA_ES8P253.prime7_tail --out "$work/PRIME7_TAIL.json"
python -m NUMERIA_ES8P253.prime7_tail --four-prime-floor --out "$work/FOUR_PRIME_RAISED_FLOOR.json"
python -m unittest NUMERIA_ES8P253.test_threshold -v
python -O -m unittest NUMERIA_ES8P253.test_threshold -v
for result in P253_THRESHOLD_504 T504_NEGATIVE_WITNESS PRIME7_TAIL FOUR_PRIME_RAISED_FLOOR; do
  cmp "NUMERIA_ES8P253/data/$result.json" "$work/$result.json"
done
```

Para regenerar os lotes com `finite_head_quadrature`, use as larguras e cortes da tabela de resultados e os seguintes parâmetros.

| Corte | Modos por setor | Painéis | Painéis por lote, salvo último | Bits GL60 | Bits Bessel |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 309 | 259 | 618 | 80 | 256 | 640 |
| 350 | 259 | 700 | 80 | 256 | 960 |
| 504 | 377 | 1.008 | 84 | 256 | 1.536 |

Use `--threads 1` em cada lote, e passe os lotes a `threshold_audit` com `--chunks`, `--matrix`, o mesmo corte, a mesma largura e os mesmos modos. O SHA-256 esperado está no certificado correspondente. Use então `low_cutoff_witness` com os mesmos parâmetros. A tentativa MPFR no corte 504 pode ser reproduzida com `python -m NUMERIA_ES8P251.gate_b.mpfr_audit --cutoff 504 --modes 377 --half-width "$width" --matrix "$matrix" --out "$work/MPFR.json"`. A fatoração encerra com código 3 no executável C ao encontrar o pivô do setor par. [Comandos da quadratura](../NUMERIA_ES8P252/README.md#reprodução)

A janela certificada pelo NUMERIA continua sendo a da P252, com meia-largura `0,95`. O novo argumento fornece uma condição suficiente para uma passagem local.

Na validação desta etapa, passaram 67 testes da P249 à P253 tanto em Python normal quanto com `-O`. O relatório do corte 504, seus vetores diádicos e os dois certificados de cauda foram reproduzidos byte a byte. Os sete hashes SHA-256 do manifesto conferiram.
