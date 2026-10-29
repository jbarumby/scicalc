# NUMERIA E8-K V15 — Sinal e robustez do bloco arquimediano 2×2

**Escopo:** análise rigorosamente racional do **sinal**, sob as hipóteses analíticas documentadas em E7/V13 e a correção dos intervalos B5. Não é um certificado global da fórmula explícita de Weil nem da Hipótese de Riemann.

## 1. Novo resultado matemático (condicional)

A E7/V14 fornece intervalos fechados para as duas entradas arquimedianas:

\[
G_{11}\in[-1.4282431751596512,\ -1.4281737477533941],
\qquad
G_{22}\in[-0.1027780844359820,\ -0.1023800061201568].
\]

Estes números são aproximações *para leitura*; o verificador utiliza **frações exatas**, não os decimais.

A representação da transformada de Fourier documentada no `fontes/TAIL_PROOF.md` do E7 dá
\(\operatorname{Re}(F_i\overline{F_j})=\cos((i-j)\pi/2)A_iA_j\).
Portanto, \(G_{12}=0\) por paridade. Este passo depende da identidade analítica de Fourier, não foi formalizado pelo presente script.

Com \(a=-\sup G_{11}>0\), \(b=-\sup G_{22}>0\), qualquer matriz real simétrica
\(H=\begin{pmatrix}G_{11}&c\\c&G_{22}\end{pmatrix}\) satisfaz
\[
 v^THv\le -a v_1^2-b v_2^2+2|c||v_1v_2|.
\]
Se \(c=0\), é estritamente negativa definida, com \(v^THv\le-\min(a,b)\Vert v\Vert^2\).
Mais forte: se \(|c|\le1/10\), então \(\det H\ge ab-1/100>0\), e \(H_{11}<0\), logo \(H\) continua estritamente negativa definida (critério de Sylvester).

A desigualdade \(ab>1/100\) é verificada com `Fraction`, e o relatório publica o valor exato de \(ab-1/100\). Margens aproximadas: \(a\approx1,42817\); \(b\approx0,1023800\); \(ab-0,01\approx0,1362164\).

## 2. Verificações realizadas

O script `VERIFICAR_E8K_V15_BLOCO.py` confere: o SHA-256 do ZIP E7 e da fonte B5; as 27 entradas do `SHA256SUMS.txt` V14; os hashes dos dois JSONL; o schema e as contagens; a correspondência exata das somas das 11 bandas; a composição exata da soma B5 com quadratura B3 e cauda Prop. 6; a situação publicada dos 10 gates V13; o sinal dos extremos; a robustez \(|c|\le1/10\); e **sete testes negativos independentes**. As condições avaliadas são *racionais e verificáveis por programa*, não uma prova de correção da fonte Python.

**Resultados locais:** 11/11 gates `PASS` e 7/7 controles negativos detectados. O JSON é determinístico e reexecutado em separado. Os relatórios originais V13/V14 são preservados com hashes para rastreabilidade.

## 3. Limites inalterados e próximo problema científico

* As cotas contínuas de segunda derivada e cauda infinita **ainda dependem das identidades analíticas explicitadas no E7/V13**; o programa apenas verifica suas consequências racionais e confronta dados.
* O sinal negativo é **somente da contribuição arquimediana**. Não substitui a matriz completa da fórmula de Weil: faltam a dedução rigorosa da normalização global, os termos de primos, polos e as restrições de projeção.
* O Lean 4 existente **não foi compilado neste ambiente**; qualquer prova formal exige `lake build` com Mathlib e ausência efetiva de `sorryAx` nas dependências.
* As verificações de dimensão finita, mesmo que se fechem, **não provam** a positividade em todas as dimensões ou a Hipótese de Riemann.

**Próxima obrigação prioritária V16:** especificar e auditar algebricamente a fórmula explícita de Weil e seus termos (primos, polos e sinal), comparando a construção da matriz total com a derivação matemática, incluindo os vínculos da projeção. Não é correto extrapolar o sinal arquimediano para a forma total.

**Status:** `CONDITIONAL_2X2_ARCHIMEDEAN_NEGATIVE_DEFINITE`.
