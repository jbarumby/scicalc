# NUMERIA–ES6 — fórmula explícita gaussiana e diferenciação (2026-11-03)

**Status:** reprodução numérica da fórmula **clássica, incondicional** de Guinand–Weil e de sua derivada em relação ao parâmetro gaussiano; **certificação aritmética apenas do termo de potências de primos** para `a=1/50`; **sem prova de RH**, sem certificação de quadratura ou de completude dos zeros, sem conexão provada ao operador ES4.

## Conteúdo

- `es6_weil_gaussian.py`: convenções explícitas, limite racional de primos, experimento gaussiano e derivada em `a`, controle adversarial ES4;
- `test_es6_weil_gaussian.py`: 12 testes, incluindo testes das séries racionais e reprodução numérica;
- `RESULTADOS_ES6.json`: duas amostras de traço, dados da derivada, certificado do termo primo, diagnósticos ES4;
- `RELATORIO_ES6.md`: fundamentação matemática, avaliação e limites;
- `MANIFESTO_SHA256.json`: hashes dos arquivos publicados;
- `TESTES_ES6.txt`: registro de execução local.

## Reprodução

```bash
python -m pip install mpmath==1.3.0 pytest==8.4.2
python -m pytest -q test_es6_weil_gaussian.py
python es6_weil_gaussian.py --dps 55 --zeros 24 --N 128 --T 128 --out RESULTADOS_ES6.json
```

Os testes são autorreferenciais apenas para propriedades elementares de aritmética racional; a verificação da fórmula global é **exploratória** porque usa a lista finita `zetazero(j)` e `mpmath.quad` sem rigor intervalar. Um resíduo pequeno **não** demonstra ausência de zeros fora da reta crítica.

## Fórmula principal

Com `h_a(t)=exp(-a t²)`, `a>0` e `z_rho=(rho-1/2)/i` possivelmente complexo, tem-se a identidade conhecida

\[
\sum_{\rho}h_a(z_\rho)=2 e^{a/4}+\frac1{2\pi}\int_{\mathbb R}h_a(t)\,[\operatorname{Re}\psi(\tfrac14+\tfrac{it}2)-\log\pi]dt-\frac1{\sqrt{\pi a}}\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}e^{-(\log n)^2/(4a)}.
\]

**A soma à esquerda inclui TODOS os zeros não triviais, com multiplicidade**; a conta exploratória inclui somente os 24 primeiros zeros conhecidos da reta crítica e seus conjugados. A forma gaussiana é admissível para o teorema clássico, mas **não** é compactamente suportada em coordenadas logarítmicas; transportar o resultado às matrizes E7/V24 exige nova análise de domínios, bases, projeções e majorantes.

## Fontes principais

- E. Carneiro, A. Chirre e M. B. Milinovich, *Bounding the log-derivative of the zeta-function*, Proposition 5: https://doi.org/10.1007/s00209-021-02820-9 (convenção de Fourier declarada no artigo).
- Paul Garrett, *Guinand's explicit formula* (notas): https://www-users.cse.umn.edu/~garrett/m/complex/notes_2020-21/12g_guinand_explicit_fml.pdf
- Projeto NUMERIA ES5 (commit `2484a37170a077089e8c7626a11edb7d92c5ee23`), V24 (`v24/README_V24.md`).
