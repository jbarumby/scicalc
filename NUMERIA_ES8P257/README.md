# NUMERIA–RIEMANN — ES8-P257 · Taper, norma e bordas de Fourier

O passo seguinte à P256 desenvolveu dois resultados matemáticos independentes da RH: (i) uma regularização exata do símbolo oscilatório de primos e uma **obstrução à aproximação em norma** do novo deslocamento de 7; (ii) uma melhoria **condicional** do orçamento de borda da P253, de `2520δ` para `2225δ`, sob as mesmas hipóteses de Guinand–Weil, Binet e cauda positiva a partir de 504. Nenhuma nova janela acima de `log7/2` está certificada.

## Como reproduzir

Python padrão, nenhuma dependência externa:

```bash
python -m unittest NUMERIA_ES8P257.test_p257 -v
python -O -m unittest NUMERIA_ES8P257.test_p257 -v
python -m NUMERIA_ES8P257.p257_strip_taper --out /tmp/P257_CERT.json
```

O auditor executa verificações **exatas em `Fraction`**, inclusive constantes para `log7`, `sqrt7`, `π`, limite de cauda, otimização de Young e lacuna de norma do primo 7. A prova das identidades funcionais e o escopo das hipóteses estão em `PROVA_ES8P257.md`. Os resultados não compilam Lean nem substituem os certificados intervalares P250–P253.

A comparação entre as cotas `2520δ` e `2225δ` só autoriza uma extensão condicional, se houver **uma margem comprovada para o modelo auxiliar acima do limiar**. A P253 relatou direções negativas em modelos auxiliares; essa dificuldade permanece.

Documentação de base: [P256](https://github.com/jbarumby/scicalc/pull/62), [P253](https://github.com/jbarumby/scicalc/pull/59), [P251 Gate C](https://github.com/jbarumby/scicalc/pull/55).
