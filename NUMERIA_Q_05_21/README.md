# NUMERIA-Q 05.21 — descoberta de circuitos reversíveis

Protótipo de síntese exata em instâncias pequenas para comparar famílias de circuitos reversíveis que calculam inverso modular em campos primos.

## Estado científico

- Benchmark: f(x) = x⁻¹ mod p para 1 ≤ x < p; f(x)=0 para x=0 e códigos inválidos.
- Transformação reversível: |x>|y>|0…0> → |x>|y XOR f(x)>|0…0>.
- Famílias: tabela-verdade por mintermos e forma normal algébrica (ANF).
- Referência clássica independente: algoritmo euclidiano estendido.
- Portas abstratas: X, CX e CCX/Toffoli.
- Avaliação: equivalência exaustiva para todos os estados da base computacional, preservação de x, saída y correta e ancillas limpas.
- Métricas: qubits, portas, profundidade lógica, Toffoli × qubits e portas × qubits.

## Execução

Requisitos: Python 3.10+, biblioteca padrão.

```bash
python run_experiment.py --out resultados_05_21
python -m unittest discover -s tests -v
```

## Limitações e cautela

Este protótipo não implementa o circuito divsteps, curvas elípticas, Berlekamp–Massey ou decodificação Reed–Solomon do trabalho da Anthropic. Não reproduz os resultados publicados, não demonstra vantagem quântica e não é evidência de descoberta inédita. As métricas são proxies lógicas e não incluem hardware, ruído, conectividade nem correção de erros.

As versões congeladas NUMERIA-Q 05.19.1 e 05.20 devem permanecer intactas. Todo o material desta etapa está isolado nesta pasta.
