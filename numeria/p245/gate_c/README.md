# NUMERIA ES8-P245 Gate C — auditável

**Resultado:** 650 entradas de quadratura em painéis 0–24, com 26 somas de 25 painéis incluídas nos blocos Arb P236; mais 156 avaliações de painéis estratificados, 52 comparadas a intervalos Arb P236; cota de erro GL60 por linha reavaliada em frações exatas sob hipóteses analíticas importadas. Não certifica a matriz inteira ou RH.

- `p245b_certified_gl60.py` e `P245B_ROOT_BRACKETS.json`: biblioteca intervalar e 60 raízes isoladas do Gate B, preservadas.
- `p245c_chunk25.py`: reconstrução e comparação de 25 painéis consecutivos, 13 pares, 2 setores.
- `p245c_stratified.py`: seis painéis internos e de borda, 13 pares, 2 setores.
- `p245c_gl60_error_budget.py`: orçamento racional independente para 1008 painéis, 754 graus.
- `P245C_*.json`: certificados dos resultados.
- `P236_CHUNK25_{even,odd}.json.gz`: entradas intervalares Arb do bloco original; hash verificado antes da comparação.
- `pilot_sources/`: quatro arquivos Arb de painéis individuais (503 e 1007).
- `test_p245c.py`: 14 testes, normal e otimizado.
- `RELATORIO_ES8P245_GATE_C.md`: metodologia, prova condicional do orçamento e limitações.
- `SHA256SUMS_P245C.txt`: manifesto de proveniência de todos os arquivos do pacote (exceto ele próprio).

**Para executar** use os comandos do relatório. O ZIP leve não inclui os arquivos de origem P236 (é necessário colocá-los nos caminhos acima); o ZIP completo os incorpora.
