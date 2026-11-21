# NUMERIA P76 — resultado do CI em FLINT/Arb

Data: 2026-11-21. Branch: `numeria-p76-cs-v2-coefficients`.

**[Workflow GitHub Actions 44949581862](https://github.com/jbarumby/scicalc/actions/runs/44949581862): concluído com sucesso.**

O trabalho `audit-b23` aprovou as etapas de instalação de `python-flint==0.9.0`, validação dos hashes das fontes, 7 testes de `unittest` em Python normal e com `-O`, comparação das 16 linhas com os decimais exatos em 192 e 256 bits, e criação do artefato `numeria-p76-b23-arb-enclosures`.

| Prova efetuada | Resultado |
| --- | --- |
| `B2_publicado > B2_recalculado` nas 16 linhas | PASSOU (Arb 192 e 256 bits) |
| `B3_publicado > B3_recalculado` nas 16 linhas | PASSOU (Arb 192 e 256 bits) |
| Menor margem rigorosa entre todas as diferenças | > 0,0000841457933772 (arredondado para baixo) |
| A3, A4, A5, C5, C6 | Derivados das fórmulas originais |
| Fonte | Chourasiya–Simonič, arXiv:2507.15184v2, Teorema 1, equações (10), (33), (34) |
| B1 | **NÃO RECONSTITUÍDO** |
| Conclusão global do Corolário 1 | **NÃO CERTIFICADA INDEPENDENTEMENTE** |

Os JSONs completos foram salvos no artefato do workflow, com SHA-256 dos bytes exportados:

- `p76_b23_192.json`: `8c013dff67f8970baf8abc8c6fd714ece687d3368b00854e54bdaf96f7fadd0a` (16.983 bytes).
- `p76_b23_256.json`: `5d798a1f0c10c20502cb478cca3811b293d06701a93b3016a020b7fabef1949b` (18.933 bytes).

A menor margem corresponde à linha `0.75000 ≤ σ ≤ 0.78125`, coeficiente B2. O audit calculou o limite inferior dyádico exato em 256 bits `2494314807213989001534656931548156787783495242521846276618043388784143200369 · 2^(-264)`.

**Cuidado metodológico:** esta auditoria verifica exclusivamente as duas colunas B2 e B3 derivadas do Teorema 1. A coluna B1 depende da quadratura rigorosa de K1 e dos demais termos K/L/A; por isso nenhuma afirmação sobre RH ou sobre nova certificação completa da P75 deve ser deduzida daqui. O PR continua em rascunho.

Próxima fase: reproduzir as 16 desigualdades de B1 com limites Arb para a integral imprópria K1, conferindo os parâmetros A0 e o arredondamento de todos os coeficientes.
