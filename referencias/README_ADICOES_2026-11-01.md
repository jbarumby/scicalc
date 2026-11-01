# Novos materiais de referência — 2026-11-01

Os três arquivos originais enviados foram conferidos com os documentos existentes no GitHub. Estes registros são **documentos de referência**, não resultados certificados automaticamente pela suíte NUMERIA.

| Anexo original | SHA-256 dos bytes enviados | Registro nesta branch |
|---|---|---|
| `TAIL_PROOF(2).md` | `f96300fc5dc0e5e88c15ae5928bdfec1eb34f0730768beedb6b3f7f9b2208cc6` | [texto integral enviado](TAIL_PROOF_ENVIADO_2026-11-01.md) |
| `TEOREMA_RESTO_DIGAMA(1).pdf` | `03c40c328cb036129dfc5ac73f288ec71b77814da8ba292f64f43e2adaaab267` | [nota científica](TEOREMA_RESTO_DIGAMA_NOTA.md) |
| `paper-v5 (12).pdf` | `19f827bee5834d61aa6dd756cdaea582492703ddbfd6bdc2058de10b93f7e814` | [referência bibliográfica e cotejo](ARTIGO_ZETA_DOIS_TERCOS.md) |

## Situação dos binários

**Os PDFs binários originais ainda não foram publicados nesta branch**: o conector de escrita disponível nesta sessão permite criar arquivos UTF-8, mas não transferir diretamente os bytes dos anexos PDF locais. O conteúdo matemático relevante e as referências estão documentados nos arquivos acima; os SHA-256 indicam exatamente quais versões originais foram recebidas.

Para anexar os PDFs originais à própria branch, abrir o editor web do GitHub e adicionar os arquivos com esses nomes em `referencias/originais/`, ou enviar os arquivos incluídos no ZIP entregue na conversa. **Não confundir as notas Markdown com cópias byte a byte dos PDFs.**

## Alertas matemáticos

1. O `TAIL_PROOF` enviado exige `R >= max(1, N*pi/L)`. Com `R=160`, `L=9/20`, os casos `N=6,8,10,20` satisfazem essa pré-condição; **`N=30` não**. Isso não invalida o lema; impede usá-lo sem aumentar `R`.
2. O documento digamma anuncia `sign(R_K)=(-1)^(K+1)`, mas sua tabela marca como **positivo** o resto em `K=8`, `x=65` e `100`: pela própria fórmula, deveria ser **negativo**. Manter a correção como alerta, sem alterar silenciosamente o PDF.
3. O artigo `paper-v5` discute limite inferior para a proporção de zeros simples na linha crítica; **não** conclui a Hipótese de Riemann nem certifica as matrizes NUMERIA. Há uma versão pública posterior no arXiv (2608.13637v2).

## Teste de precondições sem dependências

```bash
python referencias/VERIFICAR_PRECONDICOES_NOTAS.py
```

Os hashes são integridade relativa, não autenticação criptográfica da autoria.
