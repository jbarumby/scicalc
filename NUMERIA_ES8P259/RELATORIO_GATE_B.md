# ES8-P259 Gate B — direção ímpar negativa com piso em duas faixas

**Execução auditada:** [GitHub Actions #44863150729](https://github.com/jbarumby/scicalc/actions/runs/44863150729). **Modelo:** auxiliar de quatro potências primas, corte 504, peso constante 1/4 após 504 e 1/2 após 850. **Não é** a forma completa de Weil nem uma alegação de prova da RH.

## Matriz e fonte

```text
A259(L) = A_beta(L) + (1/4-beta)(I-G_504(L)) + 1/4(I-G_850(L))
L = 0.972955074527656652552676371721590
```

A quadratura ponderada original da P253 é autenticada pelo mesmo SHA fixado na P258; G504/G850 são avaliadas em Arb por identidades de extremos das Bessel esféricas. A cota de erro analítico GL60 é aplicada uma vez. O replay racional das caixas da matriz não trata pivôs inconclusivos como sinais negativos.

## Resultados reproduzidos

| Prefixo ímpar | Resultado |
| --- | --- |
| 29x29 | certificado positivo; margem > 1.0045e-23 |
| 36x36 | certificado positivo; margem > 5.3957e-25 |
| 48x48 | certificado positivo; margem > 2.2407e-25 |
| 64x64 | **direção negativa**: Rayleigh < -4.77017716e-27 |

O pivô 56 (índice zero-based) foi inconclusivo. A decisão negativa vem de uma testemunha **de 57 coordenadas diádicas**, com replay de Rayleigh por `Fraction` e todos os erros GL60. Por ser um subespaço principal do setor ímpar, essa testemunha já impede positividade dos **377 modos** na largura testada.

SHA-256 do certificado completo gerado pelo workflow: `09b5a27700fa90f9c268882dbff4459fef5e6d72bd0b528f8e6e9a52d628b69a`. Acesso ao arquivo integral no artefato da execução GitHub; o workflow confere esse digest antes da etapa de transporte.

## Transporte ao lado inferior do limiar de 7

O funcional auxiliar P259 é o da P258 acrescido de **1/4 da energia Fourier em |t|>850**. Para funções dilatadas unitariamente, Plancherel dá:

```text
||dF_L/dL||_{|t|<=T} <= 1/(2L)+T <= T+1   (0.95<=L<=1)
|d/dL [ (1/4) ||P_{|t|>850} F_L||² ]| <= (850+1)/2 < 426
```

Usando **explicitamente como premissas** as cotas analíticas de P239/P258 para o auxiliar original (resíduo <=21 na banda <=504 e derivada dos polos <=36), o orçamento P258 é 21250. Logo o orçamento da P259 é 21676. As frações racionais da série de `log(7)=2 atanh(3/4)` provam

```text
L_- = 0.9729550745276566525526763717215898648185 < log(7)/2 < L_+.
L_+-L_- = 1.351815e-34
21676*(L_+-L_-) = 2.930194194e-30
Rayleigh(L_-) < -4.76724696e-27 < 0
```

A extensão por zero preserva a mesma forma quadrática com multiplicador espectral **fixo** e termos de polos inalterados. **Condicionalmente às cotas uniformes P239/P258 citadas**, este auxiliar específico de dois pisos é indefinido em toda janela de meia-largura >=log(7)/2. A P259 Gate A apenas corrigiu a **antiga** testemunha; a nova obstrução aparece no prefixo de 64 modos.

## Reprodução e limitações

`python -m NUMERIA_ES8P259.gate_b` executa Arb; `python -m NUMERIA_ES8P259.transport` audita o sinal e as desigualdades em **frações exatas**, com testes normais e em `-O`. O workflow autentica os dados P253/P258 e publica ambos JSON.

**Não** se conclui que a forma completa de Weil é indefinida. Não se conclui RH. Para a próxima prioridade, investigar minorantes **mais finos que dois pisos** ou reter parte do símbolo oscilatório real sem substituí-lo apenas por constantes. O PR deve permanecer em rascunho até revisão independente do orçamento de dilatação e dos certificados originais.
