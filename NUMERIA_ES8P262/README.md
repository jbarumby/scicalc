# ES8-P262 — KPS: momentos theta e auditoria de consistência

Base: PR #68 (P261). Não modifica a linha matricial de Weil.

## Fonte original

Konstantopoulos, Patie, Sarkar, *Annales de l'Institut Fourier* 74 (2024), pp. 377–421, artigo https://www.numdam.org/item/10.5802/aif.3600.pdf. **Página 397, Eqs. (4.15)–(4.18)** e **página 395, Eq. (4.8)**. Foram confrontadas imagens das páginas, não apenas o resumo.

F(k) = integral de 1 a infinito de (log x)^k x^-3/4 theta_0(x) dx, com theta_0(x)=soma_{m>=1} exp(-pi*m²*x). A substituição u=log x permite calcular F(k) com a integral de u^k exp(u/4) theta_0(exp u) de zero a infinito.

## Gate A matemático

- Integrar com python-flint/Arb os termos m=1..3 sobre u de 0 a 5, em cinco segmentos. O integrando é inteiro na variável u; a quadratura é acompanhada por controle de erro automático da biblioteca.
- Limitamos o resto m>=4 e u>=5 com duas fórmulas **racionais exatas**, demonstradas em PROVA_P262.md. Os extremos dos momentos são exportados como pares diádicos exatos (mantissa, expoente).
- Calcular G(2n)=64 n(2n-1)F(2n-2)/F(2n)-1 e a fórmula de phi(n+1) como literalmente impressa na Eq. (4.18).
- **Gate de lógica:** Se G(2)>0 e G(4)>0, a expressão impressa para phi(2) é negativa. Uma função de Bernstein é não negativa. Isso determina uma obrigação de conferir a normalização, índices e sinal no artigo — **não** prova que os autores erraram, não exclui a representação corretamente formulada, e não altera a RH.

## Reprodução

Na raiz do repositório, com as dependências já pinadas:

    python -m pip install -r NUMERIA_ES8P249/requirements.txt
    python -m unittest NUMERIA_ES8P262.test_theta_moments -v
    python -O -m unittest NUMERIA_ES8P262.test_theta_moments -v
    mkdir -p replay
    python -m NUMERIA_ES8P262.theta_moments --max-n 1 --bits 192 --out replay/P262_KPS_MOMENTS.json
    sha256sum replay/P262_KPS_MOMENTS.json

O workflow da branch executa os testes e guarda os resultados. Um resultado inconclusivo **não** deve ser transformado em uma afirmação de sinal.

Próximo Gate B: confrontar a definição multiplicativa W_phi da Eq. (4.8), os coeficientes da Eq. (4.17) e a fórmula impressa (4.18) em álgebra simbólica; verificar por uma segunda fonte a normalização de Xi(t)/Xi(0), e só então testar diferenças de Bernstein para sequências maiores.
