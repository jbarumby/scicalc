# NUMERIA–RIEMANN · ES8-P261 — Lee–Yang, Fisher, Laguerre–Pólya e Bost–Connes

**Linha independente e paralela** aos PRs P259–P260 de positividade de Weil. Base: PR #67, sem modificações da P260.

## Resultado científico confirmado

O artigo original **Konstantopoulos–Patie–Sarkar (2024)**, *Annales de l'Institut Fourier* 74(1), 377–421, [DOI 10.5802/aif.3600](https://www.numdam.org/articles/10.5802/aif.3600/), Teorema 2.3, mostra que a função de Landau normalizada
`Xi(t)/Xi(0) = xi(1/2+it)/xi(1/2)`
pertence à classe de Lukacs `D_L = LP_even ∩ P_+` **se e somente se RH**. A propriedade de função característica `P_+` já vale incondicionalmente, pela densidade de Fourier positiva da representação theta (eq. 2.7). O Corolário 2.4 reformula a RH como a condição de probabilidade para o recíproco na linha imaginária.

**Nenhuma dessas reformulações é uma prova independente de RH.**

## Gate A: teste matemático efetivamente reproduzível

- Modelo finito de Ising numa cadeia aberta com peso racional `a=e^{-2K}` e fugacidade `z`: polinômio `P_N(z,a)=∑_{σ}a^{discordâncias(σ)} z^{#spins negativos}`.
- **Lee–Yang:** quatro spins ferromagnéticos `a=1/2`: prova algébrica exata que os quatro zeros de `P_4(z,a)` estão em `|z|=1` (redução palindrômica `y=z+z^{-1}`). Em `a=2` (antiferromagnético), os coeficientes permanecem positivos e palindrômicos mas a propriedade de Lee–Yang **falha**.
- **Fisher:** na fugacidade fixa `z=1`, `P_N(1,a)=2(1+a)^{N-1}`; as raízes de temperatura estão em `a=-1` e não devem ser confundidas com o teorema Lee–Yang na variável do campo.
- **Barreira para RH:** `φ(t)=3/4+(1/4)cos(t)` é função característica de uma distribuição genuína e produz matriz Toeplitz 4×4 com **15 menores principais não negativos**, mas **não** é Laguerre–Pólya, pois `(φ′)^2−φφ″` em `t=π` é `−1/8`. Ter apenas positividade probabilística/Toeplitz é insuficiente.
- **Bost–Connes:** a função de partição de equilíbrio corresponde a `ζ(β)` no regime `β>1`; o truncamento `1+2^{−s}` já possui zeros em `Re(s)=0`, não na reta crítica. Não há extensão automática de uma lei do círculo ao prolongamento analítico de `ζ`.

Todos os resultados do Gate A são **álgebra racional exata**, não extrapolação experimental nem confirmação da RH.

## Execução

Na raiz do repositório, Python 3.11+ padrão, sem dependências externas:

```bash
python -m unittest NUMERIA_ES8P261.test_p261 -v
python -O -m unittest NUMERIA_ES8P261.test_p261 -v
mkdir -p replay
python -m NUMERIA_ES8P261.leyang_fisher --out replay/P261_GATE_A.json
sha256sum replay/P261_GATE_A.json
```

O workflow da branch repete esses comandos, faz assertivas de resultados exatos e arquiva o certificado JSON com SHA-256.

## Próximo Gate B: localizar uma hipótese realmente independente

Estudar a eq. (4.18) do artigo KPS, que traduz uma possível representação `Xi=J_Ψ` em restrições para uma função de Bernstein `φ` construída com momentos theta `F(2n)`. Uma prova efetiva exigiria limites de erro **uniformes** sobre esses momentos e o controle global de todas as condições da classe, não mera verificação de alguns coeficientes.

No ramo Weil, exigir uma **identificação operacional explícita** entre um núcleo positivo de funções características e a forma completa de Weil (inclusive polos/cauda e domínios). A positividade de Toeplitz, já conhecida para Xi, **não** transfere automaticamente à forma de Weil ou aos minorantes P259.

Leia [relatório matemático](RELATORIO_P261.md) e [fontes](SOURCES.md). A P261 **não prova RH**, não altera a certificação local da P252 e não resolve a negatividade do auxiliar P259.
