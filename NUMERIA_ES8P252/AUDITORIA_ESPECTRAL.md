# P252 · Auditoria espectral da janela L = 19/20

Os dados são os 618 painéis GL60 até o corte 309, agregados nos blocos par e ímpar de 259 modos. Cada replay verifica o SHA-256 do agregado. As contas usam Arb a 1.024 bits e o orçamento analítico GL60 da P252. [Agregado](data/P252_GL60.json.gz) · [Construção e erro da quadratura](README.md) · [Hashes](SHA256SUMS)

## Menor autovalor dos blocos finitos

Seja `H` o bloco contínuo na base de Legendre e `H₀` o ponto médio da matriz intervalar montada com `β` e polos. A soma máxima dos raios por linha mais a cota analítica GL60 dá `δ ≥ ‖H−H₀‖₂`. Um vetor diádico exato `v`, escolhido por iteração inversa, fornece uma cota superior para `λ_min(H)` via quociente de Rayleigh intervalar acrescido do erro GL60. Para cada candidato `m` a limite inferior, verificam-se todos os pivôs de `H₀−(m+δ)I` por LDL intervalar. A positividade dos pivôs e a definição de `δ` implicam `λ_min(H)>m`. A iteração inversa serve para escolher o vetor. A prova usa apenas o vetor diádico obtido, o quociente intervalar e os pivôs.

| Bloco | Limite inferior estrito, arredondado para baixo | Limite superior, arredondado para cima | Pivôs positivos |
|---|---:|---:|---:|
| Par | `6,8203829349514113410285×10^-27` | `6,8203829349514113410842×10^-27` | 259/259 |
| Ímpar | `1,2379631976818498952630219×10^-23` | `1,2379631976818498952630276×10^-23` | 259/259 |

As larguras dos intervalos diádicos são inferiores a `5,562×10^-47` nos dois setores. A cota de erro GL60 domina a largura obtida. O JSON guarda os extremos diádicos e as coordenadas exatas de `v`. [Certificado dos intervalos](data/P252_SPECTRAL_BRACKET.json)

## Verificação escalar do polo ímpar

Escreva `H_ímpar = A−2ppᵀ`, com `A` o bloco contínuo sem polo. Para a margem racional exata `m=1,2379631976818498952630×10^-23`, constrói-se `M≼A−mI` usando o ponto médio do agregado ímpar, um limite inferior para `β` e a subtração dos orçamentos de erro da quadratura. Arb certifica `M≻0` e resolve `Mx=p`, com `p` avaliado por Bessel modificada. Pela congruência com `M^(-1/2)`, a condição `M−2ppᵀ≻0` equivale a `1−2pᵀM⁻¹p>0`. O limite inferior certificado para esse escalar supera `7,8412×10^-45`, o que confirma `H_ímpar≻mI`. A verificação usa os mesmos dados de quadratura e um teste diferente para o polo. [Certificado de Schur](data/P252_ODD_POLE_AUDIT.json)

A auditoria dá limites para os dois blocos finitos. A positividade de Weil sobre funções de teste na janela fixa segue pela ponte analítica e pelo complemento de Schur da P252. A extensão a uma família cofinal de janelas continua por demonstrar no NUMERIA. [Janela certificada](README.md)

## Relação com os artigos consultados

- O teorema 2.2 de Joaquim Duran compara limites espectrais de um operador de Dirac com o sinal de um problema ∂-Robin dependente de parâmetro. Isso motivou testar o sinal de uma matriz deslocada para delimitar os autovalores da P252. Aqui a justificativa é a congruência LDL de matrizes simétricas, a cota em norma e, no setor ímpar, a estrutura do polo. [Artigo de Duran](https://doi.org/10.1007/s11785-026-02074-5)
- O resumo de Sinha, Sen e Tripathy descreve o espectro de operadores triangulares inferiores de Cesàro `(p,q)` em `ℓ_s`. Aplicar resultados desse tipo aos blocos autoadjuntos de Weil exigiria demonstrar uma relação entre os operadores. A parte pública do artigo fornece a definição de Cesàro e o escopo dos resultados. [Resumo e definição](https://link.springer.com/article/10.1007/s11785-026-01960-2)
- Esfahanian e Arntzen discutem a lógica das proposições quânticas. No circuito P252, os projetores de seletor zero e de permanência nas células físicas atuam em registradores distintos e portanto comutam. Sua projeção conjunta teve probabilidade ideal de cerca de `1,90%` na entrada testada. Contagens na base computacional dão módulos quadrados de amplitudes, de modo que recuperar fases ou sinais exigiria medições adicionais. O texto de fundamentos trata da estrutura lógica dos projetores e deixa em aberto o custo computacional desse circuito. [Artigo de Esfahanian e Arntzen](https://doi.org/10.20935/AcadQuant8546) · [Simulação do circuito](data/P252_PRIME_LCU.json)

## Reprodução

Execute na raiz do repositório com as dependências clássicas da P249.

```bash
python -m pip install -r NUMERIA_ES8P249/requirements.txt
sha256sum -c NUMERIA_ES8P252/SHA256SUMS
python -m unittest NUMERIA_ES8P252.test_spectral_bracket -v
python -O -m unittest NUMERIA_ES8P252.test_spectral_bracket -v
work=$(mktemp -d)
python -m NUMERIA_ES8P252.spectral_bracket --out "$work/P252_SPECTRAL_BRACKET.json"
python -m NUMERIA_ES8P252.odd_pole_audit \
  --matrix NUMERIA_ES8P252/data/P252_GL60.json.gz \
  --margin 1.2379631976818498952630e-23 --out "$work/P252_ODD_POLE_AUDIT.json"
cmp NUMERIA_ES8P252/data/P252_SPECTRAL_BRACKET.json "$work/P252_SPECTRAL_BRACKET.json"
cmp NUMERIA_ES8P252/data/P252_ODD_POLE_AUDIT.json "$work/P252_ODD_POLE_AUDIT.json"
```
