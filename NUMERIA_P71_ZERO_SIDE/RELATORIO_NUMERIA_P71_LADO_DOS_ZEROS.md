# NUMERIA–RIEMANN — Auditoria independente após o PR #71

**Data:** 21/11/2026  
**Repositório:** https://github.com/jbarumby/scicalc/pull/71  
**Natureza:** prova computacional finita pelo lado dos zeros, com importação explícita de resultados publicados.  
**Mudanças no GitHub:** nenhuma (os arquivos deste pacote foram executados localmente).

## 1. Conclusão e alcance

Foi construído e executado um certificador **sem enumeração de primos nem potências primas**, utilizando apenas os três primeiros zeros não triviais da função zeta, intervalos inteiros e desigualdades racionais. O programa estabelece a seguinte desigualdade **para todos os reais** no intervalo examinado:

\[
\boxed{\Psi(t)>0.00012255\qquad(18\le t\le34.5).}
\]

Na combinação com o certificado do PR #71, que registra \(\Psi(t)>0\) em \(0<t\le\log(10^8)\), temos cobertura sobreposta e contínua de \(0<t\le34.5\), **desde que** se aceitem: (i) as verificações e os dados importados de Suzuki, Platt–Trudgian, FLINT e Bellotti–Wong, e (ii) o certificado anterior do PR #71. A função satisfaz \(\Psi(0)=0\).

O extremo \(t=34.5\) corresponde a \(X=e^{34.5}\approx9.619657855447764\times10^{14}\), contra o teto anterior \(10^8\). Isso significa uma ampliação por fator aproximado de \(9.62\) milhões no limite de \(X\), **sem contar novos termos aritméticos**. Não é uma demonstração da hipótese de Riemann nem um certificado para a cauda infinita.

Esta é uma reprodução independente da estratégia espectral dentro do NUMERIA; não é reivindicação de precedência científica, pois Mittermeier já relatou verificação finita pelo lado dos zeros além desta altura.

## 2. Equação de partida e dependências matemáticas

Suzuki (2023), equação (1.3), fornece a expansão incondicional

\[
\Psi(t)=\sum_{\gamma}\frac{1-\cos(\gamma t)}{\gamma^2},
\]

com \(\gamma\) percorrendo os zeros, com multiplicidades, de \(\xi(\frac12-iz)\). Os \(\gamma\) **podem ser complexos** se RH for falsa. Portanto, não se pode presumir positividade de cada parcela para todos os zeros.

Platt e Trudgian verificaram rigorosamente que todos os zeros de zeta de altura \(0<\Im\rho\le H=3\cdot10^{12}\) estão na reta crítica. Isso torna **não negativas** suas contribuições na soma de Suzuki. São retidos somente os três primeiros pares \(\pm\gamma_j\), com intervalos publicados na documentação de FLINT (`acb_dirichlet_zeta_zeros`), alargados para \(\gamma_j\in[c_j-10^{-26},c_j+10^{-26}]\):

- \(c_1=14.1347251417346937904572519836\)
- \(c_2=21.0220396387715549926284795939\)
- \(c_3=25.0108575801456887632137909926\)

Os demais zeros *abaixo* de \(H\) são ignorados como parcelas não negativas. Para os zeros *acima* de \(H\), usa-se um erro absoluto, sem pressupor que estejam na reta crítica.

A cota explícita do número de zeros de Bellotti–Wong (2025), *Mathematics of Computation*,

\[
\left|N(T)-\frac{T}{2\pi}\log\frac{T}{2\pi e}\right|
\le 0.10076\log T+0.24460\log\log T+8.08344\quad(T\ge e),
\]

permite deduzir a majorante mais simples \(N(T)\le T\log T/(2\pi)\) para \(T\ge H\). Esta passagem é segura, pois a diferença dos termos principais é linear em \(T\), enquanto o erro explícito é logarítmico. O script verifica as desigualdades numéricas auxiliares com frações exatas.

Por integração de Stieltjes:

\[
\sum_{\Im\rho>H}\frac{1}{(\Im\rho)^2}
\le 2\int_H^\infty\frac{N(u)}{u^3}\,du
\le\frac{\log H+1}{\pi H}<\frac{31}{3H}.
\]

Como \(|\Im\gamma|<1/2\) no plano espectral e \(|1-\cos(\gamma t)|\le1+\cosh(t/2)\), a soma de ambos os sinais com altura superior a \(H\) fica majorada por

\[
E_H(t)\le2(1+\cosh(t/2))\frac{31}{3H}.
\]

Para \(0\le t\le34.5\), o script prova com inteiros que \(e^{17.25}<2^{25}\), obtendo

\[
\boxed{E_H(t)<0.000231152761.}
\]

Este limite é deliberadamente conservador, mas suficiente para o sinal positivo.

## 3. Certificação racional em toda a faixa

Definimos a contribuição positiva dos três pares certificados:

\[
S_3(t)=2\sum_{j=1}^3\frac{1-\cos(\gamma_jt)}{\gamma_j^2}.
\]

O algoritmo verificou 16.501 pontos racionais \(t_k=18+k/1000\), para \(0\le k\le16.500\), usando os passos abaixo:

1. \(\pi\) foi enquadrado rigorosamente pela fórmula de Machin com séries alternadas de arctangente e resto limitado pela próxima parcela.
2. Os produtos \(\gamma_j t_k\) foram reduzidos módulo \(2\pi\) com intervalos inteiros; foi verificado, ponto a ponto, que o argumento reduzido está em \([-4,4]\).
3. O cosseno foi incluído por polinômio de Taylor até grau 80, com resto absoluto inferior a \(4^{82}/82!<10^{-50}\); todas as operações intermediárias usam arredondamentos exteriores, por divisão inteira.
4. As três parcelas foram incluídas inferiormente sem operação em ponto flutuante. Foi calculado o SHA-256 determinístico de todas as 16.501 linhas de limites inferiores.
5. A desigualdade \(|S_3'(t)|\le2(1/14+1/21+1/25)\) controla o intervalo entre pontos consecutivos: a distância até a malha é no máximo \(0.0005\).

| Componente | Valor |
|:--|--:|
| Mínimo inferior na malha de \(S_3\) | **0,000512751726189403** |
| Ponto de mínimo registrado | \(t=32,894\) |
| Perda máxima entre pontos | **0,000159047619047619** |
| Erro superior de todos os zeros acima de \(H\) | **0,000231152760666667** |
| **Margem uniforme final** | **0,000122551346475117** |

Segue \(\Psi(t)\ge S_3(t)-E_H(t)>0.00012255\) em todo \([18,34.5]\). A malha não é uma amostragem heurística: cada ponto tem um limite inferior exato, e a desigualdade de Lipschitz liga os pontos em um intervalo contínuo.

**SHA-256 das 16.501 linhas de limites inferiores (formato `t_milli,lower_scaled_integer\\n`):**  
`726d525c5b92540c290ccb2831a39547ee91280361c1599e7a5aa7d3c0cf3da9`

## 4. Execução e integridade

Executado com Python 3.13, em modo normal e com otimizações `-O`, com resultados idênticos. Não há dependências de terceiros, acesso à rede, números aleatórios, pacotes de aritmética intervalar externos ou enumeração de primos. Todas as condições críticas usam `raise`, não `assert`, e continuam ativas com `-O`.

```bash
python numeria_pr71_zero_side_audit.py
python -O numeria_pr71_zero_side_audit.py
```

A execução grava `NUMERIA_P71_ZERO_SIDE_AUDIT.json` ao lado do script. Recomenda-se comparar o `grid_sha256` e os demais campos com o JSON do pacote.

## 5. Limitações, integração e próximos gates

**Dependências externas não reproduzidas localmente:** o script importa (a) a identidade incondicional de Suzuki; (b) a verificação da reta crítica até \(3\cdot10^{12}\) de Platt–Trudgian; (c) os três enclosures de zeros da documentação FLINT; (d) a cota de contagem de Bellotti–Wong. Embora essas fontes sejam estabelecidas, uma auditoria *sem dependências externas* teria que gerar novamente os três zeros em Arb e revisar formalmente os lemas citados. O ambiente de execução aqui não dispõe de `python-flint`, então a auditoria independente não repetiu os 5.762.859 prefixos do PR #71; apenas conferiu a estrutura do PR e produziu a extensão espectral por caminho independente.

**Próxima extensão útil:** adicionar o quarto e seguintes zeros com enclosures certificados independentes, preferencialmente com FLINT/Arb, e otimizar a prova contínua com derivadas locais, em vez de usar apenas o maior limite global de Lipschitz. Depois, comparar as margens com Mittermeier e Caceres, sem usar os resultados desses preprints como axiomas do certificador.

## 6. Fontes consultadas

- Masatoshi Suzuki, *Aspects of the screw function corresponding to the Riemann zeta-function*, JLMS 2023, DOI: https://doi.org/10.1112/jlms.12785 .
- David Platt e Tim Trudgian, *The Riemann hypothesis is true up to 3·10¹²*, BLMS 2021, DOI: https://doi.org/10.1112/blms.12460 .
- Chiara Bellotti e Peng-Jie Wong, *Improved estimates for the argument and zero-counting function of the Riemann zeta-function*, *Mathematics of Computation* 2025, DOI: https://doi.org/10.1090/mcom/4133 .
- Documentação FLINT, `acb_dirichlet_zeta_zeros`: https://github.com/flintlib/flint/blob/main/doc/src/acb_dirichlet/zeta-zeros.md .
- Rainer Andreas Mittermeier, certificado pelo lado dos zeros até \(t\approx47.48\): https://zenodo.org/records/23481975 .
- Pedro Caceres, verificação finita até \(2.15\cdot10^{29}\) relatada em preprint: https://hexagonmath.org/2610.00661 .

**Estado:** RESULTADO FINITO POSITIVO, COM INSUMOS MATEMÁTICOS EXTERNOS EXPLICITADOS; RH GLOBAL E CAUDA INFINITA PERMANECEM EM ABERTO.
