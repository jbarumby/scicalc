# Diagnóstico de auditoria intervalar

- Status declarado: `certified_positive_definite_finite_matrix`
- N: `6`
- Integração: `adaptive`
- Células de quadratura: `não se aplica`
- Dimensão projetada: `4`
- Pivôs intervalares registrados: `4` de `4`

## Certificados

- `interval_ldl`: passou
- `rump_margin`: passou, delta=0.012675985577517

## Pivôs (limite inferior / superior / largura)

- Pivô 1: 0.229776478636077085417037188964238015457572 / 0.233407584202532602381234347777714577957573 / 0.003631105566455516964197158813
- Pivô 2: 0.888368308536134037969471796444125420942530 / 0.896521159570671927032830103328891045942531 / 0.008152851034537889063358306885
- Pivô 3: 0.502201242134765198308179338235133581409365 / 0.537649641914617678004453141946071081409366 / 0.03544839977985247969627380371
- Pivô 4: 1.125561870654268584957006933900198799136464 / 1.141181597243916638007048132630667549136465 / 0.01561972658964805305004119873

## Interpretação

O relatório declara certificação para a **matriz finita restrita** dos parâmetros acima; este diagnóstico não revalida a aritmética de intervalos.
Diagnóstico NÃO rigoroso (autovalores do ponto médio): 0.132604851662, 0.885612591447, 0.908189275555, 1.1421176341.

Escopo: matrizes finitas; isto não constitui demonstração da hipótese de Riemann.
