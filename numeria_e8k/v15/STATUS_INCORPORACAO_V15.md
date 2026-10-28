# Status de incorporação — NUMERIA E8-K V15

Data: 28/10/2026

## Pacotes recebidos

- `E8K_V15_RESULTADOS.zip` — SHA-256 do ZIP: a registrar no ambiente de origem.
- `E8K_V15_GITHUB_COMPLETO.zip` — pacote de integração completo.

Ambos os ZIPs passaram no teste estrutural de integridade (nenhum membro corrompido). O pacote completo contém 35 entradas, incluindo os artefatos históricos V14, os dados B5, a amostra de oráculo, a formalização Lean e o workflow V15.

## Validação local do V15

Executado a partir do pacote extraído:

```text
V15: CONDITIONAL_2X2_ARCHIMEDEAN_NEGATIVE_DEFINITE
Gates: 11 PASS: 11
G11/G22 margem: G11:1.42817374775339414933705
G22:0.102380006120156803079048
```

O processo terminou com código de saída 0. Houve uma mensagem de aquecimento de ambiente de planilhas no stderr, sem impacto no resultado do verificador.

Hashes declarados pelo manifesto V15:
- `v15/VERIFICAR_E8K_V15_BLOCO.py`: `c870d6808138f7e28f7b3cf9f0e6802a22d8e9090298f915268678455ca3f2ef`
- `v15/RESULTADO_E8K_V15.json`: `abd51fbd33aeb9dc531212a7a106517de6aaa8b037457c1f90eab5e584a72a41`
- `v15/PARECER_TECNICO_V15.md`: `82375112998e14e20efd867df5170f8abb1798541a3b4cc0eb3bb5c0ab0e4288`
- `v15/README_V15.md`: `d3bc2611921de857842aced0b00374d493dbe692b05901190e3a46e3ec2974cd`
- `.github/workflows/auditoria-e8k-v15.yml`: `982027284196f29d5c36f5c0a435f5326e5cb461a25a21e0477f98e1d477fcc4`

## Estado GitHub — importante

Este documento registra a validação local; **não afirma que os 35 artefatos do ZIP foram carregados**. Nesta sessão, a integração disponível não expõe uma operação de envio de arquivos binários locais diretamente ao repositório. O verificador V15, o JSON de resultado, os dados B5 e o workflow precisam ser publicados e conferidos no GitHub antes de declarar a incorporação completa.

## Interpretação científica

O status `CONDITIONAL_2X2_ARCHIMEDEAN_NEGATIVE_DEFINITE` é condicional às premissas analíticas documentadas. Não certifica a fórmula global de Weil nem demonstra a Hipótese de Riemann. A compilação Lean também permanece pendente de execução independente.
