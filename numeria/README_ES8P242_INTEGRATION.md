# ES8-P242 — validação local e relação com P239

## Resultado da auditoria do pacote recebido

O arquivo recebido `NUMERIA_ES8P242_REPLAY_MPFR_INDEPENDENTE.zip` foi extraído sem executar conteúdo e todos os hashes listados em `numeria/p242/SHA256SUMS_P242.txt` conferiram. O replay foi feito contra os artefatos P237 versionados, incluindo `P237_FINITE_LDL_2048.json.gz` (SHA-256 `2695768a3cca11237dd63cb40c0b6c34427b390329e5ff2ec089e31111c7898a`).

Reproduzi a exportação dos centros, compilei o C/MPFR e executei os dois setores a 1024 bits. Os novos rastros de pivôs/inversa coincidem **byte a byte** com os dois arquivos entregues (`P242_PIVOTS_even_1024.txt` e `P242_PIVOTS_odd_1024.txt`). Os seis testes fornecidos passaram. Em seguida, `p242_verify_exact.py` recontou as 142.506 entradas triangulares/radiuses com `fractions.Fraction`, verificou as duas desigualdades de Weyl e recriou `P242_AUDIT_CERT.json` byte a byte (SHA-256 `9f89c17339c9daca7d1e2e020dbfdbecc83fc47a79301b6ebde544dad0099447`).

O resultado certificado é: 377/377 pivôs positivos em cada paridade; para cada matriz simétrica pertencente às caixas intervalares P237, a cota inferior de autovalor é estritamente maior que `10^-18`. Os valores de apresentação são `δ≈3,6612×10^-47`, `ν≈4,3562×10^34` no setor par e `δ≈3,6576×10^-47`, `ν≈3,6843×10^28` no ímpar; as decisões foram feitas com as desigualdades racionais exatas, não com esses decimais.

## Interpretação e limites

P242 é uma reprodução independente do **certificado finito de matriz**: o LDL e as solves são feitos por uma implementação C/MPFR, sem usar Arb para decidir os sinais, e os raios das caixas são somados com frações racionais exatas. A sobreposição dos pivôs com os certificados Arb é apenas um controle cruzado.

P242 não regenera as 1008 integrais GL60 nem verifica independentemente a identificação da matriz com a forma de Weil. Também não deriva por si só a redução infinita. A etapa P239, em `README_ES8P239_ZHU_EXTENSION.md`, fornece separadamente as cotas analíticas de cauda/acoplamento para a janela `L=0,805`; P242 reforça o bloco finito usado nessa transferência, não substitui P239. Em conjunto, isso sustenta a conclusão na janela compacta sob os insumos analíticos e os artefatos de montagem P237; não é positividade para suportes arbitrários e não prova RH.

## Observação de compilação neste sandbox

O ambiente não tinha `libgmp-dev` (cabeçalho `gmp.h` e link de desenvolvimento). Para executar sem alterar o pacote recebido, usei um cabeçalho mínimo de declaração ABI e um symlink temporários, fora do repositório, vinculados às bibliotecas de runtime instaladas `libgmp.so.10` e `libmpfr.so.6`. A reprodução byte a byte e os testes passaram nesse modo. Para reprodução de manutenção, siga a rota recomendada no `numeria/p242/README.md` e instale `libgmp-dev`; os arquivos originais do pacote P242 foram preservados sem alterações.

Elaborado por Manus AI para o projeto NUMERIA.
