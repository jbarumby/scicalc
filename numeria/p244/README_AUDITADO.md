# P244 — reprodução auditada e escopo

A pasta `package/` preserva, sem alterações, os arquivos extraídos do ZIP recebido; seu `SHA256SUMS_P244.txt` continua sendo o manifesto original do pacote. O relatório de síntese e o wrapper integrado ficam fora dessa pasta para não alterar a proveniência recebida.

A execução local refeita em 18/11/2026 exportou os 142.506 centros dos dados P237 versionados, compilou o C/MPFR e refez os dois LDLs/inversas a 1024 bits. Os 754 pivôs foram positivos e os dois rastros recompostos coincidiram byte a byte com os originais. A auditoria por `Fraction` confirmou a margem Weyl de `1/10^18` para todas as matrizes dentro das caixas de entrada. O JSON P244 original também coincidiu byte a byte com a saída local.

O auditor P242 isolado tem uma falha de proveniência: ele lê `inverse_frobenius_squared_upper` do texto de rastreamento e exige apenas que seja positivo. Com os rastros reais e os mesmos pivôs, substituí esse campo por `1`; o auditor P242 aceitou os 754 pivôs, porque o valor adulterado torna a desigualdade de Weyl artificialmente fácil. O wrapper P244 rejeitou essa adulteração por hash. Isso invalida o uso isolado do verificador P242 com rastros não confiáveis, mas não refuta o cálculo original: o replay C/MPFR independente reproduziu os rastros originais.

O wrapper P244 recebido verifica hashes para o C e o verificador Python, mas não fixa diretamente o hash de `export_centers.py`. A cópia integrada em `p244_replay_bound.py` acrescenta pins para esse exportador e para o certificado P242. A versão integrada também distingue “rastro existente conferido” de “MPFR reexecutado nesta chamada”. Ela ainda não fixa compilador, bibliotecas, sistema operacional ou flags completas. Os hashes de saída atestam as saídas desta execução, não todo o ambiente de construção.

**Nota sobre o ambiente reproduzido:** o sandbox não tinha o cabeçalho de desenvolvimento `gmp.h`. Para compilar sem modificar nenhum arquivo recebido, usei um cabeçalho mínimo temporário em `/home/ubuntu/NUMERIA_ES8P242_INCOMING/include` e um symlink de linkedição para a biblioteca GMP de runtime já instalada. A bateria P244 incluiu e passou o caso exato 2×2 e o replay 377×2. Para reprodução independente de longo prazo, instale os pacotes de desenvolvimento GMP/MPFR oficiais e registre as versões do compilador e das bibliotecas.

A auditoria comprova apenas a positividade dos dois blocos finitos de dimensão 377, condicionada à validade das matrizes intervalares P237. Ela não reexecuta de forma independente a montagem das integrais dos 1.008 painéis nem formaliza os lemas analíticos importados. Nenhuma conclusão sobre positividade global de Weil ou RH é autorizada por este replay.

## Reprodução

A partir da raiz do repositório, com Python 3, GCC, `gmp.h`, libGMP e libMPFR disponíveis:

```sh
python3 numeria/p244/package/p244_analytic_gates.py --out /tmp/P244_ANALYTIC_CHECK.json
python3 -m unittest discover -s numeria/p244 -p 'test_p244_bound.py' -v
python3 numeria/p244/p244_replay_bound.py \
  --data numeria/p237/data \
  --traces numeria/p244/package \
  --cert numeria/p244/package/P242_AUDIT_CERT.json \
  --work /tmp/p244-replay \
  --out /tmp/P244_REPLAY_CERT.json
```

A opção `--check-existing` pode ser usada depois de preencher `--work` com os dois rastros. Ela verifica hashes/bytes e a auditoria racional, mas não compila nem executa MPFR nessa chamada.
