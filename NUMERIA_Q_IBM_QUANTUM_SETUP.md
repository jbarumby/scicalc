# NUMERIA-Q — IBM Quantum: teste de acesso e um qubit

Notebook: [NUMERIA_Q_IBM_Quantum_One_Qubit.ipynb](./NUMERIA_Q_IBM_Quantum_One_Qubit.ipynb)

## Procedimento
1. Abra o notebook no Google Colab e salve uma cópia no Drive.
2. Execute a instalação do Qiskit e do `qiskit-ibm-runtime`.
3. Crie/acesse sua conta IBM Quantum Platform e obtenha a API key e, se aplicável, o CRN da instância.
4. Insira a chave apenas no prompt oculto `getpass`; não a coloque no código ou no GitHub.
5. Execute a consulta de QPUs. Isso não submete circuitos.
6. Confira o backend/fila. A execução física exige digitar `ENVIAR` e envia 20 shots de um circuito X+medição.

O notebook não faz fallback silencioso para simulação. Se não houver QPU operacional visível, interrompe a execução. Um job físico pode consumir a cota do plano e entrar em fila.

## Fontes oficiais
- Guia de primeiro circuito físico: https://quantum.cloud.ibm.com/docs/en/guides/hello-world
- Planos de acesso: https://quantum.cloud.ibm.com/docs/en/guides/plans-overview
- Cliente Qiskit Runtime: https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/runtime-service

## Integridade NUMERIA-Q
Branch isolada: `research/ibm-quantum-one-qubit`. Nenhuma alteração intencional em `main`, NUMERIA-Q 05.19.1 ou 05.20. Nenhum job foi enviado ao criar este notebook.
