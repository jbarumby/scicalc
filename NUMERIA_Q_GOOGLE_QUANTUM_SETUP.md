# NUMERIA-Q — conexão Google Quantum AI (teste 1 qubit)

Arquivo: [NUMERIA_Q_Google_Quantum_AI_One_Qubit.ipynb](./NUMERIA_Q_Google_Quantum_AI_One_Qubit.ipynb)

## Uso
1. Abra o notebook no Colab e salve uma cópia no seu Drive.
2. Instale Cirq/cirq-google executando a primeira célula.
3. Preencha o ID do projeto Google Cloud.
4. Autentique a conta Google e execute a célula que lista processadores.
5. Se houver um processador autorizado, copie seu ID e execute a célula explicitamente marcada EXECUTAR NO HARDWARE.

O notebook não faz fallback para simulação: se não houver acesso, deve reportar o erro e parar. O circuito físico aplica X e mede um qubit, com 20 repetições.

## Requisitos de acesso
O acesso ao Google Quantum AI é restrito a uma lista aprovada. É necessário projeto Google Cloud, Quantum Engine API habilitada, IAM adequado e autorização para o serviço. Consulte:
- https://quantumai.google/cirq/google/access
- https://quantumai.google/cirq/tutorials/google/start
- https://quantumai.google/cirq/google/engine

## Integridade NUMERIA-Q
Este notebook foi adicionado numa branch separada, sem editar nem mesclar a branch main, 05.19.1 ou 05.20. Não contém credenciais. Nenhum job físico foi enviado durante a criação deste notebook.
