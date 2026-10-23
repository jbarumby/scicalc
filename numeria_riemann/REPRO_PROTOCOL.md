# NUMERIA–RIEMANN — protocolo de auditoria

Data: 2026-10-22  
Branch isolada: audit/guinand-arb-repro-2026-10-22  
Base: numeria-riemann/overlap-normalization-fix-v2  
Commit-base: 3219f6696c54c7f6e250dd7d41f88d48719b88db

## Guinand–Weil

Bombieri (2000), seção 2, usa função suave de suporte compacto no lado Mellin. Para a função de teste do protótipo, a referência adequada é a formulação geral da Proposição 5 de Chirre–Gonçalves (2022): h par, analítica na faixa |Im(s)| <= 1/2+epsilon e decaimento O((1+|s|)^-(1+delta)) para delta positivo.

O código usa uma combinação finita real de senos que se anula nos extremos. Duas integrações por partes dão F(z)=O((1+|Re z|)^-2) em faixas fixas; h(z)=F(z)F(-z) é par, inteira e O((1+|Re z|)^-4). A base estendida por zero não é necessariamente C-infinito, portanto não se afirma aplicação literal do enunciado de Bombieri; a aplicação da versão geral depende da estimativa uniforme para a base exata.

A soma inclui todos os zeros não triviais, com multiplicidades e ambos os sinais. gamma_rho=(rho-1/2)/i pode ser complexo; não se assume RH nem se troca F(gamma)F(-gamma) por |F(gamma)|^2 quando gamma é complexo. A convenção de Fourier converte o termo de primos para -2 sum Lambda(n)g(log n)/sqrt(n); os termos de polo são h(i/2)+h(-i/2)=2h(i/2). A cota de cauda e a projeção ainda requerem revisão independente.

## Arb e sinc

A versão fixada é python-flint==0.9.0. A matriz usa arb real e somas de intervalos; não chama acb.integral. A API complexa é testada separadamente. A documentação de acb.integral define integração no segmento reto e exige que a integranda certifique analiticidade da bola quando analytic=True.

Para sinc(x)=sin(x)/x, se a bola contém zero usa-se série de Taylor até x^16/17!, com resto inferior a 1e-34 quando abs_upper(x)<0.1; bola larga gera ValueError. O pré-controle run usa L*step/2<0.1, mas a condição global segura é L*step<0.1. A configuração proposta L=0.45 e step=0.01 satisfaz a condição segura. Foi adicionado test_arb_contract.py; ainda falta conectá-lo ao entry point da CI.

## Reprodução pendente

Parâmetros planejados: Python 3.11, python-flint==0.9.0, L=0.45, R=160, step=0.01, 256 termos; N=6,8,10 e 192/256 bits. O workflow existente não foi alterado nesta etapa, e não há execução CI dos seis casos nesta branch. Portanto, os seis status, códigos de saída e hashes permanecem pendentes; não afirmar que os certificados foram reproduzidos.

Certificados anteriores à correção de normalização estão substituídos. Não alterar N/L/R para forçar sucesso; caso inconclusivo não é autovalor negativo. Não alegar positividade global de Weil nem prova de RH.
