# Extensão V20 — Ze Cheng 2015 e a q-digama

**Artigo:** Ze Cheng, *Exact analytical study of ideal Bose atoms in a one-dimensional harmonic trap*, Journal of Statistical Mechanics: Theory and Experiment, 2015, P09011. DOI: https://doi.org/10.1088/1742-5468/2015/09/P09011

A referência foi identificada no arquivo bibliográfico IOPEXPORT_BIB (1).bib enviado. O resumo descreve o emprego da **q-digama** para calcular a população excitada de átomos ideais numa armadilha harmônica unidimensional. O texto integral não foi acessado: o verificador abaixo usa identidades clássicas reconstruídas independentemente, não um cotejo de todas as equações do artigo.

Se 0<q<1, a q-digama para x real positivo é

    psi_q(x) = -log(1-q) + log(q)*sum(n>=0) q^(n+x)/(1-q^(n+x)).

Para energias excitadas E_m=m*hbar*omega e potencial químico mu=0, com q=exp(-beta*hbar*omega):

    N_excited(q) = sum(m>=1) q^m/(1-q^m)
                 = [psi_q(1)+log(1-q)]/log(q).

Para q racional e x inteiro >=1, o resto de S_q(x) após M parcelas (m=x,...,x+M-1), com K=x+M, obedece à cota racional

    q^K/(1-q) <= resto <= q^K/((1-q)*(1-q^K)).

Também vale S_q(1)=sum(n>=1) d(n)q^n, onde d(n) é o número de divisores positivos; para N parcelas o resto não excede q^(N+1)*((N+1)-N*q)/(1-q)^2. O teste cruza as duas implementações, usa intervalos racionais logarítmicos do módulo V20 e verifica psi_q(x+1)-psi_q(x)=(-log q)*q^x/(1-q^x).

**Resultados esperados:** 9 inclusões (q=1/2, 3/4, 9/10; x=1,2,3), 3 conferências independentes pela função número de divisores e 6 recorrências. Executar python v20/test_q_digamma_bose.py. Evidências são geradas e publicadas no GitHub Actions.

**Limites:** Não se deve identificar q-digama real com digama de argumento complexo psi(1/4+i*r/2) na fórmula explícita de Weil. O limite q->1 pode recuperar a digama clássica sob hipóteses próprias, mas não certifica automaticamente a quadratura complexa da Fase 4. Nenhuma prova de RH.
