# Nota crítica — teorema do resto da digamma

**Fonte recebida:** `TEOREMA_RESTO_DIGAMA(1).pdf` (2 páginas).
**SHA-256 original:** `03c40c328cb036129dfc5ac73f288ec71b77814da8ba292f64f43e2adaaab267`.

## Enunciado transcrito matematicamente

Para (x>0) real e (K\ge0), escreva
[
\psi(x)=\log x-\frac1{2x}-\sum_{n=1}^K\frac{B_{2n}}{2n x^{2n}}+R_K(x).
]
Então, pela fórmula integral de Binet e expansão geométrica finita,
[
\operatorname{sgn}R_K(x)=(-1)^{K+1},\qquad
|R_K(x)|\le \frac{|B_{2K+2}|}{(2K+2)x^{2K+2}}.
]
A fórmula de Binet é um teorema clássico tomado como premissa; a dedução usa um integrando positivo. O documento relata **104/104 verificações racionais** em uma grade; essas são declarações do documento, **não reproduzidas nesta integração**.

## Correção necessária ao documento original

A tabela da segunda página afirma que (R_8(65)) e (R_8(100)) são positivos. Isso contradiz a fórmula do próprio enunciado: para (K=8), ((-1)^9=-1), portanto ambos os restos são estritamente **negativos**.

Não sobrescrever o PDF original; corrigir a tabela numa versão revisada antes de declarar a nota inteiramente consistente.

## Limitação de aplicabilidade

A prova escrita cobre (x>0) **real**. Ela **não** fornece, por si só, cotas dirigidas para a digamma de argumento complexo (¼+ir/2) da função
(Ω(r)=\Re\psi(¼+ir/2)-\log\pi). A passagem para esse argumento requer uma demonstração específica.

**Status:** `REAL_POSITIVE_DIGAMMA_REMAINDER_FORMULA_WITH_EXTERNAL_BINET_ASSUMPTION`; `PDF_TABLE_K8_SIGN_ERROR`; `COMPLEX_DIGAMMA_NOT_COVERED`.
