# ES8-P242 — Auditor MPFR independente, centro e Weyl

**Resultado obtido:** 754/754 pivôs positivos, calculados por uma nova implementação em C/MPFR 1024 bits, com enclosures dirigidos; 2/2 inversas certificadas por solves triangulares dirigidos; 142.506 entradas de raio somadas em frações racionais exatas, sem reutilizar os raios já somados em P238. A condição `(delta + 1e-18)^2 * nu < 1` passou em ambos os setores, provando positividade com margem >1e-18 para toda matriz nas caixas P237 **como serializadas**.

**Não comprovado aqui:** montagem integral desde os painéis Arb, formula explícita e transferência Weil–Zhu, modelo infinito, RH.

Conteúdo: código-fonte C da nova aritmética intervalar, exportador de centros com check SHA, auditor Python racional, testes, dois traços integrais de 377 pivôs/limite da inversa, `P242_AUDIT_CERT.json`, prova e `SHA256SUMS_P242.txt`. Não inclui os enormes arquivos P237, que precisam ser obtidos em `jbarumby/scicalc`, branch `numeria-p228-arb08-rational-audit`, pasta `numeria/p237/data/`, ou a partir dos uploads originais.

Ver `PROVA_ES8P242_CERTIFICADO_MPFR.md` para o argumento, as margens e os comandos de reprodução.
