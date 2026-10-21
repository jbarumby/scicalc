"""Verificador independente e EXATO do certificado, a partir do JSON exportado.

Não usa Arb nem ponto flutuante: lê os extremos decimais (arredondados para fora)
de ``projected_matrix_enclosures`` como racionais exatos e refaz o teste tipo Rump
em aritmética de Fractions:

  1. hull simétrico de cada par (i,j)/(j,i);  M = ponto médio,  rho = meio-raio;
  2. delta = max_i sum_j rho_ij  (cota de ||E||_2 para E simétrica, |E_ij| <= rho_ij);
  3. LDL exato de M - delta*I: se todos os pivôs são > 0, toda matriz simétrica dentro
     do enclosure é definida positiva.

Isto só valida o PASSO FINAL (enclosure exportado => positividade). A correção dos
enclosures vem do Arb e da matemática documentada em TAIL_PROOF.md e
WEIL_FORM_DERIVATION_PT.md. Opcionalmente confere se dois relatórios (ex.: 192 e 256
bits) são consistentes (todos os pares de enclosures se intersectam).

Uso:
    python -m numeria_riemann.verify_certificate RELATORIO.json [--compare OUTRO.json]
Código de saída 0 = verificado (e consistente com o status), 1 = não verificado.
"""
import argparse
import json
import sys
from decimal import Decimal
from fractions import Fraction
from pathlib import Path


def _frac(text):
    return Fraction(Decimal(text))


def load_matrix(report):
    """Matriz de intervalos [(lo, hi)] exatos, com hull simétrico."""
    rows = report["projected_matrix_enclosures"]
    n = len(rows)
    lo = [[_frac(e["lower"]) for e in row] for row in rows]
    hi = [[_frac(e["upper"]) for e in row] for row in rows]
    for i in range(n):
        for j in range(i + 1, n):
            lo[i][j] = lo[j][i] = min(lo[i][j], lo[j][i])
            hi[i][j] = hi[j][i] = max(hi[i][j], hi[j][i])
    return lo, hi


def exact_ldl_pivots(matrix):
    """LDL^T exato (Fractions). Retorna (pivôs, índice da falha ou None)."""
    n = len(matrix)
    lower = [[Fraction(0)] * n for _ in range(n)]
    diag = []
    for i in range(n):
        lower[i][i] = Fraction(1)
        pivot = matrix[i][i] - sum(lower[i][k] ** 2 * diag[k] for k in range(i))
        diag.append(pivot)
        if pivot <= 0:
            return diag, i
        for j in range(i + 1, n):
            lower[j][i] = (matrix[j][i] - sum(lower[j][k] * lower[i][k] * diag[k] for k in range(i))) / pivot
    return diag, None


def verify(report):
    lo, hi = load_matrix(report)
    n = len(lo)
    mid = [[(lo[i][j] + hi[i][j]) / 2 for j in range(n)] for i in range(n)]
    rho = [[(hi[i][j] - lo[i][j]) / 2 for j in range(n)] for i in range(n)]
    delta = max(sum(rho[i]) for i in range(n))
    shifted = [[mid[i][j] - (delta if i == j else 0) for j in range(n)] for i in range(n)]
    pivots, fail = exact_ldl_pivots(shifted)
    return {
        "verified": fail is None,
        "dimension": n,
        "delta": float(delta),
        "min_shifted_pivot": float(min(pivots)),
        "failed_at_pivot": None if fail is None else fail + 1,
    }


def compare(report_a, report_b):
    """Todos os pares de enclosures correspondentes devem se intersectar."""
    la, ha = load_matrix(report_a)
    lb, hb = load_matrix(report_b)
    if len(la) != len(lb):
        return False
    n = len(la)
    return all(la[i][j] <= hb[i][j] and lb[i][j] <= ha[i][j] for i in range(n) for j in range(n))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("report", type=Path)
    parser.add_argument("--compare", type=Path, help="outro relatório para checar consistência dos enclosures")
    args = parser.parse_args(argv)
    report = json.loads(args.report.read_text(encoding="utf-8"))
    result = verify(report)
    claimed = report.get("status") == "certified_positive_definite_finite_matrix"
    result["agrees_with_declared_status"] = result["verified"] == claimed
    ok = result["verified"] and claimed
    if args.compare:
        other = json.loads(args.compare.read_text(encoding="utf-8"))
        result["consistent_with_other_report"] = compare(report, other)
        ok = ok and result["consistent_with_other_report"]
    print(json.dumps(result, indent=2))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
