"""NUMERIA P78-S1 (prototipo) — Sensibilidade do B1 aos b4 publicos do P78.

Executa o item 3 do CODEX_NEXT_STEPS do P78 como PROTOTIPO: novo auditor
com b4 por linha, SEM editar o NUMERIA_P77_B1 original. Reusa
constants/TABLE do P76 e a mesma quadratura rigorosa acb.integral do P77.

K1 e K2 nao dependem de b4: calculados uma vez por linha; a algebra
b4-dependente e avaliada com b4 = 0.54 (baseline P77) e com o b4 por linha
do certificado P78 (fonte publica Ramare-Zuniga v1).

ATENCAO: prototipo single-run, nao substitui o certificado dyadic + CI +
SHA256 do pipeline NUMERIA. O sinal do resultado e estruturalmente robusto
(k3 > 1 em todas as linhas => canais max(1, k3^{-...}) saturados em 1 =>
B1 estritamente crescente em b4), mas os VALORES exatos devem ser
re-certificados no padrao do P77.

Uso: python3 p78_b1_sensitivity.py  (rodar de qualquer lugar do repo)
"""
import json
import sys
from pathlib import Path


def find_repo():
    """Sobe na arvore a partir deste arquivo ate achar NUMERIA_P76_COEFFICIENTS."""
    here = Path(__file__).resolve().parent
    for cand in [here] + list(here.parents):
        if (cand / "NUMERIA_P76_COEFFICIENTS").is_dir():
            return cand
    return Path("/home/z/my-project/scicalc-p78")  # fallback sandbox


REPO = find_repo()
sys.path.insert(0, str(REPO))

from flint import acb, arb, ctx  # noqa: E402
from NUMERIA_P76_COEFFICIENTS.audit_b23 import constants, q, dyadic  # noqa: E402
from NUMERIA_P76_COEFFICIENTS.table import TABLE, validate_table  # noqa: E402

T0 = 3 * 10**12


def rat(n, d):
    return arb(n) / arb(d)


def pw(x, e):
    if not x > 0:
        raise ArithmeticError("non-positive real power base")
    return (x.log() * e).exp()


def row_parts(row, cs):
    """Partes independentes de b4 (A0, a1, a2, c...)."""
    d, kap, delta = q(row.d), q(row.kappa), q(row.delta)
    L, pi, m0 = cs["logh"], cs["pi"], cs["m0"]
    A0 = (arb(10).log() * q(row.log10_a0)).exp()
    La = A0.log()
    half = (A0 / 2).log()
    a1 = A0 / (2 * pi * pi) * half**4 + q("34.842") * A0 * half**3 + arb("1.3897e16")
    r = 1 - arb.const_log2() / La
    a2 = r**4 / (2 * pi * pi) + q("34.842") / La * r**3 + arb("9.5079e15") / (A0 * La**4)
    c = 3 / (4 * kap)
    return d, kap, delta, L, pi, m0, A0, a1, a2, c


def k1_k2(row, cs, parts):
    """K1 (integral rigorosa + cauda elementar) e K2 — idem P77, sem b4."""
    d, kap, delta, L, pi, m0, A0, a1, a2, c = parts

    def f(x, analytic):
        v = x**6
        A = acb(1) + (acb(1) + acb(c) * v).log(analytic=analytic) / (2 * acb(L))
        B = acb(cs["c1"] * c.sqrt()) * x**3 + acb(cs["c2"])
        return 6 * v * (-v).exp() * A.pow(acb(4) / 3, analytic=analytic) * B.pow(acb(2) / 3, analytic=analytic)

    I = acb.integral(f, 0, 2, rel_tol=arb('1e-19'),
                     abs_tol=arb('1e-19'), eval_limit=30000, depth_limit=30)
    if not I.is_finite() or not I.imag.contains(0) or not I.real.lower() > 0:
        raise ArithmeticError("quadrature failed or inconclusive")
    U = arb(64)
    D = (1 + c / (2 * L))**2 * (1 + cs["c1"] * c.sqrt() + cs["c2"])
    tail = (D * (-U).exp() * (U**3 + 3 * U**2 + 6 * U + 6)).upper()
    lower, upper = I.real.lower(), (I.real.upper() + tail).upper()
    F = pw(2 * a2 * c.sqrt(), rat(1, 3))
    correction = (pw(2 * a1, rat(1, 3)) * pw(cs["c2"] + A0 * cs["c1"] / T0, rat(2, 3))
                  / (pw(arb(T0), rat(1, 3)) * pw(L, rat(4, 3))))
    k1lo = (F * lower + correction).lower()
    k1 = (F * upper + correction).upper()
    if not k1lo > 0:
        raise ArithmeticError("K1 lower bound inconclusive")
    two13 = pw(arb(2), rat(1, 3))
    k2 = two13 * c.sqrt() * pi.sqrt() / 2 + two13 * (4 * kap / 3) * pw(A0 / arb(T0), arb(3))
    return k1, k2


def b1_upper(row, cs, parts, k1, k2, b4):
    """Algebra b4-dependente, identica a P77 mas com b4 parametrizado."""
    d, kap, delta, L, pi, m0, A0, a1, a2, c = parts
    s1, s2 = rat(row.sigma_start, 32), rat(row.sigma_start + 1, 32)
    c3 = pi * m0 * b4 / (2 * delta) * (2 * delta * arb.const_euler() / L).exp()
    c4 = b4 / (5 * delta * delta.exp()) + q(".605") * (-2 * delta).exp() / L
    k3 = c3 + c4 * (pi.sqrt() / (2 * (2 * kap).sqrt()) + pi * m0 / T0)
    l0 = (kap * (1 + (1 / (4 * arb(T0)) * (1 + 2 * delta / L))**2)).exp() / (1 - arb(1) / T0)
    l1 = max(arb(1), pw(k1, 3 * delta / L)) * max(arb(1), pw(k3, -3 * delta / L)) * (3 * delta + 4 * delta * L.log() / L).exp()
    l2 = max(arb(1), pw(k2, 3 * delta / L)) * max(arb(1), pw(k3, -3 * delta / L)) * (3 * delta).exp()
    common = max(arb(1), pw(l0, -2 * d / L)) * max(arb(1), pw(k3, -3 * d / L))
    common *= max(pw(l0, 2 / (2 - s1)), pw(l0, 2 / (2 - s2)))
    common *= max(pw(k3, (2 * s1 - 1) / (2 - s1)), pw(k3, (2 * s2 - 1) / (2 - s2)))
    A1 = l1 * (d * (3 + 4 * L.log() / L)).exp() * common
    A1 *= max(arb(1), pw(k1, 3 * d / L))
    A1 *= max(pw(k1, 3 * (1 - s1) / (2 - s1)), pw(k1, 3 * (1 - s2) / (2 - s2)))
    A2 = l2 * (3 * d).exp() * common
    A2 *= max(arb(1), pw(k2, 3 * d / L))
    A2 *= max(pw(k2, 3 * (1 - s1) / (2 - s1)), pw(k2, 3 * (1 - s2) / (2 - s2)))
    B1 = (A1 + A2 * pw(L, -6 * (1 - s2) / (2 - s2))) / (2 * pi * d)
    return B1.upper(), k3


def main():
    cert = json.loads((REPO / "NUMERIA_P78_RAMARE_PUBLIC" / "certificate.json").read_text())
    p78_b4 = {r["row"]: r["b4"] for r in cert["rows"]}
    from fractions import Fraction as Fr
    validate_table()
    out = []
    with ctx.workprec(192):
        cs = constants()
        cs["c1"] = 6 / cs["pi"]**2 + q("1.048") / cs["logh"]
        cs["c2"] = (cs["pi"] * cs["m0"] * q(".62") / cs["logh"]
                    + 6 * cs["m0"] / (cs["pi"] * T0)
                    + cs["pi"] * cs["m0"] * q("1.048") / (arb(T0) * cs["logh"]))
        for idx in range(10, 17):  # apenas as linhas com b4 substituido no P78
            row = TABLE[idx - 1]
            parts = row_parts(row, cs)
            k1, k2 = k1_k2(row, cs, parts)
            b54 = q(".54")
            fr78 = Fr(p78_b4[idx])
            b78 = arb(fr78.numerator) / arb(fr78.denominator)
            up54, k3_54 = b1_upper(row, cs, parts, k1, k2, b54)
            up78, k3_78 = b1_upper(row, cs, parts, k1, k2, b78)
            printed = q(row.b1)
            m54 = (printed - up54).lower()
            m78 = (printed - up78).lower()
            out.append({
                "row": idx, "delta": row.delta, "b4_P78": p78_b4[idx],
                "printed_B1": row.b1,
                "B1_upper_b4_0.54": float(up54), "B1_upper_b4_P78": float(up78),
                "margin_0.54": float(m54), "margin_P78": float(m78),
                "strict_0.54": bool(m54 > 0), "strict_P78": bool(m78 > 0),
                "k3_at_0.54": float(k3_54), "k3_at_P78": float(k3_78),
                "delta_B1_pct": 100 * (float(up78) - float(up54)) / float(up54),
            })
    print(json.dumps(out, indent=1))
    fails = [o["row"] for o in out if not o["strict_P78"]]
    print("ROWS FAILING WITH P78 b4:", fails if fails else "none")
    (REPO.parent / "p78_b1_sensitivity_results.json").write_text(
        json.dumps({"schema": "NUMERIA.P78.S1.PROTOTYPE", "rows": out,
                    "rows_failing_with_public_b4": fails}, indent=1) + "\n")
    if fails:
        raise SystemExit("PROTOTYPE GATE: printed B1 not majorized by public b4 in rows " + str(fails))


if __name__ == "__main__":
    main()