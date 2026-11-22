"""NUMERIA P78-S1 (prototipo aprimorado) — Sensibilidade rigorosa do B1 aos b4 publicos do P78.

Melhorias de rigor intervalar:
1. Conversao direcional estrita para float64 (ARF_RND_CEIL / ARF_RND_FLOOR) em vez de float() midpoint.
2. Exportacao das bolas arb completas (midpoint +/- radius e extremos [lo, hi] em string).
3. Cauda da integral com expoentes exatos (4/3 e 2/3) e validacao explicita de monotonicidade/max.
4. Calculo do novo teto sugerido de B1 (2 casas decimais) para re-certificacao no P78.
"""
from fractions import Fraction as Fr
import json
import math
import sys
from pathlib import Path


def find_repo():
    """Sobe na arvore a partir deste arquivo ate achar NUMERIA_P76_COEFFICIENTS."""
    here = Path(__file__).resolve().parent
    for cand in [here] + list(here.parents):
        if (cand / "NUMERIA_P76_COEFFICIENTS").is_dir():
            return cand
    return Path("/home/z/my-project/scicalc-p78")


REPO = find_repo()
sys.path.insert(0, str(REPO))

from flint import acb, arb, ctx  # noqa: E402
from NUMERIA_P76_COEFFICIENTS.audit_b23 import constants, q  # noqa: E402
from NUMERIA_P76_COEFFICIENTS.table import TABLE, validate_table  # noqa: E402

T0 = 3 * 10**12


def rat(n, d):
    return arb(n) / arb(d)


def pw(x, e):
    if not x > 0:
        raise ArithmeticError("non-positive real power base")
    return (x.log() * e).exp()


def arb_upper_f64(x: arb) -> float:
    """Majorante float64 rigoroso (nunca arredonda para baixo)."""
    u = x.upper()
    # Garante margem de 1 ULP acima caso a conversao decimal/binaria aproxime
    val = float(u.str(17, radius=False))
    return math.nextafter(val, math.inf) if arb(val) < u else val


def arb_lower_f64(x: arb) -> float:
    """Minorante float64 rigoroso (nunca arredonda para cima)."""
    lo = x.lower()
    val = float(lo.str(17, radius=False))
    return math.nextafter(val, -math.inf) if arb(val) > lo else val


def arb_bounds_str(x: arb, digits: int = 22) -> dict:
    """Representacao textual rigorosa do intervalo [lower, upper]."""
    return {
        "lower": x.lower().str(digits, radius=False),
        "upper": x.upper().str(digits, radius=False),
        "ball": x.str(digits, radius=True),
    }


def max_upper(a: arb, b: arb) -> arb:
    """Majorante rigoroso de max(a, b) seguro contra sobreposicao de bolas."""
    return arb.max(a.upper(), b.upper())


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
    """K1 (integral rigorosa + cauda refinada) e K2 — independentes de b4."""
    d, kap, delta, L, pi, m0, A0, a1, a2, c = parts

    def f(x, analytic):
        v = x**6
        A = acb(1) + (acb(1) + acb(c) * v).log(analytic=analytic) / (2 * acb(L))
        B = acb(cs["c1"] * c.sqrt()) * x**3 + acb(cs["c2"])
        return (
            6 * v * (-v).exp()
            * A.pow(acb(4) / 3, analytic=analytic)
            * B.pow(acb(2) / 3, analytic=analytic)
        )

    I = acb.integral(
        f, 0, 2,
        rel_tol=arb("1e-22"),
        abs_tol=arb("1e-22"),
        eval_limit=40000,
        depth_limit=35,
    )
    if not I.is_finite() or not I.imag.contains(0) or not I.real.lower() > 0:
        raise ArithmeticError(f"quadrature failed or inconclusive on row {row.d}")

    # Cauda rigorosa para x >= 2 (v = x^6 >= 64) preservando expoentes 4/3 e 2/3
    U = arb(64)
    D = pw(1 + c / (2 * L), rat(4, 3)) * pw(1 + cs["c1"] * c.sqrt() + cs["c2"], rat(2, 3))
    tail = (D * (-U).exp() * (U**3 + 3 * U**2 + 6 * U + 6)).upper()

    lower, upper = I.real.lower(), (I.real.upper() + tail).upper()
    F = pw(2 * a2 * c.sqrt(), rat(1, 3))
    correction = (
        pw(2 * a1, rat(1, 3))
        * pw(cs["c2"] + A0 * cs["c1"] / T0, rat(2, 3))
        / (pw(arb(T0), rat(1, 3)) * pw(L, rat(4, 3)))
    )
    k1_ball = F * arb(lower).union(arb(upper)) + correction
    if not k1_ball.lower() > 0:
        raise ArithmeticError("K1 lower bound inconclusive")

    two13 = pw(arb(2), rat(1, 3))
    k2_ball = two13 * c.sqrt() * pi.sqrt() / 2 + two13 * (4 * kap / 3) * pw(A0 / arb(T0), arb(3))
    return k1_ball, k2_ball


def b1_upper(row, cs, parts, k1_ball, k2_ball, b4):
    """Algebra b4-dependente com propagacao rigorosa de majorantes."""
    d, kap, delta, L, pi, m0, A0, a1, a2, c = parts
    k1 = k1_ball.upper()
    k2 = k2_ball.upper()

    s1, s2 = rat(row.sigma_start, 32), rat(row.sigma_start + 1, 32)
    c3 = pi * m0 * b4 / (2 * delta) * (2 * delta * arb.const_euler() / L).exp()
    c4 = b4 / (5 * delta * delta.exp()) + q(".605") * (-2 * delta).exp() / L
    k3 = c3 + c4 * (pi.sqrt() / (2 * (2 * kap).sqrt()) + pi * m0 / T0)

    l0 = (kap * (1 + (1 / (4 * arb(T0)) * (1 + 2 * delta / L))**2)).exp() / (1 - arb(1) / T0)

    l1 = (
        max_upper(arb(1), pw(k1, 3 * delta / L))
        * max_upper(arb(1), pw(k3, -3 * delta / L))
        * (3 * delta + 4 * delta * L.log() / L).exp()
    )
    l2 = (
        max_upper(arb(1), pw(k2, 3 * delta / L))
        * max_upper(arb(1), pw(k3, -3 * delta / L))
        * (3 * delta).exp()
    )

    common = (
        max_upper(arb(1), pw(l0, -2 * d / L))
        * max_upper(arb(1), pw(k3, -3 * d / L))
        * max_upper(pw(l0, 2 / (2 - s1)), pw(l0, 2 / (2 - s2)))
        * max_upper(pw(k3, (2 * s1 - 1) / (2 - s1)), pw(k3, (2 * s2 - 1) / (2 - s2)))
    )

    A1 = (
        l1
        * (d * (3 + 4 * L.log() / L)).exp()
        * common
        * max_upper(arb(1), pw(k1, 3 * d / L))
        * max_upper(pw(k1, 3 * (1 - s1) / (2 - s1)), pw(k1, 3 * (1 - s2) / (2 - s2)))
    )
    A2 = (
        l2
        * (3 * d).exp()
        * common
        * max_upper(arb(1), pw(k2, 3 * d / L))
        * max_upper(pw(k2, 3 * (1 - s1) / (2 - s1)), pw(k2, 3 * (1 - s2) / (2 - s2)))
    )

    B1 = (A1 + A2 * pw(L, -6 * (1 - s2) / (2 - s2))) / (2 * pi * d)
    return B1, k3


def ceil_2dp_strict(x: arb) -> str:
    """Menor decimal com 2 casas estritamente maior que x.upper()."""
    hi_100 = (x.upper() * 100).upper()
    n = int(float(hi_100.str(20, radius=False)))
    while arb(n) / 100 <= x.upper():
        n += 1
    return f"{n / 100:.2f}"


def main():
    cert = json.loads((REPO / "NUMERIA_P78_RAMARE_PUBLIC" / "certificate.json").read_text())
    p78_b4 = {r["row"]: r["b4"] for r in cert["rows"]}
    validate_table()
    out = []

    with ctx.workprec(256):
        cs = constants()
        cs["c1"] = 6 / cs["pi"]**2 + q("1.048") / cs["logh"]
        cs["c2"] = (
            cs["pi"] * cs["m0"] * q(".62") / cs["logh"]
            + 6 * cs["m0"] / (cs["pi"] * T0)
            + cs["pi"] * cs["m0"] * q("1.048") / (arb(T0) * cs["logh"])
        )

        for idx in range(10, 17):
            row = TABLE[idx - 1]
            parts = row_parts(row, cs)
            k1_ball, k2_ball = k1_k2(row, cs, parts)

            b54 = q(".54")
            fr78 = Fr(p78_b4[idx])
            b78 = arb(fr78.numerator) / arb(fr78.denominator)

            b1_54, k3_54 = b1_upper(row, cs, parts, k1_ball, k2_ball, b54)
            b1_78, k3_78 = b1_upper(row, cs, parts, k1_ball, k2_ball, b78)

            printed = q(row.b1)
            m54 = printed - b1_54
            m78 = printed - b1_78

            # Decisao logica feita diretamente sobre bolas arb (zero erro float64)
            strict_54 = bool(m54.lower() > 0)
            strict_78 = bool(m78.lower() > 0)

            pct_ball = 100 * (b1_78 - b1_54) / b1_54

            out.append({
                "row": idx,
                "delta": row.delta,
                "b4_P78": p78_b4[idx],
                "printed_B1": row.b1,
                "suggested_B1_P78_2dp": ceil_2dp_strict(b1_78),
                "B1_upper_b4_0.54_f64_ceil": arb_upper_f64(b1_54),
                "B1_upper_b4_P78_f64_ceil": arb_upper_f64(b1_78),
                "margin_0.54_f64_floor": arb_lower_f64(m54),
                "margin_P78_f64_floor": arb_lower_f64(m78),
                "strict_0.54": strict_54,
                "strict_P78": strict_78,
                "k3_0.54_interval": arb_bounds_str(k3_54, 18),
                "k3_P78_interval": arb_bounds_str(k3_78, 18),
                "B1_0.54_interval": arb_bounds_str(b1_54, 22),
                "B1_P78_interval": arb_bounds_str(b1_78, 22),
                "delta_B1_pct_f64_ceil": arb_upper_f64(pct_ball),
            })

    print(json.dumps(out, indent=2))
    fails = [o["row"] for o in out if not o["strict_P78"]]
    print("ROWS FAILING WITH P78 b4:", fails if fails else "none")

    (REPO.parent / "p78_b1_sensitivity_results.json").write_text(
        json.dumps(
            {
                "schema": "NUMERIA.P78.S1.PROTOTYPE.V2",
                "precision_bits": 256,
                "rows": out,
                "rows_failing_with_public_b4": fails,
            },
            indent=2,
        )
        + "\n"
    )
    if fails:
        raise SystemExit(
            "PROTOTYPE GATE: printed B1 not majorized by public b4 in rows " + str(fails)
        )


if __name__ == "__main__":
    main()