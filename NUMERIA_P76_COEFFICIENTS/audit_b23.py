"""Rigorous check of *two* coefficients of CS v2, Table 1: B2 and B3.

Implemented from arXiv:2507.15184v2 equations (10), (33), (34),
Theorem 1 definitions of A3,A4,A5, and Lemma 5's exact m0.

This check does NOT establish B1, which additionally requires the K1
integral and A1,A2/L/K calculations. Do not label the full Table 1
independently certified when only B2 and B3 have passed.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path

from flint import arb, ctx
from .table import TABLE, validate_table

H = 3 * 10**12
SOURCE = "https://arxiv.org/html/2507.15184v2"


def q(s):
    """Enclose a decimal parsed as an exact rational, not a binary float."""
    x = Fraction(s)
    return arb(x.numerator) / arb(x.denominator)


def dyadic(x):
    """Return exact dyadic mantissa/exponent for one already-exact Arb endpoint."""
    if not x.is_exact():
        raise ArithmeticError("cannot serialize an inexact ball as an exact endpoint")
    m, e = x.man_exp()
    return [str(m), int(e)]


def interval(ball):
    return {"lo_dyad": dyadic(ball.lower()), "hi_dyad": dyadic(ball.upper())}


def constants():
    """Use the original displayed formulas; preserve every Arb rounding radius."""
    u = arb(H)
    logh = u.log()
    pi = arb.pi()
    log2 = arb.const_log2()
    gamma = arb.const_euler()
    m0 = (1 + (arb(2) / 3) * (arb(6) / 5).sqrt()).sqrt()
    if not 0 < m0 < q("1.31541"):
        raise ArithmeticError("Lemma 5 m0 check failed")

    f = 1 + (2 + gamma + arb(7) / (36 * u)) / logh
    square = 9 * f**2 * logh**2 / u
    if not 0 < square < 1:
        raise ArithmeticError("C5/C6 logarithm branch invalid")

    # Eq. (33), h=1, u=H=T0, Lemma 7.
    c5 = -(1 - square).log() * (3 + 8 * pi * m0) * (1 + q("3.31") / logh) / (9 * f**2) * logh

    # Eq. (34), u=H0=H, Lemma 8.
    c6 = (220
          + (pi / log2) * (1 + 1/(u * u.sqrt())).log()
          + (4 / log2) * (
              ((arb(5)/(2*u))**2 + (1 + 2/u)**2).sqrt() / (2*pi)
          ).log()
          - (2*pi/log2) * (1 - square).log())
    if not c5 > 0 or not c6 > 0:
        raise ArithmeticError("C5/C6 positivity lost")
    return {"c5": c5, "c6": c6, "logh": logh, "pi": pi, "log2": log2, "m0": m0}


def calc_b23(row, cs=None):
    cs = constants() if cs is None else cs
    d = q(row.d)
    sigma1 = arb(row.sigma_start)/32
    z = (3 - 2*sigma1)/(2*d) + 1/cs["logh"]
    a3 = z/(cs["pi"]*cs["log2"])
    a4 = 3*z/(2*cs["log2"])
    a5 = (cs["c6"]/(2*cs["pi"]))*z + cs["c5"]/(2*cs["pi"]*d)
    # Equation (10), with H=T0=H_RH and h=1.
    return a3 + a4, a3*cs["logh"] + a5


def audit(precision=192, table=TABLE):
    if not isinstance(precision, int) or precision < 128:
        raise ValueError("precision must be >= 128 bits")
    validate_table()
    rows=[]
    with ctx.workprec(precision):
        cs=constants()
        lower_margins=[]
        for index, row in enumerate(table, 1):
            b2, b3=calc_b23(row,cs)
            # Publications display exact decimal bounds. Compare to exact rational
            # decimal; rounding *down* would fail closed here.
            d2=(q(row.b2)-b2).lower()
            d3=(q(row.b3)-b3).lower()
            passed=bool(d2>0 and d3>0)
            lower_margins.extend((d2,d3))
            rows.append({
                "row":index,
                "sigma_interval": [f"{row.sigma_start}/32", f"{row.sigma_start+1}/32"],
                "d":row.d,
                "b2_published":row.b2,
                "b3_published":row.b3,
                "b2_recomputed":interval(b2),
                "b3_recomputed":interval(b3),
                "b2_margin_lower_dyad":dyadic(d2),
                "b3_margin_lower_dyad":dyadic(d3),
                "both_strictly_upper":passed,
            })
        if len(rows)!=16 or not all(item["both_strictly_upper"] for item in rows):
            raise ArithmeticError("B2/B3 theorem-to-table check FAILED; do not infer B1 validity")
        minimum=min(lower_margins)
        return {
            "schema":"NUMERIA.P76.CS_V2.B23.1",
            "source":SOURCE,
            "precision_bits":precision,
            "hrh":H,
            "c5_interval":interval(cs["c5"]),
            "c6_interval":interval(cs["c6"]),
            "rows":rows,
            "b2_b3_strict_upper_verified":True,
            "minimum_rounded_coefficient_margin_lower_dyad":dyadic(minimum),
            "b1_recomputed":False,
            "full_corollary_independently_audited":False,
            "note":"The B1 proof requires independently enclosing K1, K2 and related A1,A2,L1,L2; not done here.",
        }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision",type=int,default=192)
    io=parser.add_mutually_exclusive_group(required=True)
    io.add_argument("--out",type=Path)
    io.add_argument("--verify",type=Path)
    args=parser.parse_args()
    result=audit(args.precision)
    if args.verify is not None:
        actual=json.loads(args.verify.read_text(encoding="utf-8"))
        if result!=actual:
            raise SystemExit("P76 B2/B3 certificate mismatch")
        print("P76 B2/B3 verified, precision bits:",args.precision)
    else:
        args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
        print("P76 B2/B3 16 rows verified at",args.precision,"bits")

if __name__=="__main__":
    main()
