"""Arb congruence for the five hybrid weights of P93-V.39 §2.

The matrices are represented in an orthonormal parity basis. Establishing
that H is the compression of the Fourier concentration complement and
that the two Schur caps are valid are separate operator-level premises.
"""

import argparse
from fractions import Fraction
import json
from math import factorial
from pathlib import Path

from flint import arb, arb_mat, ctx

from .interval_ldl import certify_inertia, dyadic


THETAS = tuple(Fraction(i, 4) for i in range(5))


def _rational(value):
    if type(value) not in (Fraction, int, str):
        raise ValueError("expected an exact rational")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError("invalid rational") from error


def _ball(value):
    x = _rational(value)
    return arb(x.numerator) / x.denominator


def _square(matrix, size):
    return isinstance(matrix, arb_mat) and matrix.nrows() == matrix.ncols() == size


def _submatrix(matrix, first_row, stop_row, first_col, stop_col):
    return arb_mat([[matrix[i, j] for j in range(first_col, stop_col)]
                    for i in range(first_row, stop_row)])


def _shift(matrix, shift):
    result = arb_mat(matrix)
    for i in range(matrix.nrows()):
        result[i, i] -= shift
    return result


def _edit_trailing_block(base, block, m):
    for i in range(block.nrows()):
        for j in range(block.ncols()):
            base[m + i, m + j] = block[i, j]


def _symmetric(matrix):
    """Intersect independently rounded entries of an analytically symmetric product."""
    for i in range(matrix.nrows()):
        for j in range(i):
            value = matrix[i, j].intersection(matrix[j, i])
            matrix[i, j] = matrix[j, i] = value
    return matrix


def _congruence_family(a0, penalty, hq, c, d, lam, margin):
    """Assemble once the theta-independent part and the affine correction."""
    if not isinstance(a0, arb_mat) or not isinstance(penalty, arb_mat):
        raise ValueError("expected Arb matrices for A0 and B")
    n, m = a0.nrows(), penalty.nrows()
    if not 0 < m <= n or not _square(a0, n) or not _square(penalty, m):
        raise ValueError("incompatible localized head and penalty")
    if not (c > 0 and d > 0 and 0 < lam < 1 and margin >= 0):
        raise ValueError("invalid C, D, lambda or margin")
    r = n - m
    if (r and not _square(hq, r)) or (not r and hq is not None):
        raise ValueError("H_Q must describe the trailing block")

    base = _shift(a0, _ball(margin))
    for i in range(m):
        for j in range(m):
            base[i, j] -= penalty[i, j]
    if not r:
        return base, base, None

    s = lam * d / ((1 - lam) * c)
    j = _shift(hq * _ball(s), arb(-1)) * _ball(1 / (1 + s))
    trailing = _shift(_submatrix(base, m, n, m, n), _ball(c))
    localized = arb_mat(base)
    _edit_trailing_block(localized, trailing, m)
    cross_top = _submatrix(base, 0, m, m, n) * j
    cross_bottom = j * _submatrix(base, m, n, 0, m)
    for i in range(m):
        for k in range(r):
            base[i, m + k] = cross_top[i, k]
            base[m + k, i] = cross_bottom[k, i]
    _edit_trailing_block(base, j * trailing * j, m)
    correction = (j * j) * _ball(c - c / lam) + j * _ball(c / (lam * (1 + s)))
    return _symmetric(localized), _symmetric(base), _symmetric(correction)


def _at_weight(base, correction, theta):
    result = arb_mat(base)
    m = result.nrows() - correction.nrows()
    weight = _ball(theta)
    for i in range(correction.nrows()):
        for j in range(correction.ncols()):
            result[m + i, m + j] += correction[i, j] * weight
    return result


def assemble_congruence(a0, penalty, hq, c, d, lam, theta, margin=0):
    """Enclose J(M_theta - margin I)J with no interval matrix inverse.

    The theta=0 endpoint may be tested directly by ``certify_hybrid``.
    For dimensions not exceeding M, the caller passes the compression of
    the localized penalty and H_Q=None; the returned heads coincide.
    """
    c, d, lam, theta, margin = map(_rational, (c, d, lam, theta, margin))
    if theta not in THETAS:
        raise ValueError("theta must be one of 0, 1/4, 1/2, 3/4, 1")
    _, base, correction = _congruence_family(a0, penalty, hq, c, d, lam, margin)
    return base if correction is None else _at_weight(base, correction, theta)


def certify_hybrid(a0, penalty, hq, c, d, lam, margin=0, *, block_size=32):
    """Check all five weights, reusing the common congruence and correction."""
    c, d, lam, margin = map(_rational, (c, d, lam, margin))
    localized, base, correction = _congruence_family(a0, penalty, hq, c, d, lam, margin)
    signs = (1,) * a0.nrows()
    direct = certify_inertia(localized, signs, block_size=block_size)
    results = []
    for theta in THETAS:
        if correction is None or theta == 0:
            ldl = direct
        else:
            ldl = certify_inertia(_at_weight(base, correction, theta), signs, block_size=block_size)
        results.append({"theta_exact": str(theta), "ldl": ldl})
    return {
        "schema": "NUMERIA.P93V40.HYBRID.1",
        "basis": "orthonormal parity",
        "dimension": a0.nrows(),
        "localized_dimension": penalty.nrows(),
        "margin_exact": str(margin),
        "distinct_assembled_heads": 1 if correction is None else len(THETAS),
        "results": results,
    }


def conditional_schur(c, d, lam, *, rho="33/25", modes=3309):
    """Evaluate the V.39 §5 budget, conditional on the required finite heads."""
    c, d, lam, rho = map(_rational, (c, d, lam, rho))
    if not (c > 0 and d > 0 and 0 < lam < 1 and rho > 1):
        raise ValueError("invalid Schur parameters")
    if type(modes) is not int or not 600 <= modes <= 5000:
        raise ValueError("expected 600 to 5000 modes")
    a, t, h = Fraction(1387, 1000), 4500, Fraction(1, 10**60)
    m, epsilon = _ball(9 * h / 10), _ball(Fraction(1, 2**2000))
    ph = _ball(h / 10)
    s = lam * d / ((1 - lam) * c)
    f1 = (c / lam) * s / (1 + s)
    prefactor = (_ball(8 * a * t / 3) * _ball(a * t * (rho - 1 / rho)).exp()
                 / _ball((1 - 1 / rho)**2 * (1 - rho**-4)))
    rb = _ball(rho)

    rows = []
    for theta in THETAS:
        penalty = (1 - theta) * c + theta * f1
        coefficient = _ball(50 + theta * d / (1 - lam))

        def at(n):
            # The even sector gives the larger g; the odd sector gains rho^-2.
            g = prefactor / rb**(4 * n)
            reserve = _ball(2 - penalty) - 50 * g - epsilon - m
            slack = ph * reserve - 2 * coefficient**2 * g - 2 * epsilon**2
            return reserve, slack

        lo, hi = 599, 5000
        if not all(item > 0 for item in at(hi)):
            raise ArithmeticError("no certified tail budget in the search range")
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if all(item > 0 for item in at(mid)):
                hi = mid
            else:
                lo = mid
        if lo >= 600 and not any(item.upper() <= 0 for item in at(lo)):
            raise ArithmeticError("the previous dimension is numerically undecided")
        fac = factorial(2 * hi)
        if not (48 * 2**2000 < fac and 72 * 2**2000 < fac**2):
            raise ArithmeticError("pole inequalities do not hold")
        reserve, slack = at(modes)
        rows.append({
            "theta_exact": str(theta),
            "asymptotic_reserve_exact": str(2 - penalty),
            "tested_modes": modes,
            "reserve_lower_dyadic": dyadic(reserve.lower()),
            "slack_lower_dyadic": dyadic(slack.lower()),
            "tested_budget_passes": bool(reserve > 0 and slack > 0),
            "first_sufficient_modes_for_this_budget": hi,
        })
    return {
        "schema": "NUMERIA.P93V40.HYBRID.SCHUR.1",
        "rho_exact": str(rho),
        "conditional_on_full_shifted_heads": True,
        "results": rows,
    }


def run(modes=12, bits=384):
    """Rebuild both physical parity heads before checking the hybrid weights."""
    from .assembly import ARCHIVE, WIDTH, assemble

    constants = json.loads((ARCHIVE / "hybrid39.json").read_text(encoding="utf-8"))
    c, d, lam = (constants[key] for key in ("c", "d", "lambda"))
    with ctx.workprec(bits):
        sectors = []
        for parity, sector in enumerate(assemble(modes, bits)):
            m = sector["b"].nrows()
            hq = (_shift(-_submatrix(sector["g"], m, modes, m, modes), arb(-1))
                  if modes > m else None)
            sectors.append({"parity": parity, **certify_hybrid(
                sector["a0"], sector["b"], hq, c, d, lam, "1e-60")})
        return {
            "schema": "NUMERIA.P93V40.HYBRID.RUN.1",
            "precision_bits": bits,
            "width_exact": str(WIDTH),
            "sectors": sectors,
            "conditional_schur": conditional_schur(c, d, lam),
            "scope": "finite congruence and inherited Schur budget; no full crossing claim",
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modes", type=int, default=12)
    parser.add_argument("--bits", type=int, default=384)
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--out", type=Path)
    output.add_argument("--verify", type=Path)
    args = parser.parse_args()
    content = json.dumps(run(args.modes, args.bits), indent=2) + "\n"
    if args.out is not None:
        args.out.write_text(content, encoding="utf-8")
    elif args.verify.read_text(encoding="utf-8") != content:
        raise ValueError("the harmonic replay differs from the certificate")
    else:
        print(f"Verified {args.verify}")


if __name__ == "__main__":
    main()
