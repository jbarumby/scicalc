"""Assemble and certify finite fixed-beta Zhu heads from enclosed GL60 sums."""

from fractions import Fraction
from math import isqrt
from operator import index

from flint import arb, arb_mat, ctx

from NUMERIA_ES8P249.interval_audit import _dyadic
from NUMERIA_ES8P249.legendre_tail import ARCHIVED_CUTOFF
from NUMERIA_ES8P249.prime_operator import P236_HALF_WIDTH, P236_PRIME_POWERS


_PARITY = {"even": 0, "odd": 1}


def _ceil_sqrt(value):
    floor = isqrt(value.numerator // value.denominator)
    return floor + (floor * floor < value)


def gl60_error_budget(sector, modes=259, cutoff=309):
    """Return exact degreewise GL60 bounds as rational numbers.

    Uses the P237 classical Bernstein ellipse and Gauss–Legendre estimates
    with the same symbol bound <50 on the retained half-unit panels. The
    ellipse gives |Im t| <= 1/3 and |Re t| <= cutoff + 1/6. Rodrigues' bound
    and the plane-wave identity give |j_n(Lt)| <= min(2, 2 Y^n/(2n+1)!!)
    for Y >= |Lt|, since exp(|Im Lt|) < 2.
    """
    modes, cutoff = index(modes), index(cutoff)
    if sector not in _PARITY or modes < 1 or not 1 <= cutoff <= ARCHIVED_CUTOFF:
        raise ValueError("expected a parity, positive mode count and integer cutoff in [1, 504]")
    parity = _PARITY[sector]
    ellipse_radius_squared = P236_HALF_WIDTH**2 * ((cutoff + Fraction(1, 6))**2 + Fraction(1, 9))
    argument_upper = _ceil_sqrt(ellipse_radius_squared)
    envelope = []
    double_factorial = 1
    for degree in range(2 * modes - 1 + parity):
        double_factorial *= 2 * degree + 1
        if degree % 2 == parity:
            normalization = _ceil_sqrt(2 * P236_HALF_WIDTH * (2 * degree + 1))
            bessel = min(Fraction(2), Fraction(2 * argument_upper**degree, double_factorial))
            envelope.append(normalization * bessel)

    # On the rescaled ellipse the integrand is bounded by 50 a_i a_j / 12.
    # The degree-119 Chebyshev tail and GL60 error multiply this by 12 / 3^120.
    factor = Fraction(2 * cutoff * 50, 3**120)
    max_envelope = max(envelope)
    return {
        "panels": 2 * cutoff,
        "argument_upper": argument_upper,
        "panel_factor": Fraction(50, 3**120),
        "all_panels_factor": factor,
        "max_entry_bound": factor * max_envelope**2,
        "row_sum_bound": factor * max_envelope * sum(envelope),
    }


def _dimension(matrix):
    if not isinstance(matrix, arb_mat) or matrix.nrows() != matrix.ncols() or matrix.nrows() < 1:
        raise ValueError("expected a nonempty square Arb matrix")
    return matrix.nrows()


def _endpoints(value):
    return [_dyadic(value.lower()), _dyadic(value.upper())]


def assemble_finite_head(quadrature, sector, *, cutoff=309):
    """Add beta and the P237 pole term to the enclosed GL60 quadrature.

    The GL60 truncation budget is returned separately as an operator-norm
    perturbation. It is therefore added only once by the Weyl step.
    """
    dimension = _dimension(quadrature)
    budget = gl60_error_budget(sector, dimension, cutoff)
    half_width = arb(P236_HALF_WIDTH.numerator) / P236_HALF_WIDTH.denominator
    pi = arb.pi()
    mass = sum((2 * arb(p).log() / arb(q).sqrt() for q, p in P236_PRIME_POWERS), arb(0))
    beta = (arb(ARCHIVED_CUTOFF) / (2 * pi)).log() - arb(1) / ARCHIVED_CUTOFF - mass
    if not beta > 0:
        raise ArithmeticError("fixed P237 beta is not certified positive")
    parity = _PARITY[sector]
    argument = half_width / 2
    pole = [
        (2 * half_width * (2 * degree + 1)).sqrt()
        * (pi / (2 * argument)).sqrt()
        * argument.bessel_i(arb(2 * degree + 1) / 2)
        for degree in range(parity, 2 * dimension - 1 + parity, 2)
    ]
    signed_factor = -2 if parity else 2
    assembled = arb_mat(dimension, dimension)
    for row in range(dimension):
        for column in range(row, dimension):
            left, right = quadrature[row, column], quadrature[column, row]
            if not left.is_finite() or not right.is_finite() or not left.overlaps(right):
                raise ValueError(f"incompatible quadrature enclosures at ({row}, {column})")
            # Both original entries enclose the same exact symmetric integral sum.
            entry = left.intersection(right) if row != column else left
            entry += signed_factor * pole[row] * pole[column]
            if row == column:
                entry += beta
            assembled[row, column] = assembled[column, row] = entry
    return assembled, {
        "beta_endpoints_dyadic": _endpoints(beta),
        "pole_sign_and_factor": signed_factor,
        "gl60_error_budget": budget,
    }


def _midpoint_ldl_pivots(midpoint):
    """Enclose consecutive LDL pivots by symmetric Schur updates."""
    work = midpoint.tolist()
    dimension = midpoint.nrows()
    pivots = []
    for k in range(dimension):
        pivot = work[k][k]
        if not pivot.is_finite() or not pivot > 0:
            return pivots, pivot
        pivots.append(pivot)
        multipliers = [work[j][k] / pivot for j in range(k + 1, dimension)]
        for i in range(k + 1, dimension):
            factor = work[i][k]
            for j, multiplier in enumerate(multipliers[:i-k], start=k + 1):
                work[i][j] -= factor * multiplier
    return pivots, None


def _max_radius_row_sum(matrix):
    dimension = matrix.nrows()
    largest = arb(0)
    for i in range(dimension):
        radius = arb(0)
        for j in range(dimension):
            entry, mirror = matrix[i, j], matrix[j, i]
            if not entry.is_finite() or entry.mid() != mirror.mid() or entry.rad() != mirror.rad():
                raise ValueError("expected a finite and symmetric box of matrices")
            radius += entry.rad().upper()
        largest = max(largest, radius.upper())
    return largest


def certify_center_weyl(matrix, analytic_row_error=Fraction(0)):
    """Certify positivity of a symmetric interval box plus a known norm budget."""
    dimension = _dimension(matrix)
    analytic_row_error = Fraction(analytic_row_error)
    if analytic_row_error < 0:
        raise ValueError("expected a nonnegative analytic error bound")
    roundoff = _max_radius_row_sum(matrix)
    analytic_upper = (arb(analytic_row_error.numerator) / analytic_row_error.denominator).upper()
    perturbation = (roundoff + analytic_upper).upper()
    midpoint = matrix.mid()
    pivots, failed_pivot = _midpoint_ldl_pivots(midpoint)
    result = {
        "dimension": dimension,
        "status": "MIDPOINT_LDL_INCONCLUSIVE",
        "certified_positive": False,
        "positive_pivot_count": len(pivots),
        "pivot_lower_dyadic": [_dyadic(pivot.lower()) for pivot in pivots],
        "quadrature_and_assembly_radius_row_upper_dyadic": _dyadic(roundoff),
        "gl60_analytic_row_upper_dyadic": _dyadic(analytic_upper),
        "perturbation_upper_dyadic": _dyadic(perturbation),
    }
    if failed_pivot is not None:
        result["first_inconclusive_pivot_index"] = len(pivots)
        result["first_inconclusive_pivot_endpoints_dyadic"] = (
            _endpoints(failed_pivot) if failed_pivot.is_finite() else None
        )
        return result

    result["status"] = "MIDPOINT_INVERSE_INCONCLUSIVE"
    try:
        inverse = midpoint.inv()
    except ZeroDivisionError:
        return result
    if not all(entry.is_finite() for entry in inverse.entries()):
        return result
    inverse_squared_upper = sum((entry.abs_upper()**2 for entry in inverse.entries()), arb(0)).upper()
    if not inverse_squared_upper.is_finite() or not inverse_squared_upper > 0:
        return result
    lambda_lower = (1 / inverse_squared_upper.sqrt()).lower()
    margin_lower = (lambda_lower - perturbation).lower()
    certified = margin_lower > 0
    result.update({
        "status": "CERTIFIED_POSITIVE" if certified else "WEYL_INCONCLUSIVE",
        "certified_positive": bool(certified),
        "inverse_frobenius_squared_upper_dyadic": _dyadic(inverse_squared_upper),
        "midpoint_lambda_lower_dyadic": _dyadic(lambda_lower),
        "margin_lower_dyadic": _dyadic(margin_lower),
    })
    return result


def certify_finite_head(quadrature, sector, *, cutoff=309, precision_bits=1024):
    """Assemble and audit a finite head for a caller-supplied rigorous GL60 sum."""
    cutoff, precision_bits = index(cutoff), index(precision_bits)
    if precision_bits < 128:
        raise ValueError("expected at least 128 Arb bits")
    with ctx.workprec(precision_bits):
        matrix, assembly = assemble_finite_head(quadrature, sector, cutoff=cutoff)
        budget = assembly["gl60_error_budget"]
        certificate = certify_center_weyl(matrix, budget["row_sum_bound"])
        assembly["gl60_error_budget"] = {
            key: str(value) if isinstance(value, Fraction) else value for key, value in budget.items()
        }
        return {
            "schema": "NUMERIA.ES8P249.FINITE_HEAD_CENTER_WEYL.1",
            "sector": sector,
            "cutoff": cutoff,
            "fixed_beta_reference_cutoff": ARCHIVED_CUTOFF,
            "half_width": str(P236_HALF_WIDTH),
            "precision_bits": precision_bits,
            "assumes_classical_P237_GL60_E3_budget": True,
            "assumes_input_contains_retained_GL60_quadrature": True,
            "assembly": assembly,
            **certificate,
        }
