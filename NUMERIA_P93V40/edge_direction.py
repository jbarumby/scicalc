"""Bandwise Loewner upper bound for the squared low-spectrum cross block."""

from fractions import Fraction

from flint import arb

from .assembly import T, WIDTH, band_gram
from .transfer import RHO, WIDTH_CAP, ball, count_mass, low_symbol_cap, pole_audit


BAND_EDGES = (0, 500, 1000, 2000, 3000, 4000, 4250, T)
SYMBOL_CAP = 18


def _validate_split(split, parity):
    if type(split) is not int or split < 1 or type(parity) is not int or parity not in (0, 1):
        raise ValueError("the split must be positive and the parity must be 0 or 1")


def band_tail_concentration_upper(split, parity, left, right):
    """Bound the band concentration on the complementary Legendre modes.

    The Chebyshev ellipse estimate uses the squared coefficient bound with
    I_0(atd/2)^2 <= exp(atd)*min(1, pi/(4atd)), where d=rho-1/rho.
    """
    _validate_split(split, parity)
    left, right = Fraction(left), Fraction(right)
    if not 0 <= left < right <= T:
        raise ValueError("require a nonempty band in [0, T]")
    a, rho = ball(WIDTH), ball(RHO)
    delta = rho - 1 / rho
    rate = a * delta
    tail_degree = 2 * split + parity
    coefficient = 2 * rho**(-2 * tail_degree) / (delta * (1 - rho**-4))
    left_product = rate * ball(left)
    exponential_difference = left_product.exp() * (rate * ball(right - left)).expm1()
    if left:
        bessel_factor = min(arb(1), (arb.pi() / (4 * left_product)).upper())
    else:
        bessel_factor = arb(1)
    return min(arb(1), (coefficient * exponential_difference * bessel_factor).upper())


def _exact_fraction(value):
    mantissa, exponent = value.man_exp()
    return Fraction(int(mantissa)) * Fraction(2)**int(exponent)


def directional_cross_majorant(modes, split, parity):
    """Enclose a Loewner upper bound for C C^T on a prefix of the head.

    For theta=0 with localization contained in the head, the cross block
    C consists of the signed low spectrum and the pole term. In band j,
    C_j C_j^T <= 18^2 g_j G_j. Cauchy–Schwarz for the sum of all band
    operators and the pole term costs a factor equal to the number of terms.
    """
    _validate_split(split, parity)
    if type(modes) is not int or not 0 < modes <= split:
        raise ValueError("require a nonempty prefix within the head")
    low_symbol_cap(count_mass())
    width, width_cap = ball(WIDTH), ball(WIDTH_CAP)
    if not (width < width_cap and 2 * width_cap * width_cap.exp() < 15
            and 2 * width_cap.sinh() < 9):
        raise ArithmeticError("the inherited pole bounds do not apply")

    weights = {}
    for left, right in zip(BAND_EDGES, BAND_EDGES[1:]):
        g = _exact_fraction(band_tail_concentration_upper(split, parity, left, right))
        weights[right] = weights.get(right, 0) + g
        if left:
            weights[left] = weights.get(left, 0) - g
    endpoints = [(endpoint, weight) for endpoint, weight in sorted(weights.items()) if weight]
    term_count = len(BAND_EDGES)
    majorant = band_gram(endpoints, modes)[parity] * (term_count * SYMBOL_CAP**2)
    pole_cross, _ = pole_audit(split)
    diagonal = ball(term_count * pole_cross**2)
    for i in range(modes):
        majorant[i, i] += diagonal
    if not all(entry.is_finite() for entry in majorant.entries()):
        raise ArithmeticError("the directional majorant enclosure was not finite")
    return majorant
