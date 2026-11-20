"""Exact rational modal tail budgets for the P239 window and fixed beta."""

from fractions import Fraction
from math import isqrt, prod
from operator import index


HALF_WIDTH = Fraction(161, 200)
SYMBOL_BOUND = 20
ARCHIVED_CUTOFF = 504
RATIO_BOUND = Fraction(3, 40)


def _squared_ratio(degree, argument):
    # The Legendre normalization grows with n, so this ratio includes
    # sqrt((2*n+5)/(2*n+1)) from the prefactor.
    return argument**4 / ((2 * degree + 1) * (2 * degree + 3)**2 * (2 * degree + 5))


def modal_tail_bounds(cutoff, modes_per_parity, parity):
    """Bound ||B||² and ||D - beta I|| using P239/P243 analytic premises.

    Uses the spherical Bessel envelopes and the fixed-beta symbol bound 20
    on [0, 504]. The result concerns the complement and its coupling, without
    any assumption about the spectrum of the retained finite block.
    """
    cutoff = Fraction(cutoff)
    modes_per_parity, parity = index(modes_per_parity), index(parity)
    if not 0 < cutoff <= ARCHIVED_CUTOFF or modes_per_parity < 1 or parity not in (0, 1):
        raise ValueError('use 0 < cutoff <= 504, positive modes and parity 0 or 1')

    degree = 2 * modes_per_parity + parity
    x = HALF_WIDTH * cutoff
    z = HALF_WIDTH / 2
    if max(_squared_ratio(degree, x), _squared_ratio(degree, z)) >= RATIO_BOUND**2:
        raise ArithmeticError('first discarded degree does not satisfy the geometric bound')

    # Only the seed uses an integer cap. For later degrees the geometric
    # ratio is applied to the envelope with its degree-dependent prefactor.
    prefactor_squared = 2 * HALF_WIDTH * (2 * degree + 1)
    prefactor_upper = isqrt(prefactor_squared.numerator // prefactor_squared.denominator) + 1
    factorial = prod(range(1, 2 * degree + 2, 2))
    tail_transform = prefactor_upper * x**degree / factorial

    u = z**2 / 2
    cosh_upper = Fraction(1, 1) / (1 - u)
    # cosh(z) <= exp(z²/2) <= 1/(1-z²/2). Parseval then bounds
    # the head pole l2 norm by sqrt(2L) cosh(z), which is below 2.
    if not 0 < u < Fraction(1, 2) or 2 * HALF_WIDTH * cosh_upper**2 >= 4:
        raise ArithmeticError('pole bounds do not apply')
    tail_pole = prefactor_upper * cosh_upper * z**degree / factorial
    head_pole_l2_upper = 2
    head_pole_l1_upper = head_pole_l2_upper * (isqrt(modes_per_parity) + 1)
    head_transform_sum_upper = modes_per_parity * prefactor_upper

    integration = cutoff * SYMBOL_BOUND / 3  # pi > 3
    geometric = 1 / (1 - RATIO_BOUND)
    transform_sum = geometric * tail_transform
    pole_sum = geometric * tail_pole
    row = integration * tail_transform * head_transform_sum_upper + 2 * tail_pole * head_pole_l1_upper
    column = integration * transform_sum * prefactor_upper + 2 * pole_sum * head_pole_l2_upper
    deviation = integration * tail_transform * transform_sum + 2 * tail_pole * pole_sum
    return row * column, deviation


def reduced_tail_report(cutoff):
    """Replay the archived scale of tail bounds for 259 modes per parity."""
    cutoff = index(cutoff)
    modes = 259
    coupling_squared_limit = Fraction(1, 10**180)
    deviation_limit = Fraction(1, 10**190)
    for parity in (0, 1):
        coupling_squared, deviation = modal_tail_bounds(cutoff, modes, parity)
        if not (coupling_squared < coupling_squared_limit and deviation < deviation_limit):
            raise ArithmeticError('reduced-cutoff modal tail budget did not close')
    return {
        'cutoff': cutoff,
        'symbol_abs_strict_upper_assumed': SYMBOL_BOUND,
        'modes_per_parity': modes,
        'first_discarded_degrees': {'even': 2 * modes, 'odd': 2 * modes + 1},
        'coupling_norm_strict_upper': '1e-90',
        'tail_deviation_norm_strict_upper': '1e-190',
    }
