"""Exact rational modal tail budgets for the P239 window and fixed beta."""

from fractions import Fraction
from math import isqrt, prod
from operator import index


HALF_WIDTH = Fraction(161, 200)
MAX_FIXED_SYMBOL_HALF_WIDTH = Fraction(19, 20)
MAX_FOUR_PRIME_MODEL_HALF_WIDTH = Fraction(973, 1000)
SYMBOL_BOUND = 20
ARCHIVED_CUTOFF = 504
RATIO_BOUND = Fraction(3, 40)
WIDE_WINDOW_RATIO_BOUND = Fraction(1, 8)


def fixed_symbol_half_width(half_width):
    """Validate a width for which the P237 prime powers and beta are fixed.

    On this interval 5 < exp(2L) < 7. The degree-8 exponential sum is
    greater than 5 at 161/100. At 19/10, adding a geometric upper bound for
    the remaining series to that sum gives less than 7. The next prime
    power after 5 is 7.
    """
    half_width = Fraction(half_width)
    if not HALF_WIDTH <= half_width <= MAX_FIXED_SYMBOL_HALF_WIDTH:
        raise ValueError('half width must be in [161/200, 19/20] for the fixed symbol')
    return half_width


def four_prime_model_half_width(half_width):
    """Validate a width for the auxiliary model with primes frozen at P237."""
    half_width = Fraction(half_width)
    if not HALF_WIDTH <= half_width <= MAX_FOUR_PRIME_MODEL_HALF_WIDTH:
        raise ValueError("four-prime model half width must be in [161/200, 973/1000]")
    return half_width


def _squared_ratio(degree, argument):
    # The Legendre normalization grows with n, so this ratio includes
    # sqrt((2*n+5)/(2*n+1)) from the prefactor.
    return argument**4 / ((2 * degree + 1) * (2 * degree + 3)**2 * (2 * degree + 5))


def modal_tail_bounds(cutoff, modes_per_parity, parity, *, half_width=HALF_WIDTH):
    """Bound ||B||² and ||D - beta I|| using P239/P243 analytic premises.

    Uses the spherical Bessel envelopes and the fixed-beta symbol bound 20
    on [0, 504] for the auxiliary four-prime model. The result concerns the
    complement and its coupling without any assumption about the head spectrum.
    """
    cutoff = Fraction(cutoff)
    modes_per_parity, parity = index(modes_per_parity), index(parity)
    half_width = four_prime_model_half_width(half_width)
    if not 0 < cutoff <= ARCHIVED_CUTOFF or modes_per_parity < 1 or parity not in (0, 1):
        raise ValueError('use 0 < cutoff <= 504, positive modes and parity 0 or 1')

    degree = 2 * modes_per_parity + parity
    x = half_width * cutoff
    z = half_width / 2
    squared_ratio = max(_squared_ratio(degree, x), _squared_ratio(degree, z))
    ratio_bound = RATIO_BOUND
    if half_width > HALF_WIDTH and squared_ratio >= ratio_bound**2:
        ratio_bound = WIDE_WINDOW_RATIO_BOUND
    if squared_ratio >= ratio_bound**2:
        raise ArithmeticError('first discarded degree does not satisfy the geometric bound')

    # Only the seed uses an integer cap. For later degrees the geometric
    # ratio is applied to the envelope with its degree-dependent prefactor.
    prefactor_squared = 2 * half_width * (2 * degree + 1)
    prefactor_upper = isqrt(prefactor_squared.numerator // prefactor_squared.denominator) + 1
    factorial = prod(range(1, 2 * degree + 2, 2))
    tail_transform = prefactor_upper * x**degree / factorial

    u = z**2 / 2
    cosh_upper = Fraction(1, 1) / (1 - u)
    # cosh(z) <= exp(z²/2) <= 1/(1-z²/2). Parseval then bounds
    # the head pole l2 norm by sqrt(2L) cosh(z), which is below 2.
    if not 0 < u < Fraction(1, 2) or 2 * half_width * cosh_upper**2 >= 4:
        raise ArithmeticError('pole bounds do not apply')
    tail_pole = prefactor_upper * cosh_upper * z**degree / factorial
    head_pole_l2_upper = 2
    head_pole_l1_upper = head_pole_l2_upper * (isqrt(modes_per_parity) + 1)
    head_transform_sum_upper = modes_per_parity * prefactor_upper

    integration = cutoff * SYMBOL_BOUND / 3  # pi > 3
    geometric = 1 / (1 - ratio_bound)
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
