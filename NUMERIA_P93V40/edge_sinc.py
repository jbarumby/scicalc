"""Directed Legendre coefficients and norms for bandlimited edge evaluations."""

from fractions import Fraction

from flint import arb, ctx

from .assembly import WIDTH


EDGE_CUTOFF = 1_600_000


def _positive_ball(value):
    exact = Fraction(value)
    if exact <= 0:
        raise ValueError("edge width and frequency cutoff must be positive")
    return arb(exact.numerator) / exact.denominator


def _endpoint_polynomial_integral(n, h, sine, cosine, tolerance):
    """Integrate sin(h*u)*(P_n(1-2*u)-1)/u for 0 <= u <= 1."""
    total = arb(0)
    inverse_h = 1 / h
    n_times_n_plus_one = n * (n + 1)
    scale = inverse_h
    binomial_product = 1
    alternating_sum = 1
    absolute_sum = 1
    parity = (-1)**n

    for j in range(n):
        index = j + 1
        next_product = binomial_product * (n - j) * (n + index) // index**2
        if j + 1 < n:
            # Every later ratio of successive positive term bounds is at
            # most (n+n*(n+1)/(j+2))/h. Use a geometric remainder only
            # when the whole unresolved tail is small at working precision.
            bound = scale * (1 + absolute_sum + next_product)
            if bound < tolerance:
                ratio = (n + arb(n_times_n_plus_one) / (j + 2)) * inverse_h
                if ratio < arb(1) / 2:
                    return total + arb(0, (2 * bound).upper())
        at_one = parity * alternating_sum - 1
        if j % 2:
            term = at_one * sine
        else:
            term = at_one * cosine + next_product
        total += (-1)**(j // 2 + 1) * scale * term
        alternating_sum += (-1)**index * next_product
        absolute_sum += next_product
        binomial_product = next_product
        scale *= index * inverse_h
    return total


def endpoint_coefficients(modes, *, width=WIDTH, omega=EDGE_CUTOFF):
    """Enclose the first `modes` coefficients of the right endpoint in each parity.

    The basis is sqrt((2n+1)/(2a))*P_n(x/a) on [-a,a]. The evaluation
    kernel is sin(omega*(a-x))/(pi*(a-x)), and the left endpoint
    coefficients are (-1)**n times those returned here.
    """
    if type(modes) is not int or modes <= 0:
        raise ValueError("modes must be a positive integer")
    bits = ctx.prec
    degree = 2 * modes - 1
    h = 2 * _positive_ball(width) * _positive_ball(omega)
    guard = 64
    if not (h > 2 * degree and arb(degree * (degree + 1)) < 16 * h):
        # The general exact polynomial formula can cancel when h is
        # small. This also covers tests at modest cutoffs.
        inverse_mantissa, inverse_exponent = (1 / h).upper().man_exp()
        inverse_log = max(0, int(inverse_exponent) + int(inverse_mantissa).bit_length())
        guard += degree * (3 + degree.bit_length() + inverse_log)

    with ctx.workprec(bits + guard):
        a = _positive_ball(width)
        h = 2 * a * _positive_ball(omega)
        sine, cosine = h.sin_cos()
        base = h.si()
        tolerance = arb(2)**(-bits - 32)
        pi = arb.pi()
        values = [
            ((2 * n + 1) / (2 * a)).sqrt()
            * (base + _endpoint_polynomial_integral(n, h, sine, cosine, tolerance)) / pi
            for n in range(degree + 1)
        ]
    return [value + arb(0) for value in values[::2]], [value + arb(0) for value in values[1::2]]


def endpoint_norm_squares(*, width=WIDTH, omega=EDGE_CUTOFF):
    """Return the squared norms of the even and odd projected endpoint kernels."""
    cutoff = _positive_ball(omega)
    h = 2 * _positive_ball(width) * cutoff
    sine, cosine = h.sin_cos()
    sine_integral = (2 * h).si()
    cosine_integral = (2 * h).ci() - (2 * h).log() - arb.const_euler()
    factor = cutoff / arb.pi()**2
    norm = factor * (sine_integral - sine**2 / h)
    cross = factor * (cosine * cosine_integral + sine * sine_integral) / h
    return (norm + cross) / 2, (norm - cross) / 2


def endpoint_tail_norm_squares(coefficients, *, width=WIDTH, omega=EDGE_CUTOFF):
    """Enclose the squared complementary norms after the supplied prefixes."""
    return tuple(
        full_norm - sum((value**2 for value in sector), arb(0))
        for full_norm, sector in zip(endpoint_norm_squares(width=width, omega=omega),
                                     coefficients, strict=True)
    )
