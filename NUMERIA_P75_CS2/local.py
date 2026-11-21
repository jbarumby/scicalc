"""Provable per-cell curvature upper bounds, independent of sampled searches."""
from flint import arb


def local_curvature_upper(gammas, midpoint, width):
    """S''(t) = 2*sum cos(gamma*t), |S'''(t)| <= 2*sum gamma.

    max S'' on [m-h/2,m+h/2] <= S''(m).upper + h*sum(gamma.upper).
    Return an upper endpoint and cap it by the global 2N upper bound.
    """
    if not width > 0:
        raise ValueError('width must be positive')
    gammas=tuple(gammas)
    if not gammas:
        raise ValueError('at least one zero required')
    for gamma in gammas:
        if not gamma > 0:
            raise ArithmeticError('zero ordinate must be positive')
    second=(2*sum(((g*midpoint).cos() for g in gammas),arb(0))).upper()
    third=(sum((g.upper() for g in gammas),arb(0)) * width.upper()).upper()
    return min(arb(2*len(gammas)), (second+third).upper()).upper()


def lower_secant_local(left_value, right_value, width, curvature_upper):
    """Whole-cell S(t) lower from endpoint lower bounds and upper on S''.

    S >= secant - M*h^2*u(1-u)/2. If the endpoint separation is d,
    optimize the parabola exactly with upper interval rounding.
    """
    left=left_value.lower(); right=right_value.lower()
    endpoint_min=min(left,right)
    m=max(arb(0),curvature_upper.upper())
    if not m>0:
        return endpoint_min
    h2=(width.upper()**2).upper()
    c=(m*h2/2).upper()
    d=abs(right-left).upper()
    if d >= c:
        return endpoint_min
    penalty=((c-d)**2/(4*c)).upper()
    return (endpoint_min-penalty).lower()
