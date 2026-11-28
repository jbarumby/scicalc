"""Certified spherical Bessel sequences from a recurrence on exact Arb points."""

from flint import arb, ctx


def _oscillatory(x, inv, last):
    sine, cosine = x.sin_cos()
    first = sine * inv
    if last == 0:
        return [first]
    second = (first - cosine) * inv
    previous, current = first.mid(), second.mid()
    error = first.rad() + second.rad()
    centers = []
    for n in range(1, last):
        estimate = ((2 * n + 1) * inv) * current - previous
        previous, current = current, estimate.mid()
        error += estimate.rad()
        centers.append(current)

    # For a=(2*n+1)/x, the map (u,v)->(a*u-v,u) preserves the square
    # (u*u+v*v-a*u*v)/(1+a/2). The resulting norm decreases as a increases.
    # Each local rounding error contributes at most its radius to this norm.
    # Its square is >= (1-a/2)*u*u, and 1-a/2 >= 1/(2*x) up to last.
    radius = (error * (2 * x).sqrt()).upper()
    return [first, second, *(arb(center, radius) for center in centers)]


def _ratios(inv, last, degree, twice_x_ceiling):
    # The j_k power series gives j_k(x) ~ x**k/(2*k+1)!! for fixed x
    # as k->infinity. Thus 0<r_k=j_(k+1)/j_k<1 eventually. The identity
    # r_(k-1)=1/(a_k-r_k) maps [0,1] into [0,1] for all k>floor(x),
    # so descending induction establishes that enclosure even at cutoff.
    cutoff = max(degree, twice_x_ceiling) + (ctx.prec + 2) // 3
    unit = arb(arb(1) / 2, arb(1) / 2)
    ratio = unit
    reverse = []
    for k in range(cutoff, last, -1):
        denominator = (2 * k + 1) * inv - ratio
        if not denominator > 0:
            raise ArithmeticError("unable to enclose a positive Bessel ratio denominator")
        ratio = (1 / denominator).intersection(unit)
        if k <= degree:
            reverse.append(ratio)
    # At k>=ceil(2*x), a_k>=4, so the tail ambiguity contracts by <=1/9
    # per step before reaching a retained ratio. This gives at least ctx.prec
    # bits of analytic contraction. Arb separately encloses rounding errors.
    return reversed(reverse)


def bessels_at_point(x, degree):
    """Return enclosures of j_0(x),...,j_degree(x) for an exact positive x."""
    if not isinstance(x, arb) or not x.is_exact() or not x.is_finite() or not x > 0:
        raise ValueError("expected an exact positive finite Arb point")
    if type(degree) is not int or degree < 0:
        raise ValueError("expected a nonnegative integer Bessel degree")
    inv = 1 / x
    last = min(degree, int(x.floor().unique_fmpz()))
    values = _oscillatory(x, inv, last)
    if last < degree:
        twice_x_ceiling = int((2 * x).ceil().unique_fmpz())
        for ratio in _ratios(inv, last, degree, twice_x_ceiling):
            values.append(values[-1] * ratio)
    if not all(value.is_finite() for value in values):
        raise ArithmeticError("the Bessel recurrence did not produce finite enclosures")
    return values
