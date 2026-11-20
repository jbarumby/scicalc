"""P254: exact-rational floor-upgrade and block-Schur checks.

Analytic premise imported from P253: Psi_4(t)-1/4 >= 0 for |t|>=504.
Source data, if supplied, must be rigorous continuous-integral enclosures.
This module does not reexecute P253 Arb quadrature.
"""
from fractions import Fraction as Q


def iv(x):
    if isinstance(x, tuple):
        if len(x) != 2:
            raise ValueError("interval must have two endpoints")
        lo, hi = map(Q, x)
    else:
        lo = hi = Q(x)
    if lo > hi:
        raise ValueError("reversed interval")
    return lo, hi


def plus(a, b):
    x, y = iv(a), iv(b)
    return x[0] + y[0], x[1] + y[1]


def minus(a, b):
    x, y = iv(a), iv(b)
    return x[0] - y[1], x[1] - y[0]


def times(a, b):
    x, y = iv(a), iv(b)
    products = (x[0]*y[0], x[0]*y[1], x[1]*y[0], x[1]*y[1])
    return min(products), max(products)


def maxabs(a):
    lo, hi = iv(a)
    return max(abs(lo), abs(hi))


def dyadic(pair):
    if len(pair) != 2:
        raise ValueError("dyadic point must be (mantissa, exponent)")
    m, e = map(int, pair)
    return Q(m * (1 << e)) if e >= 0 else Q(m, 1 << -e)


def delta_from_beta(beta, *, new_floor=Q(1, 4), cutoff=504):
    """1/4 is certified for four-prime residual only at cutoff >=504.

    Caller must independently replay/authenticate the P253 tail certificate.
    """
    if cutoff < 504:
        raise ValueError("raised floor NOT certified on [309,504); use cutoff>=504")
    delta = minus(iv(new_floor), iv(beta))
    if delta[0] <= 0:
        raise ArithmeticError("strict raised floor not established")
    return delta


def floor_upgraded_intervals(old, gram, delta):
    """Enclose M_new=M_old+delta*(I-G_T) with exact rational endpoints.

    G_T must be the CONTINUOUS Fourier concentration Gram, including
    quadrature error, in the same orthonormal basis and normalization.
    """
    d = iv(delta)
    if d[0] < 0:
        raise ValueError("floor increment must be nonnegative")
    n = len(old)
    if n == 0 or len(gram) != n or any(len(row) != n for row in old) or any(len(row) != n for row in gram):
        raise ValueError("old and gram must be equally sized square matrices")
    result = []
    for i in range(n):
        row = []
        for j in range(n):
            a, g = iv(old[i][j]), iv(gram[i][j])
            if iv(old[j][i]) != a or iv(gram[j][i]) != g:
                raise ValueError("matrix intervals must be symmetric")
            row.append(plus(a, times(d, minus(int(i == j), g))))
        result.append(row)
    return result


def gershgorin_positive_lower(matrix, target=Q(0)):
    """Certify a lower eigenvalue bound using rational interval Gershgorin."""
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("square matrix required")
    for i in range(n):
        for j in range(n):
            if iv(matrix[i][j]) != iv(matrix[j][i]):
                raise ValueError("must be symmetric")
    lower = min(iv(matrix[i][i])[0] - sum(maxabs(matrix[i][j]) for j in range(n) if i != j) for i in range(n))
    if lower <= Q(target):
        raise ArithmeticError("Gershgorin positivity inconclusive")
    return lower


def block_schur_gate(a_lower, d_lower, cross, margin=Q(0)):
    """Strict block criterion conditional on certified A>=aI and D>=dI.

    If m<a,d and (a-m)(d-m) > ||B||_F^2, then [[A,B],[B*,D]]>=mI.
    """
    a, d, m = Q(a_lower), Q(d_lower), Q(margin)
    if a <= m or d <= m or not cross:
        raise ArithmeticError("positive diagonal gaps and nonempty coupling required")
    width = len(cross[0])
    if width == 0 or any(len(row) != width for row in cross):
        raise ValueError("rectangular nonempty coupling required")
    cross_sq = sum(maxabs(entry)**2 for row in cross for entry in row)
    slack = (a-m)*(d-m) - cross_sq
    if slack <= 0:
        raise ArithmeticError("block-Schur criterion inconclusive")
    return {"a_lower":str(a), "d_lower":str(d), "margin":str(m),
            "frobenius_cross_squared_upper":str(cross_sq),
            "strict_schur_slack":str(slack), "certified_for_supplied_bounds":True}
