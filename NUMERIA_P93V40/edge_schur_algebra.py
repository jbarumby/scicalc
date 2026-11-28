"""Directed Schur reduction for a rank-one correction to a finite block."""

from flint import arb, arb_mat


class TailLoadingError(ArithmeticError):
    """The interval proof of positivity for the corrected tail did not close."""

    def __init__(self, loading, denominator):
        super().__init__("the tail loading is not certified below one")
        self.loading = loading
        self.denominator = denominator


def _symmetric_enclosure(matrix):
    """Intersect the two enclosures for each off-diagonal entry."""
    size = matrix.nrows()
    result = arb_mat(matrix)
    for i in range(size):
        for j in range(i):
            entry = matrix[i, j].intersection(matrix[j, i])
            result[i, j] = result[j, i] = entry
    return result


def finite_rank_schur(a, c, d, u, v, coefficient):
    """Enclose the Schur complement after subtracting gamma * (u,v)(u,v)^T.

    The exact symmetric tail D must be known positive by the caller.
    With z=D^-1 v and E=D^-1 C^T, the corrected tail is positive if
    1-gamma*v^T z > 0. The remaining Schur correction is
    gamma/(1-gamma*v^T z) times the outer product of u-Cz.
    """
    if not all(isinstance(matrix, arb_mat) for matrix in (a, c, d, u, v)):
        raise TypeError("the blocks and endpoint columns must be Arb matrices")
    head, tail = a.nrows(), d.nrows()
    if not (head > 0 and tail > 0 and a.ncols() == head and d.ncols() == tail
            and c.nrows() == head and c.ncols() == tail
            and u.nrows() == head and u.ncols() == 1
            and v.nrows() == tail and v.ncols() == 1):
        raise ValueError("incompatible head, cross, tail, or endpoint dimensions")
    if not isinstance(coefficient, arb) or not coefficient.is_finite() or not coefficient >= 0:
        raise ValueError("the rank-one coefficient must be a nonnegative finite Arb interval")
    if not all(entry.is_finite() for matrix in (a, c, d, u, v) for entry in matrix.entries()):
        raise ValueError("all matrix entries must be finite")
    a, d = _symmetric_enclosure(a), _symmetric_enclosure(d)
    cross_transpose = c.transpose()
    rhs = arb_mat([[*(cross_transpose[i, j] for j in range(head)), v[i, 0]]
                   for i in range(tail)])
    try:
        solved = d.solve(rhs)
    except ZeroDivisionError as exc:
        raise ArithmeticError("the Arb solve against the tail did not close") from exc
    if not all(entry.is_finite() for entry in solved.entries()):
        raise ArithmeticError("the Arb solve against the tail was not finite")
    cross_solution = arb_mat([[solved[i, j] for j in range(head)] for i in range(tail)])
    endpoint_solution = arb_mat(tail, 1, [solved[i, head] for i in range(tail)])
    tail_quadratic = (v.transpose() * endpoint_solution)[0, 0]
    tail_loading = coefficient * tail_quadratic
    denominator = 1 - tail_loading
    if not denominator > 0:
        raise TailLoadingError(tail_loading, denominator)

    base_schur = _symmetric_enclosure(a - c * cross_solution)
    effective_endpoint = u - c * endpoint_solution
    multiplier = coefficient / denominator
    correction = (effective_endpoint * effective_endpoint.transpose()) * multiplier
    schur = _symmetric_enclosure(base_schur - correction)
    corrected_elimination = cross_solution - (endpoint_solution * effective_endpoint.transpose()) * multiplier
    if not all(entry.is_finite() for entry in schur.entries()):
        raise ArithmeticError("the Schur complement enclosure was not finite")
    return {
        "schur": schur,
        "base_schur": base_schur,
        "tail_quadratic": tail_quadratic,
        "tail_loading": tail_loading,
        "denominator": denominator,
        "multiplier": multiplier,
        "effective_endpoint": effective_endpoint,
        "corrected_elimination": corrected_elimination,
    }
