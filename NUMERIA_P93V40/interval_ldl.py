"""Directed Arb LDL proof for leading principal prefixes or a signed inertia."""

from fractions import Fraction

from flint import arb, arb_mat


def dyadic(value):
    mantissa, exponent = value.man_exp()
    return [str(mantissa), int(exponent)]


def _interval(value):
    if not value.is_finite():
        return None
    return [dyadic(value.lower()), dyadic(value.upper())]


def _validate_matrix(matrix):
    if not isinstance(matrix, arb_mat) or not 0 < matrix.nrows() == matrix.ncols():
        raise ValueError("expected a nonempty square Arb matrix")
    return matrix.nrows()


def _validate_block_size(block_size):
    if type(block_size) is not int or block_size < 1:
        raise ValueError("block_size must be a positive integer")


def _shifted_lower_triangle(matrix, dimension, shift):
    result = []
    for i in range(dimension):
        row = []
        for j in range(i + 1):
            entry, mirror = matrix[i, j], matrix[j, i]
            if not entry.is_finite() or not mirror.is_finite():
                raise ValueError(f"nonfinite interval at ({i}, {j})")
            if entry.mid() != mirror.mid() or entry.rad() != mirror.rad():
                raise ValueError(f"asymmetric interval at ({i}, {j})")
            row.append(entry - shift if i == j else entry)
        result.append(row)
    return result


def _update_trailing(work, start, stop, left, right):
    size = len(work) - stop
    if not size:
        return
    width = stop - start
    lhs = arb_mat(size, width, [value for column in zip(*left, strict=True) for value in column])
    rhs = arb_mat(width, size, [value for row in right for value in row])
    update = lhs * rhs
    for i in range(size):
        row = work[stop + i]
        for j in range(i + 1):
            row[stop + j] -= update[i, j]


def _eliminate(work, expected_signs, block_size):
    size = len(work)
    pivots = []
    for start in range(0, size, block_size):
        stop = min(size, start + block_size)
        trailing_left, trailing_right = [], []
        for k in range(start, stop):
            pivot = work[k][k]
            sign = 1 if pivot > 0 else -1 if pivot < 0 else 0
            interval = _interval(pivot)
            if not sign or sign != expected_signs[k]:
                return pivots, {
                    "pivot_index": k + 1,
                    "pivot_interval_dyadic": interval,
                    "expected_sign": expected_signs[k],
                    "proven_sign": sign,
                    "status": "unexpected_pivot_sign" if sign else "indeterminate_pivot",
                }
            pivots.append(interval)
            multipliers = [work[i][k] / pivot for i in range(k + 1, size)]
            for i in range(k + 1, size):
                row = work[i]
                factor = row[k]
                for j in range(k + 1, min(i + 1, stop)):
                    row[j] -= factor * multipliers[j - k - 1]
            # All updates involving at least one panel column are now
            # applied. The trailing square receives the block product once.
            if stop < size:
                trailing_left.append([work[i][k] for i in range(stop, size)])
                trailing_right.append(multipliers[stop - k - 1:])
        _update_trailing(work, start, stop, trailing_left, trailing_right)
    return pivots, None


def _report(pivots, obstruction):
    return {
        "status": obstruction["status"] if obstruction else "certified",
        "certified": obstruction is None,
        "pivots_proved": len(pivots),
        "pivot_intervals_dyadic": pivots,
        "first_obstruction": obstruction,
    }


def certify_inertia(matrix, expected_signs, *, block_size=32):
    """Certify a prescribed sequence of strict LDL pivot signs.

    Sylvester inertia then follows by congruence for every symmetric matrix
    contained in the input ball matrix. The routine stops if a pivot fails
    its requested sign or if its interval contains zero.
    """
    dimension = _validate_matrix(matrix)
    _validate_block_size(block_size)
    signs = tuple(expected_signs)
    if len(signs) != dimension or any(type(sign) is not int or abs(sign) != 1 for sign in signs):
        raise ValueError("one expected strict sign (+1 or -1) is required for each pivot")
    work = _shifted_lower_triangle(matrix, dimension, arb(0))
    pivots, obstruction = _eliminate(work, signs, block_size)
    result = _report(pivots, obstruction)
    result["inertia"] = None if obstruction else {
        "positive": signs.count(1), "negative": signs.count(-1), "zero": 0,
    }
    return result


def certify_prefixes(matrix, sizes, margin=Fraction(0), *, block_size=32):
    """Certify A_N > margin*I at requested prefixes using one block LDL.

    A strictly negative shifted pivot proves failure of the requested
    margin, which implies a negative direction of A_N itself only at zero
    margin. A pivot containing zero leaves the sign inconclusive.
    """
    dimension = _validate_matrix(matrix)
    _validate_block_size(block_size)
    sizes = tuple(sizes)
    if (not sizes or any(type(n) is not int or not 1 <= n <= dimension for n in sizes)
            or sizes != tuple(sorted(set(sizes)))):
        raise ValueError("sizes must be strictly increasing positive integers within the matrix")
    if type(margin) not in (int, str, Fraction):
        raise ValueError("margin must be an exact nonnegative rational")
    margin = Fraction(margin)
    if margin < 0:
        raise ValueError("margin must be nonnegative")
    shift = arb(margin.numerator) / margin.denominator
    work = _shifted_lower_triangle(matrix, sizes[-1], shift)
    pivots, obstruction = _eliminate(work, (1,) * sizes[-1], block_size)
    result = _report(pivots, obstruction)
    negative = obstruction is not None and obstruction["proven_sign"] < 0
    result["margin_exact"] = str(margin)
    result["prefixes"] = [
        {
            "size": size,
            "status": ("certified" if size <= len(pivots) else
                       "negative_shifted_direction" if negative else "inconclusive"),
            "negative_direction_of_unshifted_matrix": bool(negative and margin == 0 and size > len(pivots)),
        }
        for size in sizes
    ]
    return result
