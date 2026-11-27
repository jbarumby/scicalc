"""Propose a triangular congruence numerically and verify it exactly in Arb."""

from fractions import Fraction

from flint import arb, arb_mat

from .interval_ldl import dyadic


def _margin(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def prove_negative(matrix, v, margin=0):
    """Test both signs of interest for an exact dyadic Arb vector."""
    if v.ncols() != 1 or v.nrows() != matrix.nrows() or not all(x.is_exact() for x in v.entries()):
        raise ValueError("a matching exact vector is required")
    q = (v.transpose() * (matrix * v))[0, 0]
    shifted = q - _margin(margin) * (v.transpose() * v)[0, 0]
    report = {
        "proven_negative": bool(q < 0),
        "proven_below_requested_margin": bool(shifted < 0),
        "rayleigh_numerator_dyadic": [dyadic(q.lower()), dyadic(q.upper())],
        "shifted_numerator_dyadic": [dyadic(shifted.lower()), dyadic(shifted.upper())],
    }
    if shifted < 0 or q < 0:
        report["witness_dyadic"] = [dyadic(x) for x in v.entries()]
    return report


def verify_triangular_congruence(matrix, r, sizes, margin):
    """Check positively dominant leading blocks of R^T(A-margin I)R.

    An upper triangular R with nonzero exact diagonal preserves the
    leading subspaces. Therefore each leading block uses the matching
    principal block of A, so a single product serves all requested sizes.
    """
    n = matrix.nrows()
    if r.nrows() != n or r.ncols() != n or matrix.ncols() != n:
        raise ValueError("incompatible square matrices")
    if not all(r[i, i].is_exact() and not r[i, i].is_zero() for i in range(n)):
        raise ValueError("the triangular diagonal must be exact and nonzero")
    if not all(r[i,j].is_exact() and (i<=j or r[i,j].is_zero())
               for i in range(n) for j in range(n)):
        raise ValueError("the congruence must be exact upper triangular")
    shifted = arb_mat(matrix)
    h = _margin(margin)
    for i in range(n):
        shifted[i, i] -= h
    transformed = r.transpose() * (shifted * r)
    if not all(entry.is_finite() for entry in transformed.entries()):
        raise ArithmeticError("nonfinite congruence enclosure")
    accum = [arb(0) for _ in range(n)]
    results = []
    pending = set(sizes)
    for i in range(n):
        for j in range(i):
            contribution = abs(transformed[i,j].intersection(transformed[j,i])).upper()
            accum[i] += contribution
            accum[j] += contribution
        if i+1 in pending:
            bound = min((transformed[j,j]-accum[j]).lower() for j in range(i+1))
            results.append({
                "size": i+1,
                "certified": bool(bound > 0),
                "row_margin_lower_dyadic": dyadic(bound),
            })
    return results


def _point_factor(matrix, h, panel_width=32):
    """A heuristic point LDL. Discarding radii here only proposes a basis."""
    n = matrix.nrows()
    work = [[(matrix[i,j]-h if i==j else matrix[i,j]).mid() for j in range(i+1)]
            for i in range(n)]
    lower = arb_mat(n,n)
    scaling = []
    for start in range(0,n,panel_width):
        stop = min(n,start+panel_width)
        raw, multipliers = [], []
        for k in range(start,stop):
            pivot = work[k][k]
            lower[k,k] = 1
            if not pivot > 0:
                return lower, scaling, k
            scaling.append((1/pivot.sqrt()).mid())
            column = [(work[i][k]/pivot).mid() for i in range(k+1,n)]
            for i in range(k+1,n):
                lower[i,k] = column[i-k-1]
                entry = work[i][k]
                for j in range(k+1,min(i+1,stop)):
                    work[i][j] = (work[i][j]-entry*column[j-k-1]).mid()
            if stop < n:
                raw.append([work[i][k] for i in range(stop,n)])
                multipliers.append(column[stop-k-1:])
        if stop < n:
            trailing = n-stop
            product = (arb_mat(stop-start,trailing,[x for col in raw for x in col]).transpose()
                       * arb_mat(stop-start,trailing,[x for col in multipliers for x in col]))
            for i in range(trailing):
                for j in range(i+1):
                    work[stop+i][stop+j] = (work[stop+i][stop+j]-product[i,j]).mid()
    return lower,scaling,None


def certify_by_congruence(matrix, sizes, margin):
    """Use an Arb midpoint proposal, then recompute the congruence in Arb."""
    n = matrix.nrows()
    lower, scaling, obstruction = _point_factor(matrix, _margin(margin))
    length = n if obstruction is None else obstruction+1
    try:
        tri = arb_mat([[lower[i,j] for j in range(length)] for i in range(length)])
        inverse = tri.transpose().inv().mid()
    except (ArithmeticError, ZeroDivisionError, ValueError):
        return {"status": "inconclusive_proposal", "prefixes": []}
    if obstruction is not None:
        v = arb_mat(n,1,[inverse[i,obstruction] if i<length else arb(0) for i in range(n)])
        result = prove_negative(matrix,v,margin)
        prefixes = []
        if obstruction:
            leading = arb_mat([[matrix[i,j] for j in range(obstruction)] for i in range(obstruction)])
            r = arb_mat([[ (inverse[i,j]*scaling[j]).mid() if i<=j else arb(0)
                           for j in range(obstruction)] for i in range(obstruction)])
            targets = sorted({obstruction, *(k for k in sizes if k<=obstruction)})
            prefixes = verify_triangular_congruence(leading,r,targets,margin)
        if result["proven_negative"]:
            status = "negative_direction"
        elif result["proven_below_requested_margin"]:
            status = "strictly_below_requested_margin"
        else:
            status = "inconclusive_proposal"
        return {
            "status": status,
            "unproved_point_pivot": obstruction+1,
            "rayleigh_check": result,
            "prefixes": prefixes,
        }
    r = arb_mat(n,n,[(inverse[i,j]*scaling[j]).mid() if i<=j else arb(0)
                     for i in range(n) for j in range(n)])
    results = verify_triangular_congruence(matrix,r,sizes,margin)
    return {"status": "certified" if all(row["certified"] for row in results)
            else "inconclusive_congruence", "prefixes": results}
