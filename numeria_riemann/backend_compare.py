"""Cross-check zeta evaluations with mpmath and optional python-flint.

This is a numerical cross-check, not a proof or a zero certificate.
"""
from __future__ import annotations

from typing import Iterable

import mpmath as mp


def compare_zeta_values(
    points: Iterable[complex | tuple[float, float]],
    *,
    dps: int = 30,
) -> list[dict[str, float | str | bool]]:
    """Compare mpmath.zeta(s) and python-flint acb.zeta(s), if installed.

    Each point may be a Python complex or a (real, imaginary) pair.
    The absolute difference is a diagnostic only: agreement between
    implementations does not by itself certify either result.
    """
    try:
        from flint import acb, ctx
    except ImportError as exc:
        raise ImportError(
            "python-flint is optional; install it with "
            "'python -m pip install python-flint'"
        ) from exc

    if not isinstance(dps, int) or dps < 15:
        raise ValueError("dps must be an integer >= 15")

    rows: list[dict[str, float | str | bool]] = []
    old_flint_dps = ctx.dps
    try:
        ctx.dps = dps
        with mp.workdps(dps):
            for point in points:
                if isinstance(point, tuple):
                    if len(point) != 2:
                        raise ValueError("tuple points must be (real, imaginary)")
                    re, im = point
                else:
                    re, im = point.real, point.imag
                s_mp = mp.mpc(re, im)
                z_mp = mp.zeta(s_mp)
                z_flint = acb(str(re), str(im)).zeta()
                # Use midpoint values for a transparent diagnostic. python-flint's
                # enclosure/radius should be inspected separately for certification.
                flint_mid = z_flint.mid()
                z_flint_mp = mp.mpc(str(flint_mid.real), str(flint_mid.imag))
                difference = abs(z_mp - z_flint_mp)
                rows.append({
                    "real": float(re),
                    "imag": float(im),
                    "mpmath_abs": float(abs(z_mp)),
                    "python_flint_abs": float(abs(z_flint_mp)),
                    "absolute_difference": float(difference),
                    "status": "crosscheck_only_not_certificate",
                })
    finally:
        ctx.dps = old_flint_dps
    return rows


def zero_residual_report(count: int = 10, *, dps: int = 40) -> list[dict[str, float | int | str]]:
    """Evaluate |zeta(rho_n)| at mpmath's first zero approximations."""
    if not isinstance(count, int) or count < 1:
        raise ValueError("count must be a positive integer")
    if not isinstance(dps, int) or dps < 15:
        raise ValueError("dps must be an integer >= 15")
    rows = []
    with mp.workdps(dps):
        for n in range(1, count + 1):
            rho = mp.zetazero(n)
            residual = abs(mp.zeta(rho))
            rows.append({
                "index": n,
                "real_part": float(mp.re(rho)),
                "ordinate": float(mp.im(rho)),
                "zeta_abs_residual": float(residual),
                "status": "numerical_residual_not_zero_certificate",
            })
    return rows
