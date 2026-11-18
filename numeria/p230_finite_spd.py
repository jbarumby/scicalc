#!/usr/bin/env python3
"""Certify finite-block SPD by Arb LDL + inverse norm + exact Weyl bound.

The target matrix is enclosed by the exact rational intervals recomposed in
P228. A symmetric rational P193 center matrix is certified SPD independently;
its inverse Frobenius norm gives a coarse eigenvalue lower bound. The exact
P228 infinity-norm perturbation then transfers positivity to the P227 target.
"""
from __future__ import annotations
import argparse, gzip, json, time
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path

import flint
from flint import arb, arb_mat, ctx
from p228_recompose_audit import audit_sector, dyad

PREC_BITS = 512
LAMBDA_CENTER = F(1, 10**20)
INV_FROB_SQUARED_CEILING = 10**40


def frac_decimal(x: F, digits=12) -> str:
    with localcontext() as c:
        c.prec = digits + 8
        y = Decimal(x.numerator) / Decimal(x.denominator)
        return f'{y:.{digits}E}'


def load_center_matrix(p193_dir: Path, sector: str):
    path = p193_dir / f'P193_matrix_{sector}.json.gz'
    with gzip.open(path, 'rt') as f:
        data = json.load(f)
    n = data.get('N')
    if n != 200 or data.get('sector') != sector:
        raise ValueError(f'P193 matrix metadata mismatch: {path}')
    raw = data['mid']
    exact = [[dyad(x) for x in row] for row in raw]
    for i in range(n):
        for j in range(i):
            if exact[i][j] != exact[j][i]:
                raise ValueError(f'P193 center matrix is not exactly symmetric at {(i,j)}')
    # Arb balls contain the exact dyadic center entries after rounding at 512 bits.
    rows = [[arb(int(m)) * (arb(2) ** int(e)) for m, e in row] for row in raw]
    return data, exact, arb_mat(rows)


def certify_center_spd(A: arb_mat):
    """Interval LDL proves SPD; interval inverse gives a norm-based lambda bound."""
    n = A.nrows()
    if n != 200 or A.ncols() != n:
        raise ValueError('expected a 200x200 matrix')
    # Symmetric Schur complements; strict positivity of every interval pivot
    # proves positivity of every leading principal minor (Sylvester criterion).
    schur = [[A[i, j] for j in range(i + 1)] for i in range(n)]
    for k in range(n):
        pivot = schur[k][k]
        if not (pivot > 0):
            raise ArithmeticError(f'LDL pivot {k} not proved positive: {pivot!r}')
        for i in range(k + 1, n):
            aik = schur[i][k]
            for j in range(k + 1, i + 1):
                schur[i][j] = schur[i][j] - aik * schur[j][k] / pivot

    inv = A.inv()
    sum_sq = arb(0)
    for i in range(n):
        for j in range(n):
            x = inv[i, j]
            if not x.is_finite():
                raise ArithmeticError(f'nonfinite inverse enclosure at {(i,j)}')
            u = x.abs_upper()
            sum_sq += u * u
    threshold = arb(INV_FROB_SQUARED_CEILING)
    if not (sum_sq < threshold):
        raise ArithmeticError(f'could not prove Frobenius norm squared < 1e40: {sum_sq!r}')
    # For SPD A, lambda_min(A)=1/||A^-1||_2 >= 1/||A^-1||_F > 10^-20.
    return {'ldl_positive_pivots': n,
            'inverse_frobenius_squared_upper_less_than': '10^40',
            'lambda_min_center_lower_exact': str(LAMBDA_CENTER),
            'inverse_frobenius_squared_arb_enclosure': repr(sum_sq)}


def certify_sector(chunks: Path, p193_dir: Path, sector: str):
    started = time.time()
    ctx.prec = PREC_BITS
    audit = audit_sector(chunks, p193_dir, sector)
    data, exact, A = load_center_matrix(p193_dir, sector)
    center_cert = certify_center_spd(A)
    row_error = F(audit['P227_to_P193center_infinity_norm_upper_exact'])
    lambda_target = LAMBDA_CENTER - row_error
    if lambda_target <= 0:
        raise ArithmeticError(f'Weyl margin not positive in {sector}: {lambda_target}')
    half_margin = F(1, 2) * LAMBDA_CENTER
    if not row_error < half_margin:
        raise ArithmeticError(f'perturbation is not below half the center margin in {sector}')
    # Exact symmetry plus row sums gives ||B-A||_2 <= ||B-A||_infinity.
    # Weyl then proves lambda_min(B) >= lambda_min(A)-||B-A||_2.
    return {'sector': sector,
            'dimension': 200,
            'P193_center_matrix_gzip_sha256': audit['P193_source_sha256'],
            'P193_center_exactly_symmetric': True,
            'center_spd_certificate': center_cert,
            'P227_interval_entries_compared': audit['upper_triangle_entries'],
            'P227_to_P193center_infinity_norm_upper_exact': str(row_error),
            'P227_to_P193center_infinity_norm_upper_decimal': frac_decimal(row_error),
            'lambda_min_P227_finite_lower_exact': str(lambda_target),
            'lambda_min_P227_finite_lower_decimal': frac_decimal(lambda_target),
            'simple_uniform_lower_bound_exact': str(half_margin),
            'simple_uniform_lower_bound_holds': True,
            'P227_target_symmetry_basis': 'bilinear entries use the symmetric product u_n(t)u_m(t); beta is diagonal and pole terms are symmetric, so the stored upper triangle defines a real symmetric block',
            'method': 'exact P228 rational row-sum perturbation + Weyl inequality; P193 center SPD proved by 512-bit Arb interval LDL and inverse-Frobenius bound',
            'elapsed_seconds': round(time.time() - started, 2)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--chunks', type=Path, required=True)
    p.add_argument('--p193', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    result = {'stage': 'NUMERIA ES8-P230 finite 200x200 SPD perturbation certificate',
              'arb_version': getattr(flint, '__version__', 'unknown'),
              'arb_precision_bits': PREC_BITS,
              'sectors': [certify_sector(a.chunks, a.p193, s) for s in ('even', 'odd')],
              'scope': 'finite reduced 200x200 matrices only; no P207 tail/coupling proof, global Weil positivity, or RH proof'}
    a.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
