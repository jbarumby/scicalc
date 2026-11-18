#!/usr/bin/env python3
"""Extract a rigorously boxed even-sector 4x4 prefix for P231/Lean.

Uses the exact P228 recomposition (all 14 Arb chunks, beta, poles and P229
quadrature error), rounds the P193 center to 40 decimal places, and checks
entrywise coverage plus an exact shifted LDL factorization.
"""
from __future__ import annotations
import argparse, hashlib, json
from fractions import Fraction as F
from pathlib import Path
from p228_recompose_audit import audit_sector, beta_interval, dyad, pole_component, read_chunks, read_p193
from p229_quadrature_bound import total_error_bound

N = 4
T = 200
SCALE = 10**40
EPSILON = F(1, SCALE)
SHIFT = F(1, 10**20)
RHO = 4 * EPSILON


def round_fixed(x: F, scale: int) -> F:
    """Nearest multiple of 1/scale, ties to even, using integer arithmetic."""
    y = x * scale
    sign = -1 if y < 0 else 1
    q, r = divmod(abs(y.numerator), y.denominator)
    if 2*r > y.denominator or (2*r == y.denominator and q % 2):
        q += 1
    return F(sign*q, scale)


def ldl_exact(matrix):
    """Return unit-lower L and diagonal d such that matrix = L*diag(d)*L^T."""
    n = len(matrix)
    L = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    d = []
    for k in range(n):
        pivot = matrix[k][k] - sum(L[k][a]**2 * d[a] for a in range(k))
        if pivot <= 0:
            raise ArithmeticError(f'nonpositive exact LDL pivot {k}: {pivot}')
        d.append(pivot)
        for i in range(k+1, n):
            numerator = matrix[i][k] - sum(L[i][a]*L[k][a]*d[a] for a in range(k))
            L[i][k] = numerator / pivot
    # Independent exact multiplication check for the factorization.
    for i in range(n):
        for j in range(n):
            got = sum(L[i][k]*d[k]*L[j][k] for k in range(n))
            if got != matrix[i][j]:
                raise ArithmeticError(f'LDL identity failed at {(i,j)}')
    return L, d


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--chunks', type=Path, required=True)
    p.add_argument('--p193', type=Path, required=True)
    p.add_argument('--p230', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()

    sums, manifest = read_chunks(a.chunks, 'even', T)
    _ref, _E193, cert, centers = read_p193(a.p193, 'even')
    p230 = json.loads(a.p230.read_text())
    p230_sector = next(s for s in p230['sectors'] if s['sector'] == 'even')
    if p230_sector['P193_center_matrix_gzip_sha256'] != cert['matrix_gz_sha256']:
        raise ValueError('P230 and P231 do not refer to the same P193 center matrix')
    E227 = total_error_bound(T)
    beta = beta_interval(T)
    poles = [pole_component(2*k) for k in range(N)]
    target_boxes, rounded, raw_centers = [], [], []
    max_center_rounding = F(0)
    max_box_deviation = F(0)
    box_rows = []
    for i in range(N):
        row=[]; crow=[]; rawrow=[]
        for j in range(N):
            if i <= j:
                lo, hi = sums[(i,j)]
                if i == j:
                    lo += beta[0]; hi += beta[1]
                pl, ph = poles[i], poles[j]
                lo += 2*pl[0]*ph[0] - E227
                hi += 2*pl[1]*ph[1] + E227
                c0 = centers[(i,j)]
            else:
                lo, hi = target_boxes[j][i]
                c0 = centers[(j,i)]
            cp = round_fixed(c0, SCALE)
            dev = max(abs(lo-cp), abs(hi-cp))
            if dev > EPSILON:
                raise ArithmeticError(f'entry box outside epsilon at {(i,j)}: {dev}')
            rounding = abs(c0-cp)
            max_center_rounding=max(max_center_rounding,rounding)
            max_box_deviation=max(max_box_deviation,dev)
            row.append((lo,hi));crow.append(cp);rawrow.append(c0)
            box_rows.append([i,j,str(lo),str(hi),str(cp),str(dev)])
        target_boxes.append(row);rounded.append(crow);raw_centers.append(rawrow)

    if max_center_rounding > F(1,2*SCALE):
        raise ArithmeticError('round-to-nearest error exceeded half a unit')
    shifted = [[rounded[i][j] - (SHIFT if i == j else 0) for j in range(N)] for i in range(N)]
    L, d = ldl_exact(shifted)
    full_audit = audit_sector(a.chunks, a.p193, 'even')
    full_row_bound = F(full_audit['P227_to_P193center_infinity_norm_upper_exact'])
    if full_row_bound != F(p230_sector['P227_to_P193center_infinity_norm_upper_exact']):
        raise ValueError('P228 recomputation disagrees with the P230 row bound')
    center_lower = F(p230_sector['center_spd_certificate']['lambda_min_center_lower_exact'])
    target_lower = F(p230_sector['lambda_min_P227_finite_lower_exact'])
    if center_lower <= 0 or target_lower <= 0:
        raise ValueError('P230 does not provide a positive finite lower bound')

    def matstr(m): return [[str(x) for x in row] for row in m]
    out = {
        'stage':'NUMERIA ES8-P231 exact even-sector prefix4 extraction',
        'sector':'even','T':T,'prefix_modes':[0,1,2,3],'matrix_dimension':N,
        'P193_center_gzip_sha256':cert['matrix_gz_sha256'],
        'P230_center_spd_lower_bound_exact':str(center_lower),
        'P230_200_mode_target_lambda_lower_bound_exact':str(target_lower),
        'P227_chunk_manifest':manifest,
        'rounding_scale':str(SCALE),'entrywise_epsilon_exact':str(EPSILON),
        'quadratic_error_rho_exact':str(RHO),'center_shift_exact':str(SHIFT),
        'P228_full_200_mode_row_bound_exact':str(full_row_bound),
        'max_center_rounding_error_exact':str(max_center_rounding),
        'max_P227_box_to_rounded_center_exact':str(max_box_deviation),
        'rounded_center_matrix_exact':matstr(rounded),
        'source_P193_center_prefix_exact':matstr(raw_centers),
        'P227_target_boxes_exact_upper_triangle_and_reflection':box_rows,
        'shifted_LDL_unit_lower_exact':matstr(L),
        'shifted_LDL_diagonal_exact':[str(x) for x in d],
        'checks':{'all_16_reflected_boxes_within_epsilon':True,
                  'rounding_error_at_most_half_unit':True,
                  'all_4_shifted_LDL_pivots_strictly_positive':True,
                  'shifted_LDL_identity_exact':True,
                  'rho_strictly_less_than_shift':RHO<SHIFT},
        'interpretation':'Exact rational data for a four-mode finite prefix only. Inclusion of the mathematical Weil entries in P227 boxes remains conditional on P227/P229 formulas; no global or RH conclusion.'
    }
    a.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='P227_target_boxes_exact_upper_triangle_and_reflection'},indent=2))

if __name__=='__main__':
    main()
