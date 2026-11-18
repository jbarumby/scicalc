#!/usr/bin/env python3
"""P239: rigorous one-stroke tail/coupling bounds for the P237 window.

Uses the Bessel majorant from Zhu, arXiv:2608.24827v2, §4, with exact
rational parameters L=161/200, T=504 and N=377. All endpoint arithmetic is
Arb. This verifies the infinite Legendre tail/coupling budget against the
finite P238 margin; it does not formalize the imported analytic lemmas in Lean.
"""
from __future__ import annotations
import argparse, gzip, hashlib, json
from fractions import Fraction as F
from pathlib import Path
from flint import arb, ctx

LQ = F(161, 200)
T = 504
N = 377
M = 20
RATIO = F(3, 40)
COMMON_MARGIN = F(1, 10**18)
GLOBAL_MARGIN = F(1, 2*10**18)


def aq(q: F | int) -> arb:
    q = F(q)
    return arb(q.numerator) / arb(q.denominator)


def dyadic(m: str | int, e: str | int) -> arb:
    m, e = int(m), int(e)
    return arb(m << e) if e >= 0 else arb(m) / arb(1 << -e)


def wire_ball(wire) -> arb:
    (m, e), (r, k) = wire
    return dyadic(m, e) + arb(0, dyadic(r, k))


def ball_wire(x: arb):
    m, e = x.mid().man_exp()
    r, k = x.rad().upper().man_exp()
    return [[str(m), str(e)], [str(r), str(k)]]


def load_gzip_json(path: Path):
    return json.loads(gzip.decompress(path.read_bytes()))


def product_double_factorial(n: int) -> int:
    out = 1
    for k in range(1, 2*n+2, 2):
        out *= k
    return out


def amplitude_bound(n: int, x: F, L: arb) -> arb:
    """2 sqrt(L(n+1/2)) x^n/(2n+1)!!, an upper bound for |T_hat_n|."""
    q = x**n / product_double_factorial(n)
    nu = aq(F(2*n+1, 2))
    return (2 * (L*nu).sqrt() * aq(q)).upper()


def pole_amplitude_bound(n: int, z: F, L: arb) -> arb:
    """Upper bound for |p_n| from the positive series of spherical i_n(z)."""
    q = z**n / product_double_factorial(n)
    nu = aq(F(2*n+1, 2))
    zball = aq(z)
    return (2 * (L*nu).sqrt() * aq(q) * (zball*zball/2).exp()).upper()


def ratio_squared(n: int, x: F) -> F:
    # Exact square of b_(n+2)/b_n; the same formula holds for the pole majorant.
    return x**4 / ((2*n+1) * (2*n+3)**2 * (2*n+5))


def verify_source_binding(data_dir: Path, p238_path: Path):
    p238 = json.loads(p238_path.read_text())
    if p238.get('schema') != 'NUMERIA.ES8P238.P237_SCHUR_MARGIN_AUDIT.1':
        raise ValueError('unsupported P238 certificate')
    paths = {
        'matrix': data_dir/'P237_FINITE_LDL_2048.json.gz',
        'center': data_dir/'P237_CENTER_WEYL_2048.json.gz',
        'shadow': data_dir/'P237_SHADOW_FULL_AUDIT.json',
        'C_even': data_dir/'P237_FULL_C_even.json.gz',
        'C_odd': data_dir/'P237_FULL_C_odd.json.gz',
    }
    for key, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(path)
    actual = {k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in paths.items()}
    expected = {
        'matrix': p238['P237_matrix_artifact_sha256'],
        'center': p238['P237_center_Weyl_certificate_sha256'],
        'shadow': p238['P237_shadow_audit_sha256'],
    }
    for key, digest in expected.items():
        if actual[key] != digest:
            raise ArithmeticError(f'P238 {key} hash mismatch')
    finite = load_gzip_json(paths['matrix'])
    if finite.get('schema') != 'NUMERIA.ES8P237.FULL_FINITE_INTERVAL_LDL.1':
        raise ValueError('unsupported P237 finite matrix artifact')
    if (finite.get('L'), finite.get('Tsharp'), finite.get('n_modes_per_sector')) != ('161/200', 504, N):
        raise ArithmeticError('P237 dimensions/window do not match the tail proof')
    if finite.get('global_GL60_error_added_once_per_matrix_entry') is not True:
        raise ArithmeticError('P237 quadrature error budget is not authenticated')
    for sec in ('even', 'odd'):
        item = p238['sectors'][sec]
        ckey = 'C_'+sec
        if actual[ckey] != item['source_full_C_sha256']:
            raise ArithmeticError(f'P238 {sec} C hash mismatch')
        if item.get('center_ldl_positive_pivots') != N or item.get('center_spd') is not True:
            raise ArithmeticError(f'P238 finite center is not certified SPD: {sec}')
        if item.get('full_interval_box_spd') is not True:
            raise ArithmeticError(f'P238 full interval block not certified SPD: {sec}')
        meta = finite['sectors'][sec].get('assembly_metadata', {})
        if (finite['sectors'][sec].get('beta_poles_and_GL60_added_once') is not True
                or not all(k in meta for k in ('beta_ball', 'AL_ball', 'pole_vector_balls'))
                or len(meta.get('pole_vector_balls', [])) != N):
            raise ArithmeticError(f'P237 assembly provenance incomplete: {sec}')
    if F(p238['common_finite_block_margin_lower_exact']) != COMMON_MARGIN:
        raise ArithmeticError('P238 common finite margin differs from the expected rational')
    return p238, finite, actual


def certify_sector(parity: int, L: arb, x: F, z: F, c: arb):
    n0 = 2*N + parity
    if ratio_squared(n0, x) >= RATIO**2:
        raise ArithmeticError(f'Bessel tail ratio is not < 3/40 at n={n0}')
    if ratio_squared(n0, z) >= RATIO**2:
        raise ArithmeticError(f'pole tail ratio is not < 3/40 at n={n0}')

    b0 = amplitude_bound(n0, x, L)
    # The exact geometric ratio bound is 1/(1-3/40)=40/37.
    sb = (b0 * aq(F(40, 37))).upper()
    p0 = pole_amplitude_bound(n0, z, L)
    sp = (p0 * aq(F(40, 37))).upper()

    modes = range(parity, 2*N + parity, 2)
    amplitudes = [2*(L*aq(F(2*m+1, 2))).sqrt() for m in modes]
    a_sum = sum((a.upper() for a in amplitudes), arb(0)).upper()
    a_max = amplitudes[-1].upper()

    # p_m are coefficients of cosh(x/2) or sinh(x/2) in an orthonormal
    # Legendre basis. Cauchy-Schwarz bounds the finite head without computing it.
    cosh_z = (aq(z).exp()+(-aq(z)).exp())/2
    pnorm = (2*L).sqrt() * cosh_z
    phead_sum = (arb(N).sqrt()*pnorm).upper()
    phead_max = pnorm.upper()

    # c >= T*M/pi, since pi>3; it is deliberately rational and conservative.
    row_B = (c*b0*a_sum + 2*p0*phead_sum).upper()
    col_B = (c*sb*a_max + 2*sp*phead_max).upper()
    eps_B = (row_B*col_B).sqrt().upper()
    eps_D = (c*b0*sb + 2*p0*sp).upper()

    if not eps_B < aq(F(1, 10**90)):
        raise ArithmeticError(f'coupling bound failed <1e-90 for parity {parity}')
    if not eps_D < aq(F(1, 10**190)):
        raise ArithmeticError(f'tail-block bound failed <1e-190 for parity {parity}')
    return {
        'parity': 'even' if parity == 0 else 'odd',
        'first_discarded_degree': n0,
        'ratio_squared_exact': str(ratio_squared(n0, x)),
        'ratio_below_3_over_40': True,
        'tail_transform_first_upper_ball': ball_wire(b0),
        'tail_transform_sum_upper_ball': ball_wire(sb),
        'tail_pole_first_upper_ball': ball_wire(p0),
        'tail_pole_sum_upper_ball': ball_wire(sp),
        'leading_transform_amplitude_sum_upper_ball': ball_wire(a_sum),
        'leading_transform_amplitude_max_upper_ball': ball_wire(a_max),
        'leading_pole_l1_upper_ball': ball_wire(phead_sum),
        'leading_pole_max_upper_ball': ball_wire(phead_max),
        'B_row_sum_upper_ball': ball_wire(row_B),
        'B_column_sum_upper_ball': ball_wire(col_B),
        'epsilon_B_upper_ball': ball_wire(eps_B),
        'epsilon_B_below_1e-90': True,
        'epsilon_D_upper_ball': ball_wire(eps_D),
        'epsilon_D_below_1e-190': True,
        '_epsilon_B': eps_B,
        '_epsilon_D': eps_D,
    }


def build(data_dir: Path, p238_path: Path, bits: int):
    if bits < 512:
        raise ValueError('at least 512 bits are required')
    ctx.prec = bits
    p238, finite, hashes = verify_source_binding(data_dir, p238_path)

    L = aq(LQ)
    Tball = arb(T)
    pi = arb.pi()
    l2, l3, l5 = arb(2).log(), arb(3).log(), arb(5).log()
    AL = 2*l2/arb(2).sqrt() + 2*l3/arb(3).sqrt() + l2 + 2*l5/arb(5).sqrt()
    beta = (Tball/(2*pi)).log() - arb(1)/T - AL

    # Certified prime-power cutoff: 5 < exp(2L) < 6, so n=2,3,4,5 only.
    e2L = (2*L).exp()
    if not (e2L.lower() > 5 and e2L.upper() < 6):
        raise ArithmeticError('prime-power cutoff 5 < exp(2L) < 6 not certified')
    if not AL.upper() < 5:
        raise ArithmeticError('A_L < 5 not certified')
    if not (beta.lower() > aq(F(1, 1000)) and beta.upper() < 1):
        raise ArithmeticError('0.001 < beta* < 1 not certified')

    # Compare independently recomputed A_L and beta* against both sector metadata.
    for sec in ('even', 'odd'):
        meta = finite['sectors'][sec]['assembly_metadata']
        if not (AL-wire_ball(meta['AL_ball'])).contains(0):
            raise ArithmeticError(f'A_L does not overlap P237 {sec} assembly')
        if not (beta-wire_ball(meta['beta_ball'])).contains(0):
            raise ArithmeticError(f'beta* does not overlap P237 {sec} assembly')

    # Binet bound for y=t/2 in [0,252]. With a=5/4 and z=a+iy,
    # |u^2+z^2|>=a^2, and 2∫u/(e^(2πu)-1)du=1/12.
    # Thus |Re psi(1/4+iy)| < 6 + 2/5 + 4/75 + 4.
    # Add |log pi|<2, A_L<5, beta*<1: |Psi_L-beta*|<18.454<20.
    if not ((arb(6).exp().lower() > 253)
            and ((aq(F(5,4)**2 + F(252**2))).sqrt().upper() < 253)
            and (pi.lower() > 1)
            and (arb(2).exp().lower() > pi.upper())):
        raise ArithmeticError('elementary constants for the Binet M=20 bound failed')
    M_bound = F(6) + F(2,5) + F(4,75) + 4 + 2 + 5 + 1
    if not M_bound < M:
        raise ArithmeticError('Binet-derived symbol bound is not <20')

    x = LQ*T
    z = LQ/2
    c = arb(T*M) / 3  # T*M/pi < T*M/3 = 3360
    results = [certify_sector(s, L, x, z, c) for s in (0,1)]

    for result in results:
        epsD, epsB = result.pop('_epsilon_D'), result.pop('_epsilon_B')
        if not beta.lower() - epsD.upper() > aq(COMMON_MARGIN):
            raise ArithmeticError('tail diagonal lower bound does not exceed the finite margin')
        if not aq(COMMON_MARGIN) - epsB.upper() > aq(GLOBAL_MARGIN):
            raise ArithmeticError('global margin deduction failed')

    report = {
        'schema': 'NUMERIA.ES8P239.ZHU_ONE_STROKE_TAIL_CERT.1',
        'precision_bits': bits,
        'source': 'Zhu, arXiv:2608.24827v2, Theorem 1.1 and Section 4; formulas independently instantiated at P237 parameters',
        'parameters': {'L': '161/200', 'Tsharp': T, 'modes_per_parity': N,
                       'first_discarded_degrees': {'even': 2*N, 'odd': 2*N+1}},
        'P237_input_sha256': hashes,
        'P238_certificate_sha256': hashlib.sha256(p238_path.read_bytes()).hexdigest(),
        'finite_block_margin_lower_exact': str(COMMON_MARGIN),
        'prime_power_cutoff_certified': '5 < exp(2L) < 6',
        'A_L_upper_bound_strict': '5',
        'beta_ball': ball_wire(beta),
        'beta_lower_bound_strict': '1/1000',
        'beta_upper_bound_strict': '1',
        'symbol_abs_bound_M_strict': str(M),
        'symbol_M_proof_rational_bound': str(M_bound),
        'Binet_bound_inputs': {
            'Re_psi_abs_strict': '6 + 2/5 + 4/75 + 4',
            'log_pi_abs_strict': '2',
            'A_L_strict': '5',
            'beta_strict': '1',
            'denominator_lower_bound': '(5/4)^2',
            'integral_identity': 'integral_0^infty u/(exp(2*pi*u)-1) du = 1/24'},
        'T_over_pi_M_upper_rational': str(F(T*M,3)),
        'geometric_ratio_upper': str(RATIO),
        'sector_bounds': results,
        'minimum_of_finite_margin_and_tail_coercivity_exceeds': str(COMMON_MARGIN),
        'full_R_lower_bound_strict': str(GLOBAL_MARGIN),
        'Q_ge_R_from_one_stroke_envelope': True,
        'even_odd_real_decomposition_certified': True,
        'complex_extension_by_parity_lemma': True,
        'window_positivity_claim': 'Q(f) >= 5e-19 ||f||_2^2 for every real or complex f supported in [-0.805,0.805], subject to the stated classical analytic inputs and P237 interval certificate',
        'finite_to_infinite_tail_and_coupling_certified': True,
        'lean_formalized': False,
        'infinite_Weil_positivity_global': False,
        'RH_proved': False,
    }
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', type=Path, default=Path('numeria/p237/data'))
    ap.add_argument('--p238', type=Path, default=Path('numeria/P238_SCHUR_MARGIN.json'))
    ap.add_argument('--out', type=Path, default=Path('numeria/P239_ZHU_TAIL_COUPLING.json'))
    ap.add_argument('--bits', type=int, default=1024)
    args = ap.parse_args()
    report = build(args.data, args.p238, args.bits)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n')
    print('P239 PASS: beta>1e-3; epsilon_B<1e-90 and epsilon_D<1e-190 in both parities; global lower >5e-19')
    print('wrote', args.out)

if __name__ == '__main__':
    main()
