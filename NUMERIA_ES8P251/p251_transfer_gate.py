"""P251: proven obstruction to naive Weyl transfer T=504 -> T=309.

This is NOT a proof that the reduced matrix has a negative eigenvalue.
The analytic Binet lower bound and the fixed-beta convention are inherited
from P249/P250; all interval sign decisions use P245B directed integers.
"""
import hashlib
import json
from fractions import Fraction as F

from NUMERIA_ES8P250.adaptive_signed_tail import load_p245b, run

L = F(161, 200)
OLD_MARGIN = F(1, 10**18)
PINNED_P250_SHA256 = '2bde93ccce25b8d2a1c6e59b49bcce43475fb2c7bc83366edaa1a19b7d4efe38'


def need(ok, explanation):
    if not ok:
        raise ArithmeticError(explanation)


def certify():
    p250 = run()  # Recompute the signed-tail certificate, not just check archived flags.
    wire = json.dumps(p250, indent=2, ensure_ascii=False) + '\n'
    need(hashlib.sha256(wire.encode()).hexdigest() == PINNED_P250_SHA256,
         'P250 certificate bytes changed')
    g = load_p245b()
    need(g.PREC == 256 and g.PI.u < 4 * g.S, 'expected directed arithmetic and pi < 4')
    delta = F(int(p250['min_lower_dyadic_numerator']), g.S)
    need(delta > 0 and p250['count_cells'] == 230, 'signed tail certificate failed')

    # Four cells partition [309,310]. A 1-Lipschitz sine estimate bounds
    # sin(L*t) by its centre upper plus L/8 on each quarter-cell.
    sine_bounds = []
    for k in range(4):
        t_centre = F(309) + F(2*k + 1, 8)
        sine, _ = g.sin_cos(g.iv(L * t_centre))
        upper = F(sine.u, g.S) + L/8
        need(upper < -F(1, 2), 'sine magnitude not certified')
        sine_bounds.append(str(upper))

    # For 309<=t<=310: |u_0(t)|² > 2L/(4*250²), since Lt<250.
    # H00 = (1/pi) int_(309)^(504) (Phi-beta)|u0|² dt.
    # P250+Binet give Phi-beta >= delta, and 1/pi>1/4.
    h00_lower = delta * F(1, 4) * (2*L) * F(1, 4) / (250**2)
    need(h00_lower > F(288, 10**12) and h00_lower > OLD_MARGIN,
         'nontransfer lower bound too small')
    return {
        'schema': 'NUMERIA.ES8P251.NONTRANSFER_GATE.1',
        'p250_replayed_sha256': PINNED_P250_SHA256,
        'interval_delta_lower_exact': str(delta),
        'four_sine_cell_upper_bounds_exact': sine_bounds,
        'H00_lower_exact': str(h00_lower),
        'old_margin_exact': str(OLD_MARGIN),
        'H00_to_old_margin_ratio_exact': str(h00_lower/OLD_MARGIN),
        'conclusion': '||H|| >= H00 > old margin; naive Weyl m-||H|| cannot certify T=309',
        'reduced_259_mode_matrix_certified': False,
        'reduced_matrix_indefinite_proved': False,
        'RH_proved': False,
    }


if __name__ == '__main__':
    print(json.dumps(certify(), indent=2, ensure_ascii=False))
