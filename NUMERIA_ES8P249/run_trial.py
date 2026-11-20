"""Generate the finite ES8-P249 experiment and its reproducibility record."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

from flint import arb, ctx
import numpy as np

from .interval_audit import run_audit
from .interval_fft import toeplitz_matvec_arb
from .legendre_tail import reduced_tail_report
from .p236_replay import panel_cell_lags, replay_p236_archive, replay_panel
from .prime_operator import (
    P236_HALF_WIDTH,
    P236_PRIME_POWERS,
    ToeplitzFFT,
    dyadic_prime_components,
    fourier_cell_lags,
    prime_symbol,
)


P247_BASE = '6adabcca7cc09f60df740934114d49966893d3d7'
PILOT_DIR = Path(__file__).resolve().parents[1] / 'numeria/p245/package/pilot_sources'
PANELS = {0: '0', 503: '503_900', 1007: '1007_1152'}
CELL_COUNT = 128
REPLAY_PRECISION = 160
ARCHIVE_PRECISION = 900
MAX_LEVEL = 5


def _dyadic(point):
    mantissa, exponent = point.man_exp()
    return [str(mantissa), int(exponent)]


def _interval(ball):
    return [_dyadic(ball.lower()), _dyadic(ball.upper())]


def _panel_fft_report(panel, cells):
    roots, gl_weights = np.polynomial.legendre.leggauss(60)
    frequencies = (2 * panel + 1) / 4 + roots / 4
    weights = gl_weights / (4 * np.pi)
    vector = np.array([((i * i + 5 * i) % 19 - 9) / 8 for i in range(cells)])
    lags = fourier_cell_lags(cells, P236_HALF_WIDTH, frequencies, weights * prime_symbol(frequencies))
    computed = ToeplitzFFT(lags).matvec(vector)
    exact_lags = panel_cell_lags(panel, cells, bits=REPLAY_PRECISION)
    enclosure = toeplitz_matvec_arb(exact_lags, vector.tolist(), precision_bits=REPLAY_PRECISION)

    terms = dyadic_prime_components(frequencies, MAX_LEVEL)
    components = [
        ToeplitzFFT(fourier_cell_lags(cells, P236_HALF_WIDTH, frequencies, weights * term)).matvec(vector)
        for term in terms
    ]
    reconstruction = np.sum(components, axis=0)

    with ctx.workprec(REPLAY_PRECISION):
        errors = [(arb(float(value)) - exact.real).abs_upper() for value, exact in zip(computed, enclosure)]
        curvature = sum(
            (2 * arb(p).log() * arb(q).log() ** 2 / arb(q).sqrt() for q, p in P236_PRIME_POWERS),
            arb(0),
        )
        # Sum of positive GL weights divided by pi is 1/(2*pi) per panel.
        half_width = arb(P236_HALF_WIDTH.numerator) / P236_HALF_WIDTH.denominator
        residual_norm_bound = (half_width * curvature / (arb.pi() * 8 * 2 ** (2 * MAX_LEVEL))).upper()
        return {
            'input_vector_dyadic': [_dyadic(arb(float(v))) for v in vector],
            'binary64_max_abs_error_upper_dyadic': _dyadic(max(errors)),
            'analytic_operator_norm_bound_if_residual_omitted': _dyadic(residual_norm_bound),
            'levels_through': MAX_LEVEL,
            'component_action_inf_norms_binary64': [float(np.linalg.norm(y, ord=np.inf)) for y in components],
            'reconstruction_max_abs_error_observed': float(np.max(np.abs(reconstruction - computed))),
        }


def _replay_report(panel, degrees, archive):
    replay = replay_panel(panel, CELL_COUNT, degrees, bits=REPLAY_PRECISION)
    comparison = replay_p236_archive(archive, bits=ARCHIVE_PRECISION)
    samples = comparison.pop('sample_checks')
    if len(samples) != 534 or any(not item['archive_contained_in_replay'] for item in samples):
        raise ArithmeticError(f'archived P236 interval not contained in replay: {archive.name}')
    mantissa, exponent = map(int, comparison['max_replay_width'])
    if not mantissa * Fraction(2) ** exponent < Fraction(2) ** -320:
        raise ArithmeticError(f'P236 replay interval too wide: {archive.name}')
    comparison['archive_intervals_contained'] = len(samples)
    intervals = []
    with ctx.workprec(REPLAY_PRECISION):
        for entry in replay['entries']:
            residual_upper = entry['grid_residual'].abs_upper()
            # A zero bound is exact for the constant polynomial in both slots.
            if entry['grid_bound'] != 0 and not residual_upper < entry['grid_bound']:
                raise ArithmeticError(f'grid bridge bound did not close: {entry["degrees"]}')
            intervals.append({
                'degrees': entry['degrees'],
                'legendre_prime': _interval(entry['legendre_prime']),
                'cell_projected_prime': _interval(entry['cell_projected_prime']),
                'grid_residual': _interval(entry['grid_residual']),
                'analytic_grid_error_upper': _dyadic(entry['grid_bound']),
                'p236_full_quadrature': _interval(entry['p236_full_quadrature']),
            })
    return {'source_archive': archive.name, 'comparison': comparison, 'entries': intervals}


def run_trial():
    audit = run_audit()
    tail_cutoff = audit['p236_tail_start']['certified_range'][0]
    panels = []
    for panel, suffix in PANELS.items():
        replays = [
            _replay_report(panel, tuple(2 * n + offset for n in range(3)), PILOT_DIR / f'P236_{sector}_{suffix}.json')
            for offset, sector in enumerate(('even', 'odd'))
        ]
        panels.append({'panel': panel, 'fft': _panel_fft_report(panel, CELL_COUNT), 'replays': replays})
    return {
        'schema': 'NUMERIA.ES8P249.FINITE_PRIME_OPERATOR.1',
        'source_base': P247_BASE,
        'half_width_exact': str(P236_HALF_WIDTH),
        'prime_powers': P236_PRIME_POWERS,
        'panel_cells': CELL_COUNT,
        'panel_precision_bits': REPLAY_PRECISION,
        'audit': audit,
        'reduced_cutoff_modal_tail': reduced_tail_report(tail_cutoff),
        'p236_panels': panels,
        'scope': {
            'signed_tail_extension_uses_archimedean_lower_bound': True,
            'p237_low_block_recertified_for_earlier_cutoff': False,
            'full_continuous_weil_operator_certified': False,
            'cofinal_positivity_established': False,
            'riemann_hypothesis_proved': False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    output = json.dumps(run_trial(), indent=2, allow_nan=False) + '\n'
    if args.out is None:
        print(output, end='')
    else:
        args.out.write_text(output, encoding='utf-8')


if __name__ == '__main__':
    main()
