"""Interval enclosure audit for a finite NUMERIA-RIEMANN Weil matrix (v2).

Mudanças em relação à versão PRONTO_GITHUB (mesma matemática, enclosures mais justos):

1. Integral arquimediana: por padrão usa ``acb.integral`` (integração adaptativa
   certificada do Arb) em vez de somas de intervalo em passo fixo. A soma por
   células sofre do problema de dependência (cada fator é enclosed
   separadamente) e dava raios de ordem 0.1-1 nas entradas diagonais; a
   integração adaptativa dá raios ~1e-12. O caminho antigo continua disponível
   com ``--archimedean rangesum`` para comparação cruzada.
2. Omega(r) = Re psi(1/4 + i r/2) - log(pi) é escrita como função analítica de r,
   (psi(1/4+ir/2) + psi(1/4-ir/2))/2 - log(pi), avaliada com digamma do Arb.
3. Entradas com i-j ímpar são zero estrutural (Re(F_i conj F_j) = 0): ficam
   exatamente 0, inclusive na cauda (antes recebiam raio de cauda à toa).
4. Certificação: além do LDL intervalar, há um teste tipo Rump
   (M - ||rad||_inf * I definida positiva => toda matriz no enclosure é SPD),
   que não sofre do efeito de empacotamento do LDL intervalar.
5. Relatório: bolas salvas com ponto médio e raio superior em decimal (sem o
   arredondamento de 3 dígitos), limites superiores dos pivôs e o índice onde
   a fatoração parou. Autovalores do ponto médio entram só como diagnóstico
   NÃO rigoroso, claramente rotulado.

Escopo: certifica no máximo a matriz finita restrita para a base/parâmetros
dados. Não é prova do critério de Weil global nem da Hipótese de Riemann.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import flint
from flint import acb, arb, ctx

try:
    from .exact_mesh import exact_cell_count, validate_decimal
except ImportError:  # standalone audit execution
    from exact_mesh import exact_cell_count, validate_decimal

SCHEMA_VERSION = 3
CERTIFIED = "certified_positive_definite_finite_matrix"
INCONCLUSIVE = "inconclusive_interval_ldl"


# --------------------------------------------------------------------------
# utilidades
# --------------------------------------------------------------------------
def _assert_ball_radius_semantics():
    """Fail closed if installed python-flint ball radius semantics differ."""
    for radius in (arb("0.1"), arb("1e-30"), arb("123456789")):
        ball = arb(0, radius)
        if not (ball.contains(radius) and ball.contains(-radius)):
            raise ArithmeticError("Arb midpoint/radius constructor failed containment test")


def _a(value):
    """Parse a user decimal as an Arb ball, avoiding binary-float constants."""
    return arb(validate_decimal(value))


def _digits(precision):
    return int(precision * 0.30103) + 8


def _exact_fraction(x):
    """Valor exato de um arb SEM raio (ponto) como Fraction."""
    from fractions import Fraction

    man, exp = x.man_exp()
    man, exp = int(man), int(exp)
    return Fraction(man) * (Fraction(2) ** exp)


def _decimal_directed(point, digits, direction):
    """Decimal com `digits` algarismos significativos, arredondado com direção.

    direction = "floor" (limite inferior, nunca acima do valor) ou "ceil"
    (limite superior, nunca abaixo). `point` é um arb exato (raio 0). Usa
    Fractions: o único arredondamento é o direcionado.
    """
    from decimal import Decimal
    from fractions import Fraction

    if not point.is_finite():
        raise ArithmeticError("não é possível exportar um enclosure não finito")
    v = _exact_fraction(point)
    if v == 0:
        return "0"
    man, exp = point.man_exp()
    e10 = int((int(man).bit_length() + int(exp)) * 0.30103) - digits - 2
    scaled = v / (Fraction(10) ** e10)
    q = scaled.__floor__() if direction == "floor" else scaled.__ceil__()
    return str(Decimal(f"{q}E{e10}"))  # parse exato; scaleb arredondaria ao contexto (28 dígitos)


def _ball_json(x, digits):
    """Serializa uma bola de forma rigorosa.

    - ``lower``/``upper``: extremos decimais arredondados PARA FORA (o intervalo
      exportado contém a bola original);
    - ``mid``/``rad``: ponto médio e raio exatos como [mantissa, expoente]
      (valor = mantissa * 2**expoente), sem perda.
    """
    if not x.is_finite():
        raise ArithmeticError("não é possível exportar um enclosure não finito")
    mid_man, mid_exp = x.mid().man_exp()
    rad_man, rad_exp = x.rad().man_exp()
    return {
        "lower": _decimal_directed(x.lower(), digits, "floor"),
        "upper": _decimal_directed(x.upper(), digits, "ceil"),
        "mid_exact": [str(int(mid_man)), str(int(mid_exp))],
        "rad_exact": [str(int(rad_man)), str(int(rad_exp))],
    }


def _lower_decimal(x, digits=40):
    return _decimal_directed(x.lower(), digits, "floor")


def _upper_decimal(x, digits=40):
    return _decimal_directed(x.upper(), digits, "ceil")


def _sinc(x):
    """Enclose sin(x)/x for real balls, including balls containing zero."""
    return x.sinc()


# --------------------------------------------------------------------------
# parte arquimediana
# --------------------------------------------------------------------------
def _omega_series(r, terms):
    """Enclose Re digamma(1/4+i*r/2)-log(pi) for r>=0 via série (caminho legado).

    Re psi(1/4+ir/2) = Re psi(1/2+ir) - log(2) - pi/(2 cosh(pi r)),
    Re psi(1/2+ir) = -gamma - 2 log(2) + sum_{n>=0} r^2/((n+1/2)((n+1/2)^2+r^2)).
    O resto da série fica entre a integral e a integral mais o primeiro termo
    omitido; com r ~ R é preciso terms >> R, senão o raio explode.
    """
    pi = arb.pi()
    r2 = r * r
    total = -arb.const_euler() - 3 * arb(2).log() - pi.log()
    for n in range(terms):
        y = arb(n) + arb("0.5")
        total += r2 / (y * (y * y + r2))
    y = arb(terms) + arb("0.5")
    integral_tail = arb("0.5") * (1 + r2 / (y * y)).log()
    first_omitted = r2 / (y * (y * y + r2))
    sum_tail = integral_tail + arb(first_omitted / 2, first_omitted / 2)
    q = (-pi * r).exp()
    sech_term = pi * q / (1 + q * q)
    return total + sum_tail - sech_term


def _omega_analytic(r):
    """Omega(r) como função analítica de r (faixa |Im r| < 1/2), r acb.

    Para r real: Re psi(1/4+ir/2) = (psi(1/4+ir/2) + psi(1/4-ir/2))/2.
    """
    z = acb(0, 1) * r / 2
    quarter = acb("0.25")
    return (quarter + z).digamma() / 2 + (quarter - z).digamma() / 2 - arb.pi().log()


def _fourier_amplitude(r, L, k):
    """Real amplitude in F_k(r)=i*exp(-i*k*pi/2)*amplitude (r arb ou acb)."""
    w = arb(k) * arb.pi() / (2 * L)
    return 2 * L.sqrt() * w / (w + r) * _sinc(L * (r - w))


def _phase(i, j):
    """Fase de Re(F_i conj F_j)/(amp_i amp_j): 0 se i-j ímpar, senão +-1 (inteiro)."""
    d = i - j
    if d % 2:
        return 0
    return -1 if (abs(d) // 2) % 2 else 1


def _tail_bound(L, N, R):
    """Raio da cauda ∫_R^∞ para cada entrada (i, j) não nula por estrutura.

    Para r >= R >= 2*w_max: |F_k(r)| <= 4 w_k/(sqrt(L) r^2) e
    |Omega(r)| <= 8 + log(3r) para r >= 1 (prova em TAIL_PROOF.md: Lemas 1-4 e Proposição 5).
    """
    pi = arb.pi()
    wmax = arb(N) * pi / (2 * L)
    if not (R.lower() >= (2 * wmax).upper()):
        raise ValueError("R must be rigorously >= 2*w_max for the tail bound")
    if not (R.lower() >= 1):
        raise ValueError("R must be rigorously >= 1 for elementary tail bound")
    tail_factor = (8 + (3 * R).log() + arb(1) / 3) / (3 * R ** 3)
    tails = [[arb(0) for _ in range(N)] for _ in range(N)]
    for i in range(N):
        wi = arb(i + 1) * pi / (2 * L)
        for j in range(i, N):
            if _phase(i + 1, j + 1) == 0:
                continue  # zero estrutural: a cauda também é exatamente 0
            wj = arb(j + 1) * pi / (2 * L)
            tails[i][j] = tails[j][i] = 16 * wi * wj / (pi * L) * tail_factor
    return tails


def _archimedean_matrix_adaptive(L, N, R, tol_bits):
    """(1/pi) ∫_0^R Re(F_i conj F_j) Omega dr por integração adaptativa certificada."""
    pi = arb.pi()
    tol = arb(2) ** (-tol_bits)
    G = [[arb(0) for _ in range(N)] for _ in range(N)]
    for i in range(1, N + 1):
        for j in range(i, N + 1):
            ph = _phase(i, j)
            if ph == 0:
                continue

            def integrand(r, analytic, i=i, j=j, ph=ph):
                return ph * _fourier_amplitude(r, L, i) * _fourier_amplitude(r, L, j) \
                    * _omega_analytic(r) / pi

            value = acb.integral(integrand, 0, R, rel_tol=tol, abs_tol=tol).real
            if not value.is_finite():
                raise ArithmeticError(f"integral não convergiu para a entrada ({i},{j})")
            G[i - 1][j - 1] = G[j - 1][i - 1] = value
    return G


def _archimedean_matrix_rangesum(L, N, R, step_input, R_input, terms):
    """Caminho legado: somas de intervalo por células (para comparação cruzada)."""
    count = exact_cell_count(R_input, step_input)
    width = R / count
    pi = arb.pi()
    G = [[arb(0) for _ in range(N)] for _ in range(N)]
    for cell in range(count):
        center = (arb(2 * cell + 1) * R) / (2 * count)
        r = center + arb(0, width / 2)
        omega = _omega_series(r, terms)
        amplitudes = [_fourier_amplitude(r, L, k) for k in range(1, N + 1)]
        for i in range(N):
            for j in range(i, N):
                phase = _phase(i + 1, j + 1)
                if phase == 0:
                    continue
                contribution = width * (phase * amplitudes[i] * amplitudes[j] * omega / pi)
                G[i][j] += contribution
                if i != j:
                    G[j][i] += contribution
    return G, count


# --------------------------------------------------------------------------
# parte aritmética (primos) e restrições
# --------------------------------------------------------------------------
def _prime_power_base(n):
    """Retorna p se n = p^k (k>=1), senão None; Lambda(n) = log p nesse caso."""
    m, p, d = n, None, 2
    while d * d <= m:
        if m % d == 0:
            p = d
            while m % d == 0:
                m //= d
            break
        d += 1
    if p is None:
        return n
    return p if m == 1 else None


def _overlap(a, L, i, j):
    """Exact analytic overlap enclosure, symmetrized as in the prototype."""
    pi = arb.pi()
    b = 2 * L - a

    def directed(k, m):
        wk = arb(k) * pi / (2 * L)
        wm = arb(m) * pi / (2 * L)
        first_phase = wk * a
        if k == m:
            first = b * first_phase.cos()
        else:
            delta = wk - wm
            first = ((delta * b + first_phase).sin() - first_phase.sin()) / delta
        total = wk + wm
        second = ((total * b + first_phase).sin() - first_phase.sin()) / total
        return (first - second) / (2 * L)

    return (directed(i, j) + directed(j, i)) / 2


def _prime_matrix(L, N):
    """Construct the prime-power matrix, with rigorous cutoff decisions."""
    twol = 2 * L
    P = [[arb(0) for _ in range(N)] for _ in range(N)]
    active = []
    n = 2
    while True:
        logn = arb(n).log()
        if logn.lower() > twol.upper():
            break
        if logn.upper() >= twol.lower():
            raise ArithmeticError(f"Cannot separate prime-power cutoff at n={n}; increase precision")
        p = _prime_power_base(n)
        if p is not None:
            active.append(n)
            coefficient = -2 * arb(p).log() / arb(n).sqrt()
            for i in range(1, N + 1):
                for j in range(i, N + 1):
                    value = coefficient * _overlap(logn, L, i, j)
                    P[i - 1][j - 1] += value
                    if i != j:
                        P[j - 1][i - 1] += value
        n += 1
    return P, active


def _constraint_nullspace(L, N):
    """Interval-enclosed basis B for the exact two pole constraints C*B=0."""
    pi = arb.pi()
    rootL = L.sqrt()
    plus, minus = [], []
    for k in range(1, N + 1):
        w = arb(k) * pi / (2 * L)
        parity = (-1) ** k
        denom = rootL * (w * w + arb("0.25"))
        plus.append(w * ((-L / 2).exp() - parity * (L / 2).exp()) / denom)
        minus.append(w * ((L / 2).exp() - parity * (-L / 2).exp()) / denom)
    c00, c01 = plus[0], plus[1]
    c10, c11 = minus[0], minus[1]
    determinant = c00 * c11 - c01 * c10
    if determinant.contains(0):
        raise ArithmeticError("2x2 pole-constraint pivot is not separated from zero")
    dimension = N - 2
    B = [[arb(0) for _ in range(dimension)] for _ in range(N)]
    for col, j in enumerate(range(2, N)):
        B[0][col] = (-plus[j] * c11 + c01 * minus[j]) / determinant
        B[1][col] = (-c00 * minus[j] + plus[j] * c10) / determinant
        B[j][col] = arb(1)
    return B


def _project(B, W):
    """A = B^T W B, com W simétrica (O(N^3) em vez de O(N^4))."""
    N, dim = len(B), len(B[0])
    WB = [[sum((W[k][e] * B[e][j] for e in range(N)), arb(0)) for j in range(dim)] for k in range(N)]
    A = [[arb(0) for _ in range(dim)] for _ in range(dim)]
    for i in range(dim):
        for j in range(i, dim):
            value = sum((B[k][i] * WB[k][j] for k in range(N)), arb(0))
            A[i][j] = A[j][i] = value
    return A


# --------------------------------------------------------------------------
# certificação
# --------------------------------------------------------------------------
def _ldl_pivots(A):
    """LDL^T em aritmética de bolas. Retorna (pivôs, certificado, índice_da_falha)."""
    n = len(A)
    lower = [[arb(0) for _ in range(n)] for _ in range(n)]
    diagonal = []
    for i in range(n):
        lower[i][i] = arb(1)
        pivot = A[i][i]
        for k in range(i):
            pivot -= lower[i][k] * lower[i][k] * diagonal[k]
        diagonal.append(pivot)
        if not (pivot.lower() > 0):
            return diagonal, False, i
        for j in range(i + 1, n):
            value = A[j][i]
            for k in range(i):
                value -= lower[j][k] * lower[i][k] * diagonal[k]
            lower[j][i] = value / pivot
    return diagonal, True, None


def _symmetrize(A):
    n = len(A)
    S = [[arb(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            S[i][j] = S[j][i] = A[i][j].union(A[j][i])
    return S


def _rump_certificate(A):
    """Teste de positividade definida para TODA matriz simétrica dentro do enclosure.

    Sejam M os pontos médios (exatos) e rho os raios. Para E simétrica com
    |E_ij| <= rho_ij vale ||E||_2 <= ||rho||_inf =: delta. Se M - delta*I é
    definida positiva (LDL com pivôs certificados > 0), então lambda_min(M+E) > 0.
    Retorna (certificado, delta, pivôs de M - delta*I, índice da falha).
    """
    S = _symmetrize(A)
    n = len(S)
    M = [[S[i][j].mid() for j in range(n)] for i in range(n)]
    rho = [[S[i][j].rad() for j in range(n)] for i in range(n)]
    delta = max((sum((rho[i][j] for j in range(n)), arb(0)) for i in range(n)),
                key=lambda x: float(x.upper())).upper()
    shifted = [[M[i][j] - (delta if i == j else 0) for j in range(n)] for i in range(n)]
    pivots, ok, fail = _ldl_pivots(shifted)
    return ok, delta, pivots, fail


def _midpoint_eigenvalues(A, digits):
    """DIAGNÓSTICO NÃO RIGOROSO: autovalores do ponto médio (mpmath)."""
    try:
        import mpmath as mp
    except ImportError:
        return None
    mp.mp.dps = digits
    n = len(A)
    S = _symmetrize(A)
    M = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            M[i, j] = mp.mpf(S[i][j].mid().str(digits, radius=False))
    evals = mp.eigsy(M, eigvals_only=True)
    return sorted(mp.nstr(e, 12) for e in evals)


# Apelidos dos nomes da versão anterior (testes e scripts externos).
_omega = _omega_series
_von_mangoldt_prime = _prime_power_base


def _interval_ldl_pivots(A):
    """Compatibilidade: retorna (pivôs, certificado) como na versão anterior."""
    pivots, ok, _ = _ldl_pivots(A)
    return pivots, ok


# --------------------------------------------------------------------------
# execução
# --------------------------------------------------------------------------
def run(L_value="0.45", N=6, R_value="160", step_value="0.02", terms=128,
        precision=192, archimedean="adaptive", integral_tol_bits=50):
    if N < 4 or precision < 64:
        raise ValueError("Require N>=4 and precision>=64 bits")
    if archimedean not in ("adaptive", "rangesum"):
        raise ValueError("archimedean must be 'adaptive' or 'rangesum'")
    if archimedean == "rangesum" and terms < 16:
        raise ValueError("rangesum requires terms>=16")
    old_precision = ctx.prec
    ctx.prec = precision
    try:
        L_in, R_in, step_in = (validate_decimal(v) for v in (L_value, R_value, step_value))
        L, R, step = _a(L_in), _a(R_in), _a(step_in)
        if not (L.lower() > 0 and R.lower() > 0 and step.lower() > 0):
            raise ValueError("L, R and step must be strictly positive")
        _assert_ball_radius_semantics()

        tails = _tail_bound(L, N, R)
        if archimedean == "adaptive":
            G = _archimedean_matrix_adaptive(L, N, R, integral_tol_bits)
            cells = None
        else:
            if not ((L * step / 2).upper() < arb("0.1")):
                raise ValueError("step too wide for sinc enclosure")
            G, cells = _archimedean_matrix_rangesum(L, N, R, step_in, R_in, terms)
        for i in range(N):
            for j in range(N):
                if tails[i][j].upper() > 0:
                    G[i][j] += arb(0, tails[i][j])

        P, active = _prime_matrix(L, N)
        W = [[G[i][j] + P[i][j] for j in range(N)] for i in range(N)]
        B = _constraint_nullspace(L, N)
        projected = _project(B, W)

        pivots, ldl_ok, ldl_fail = _ldl_pivots(projected)
        rump_ok, delta, rump_pivots, rump_fail = _rump_certificate(projected)
        certified = ldl_ok or rump_ok
        digits = _digits(precision)
        dim = N - 2

        report = {
            "schema_version": SCHEMA_VERSION,
            "project": "NUMERIA-RIEMANN",
            "method": ("Arb adaptive certified quadrature" if archimedean == "adaptive"
                       else "Arb interval range-sum quadrature")
                      + " + analytic tail + interval LDL + Rump-type margin test",
            "status": CERTIFIED if certified else INCONCLUSIVE,
            "claim_scope": "finite constrained sine-basis matrix only; not global Weil positivity or RH",
            "environment": {"python_flint": flint.__version__},
            "parameters": {
                "L": L_in, "N": N, "R": R_in, "archimedean": archimedean,
                "integral_tol_bits": integral_tol_bits if archimedean == "adaptive" else None,
                "requested_step": step_in if archimedean == "rangesum" else None,
                "quadrature_cells": cells,
                "digamma_series_terms": terms if archimedean == "rangesum" else None,
                "arb_precision_bits": precision,
            },
            "active_prime_powers": active,
            "projected_dimension": dim,
            "certificates": {
                "interval_ldl": {"passed": ldl_ok, "failed_at_pivot": None if ldl_fail is None else ldl_fail + 1},
                "rump_margin": {"passed": rump_ok, "delta_upper": _upper_decimal(delta, 12),
                                "failed_at_pivot": None if rump_fail is None else rump_fail + 1},
            },
            "interval_ldl_pivots": [str(p) for p in pivots],  # só exibição (arredonda); rigoroso: *_enclosures e *_bounds
            "interval_ldl_pivot_enclosures": [_ball_json(p, digits) for p in pivots],
            "interval_ldl_pivot_lower_bounds": [_lower_decimal(p) for p in pivots],
            "interval_ldl_pivot_upper_bounds": [_upper_decimal(p) for p in pivots],
            "rump_shifted_pivot_lower_bounds": [_lower_decimal(p) for p in rump_pivots],
            "archimedean_entry_enclosures": [[_ball_json(G[i][j], digits) for j in range(N)] for i in range(N)],
            "projected_matrix_enclosures": [[_ball_json(projected[i][j], digits) for j in range(dim)] for i in range(dim)],
            "diagnostic_only_midpoint_eigenvalues": {
                "rigorous": False,
                "values": _midpoint_eigenvalues(projected, 30),
            },
            "display_only_keys": ["interval_ldl_pivots"],
            "limitations": [
                "The archimedean integral over [0,R] is enclosed by Arb certified integration (or range sums in legacy mode).",
                "The infinite archimedean tail is enclosed using the documented analytic bound.",
                "The pole constraints and prime-overlap integrals are enclosed with Arb arithmetic.",
                "A positive certificate covers only this finite-dimensional constrained matrix.",
                "A failed test is inconclusive and is not evidence of a negative eigenvalue.",
                "Midpoint eigenvalues are a non-rigorous diagnostic and certify nothing.",
                "Exported decimal bounds are rounded outward (lower: floor, upper: ceil); keys listed in display_only_keys are not.",
                "No result here proves the global Weil criterion or the Riemann Hypothesis.",
            ],
        }
        return report
    finally:
        ctx.prec = old_precision


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--L", type=str, default="0.45")
    parser.add_argument("--N", type=int, default=6)
    parser.add_argument("--R", type=str, default="160")
    parser.add_argument("--archimedean", choices=("adaptive", "rangesum"), default="adaptive")
    parser.add_argument("--integral-tol-bits", type=int, default=50,
                        help="tolerância absoluta/relativa 2^-bits do integrador adaptativo")
    parser.add_argument("--step", type=str, default="0.02", help="só para --archimedean rangesum")
    parser.add_argument("--series-terms", type=int, default=128, help="só para --archimedean rangesum")
    parser.add_argument("--precision", type=int, default=192)
    parser.add_argument("--out", type=Path, default=Path("weil-interval-audit.json"))
    parser.add_argument("--require-positive", action="store_true",
                        help="exit nonzero unless a certificate passes")
    args = parser.parse_args(argv)
    report = run(args.L, args.N, args.R, args.step, args.series_terms, args.precision,
                 args.archimedean, args.integral_tol_bits)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in (
        "status", "parameters", "active_prime_powers", "projected_dimension", "certificates",
        "interval_ldl_pivot_lower_bounds", "diagnostic_only_midpoint_eigenvalues", "limitations")}, indent=2))
    if args.require_positive and report["status"] != CERTIFIED:
        raise SystemExit("no certificate passed for the reference finite matrix")


if __name__ == "__main__":
    main()
