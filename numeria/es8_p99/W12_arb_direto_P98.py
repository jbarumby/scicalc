#!/usr/bin/env python3
"""ES8-P99: robustified direct P92 W12 engine (inherits P98 experimental reformulation).

Usa python-flint Acb.integral para integrais validadas e calcula
  W12 = P12 - S12 - (gamma+log(4*pi))*R(delta) - I_delta
na normalização b(x)=exp(-1/(1-(x/h)^2)), h=2/5, delta=log(2).
A correlação R(s) é integrada diretamente, sem carregar intervalos de outra fase.

O script não lê nem importa arquivos P89/P91/P92. Só declara ENCLOSED se todas
as integrais Arb retornarem bolas finitas e se as hipóteses de analiticidade
por painel passarem pelas verificações implementadas. Qualquer falha retorna
INCONCLUSIVE. A aritmética é validada pelo Arb; a decomposição em painéis e as
regras de extensão do integrando devem ser auditadas matematicamente antes de
usar o resultado como certificado publicado.

Dependência: python -m pip install python-flint
Exemplo exploratório: python W12_arb_direto_P92.py --prec 128 --tol 1e-12
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import platform
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Callable

try:
    from flint import acb, arb, ctx, fmpq
except ImportError as exc:
    raise SystemExit("Instale python-flint: python -m pip install python-flint") from exc


def Q(p: int, q: int = 1) -> arb:
    return arb(fmpq(p, q))


def finite(x) -> bool:
    return bool(x.is_finite())


def acb_nan() -> acb:
    return acb(arb.nan(), arb.nan())


def real_bounds(z: acb):
    return z.real.lower(), z.real.upper()


def imag_is_zero(z: acb) -> bool:
    return bool(z.imag.is_zero())


@dataclass
class W12P92:
    h: arb
    tol: arb
    eps: arb
    eval_limit: int
    depth_limit: int
    calls_R: int = 0
    integral_count: int = 0
    outer_analytic_denials: int = 0
    inner_support_trims: int = 0
    outer_panel_count: int = 0
    derivative_intersections: int = 0
    eps_denominator: int = 2048

    @classmethod
    def make(cls, tol_text: str, eval_limit: int, depth_limit: int, eps_denominator: int = 2048):
        h = Q(2, 5)
        tol = arb(tol_text)
        if eps_denominator < 16 or eps_denominator & (eps_denominator - 1):
            raise ValueError("eps_denominator must be a power of two >= 16")
        eps = Q(1, eps_denominator)
        return cls(h, tol, eps, eval_limit, depth_limit, eps_denominator=eps_denominator)

    def bump(self, z: acb, analytic: bool) -> acb:
        """Bump P92, with real-axis enclosures for nonanalytic requests.

        The genuine integration path is real. For analytic=False, only z.real
        is relevant; the imaginary coordinate of a complex *query* does not
        represent a complex value of the non-analytic compactly supported bump.
        All real values lie in [0, exp(-1)], including support-crossing cells.
        Analytic=True is accepted ONLY for disks wholly in one analytic branch,
        i.e. strictly inside support or strictly outside it.
        """
        rl, ru = real_bounds(z)
        h_lo, h_hi = self.h.lower(), self.h.upper()
        entirely_inside = bool(rl > -h_lo and ru < h_lo)
        entirely_right = bool(rl >= h_hi)
        entirely_left = bool(ru <= -h_hi)
        if entirely_right or entirely_left:
            return acb(0)
        if entirely_inside:
            x = z if analytic else acb(z.real)
            den = 1 - (x / self.h) ** 2
            if den.contains(0):
                return acb_nan() if analytic else self.bump_range()
            return (-1 / den).exp()
        return acb_nan() if analytic else self.bump_range()

    @staticmethod
    def bump_range() -> acb:
        # b >= 0; b <= e^-1 for all real x. Real-only interval enclosure.
        e_inv = 1 / arb.const_e()
        return acb(e_inv / 2 + arb(0, (e_inv / 2).upper()))

    def integrate(self, f: Callable, a, b, label: str) -> acb:
        self.integral_count += 1
        val = acb.integral(
            f, a, b,
            abs_tol=self.tol,
            rel_tol=self.tol,
            eval_limit=self.eval_limit,
            depth_limit=self.depth_limit,
        )
        if not finite(val):
            raise ArithmeticError(f"Acb.integral retornou bola não finita em {label}: {val}")
        return val

    def moment_A(self) -> acb:
        def f(y, analytic):
            by = self.bump(y, analytic)
            return (y / 2).exp() * by
        return self.integrate(f, -self.h, self.h, "momento A")

    def correlation(self, s: acb) -> acb:
        """Integrate over a provable superset of the overlap of two supports.

        For s in [slo,shi], y must lie in
        [-h,h] intersect [-h-shi,h-slo]. This reduces the inner integration
        domain without relying on prior P89/P95 values or float point geometry.
        """
        self.calls_R += 1
        sl, su = real_bounds(s)
        h = self.h
        if sl >= (2*h).upper() or su <= -(2*h).lower():
            return acb(0)
        low = max((-h).lower(), (-h-su).lower())
        high = min(h.upper(), (h-sl).upper())
        if not bool(low < high):
            return acb(0)
        if bool(low > (-h).lower() or high < h.upper()):
            self.inner_support_trims += 1

        # The integration argument y is real. For analytic=False, do NOT
        # claim holomorphy of the dependence s -> R(s).
        def f(y, analytic):
            if not analytic:
                y = acb(y.real)
            return self.bump(y+s, analytic) * self.bump(y, analytic)
        return self.integrate(f, low, high, f"R(s={s})")

    def correlation_real_log(self, n: int) -> acb:
        return self.correlation(acb(arb(n).log()))

    def R_delta(self, delta: acb) -> acb:
        return self.correlation(delta)

    def small_kernel_integral(self) -> acb:
        """Integral de K(t)=(e^(t/2)-1)/sinh(t), forma sem singularidade."""
        def kfun(t, analytic):
            y = (-t / 2).exp()
            return acb(2) * y / ((1 + y) * (1 + y*y))
        return self.integrate(kfun, arb(0), self.eps, "kernel K em [0,eps]")

    @staticmethod
    def frac_endpoint(x: arb, upper: bool = False) -> Fraction:
        endpoint = x.upper() if upper else x.lower()
        q = endpoint.fmpq()
        return Fraction(int(q.p), int(q.q))

    @staticmethod
    def ball_from_rational_interval(lo: Fraction, hi: Fraction) -> arb:
        if lo > hi:
            raise ArithmeticError("empty intersection")
        mid = (lo+hi)/2
        rad = (hi-lo)/2
        return arb(fmpq(mid.numerator, mid.denominator),
                   fmpq(rad.numerator, rad.denominator))

    @classmethod
    def intersect_real_enclosures(cls, x: arb, y: arb) -> arb:
        lo = max(cls.frac_endpoint(x), cls.frac_endpoint(y))
        hi = min(cls.frac_endpoint(x, True), cls.frac_endpoint(y, True))
        return cls.ball_from_rational_interval(lo, hi)

    def outer_integrand(self, delta: acb, rdelta: acb, panel_lo: arb, panel_hi: arb):
        def f(t, analytic):
            # Holomorphy of correlation R(delta +/- t) on COMPLEX disks has
            # not been proven. Never supply a potentially false ellipse bound.
            if analytic:
                self.outer_analytic_denials += 1
                return acb_nan()
            t = acb(t.real)
            rplus = self.correlation(t + delta)
            rminus = self.correlation(t - delta)
            F = (rplus + rminus) / 2
            den = t.sinh()
            if den.contains(0):
                return acb_nan()
            # Independent C^2 bound for the real difference; intersect only
            # enclosures that both have a mathematical justification.
            lo_t = self.frac_endpoint(t.real)
            hi_t = self.frac_endpoint(t.real, True)
            t_sup = max(abs(lo_t), abs(hi_t))
            radius = t_sup*t_sup/Fraction(2,5)
            radius_ball = arb(0, fmpq(radius.numerator, radius.denominator))
            if not (F.imag.contains(0) and rdelta.imag.contains(0)):
                raise ArithmeticError("correlations unexpectedly nonreal")
            d_raw = F.real-rdelta.real
            d = self.intersect_real_enclosures(d_raw, radius_ball)
            self.derivative_intersections += 1
            y = (-t/2).exp()
            K = acb(2)*y/((1+y)*(1+y*y))
            return acb(rdelta.real*K.real + (t/2).exp().real*d/den.real)
        return f

    @staticmethod
    def geometric_panels(eps: arb, stop: arb) -> list[arb]:
        """Dyadic near-zero panels, so interval divisions avoid giant widths."""
        pts = [eps]
        while bool((pts[-1]*2).upper() < stop.lower()):
            pts.append(pts[-1]*2)
            if len(pts) > 64:
                raise RuntimeError("excessive number of dyadic panels")
        pts.append(stop)
        assert all(bool(pts[i].upper() < pts[i+1].lower())
                   for i in range(len(pts)-1))
        return pts

    def compute(self) -> dict:
        h, tol, eps = self.h, self.tol, self.eps
        delta = acb(arb(2).log())
        T = 2*h + delta.real
        a = self.moment_A()
        pole = a*a*(arb(2).sqrt()+1/arb(2).sqrt())
        rdelta = self.R_delta(delta)

        # S12: k=2,3,4 (inclui 4=2^2), com Lambda(k)=log(prime(k)).
        prime = acb(0)
        for k,p in ((2,2),(3,3),(4,2)):
            logk = arb(k).log()
            coeff = arb(p).log()/arb(k).sqrt()
            prime += coeff*(self.correlation(acb(logk)+delta)
                            + self.correlation(acb(logk)-delta))

        gamma_log = arb.const_euler() + (4*arb.pi()).log()

        # P96 lemma: ||b'||_infty <= 8/(e^2 h), ||b'||_1 = 2/e,
        # ||R''||_infty <= 16/(e^3 h) < 2/h. Symmetric second difference gives
        # |F_delta(t)-R(delta)| <= t^2/h. For 0<t<eps<1, exp(t/2)/sinh(t)
        # <= 2/t. Hence the omitted corrective integral is <= eps^2/h.
        small_main = rdelta * self.small_kernel_integral()
        correction_bound = eps*eps/h
        small_correction = acb(arb(0, correction_bound), arb(0))

        # Divisão nos pontos de contato: 2h-delta, t=delta (s=0), T=2h+delta.
        p1 = 2*h-delta.real
        p2 = delta.real
        points = self.geometric_panels(eps, p1) + [p2, T]
        if not all(points[i].upper() < points[i+1].lower()
                   for i in range(len(points)-1)):
            raise ArithmeticError("painéis não ficaram ordenados/disjuntos em Arb")
        I_mid = acb(0)
        self.outer_panel_count = len(points)-1
        for lo,hi in zip(points[:-1],points[1:]):
            I_mid += self.integrate(self.outer_integrand(delta,rdelta,lo,hi),
                                    lo,hi,f"I painel [{lo},{hi}]")

        # Para t>=T, F_delta(t)=0; cauda exata R(delta)*log(tanh(T/2)).
        tail = rdelta * ((T/2).tanh()).log()
        I = small_main + small_correction + I_mid + tail
        W = pole - prime - gamma_log*rdelta - I
        for name,val in (("A",a),("R_delta",rdelta),("pole",pole),
                         ("prime_sum",prime),("I",I),("W12_P92",W)):
            if not finite(val):
                raise ArithmeticError(f"resultado não finito em {name}: {val}")
        sign_certified = bool(W.imag.contains(0) and W.real.upper() < 0)
        status = ("DIRECT_P92_W12_STRICTLY_NEGATIVE_ENCLOSURE"
                  if sign_certified else "DIRECT_P92_FINITE_ENCLOSURE_ZERO_NOT_EXCLUDED")
        def exact_bounds(z):
            def fstr(q):
                qq = q.fmpq()
                return f"{int(qq.p)}/{int(qq.q)}"
            return {"lo": fstr(z.real.lower()),
                    "hi": fstr(z.real.upper()),
                    "im_lo": fstr(z.imag.lower()),
                    "im_hi": fstr(z.imag.upper())}
        return {
            "status": status,
            "version": "ES8-P99_DIRECT_P98_METHOD",
            "derivative_intersections": self.derivative_intersections,
            "eps_denominator": self.eps_denominator,
            "precision_bits": ctx.prec,
            "rigor_scope": "interval arithmetic, strict outer analytic denial; full analytic justification of nested parameter remains subject to independent review",
            "exact_real_bounds": {name:exact_bounds(val) for name,val in
                (("moment_A",a),("R_delta",rdelta), ("pole",pole),
                 ("prime_sum",prime),("I_mid",I_mid), ("I_tail",tail),
                 ("I_total",I),("W12",W))},
            "outer_analytic_denials": self.outer_analytic_denials,
            "inner_support_trims": self.inner_support_trims,
            "outer_panel_count": self.outer_panel_count,
            "min_support_span": "computed from exact Arb endpoint bounds",

            "normalization": "b(x)=exp(-1/(1-(x/h)^2)); h=2/5",
            "centers": [1,2], "delta": str(delta.real), "T": str(T),
            "enclosures_arb": {"A":str(a),"R_delta":str(rdelta),
                               "pole":str(pole),"prime_sum":str(prime),
                               "small_I_main":str(small_main),
                               "small_I_correction_bound":str(correction_bound),
                               "I_mid":str(I_mid),"I_tail":str(tail),
                               "I_total":str(I),"W12":str(W)},
            "W12_real_lower": str(W.real.lower()),
            "W12_real_upper": str(W.real.upper()),
            "is_strictly_negative": bool(W.real.upper() < 0),
            "integration_calls": self.integral_count,
            "correlation_calls": self.calls_R,
            "quadrature": "python-flint Acb.integral, interval ball arithmetic; nested R(s) integrated directly",
            "provenance": "No P89/P91/P92 JSON/archive is read by this script.",
            "warning": "Outer analytic=True queries are always denied. Interval calculation is separate from full mathematical proof of the Bombieri normalization and small-t remainder. Zero sign certification unless strictly bounded. Not a RH proof."
        }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prec',type=int,default=256)
    ap.add_argument('--eps-den',type=int,default=2048)
    ap.add_argument('--tol',default='1e-13')
    ap.add_argument('--eval-limit',type=int,default=100000)
    ap.add_argument('--depth-limit',type=int,default=30)
    ap.add_argument('--output',type=Path,default=Path('W12_ARB_P92_DIRECT.json'))
    args=ap.parse_args()
    ctx.prec=args.prec
    start=time.time()
    if args.prec < 128 or args.eval_limit < 100 or args.depth_limit < 5:
        raise SystemExit("Preconditions: prec>=128, eval_limit>=100, depth_limit>=5")
    try:
        calc=W12P92.make(args.tol,args.eval_limit,args.depth_limit,eps_denominator=args.eps_den)
        result=calc.compute()
        result['elapsed_seconds']=time.time()-start
    except Exception as exc:
        result={"status":"INCONCLUSIVE","reason":f"{type(exc).__name__}: {exc}",
                "elapsed_seconds":time.time()-start,
                "warning":"Não foi produzida inclusão finita de W12; nenhum certificado deve ser inferido."}
    txt=json.dumps(result,indent=2,ensure_ascii=False)
    print(txt)
    args.output.write_text(txt+'\n',encoding='utf-8')
    print(f"Relatório salvo: {args.output}",file=sys.stderr)
    if not result['status'].startswith('DIRECT_P92_'):
        raise SystemExit(2)

if __name__=='__main__':
    main()