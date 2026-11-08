#!/usr/bin/env python3
"""Nova tentativa, independente da P89, de enclausurar W12 na normalização P92.

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
import json
import math
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

    @classmethod
    def make(cls, tol_text: str, eval_limit: int, depth_limit: int):
        h = Q(2, 5)
        tol = arb(tol_text)
        # eps=1/65536 conforme a cota racional P93 para a parte corretiva.
        eps = Q(1, 65536)
        return cls(h, tol, eps, eval_limit, depth_limit)

    def bump(self, z: acb, analytic: bool) -> acb:
        """Extensão intervalar segura do bump P92, real e compactamente apoiado.

        No interior, avalia a fórmula por Arb/Acb. Fora do suporte, o valor
        real é zero. Em células que cruzam a borda, devolve uma caixa complexa
        ampla contendo a imagem; se o integrador exige analiticidade nessa
        célula, devolve não-finito para forçar subdivisão/falha fechada.
        """
        rl, ru = real_bounds(z)
        h_lo, h_hi = self.h.lower(), self.h.upper()
        entirely_inside = bool(rl > -h_lo and ru < h_lo)
        entirely_right = bool(rl >= h_hi)
        entirely_left = bool(ru <= -h_hi)
        if entirely_inside:
            den = 1 - (z / self.h) ** 2
            if den.contains(0):
                return acb_nan() if analytic else acb(arb(0, 1), arb(0, 1))
            return (-1 / den).exp()
        if entirely_right or entirely_left:
            # A integração é ao longo da reta real; nessa região a função
            # compactamente apoiada é identicamente zero e analítica.
            if imag_is_zero(z) or analytic:
                return acb(0)
            return acb(arb(0, 1), arb(0, 1))
        if analytic:
            return acb_nan()
        # |b_P92(x)| <= e^{-1} em x real, e a caixa quadrada contém essa faixa.
        rad = (1 / arb.const_e()).upper()
        return acb(arb(0, rad), arb(0, rad))

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
        self.calls_R += 1
        # Se todo o enclosure de s fica fora de [-2h,2h], R(s)=0 por suporte.
        sl, su = real_bounds(s)
        support = 2 * self.h
        if sl >= support.upper() or su <= -support.lower():
            return acb(0)

        def f(y, analytic):
            return self.bump(y + s, analytic) * self.bump(y, analytic)
        return self.integrate(f, -self.h, self.h, f"R(s={s})")

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

    def outer_integrand(self, delta: acb, rdelta: acb, panel_lo: arb, panel_hi: arb):
        def f(t, analytic):
            if analytic:
                # Acb.integral requer a verificação explícita de que o disco
                # de t não cruza origem nem os contatos do suporte (ou t=delta).
                tl, tu = real_bounds(t)
                if not (tl > panel_lo and tu < panel_hi):
                    return acb_nan()
                critical = (arb(0), 2*self.h-delta.real, delta.real,
                            2*self.h+delta.real)
                for c in critical:
                    if tl <= c.upper() and tu >= c.lower():
                        return acb_nan()
            rplus = self.correlation(t + delta)
            rminus = self.correlation(t - delta)
            F = (rplus + rminus) / 2
            den = t.sinh()
            if den.contains(0):
                return acb_nan()
            return ((t/2).exp()*F-rdelta)/den
        return f

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

        # Para [0,eps], usa identidade K + termo corretivo. A cota analítica
        # rigorosa do termo corretivo é 12 eps^2/h, válida também para bump P92.
        small_main = rdelta * self.small_kernel_integral()
        correction_bound = 12*eps*eps/h
        small_correction = acb(arb(0, correction_bound), arb(0))

        # Divisão nos pontos de contato: 2h-delta, t=delta (s=0), T=2h+delta.
        p1 = 2*h-delta.real
        p2 = delta.real
        points = [eps, p1, p2, T]
        if not all(points[i].upper() < points[i+1].lower()
                   for i in range(len(points)-1)):
            raise ArithmeticError("painéis não ficaram ordenados/disjuntos em Arb")
        I_mid = acb(0)
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
        return {
            "status": status,
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
            "warning": "No P89 data were used. A finite Arb ball is a rigorous enclosure only if the piecewise bump enclosures, Acb.integral analyticity flags, and small-t remainder lemma are audited as stated. A negative sign is certified only if the final real upper bound is <0. This is one matrix entry, not positivity of a full quadratic form or RH."
        }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prec',type=int,default=128)
    ap.add_argument('--tol',default='1e-12')
    ap.add_argument('--eval-limit',type=int,default=100000)
    ap.add_argument('--depth-limit',type=int,default=30)
    ap.add_argument('--output',type=Path,default=Path('W12_ARB_P92_DIRECT.json'))
    args=ap.parse_args()
    ctx.prec=args.prec
    start=time.time()
    try:
        calc=W12P92.make(args.tol,args.eval_limit,args.depth_limit)
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