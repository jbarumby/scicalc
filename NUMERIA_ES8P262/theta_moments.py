"""P262: directed Arb integrals for KPS theta moments and printed Eq. (4.18).

Source: Konstantopoulos–Patie–Sarkar (2024), AIF 74, p. 397,
https://www.numdam.org/item/10.5802/aif.3600.pdf
Numerical claims below apply to the formula AS PRINTED, NOT to RH.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from math import factorial, comb
from pathlib import Path
import json

from flint import acb, arb, ctx

SOURCE = "https://www.numdam.org/item/10.5802/aif.3600.pdf"
FORMULA_PAGE = 397
U_END = 5
M_MAX = 3

def need(value, message):
    if not value:
        raise ArithmeticError(message)

def tail_rational(k):
    """Upper bound for missing theta m>=4 on [0,5] and u>=5 all m.

    u = log x; integrand u^k exp(u/4) sum_m exp(-pi*m^2 exp(u)).
    For m>=4: theta <= 2exp(-48)exp(-48u), using pi>3,
    exp(u)>=1+u, exp(-48)<2^-48.
    For u>=5: exp(u)>100(1+u-5), theta <=
    2exp(-300)exp(-300(u-5)), exp(5/4)<4, and e>2.
    """
    if type(k) is not int or k < 0 or k > 30:
        raise ValueError("even nonnegative k up to 30 expected")
    theta_m_tail = Fraction(2 * factorial(k), 2**48 * 47**(k+1))
    u_tail = Fraction(8, 2**300) * sum(
        (Fraction(comb(k,j)*5**(k-j)*factorial(j),299**(j+1))
         for j in range(k+1)), Fraction(0))
    return theta_m_tail + u_tail

def rational_ball(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator

def endpoints(value):
    """Exact outward binary endpoints; [mantissa, exponent] pairs."""
    lower = value.lower().man_exp()
    upper = value.upper().man_exp()
    return [[str(lower[0]), int(lower[1])],
            [str(upper[0]), int(upper[1])]]

def moment(k, *, precision_bits=192):
    """Enclose F(k)=int_1^inf log(x)^k x^-3/4 theta0(x)dx.

    Complete finite segment uses acb.integral with analytic entire callback.
    Explicit rational error bounds cover theta m>=4 and u>=5.
    """
    if type(k) is not int or k < 0 or k%2:
        raise ValueError("nonnegative even k required")
    if precision_bits < 128:
        raise ValueError("precision must be at least 128 Arb bits")
    with ctx.workprec(precision_bits):
        pi=arb.pi()
        total=arb(0)
        segments=((0,1),(1,2),(2,3),(3,4),(4,5))
        for m in range(1, M_MAX+1):
            for start,stop in segments:
                def integrand(u, analytic):
                    # All factors are entire in u, no branch cut.
                    return (u**k) * (u/4).exp() * (-pi*m*m*u.exp()).exp()
                value=acb.integral(integrand, start, stop, rel_tol=arb(2)**(-96),
                                   abs_tol=arb(2)**(-120),
                                   eval_limit=12000, depth_limit=60)
                need(value.real.is_finite() and value.imag.is_finite(),
                     f"nonfinite Arb integral at k={k}, m={m}, segment={start,stop}")
                need(value.imag.overlaps(arb(0)),
                     f"unexpected imaginary integral at k={k},m={m}")
                total += value.real
        error=tail_rational(k)
        total += arb(0, rational_ball(error).upper())
        need(total > 0, f"moment positivity inconclusive: F({k})")
        return total, error

def paper_G(n, moments):
    """Equation (4.18) G(2n), as printed in journal PDF."""
    if n < 1:
        raise ValueError("n>=1")
    return 64*n*(2*n-1)*moments[2*n-2]/moments[2*n]-1

def paper_phi(n, G):
    """Printed Eq. 4.18: phi(n+1)=-G(2n)/(8(n+1)G(2n+2))."""
    return -G[n]/(8*(n+1)*G[n+1])

def audit(*, max_n=3, precision_bits=192):
    """Compute F(0),...,F(2(max_n+1)), and gates for n<=max_n.

    A negative phi is ONLY an incompatibility between the PRINTED formula
    and Bernstein nonnegativity, subject to equation normalization audit.
    It is NOT a refutation of RH, nor the Bernstein representation question.
    """
    if not 1 <= max_n <= 8:
        raise ValueError("max_n must be 1..8")
    with ctx.workprec(precision_bits):
        F={}
        tail={}
        for k in range(0,2*(max_n+2),2):
            F[k],tail[k]=moment(k,precision_bits=precision_bits)
        G={n:paper_G(n,F) for n in range(1,max_n+2)}
        phi={n:paper_phi(n,G) for n in range(1,max_n+1)}
        result={
            "schema":"NUMERIA.ES8P262.KPS_THETA_BERNSTEIN_SIGN_GATE.1",
            "paper":{"url":SOURCE,"page_printed":FORMULA_PAGE,"equations":["4.15","4.16","4.17","4.18"]},
            "moment_formula":"F(k)=int_0^inf u^k exp(u/4) sum_{m>=1} exp(-pi*m*m*exp(u)) du",
            "interval_math":"python-flint acb.integral, directed Arb enclosures",
            "truncation":"u<=5,m<=3; rational error for both omitted tails",
            "precision_bits":precision_bits,
            "moments":{
                 str(k):{"arb_enclosure_dyadic":endpoints(F[k]),
                         "omitted_tail_upper_exact":str(tail[k]),
                         "positive_proved":bool(F[k]>0)}
                 for k in F},
            "G_as_printed":{
                 str(2*n):{"arb_enclosure_dyadic":endpoints(G[n]),
                            "sign":"positive" if G[n]>0 else "negative" if G[n]<0 else "inconclusive"}
                 for n in G},
            "phi_as_printed":{
                 str(n+1):{"arb_enclosure_dyadic":endpoints(phi[n]),
                           "sign":"negative" if phi[n]<0 else "positive" if phi[n]>0 else "inconclusive"}
                 for n in phi},
            "bernstein_necessary":"phi(x)>=0 for x>=0",
            "negative_printed_phi_detected":any(p<0 for p in phi.values()),
            "printed_formula_normalization_requires_independent_check":True,
            "equivalence_RH_not_used":True,
            "original_article_error_claimed":False,
            "Bernstein_candidate_disproved":False,
            "Riemann_hypothesis_proved_or_disproved":False,
            "Weil_positivity_claimed":False,
        }
        need(result["negative_printed_phi_detected"],
             "expected sign obstruction not reproduced")
        return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-n",type=int,default=3)
    parser.add_argument("--bits",type=int,default=192)
    args=parser.parse_args()
    result=audit(max_n=args.max_n,precision_bits=args.bits)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("P262 KPS sign pilot: ", {k:v["sign"] for k,v in result["phi_as_printed"].items()})
if __name__=="__main__":
    main()
