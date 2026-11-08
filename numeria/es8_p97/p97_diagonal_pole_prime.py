#!/usr/bin/env python3
"""NUMERIA ES8-P97: calculate diagonal pole/prime terms with independent Arb integrals.

Normalization P89: b_h(x)=exp(1-1/(1-(x/h)^2)) on |x|<h; h=2/5.
For both centers log(1) and log(2), translation difference delta=0:
P_ii=2*A^2, S_ii=-2*(log 2)/sqrt(2)*R(log 2),
A=integral exp(y/2)*b(y) dy; R(s)=integral b(y+s)b(y)dy.
This is NOT the archimedean term or a complete Weil certificate.
"""
import argparse,json,sys
from fractions import Fraction
from flint import arb, acb, ctx, fmpq
def binary(x):
    return [str(y) for v in (x.mid(),x.rad()) for y in v.man_exp()]
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--bits",type=int,default=256)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    if a.bits<192:raise SystemExit("at least 192 bits")
    ctx.prec=a.bits
    h=arb(2)/5
    log2=arb(2).log()
    e=arb.const_e()
    def bump(z,analytic):
        re=z.real
        if bool(re.lower()>-h) and bool(re.upper()<h):
            return (1-acb(1)/(1-(z/h)**2)).exp()
        if bool(re.upper()<-h) or bool(re.lower()>h):
            return acb(0)
        if analytic:return acb(arb("nan"))
        if not bool(z.imag==0):return acb(arb("nan"))
        return acb(arb("0.5 +/- 0.5"))
    def integrate(f,a,b,label):
        v=acb.integral(f,a,b,abs_tol=arb("1e-13"),rel_tol=arb("1e-24"),
                       deg_limit=180,eval_limit=300000,depth_limit=90)
        if not v.is_finite():raise ArithmeticError(label+" nonfinite")
        if not bool(v.imag.contains(0)):raise ArithmeticError(label+" imaginary part nonzero")
        return v.real
    moment=integrate(lambda y,analytic:(y/2).exp()*bump(y,analytic),-h,h,"moment")
    overlap_high = h-log2
    if not bool(overlap_high> -h):raise ArithmeticError("empty overlap")
    correlation=integrate(lambda y,analytic:bump(y+acb(log2),analytic)*bump(y,analytic),
                          -h,overlap_high,"R(log2)")
    pole=2*moment**2
    prime=-(2*log2/arb(2).sqrt())*correlation
    assert bool(moment>0 and correlation>0)
    report={"status":"ARB_DIAGONAL_POLE_PRIME_COMPONENTS",
            "normalization":"P89", "bits":a.bits,
            "centers":[1,2], "diagonal_translation_invariance":True,
            "moment_A":binary(moment), "R_log2":binary(correlation),
            "pole_11_22":binary(pole),"prime_11_22":binary(prime),
            "warnings":["Arb internal integral: not mathematically independent from P89 moment strategy",
                        "independent implementation and invocation",
                        "archimedean integral and rigorously bounded tail needed for W_ii"]}
    with open(a.out,"w") as f:json.dump(report,f,indent=2)
    print(report["status"],str(pole),str(prime),flush=True)
if __name__=="__main__":main()
