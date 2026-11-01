#!/usr/bin/env python3
"""V23 rigorous *Fourier-norm* error for C-infinity tapered NUMERIA modes.

This is not yet a Weil-form or archimedean integral error bound.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,json
from AUDITAR_V23_DERIVADAS import L,PI_HI,SQLO,SQHI,decimal_out

def gaussian_mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def gaussian_add(a,b):return (a[0]+b[0],a[1]+b[1])

def symmetry_adversary():
    # P(1/2+t)=(t^2+17/16)^2-t^2/4=t^4+(15/8)t²+289/256.
    # Has off-line zeros t=+/-1/4 +/- i, despite P(s)=P(1-s).
    roots=[]
    for a in (-1,1):
        for b in (-1,1):
            t=(F(a,4),F(b))
            t2=gaussian_mul(t,t);t4=gaussian_mul(t2,t2)
            p=gaussian_add(gaussian_add(t4,(F(15,8)*t2[0],F(15,8)*t2[1])),(F(289,256),F(0)))
            assert p==(0,0)
            roots.append([str(F(1,2)+t[0]),str(t[1])])
    return {"polynomial_in_t":"t^4+(15/8)t^2+289/256, t=s-1/2",
            "functional_symmetry":"P(s)=P(1-s), real coefficients",
            "off_line_zeros":roots,"off_line_symmetry_does_not_imply_RH":True}

def taper_bounds(delta):
    d=F(delta)
    if not F(0)<d<L:raise ValueError("require 0<delta<L")
    # phi_k(x)=cos(pi*x/(2L))/sqrt(L) (k=1), sin(pi*x/L)/sqrt(L) (k=2);
    # each has zero boundary values and |phi_k'| <= k*pi/(2L sqrt L).
    # Define smooth cutoff chi_delta=eta((L-|x|)/delta),
    # eta(t)=b(t)/(b(t)+b(1-t)) where b(t)=exp(-1/t) for t>0, else 0.
    # chi=1 for |x|<=L-delta and 0 for |x|>=L; no cusp near x=0.
    # Integral |phi-chi phi| <= 2*delta^2*sup|phi'| (conservative).
    lip1=PI_HI/(2*L*SQLO)
    lip2=PI_HI/(L*SQLO)
    e1=2*d*d*lip1;e2=2*d*d*lip2
    supFourier=2*SQHI
    # For Fourier Fi-F_i,taper, triangle inequality in L1,
    # |FiFjbar-F'iF'jbar| <= 2B(ei+ej)+ei*ej (conservative).
    errors=[[2*supFourier*(a+b)+a*b for b in (e1,e2)] for a in (e1,e2)]
    return {"delta":str(d),"mode1_Fourier_sup_error_upper":decimal_out(e1,True),
            "mode2_Fourier_sup_error_upper":decimal_out(e2,True),
            "product_pointwise_fourier_error_upper":[[decimal_out(x,True) for x in row] for row in errors],
            "normalization":"Fourier: integral_{-L}^L phi_k(x)exp(irx) dx, no additional 2*pi factor",
            "limitation":"Pointwise Fourier error does not imply global weighted Weil integral error without tail analysis."}

def audit():
    a=taper_bounds(F(1,100));b=taper_bounds(F(1,1000))
    # For delta ten times smaller: L1 Fourier error 100 times smaller.
    assert F(a["mode1_Fourier_sup_error_upper"])>=99*F(b["mode1_Fourier_sup_error_upper"])
    assert F(a["mode2_Fourier_sup_error_upper"])>=99*F(b["mode2_Fourier_sup_error_upper"])
    return {"schema":"NUMERIA.V23.SMOOTH_CUTOFF_AND_SYMMETRY_ADVERSARY.1",
            "status":"C_INFINITY_TAPER_FOURIER_ERROR_BOUNDED_SYMMETRY_COUNTEREXAMPLE_PROVED",
            "tapers":[a,b],"symmetry_adversary":symmetry_adversary(),
            "analytic_inputs_not_formalized":["Existence of C-infinity flat cutoff with standard exp(-1/t) construction",
                         "Mean-value theorem and Fourier L1 to L-infinity bound"],
            "not_certified":["Error in the global Weil quadratic form after taper",
                             "Applicability of Burnol short-support positivity for L=9/20",
                             "Positivity for arbitrary test functions or RH"]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
    try:result=audit()
    except Exception as e:print("FAIL",type(e).__name__,e);return 2
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("PASS: quantitative C-infinity taper Fourier bounds and exact functional-symmetry counterexample")
    return 0
if __name__=="__main__":raise SystemExit(main())
