"""Rigorous Arb audit of Chourasiya--Simonic v2 Table 1 B1.

Integral K1: substitute u=x^6; integrate x in [0,2] using acb.integral
and an explicit elementary bound for u>=64. No numerical point samples
are accepted as a certificate.
"""
import argparse
import json
from pathlib import Path
from flint import acb, arb, ctx
from NUMERIA_P76_COEFFICIENTS.audit_b23 import constants, q, dyadic
from NUMERIA_P76_COEFFICIENTS.table import TABLE, validate_table

T0=3*10**12
CUT=64
SOURCE="https://arxiv.org/html/2507.15184v2"

def rat(n,d):
    return arb(n)/arb(d)

def pw(x,e):
    if not x>0:
        raise ArithmeticError("non-positive real power base")
    return (x.log()*e).exp()

def tail_upper(c,c1,c2,L):
    """For u>=64: integrand <=D*exp(-u)*u^3, integrate exactly."""
    if not all(z>0 for z in (c,c1,c2,L)):
        raise ArithmeticError("tail requires positive constants")
    U=arb(CUT)
    D=(1+c/(2*L))**2*(1+c1*c.sqrt()+c2)
    t=(D*(-U).exp()*(U**3+3*U**2+6*U+6)).upper()
    if not t>0 or not t.is_finite():
        raise ArithmeticError("invalid tail")
    return t

def k1_integral(c,c1,c2,L):
    """Rigorous acb.integral enclosure on [0,64] plus strict tail bound."""
    def f(x,analytic):
        v=x**6
        A=acb(1)+(acb(1)+acb(c)*v).log(analytic=analytic)/(2*acb(L))
        B=acb(c1*c.sqrt())*x**3+acb(c2)
        return 6*v*(-v).exp()*A.pow(acb(4)/3,analytic=analytic)*B.pow(acb(2)/3,analytic=analytic)
    I=acb.integral(f,0,2,rel_tol=arb('1e-19'),
                   abs_tol=arb('1e-19'),eval_limit=30000,depth_limit=30)
    if not I.is_finite() or not I.imag.contains(0) or not I.real.lower()>0:
        raise ArithmeticError("quadrature failed or inconclusive")
    tail=tail_upper(c,c1,c2,L)
    return I.real.lower(),(I.real.upper()+tail).upper(),tail

def row_audit(row,cs):
    """Reconstruct every factor in Theorem 1 and (10), (35)--(40)."""
    d,kap,delta=q(row.d),q(row.kappa),q(row.delta)
    s1,s2=rat(row.sigma_start,32),rat(row.sigma_start+1,32)
    L,pi,m0=cs["logh"],cs["pi"],cs["m0"]
    # Choose EXACT real A0=10^printed log10 A0; no unknown optimizer digits.
    A0=(arb(10).log()*q(row.log10_a0)).exp()
    if not arb(10)**10<A0<arb(T0):
        raise ArithmeticError("invalid A0")
    La=A0.log()
    half=(A0/2).log()
    a1=A0/(2*pi*pi)*half**4+q("34.842")*A0*half**3+arb("1.3897e16")
    r=1-arb.const_log2()/La
    a2=r**4/(2*pi*pi)+q("34.842")/La*r**3+arb("9.5079e15")/(A0*La**4)
    if not a1>0 or not a2>0:raise ArithmeticError("invalid a1/a2")
    c=3/(4*kap)
    lower,upper,tail=k1_integral(c,cs["c1"],cs["c2"],L)
    F=pw(2*a2*c.sqrt(),rat(1,3))
    correction=pw(2*a1,rat(1,3))*pw(cs["c2"]+A0*cs["c1"]/T0,rat(2,3))/(pw(arb(T0),rat(1,3))*pw(L,rat(4,3)))
    k1lo=(F*lower+correction).lower()
    k1=(F*upper+correction).upper()
    if not k1lo>0:raise ArithmeticError("K1 lower bound inconclusive")
    two13=pw(arb(2),rat(1,3))
    k2=two13*c.sqrt()*pi.sqrt()/2+two13*(4*kap/3)*pw(A0/arb(T0),arb(3))
    c3=pi*m0*q(".54")/(2*delta)*(2*delta*arb.const_euler()/L).exp()
    c4=q(".54")/(5*delta*delta.exp())+q(".605")*(-2*delta).exp()/L
    k3=c3+c4*(pi.sqrt()/(2*(2*kap).sqrt())+pi*m0/T0)
    l0=(kap*(1+(1/(4*arb(T0))*(1+2*delta/L))**2)).exp()/(1-arb(1)/T0)
    l1=max(arb(1),pw(k1,3*delta/L))*max(arb(1),pw(k3,-3*delta/L))*(3*delta+4*delta*L.log()/L).exp()
    l2=max(arb(1),pw(k2,3*delta/L))*max(arb(1),pw(k3,-3*delta/L))*(3*delta).exp()
    common=max(arb(1),pw(l0,-2*d/L))*max(arb(1),pw(k3,-3*d/L))
    common*=max(pw(l0,2/(2-s1)),pw(l0,2/(2-s2)))
    common*=max(pw(k3,(2*s1-1)/(2-s1)),pw(k3,(2*s2-1)/(2-s2)))
    A1=l1*(d*(3+4*L.log()/L)).exp()*common
    A1*=max(arb(1),pw(k1,3*d/L))
    A1*=max(pw(k1,3*(1-s1)/(2-s1)),pw(k1,3*(1-s2)/(2-s2)))
    A2=l2*(3*d).exp()*common
    A2*=max(arb(1),pw(k2,3*d/L))
    A2*=max(pw(k2,3*(1-s1)/(2-s1)),pw(k2,3*(1-s2)/(2-s2)))
    B1=(A1+A2*pw(L,-6*(1-s2)/(2-s2)))/(2*pi*d)
    upper_B1=B1.upper()
    margin=(q(row.b1)-upper_B1).lower()
    return {
        "sigma":[f"{row.sigma_start}/32",f"{row.sigma_start+1}/32"],
        "printed_B1":row.b1,
        "exact_log10_A0_choice":row.log10_a0,
        "k1_integral_lower_dyadic":dyadic(lower),
        "k1_integral_upper_dyadic":dyadic(upper),
        "k1_tail_upper_dyadic":dyadic(tail),
        "k1_lower_dyadic":dyadic(k1lo),
        "k1_upper_dyadic":dyadic(k1),
        "B1_upper_dyadic":dyadic(upper_B1),
        "margin_lower_dyadic":dyadic(margin),
        "strict":bool(margin>0)
    }

def audit(precision=192,one=None):
    if precision not in (192,256,320):raise ValueError("precision 192,256,320 only")
    if one is not None and (type(one) is not int or not 1<=one<=16):
        raise ValueError("row index must be 1..16")
    validate_table()
    with ctx.workprec(precision):
        cs=constants()
        L,pi,m0=cs["logh"],cs["pi"],cs["m0"]
        cs["c1"]=6/pi**2+q("1.048")/L
        cs["c2"]=pi*m0*q(".62")/L+6*m0/(pi*T0)+pi*m0*q("1.048")/(arb(T0)*L)
        selected=[TABLE[one-1]] if one else TABLE
        rows=[]
        for index,row in enumerate(selected,one or 1):
            r=row_audit(row,cs);r["index"]=index;rows.append(r)
        passed=all(r["strict"] for r in rows)
        return {
            "schema":"NUMERIA.P77.CS_V2.B1.1",
            "source":SOURCE,"precision_bits":precision,
            "quadrature":"analytic u=x^6 over x=0..2; rigorous elementary u>=64 tail",
            "rows":rows,"tested_rows_passed":passed,
            "full_16_B1_verified":passed and len(rows)==16,
            "global_RH_proved":False
        }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--precision",type=int,default=192)
    p.add_argument("--row",type=int)
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    data=audit(a.precision,a.row)
    if a.out:a.out.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"rows":len(data["rows"]),"passed":data["tested_rows_passed"],"precision":a.precision}))
    if not data["tested_rows_passed"]:raise SystemExit("P77 proof gate FAILED")

if __name__=="__main__":main()
