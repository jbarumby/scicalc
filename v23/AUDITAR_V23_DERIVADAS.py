#!/usr/bin/env python3
"""V23: independently reproduce analytic upper bounds M>=sup|f''| of B3.

Uses only exact rational arithmetic, elementary calculus inequalities,
integral comparisons, Euler-Maclaurin remainder (stated, not Lean-proven).
Separate implementation from E7 sources; rejects if any independent bound
exceeds the historical fixed-point bound. Does NOT claim RH.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, math, zipfile

E7_SHA="4acf99a2625a262b670ee73ed104b970da0951ffd98250cdd254f387110b02a7"
L=F(9,20); A=F(1,4); TERMS=128
SOURCES={"G11":"corrigidos/RESULTADO_B3_G11_E2_CORRIGIDO.json",
         "G22":"corrigidos/RESULTADO_B3_G22_E2_CORRIGIDO.json"}

def atan_unit(d,n=70):
    z=F(1,d)
    val=sum(((-1)**k*z**(2*k+1)/F(2*k+1) for k in range(n)),F(0))
    rem=z**(2*n+1)/F(2*n+1)
    return (val,val+rem) if n%2==0 else (val-rem,val)

def pi_enclosure():
    a,b=atan_unit(5);c,d=atan_unit(239)
    return 16*a-4*d,16*b-4*c

PI_LO,PI_HI=pi_enclosure()
Q=F(1,3)
LOG2_LO=sum((2*Q**(2*k+1)/F(2*k+1) for k in range(100)),F(0))
LOG2_HI=LOG2_LO+2*Q**201/(F(201)*(1-Q*Q))

def log_enclosure(x,terms=100):
    x=F(x)
    if x<=0:raise ValueError("log argument must be positive")
    power=0
    while x>=2:x/=2;power+=1
    while x<1:x*=2;power-=1
    q=(x-1)/(x+1)
    s=sum((2*q**(2*k+1)/F(2*k+1) for k in range(terms)),F(0))
    rem=2*q**(2*terms+1)/(F(2*terms+1)*(1-q*q))
    if power>=0:return s+power*LOG2_LO,s+rem+power*LOG2_HI
    return s+power*LOG2_HI,s+rem+power*LOG2_LO

def gamma_enclosure(n=1024):
    # Euler-Maclaurin alternating remainder for H_n-log n.
    B=[F(1,6),-F(1,30),F(1,42),-F(1,30),F(5,66),
       -F(691,2730),F(7,6),-F(3617,510)]
    h=sum((F(1,j) for j in range(1,n+1)),F(0))-F(1,2*n)
    t=[B[k]/F(2*(k+1)*n**(2*(k+1))) for k in range(8)]
    a,b=log_enclosure(n)
    return h+sum(t,F(0))-b, h+sum(t[:7],F(0))-a

GAM_LO,GAM_HI=gamma_enclosure()
LOGPI_LO=log_enclosure(PI_LO)[0]
LOGPI_HI=log_enclosure(PI_HI)[1]
OMEGA0_LO=-GAM_HI-PI_HI/2-3*LOG2_HI-LOGPI_HI
OMEGA0_HI=-GAM_LO-PI_LO/2-3*LOG2_LO-LOGPI_LO
assert F(-6)<OMEGA0_LO<OMEGA0_HI<0

def omega_enclosure(r):
    r=F(r);b2=r*r/4
    subtotal=sum((b2/(u*(u*u+b2)) for u in (F(j)+A for j in range(TERMS))),F(0))
    u=F(TERMS)+A
    low,high=log_enclosure(1+b2/(u*u))
    return (OMEGA0_LO+subtotal+low/2,
            OMEGA0_HI+subtotal+high/2+b2/(u*(u*u+b2)))

def omega_derivatives_bound(left,right):
    a,b=F(left),F(right)
    first=F(0);second=F(0)
    for j in range(TERMS):
        u=F(j)+A;u2=u*u
        first+=min(u*b/(2*(u2+a*a/4)**2),F(1,3*u2))
        second+=min(u*max(abs(u2-3*a*a/4),abs(u2-3*b*b/4))/
                    (2*(u2+a*a/4)**3),F(1,2*u**3))
    u=F(TERMS)+A
    # Integral comparison of the positive tail sums.
    first+=(1/u**2+1/u)/3
    second+=F(1,2*u**3)+F(1,4*u*u)
    return first,second

def sqrt_l():
    sc=10**65
    n=math.isqrt(L.numerator*sc*sc//L.denominator)
    return F(n,sc),F(n+1,sc)

SQLO,SQHI=sqrt_l()

def m2_independent(left,right,kind):
    l,h=F(left),F(right)
    assert 0<=l<h<=160 and kind in SOURCES
    omega=max(abs(z) for r in (l,h) for z in omega_enclosure(r))
    d1,d2=omega_derivatives_bound(l,h)
    if (kind=="G11" and h<=6) or (kind=="G22" and h<=12):
        a0=2*SQHI
        a1=SQHI*L
        a2=F(2,3)*SQHI*L**2
    else:
        w2=PI_HI**2/(4*L**2) if kind=="G11" else PI_HI**2/L**2
        D=l*l-w2
        if D<=0:raise ValueError("cannot bound Fourier denominator on this band")
        C=(PI_HI if kind=="G11" else 2*PI_HI)/(L*SQLO)
        a0=C/D
        a1=C*(L/D+2*h/D**2)
        a2=C*(L**2/D+4*L*h/D**2+2/D**2+8*h*h/D**3)
    # Derivative bound for f=A(r)^2*Omega(r)/pi by product rule.
    M=((2*a1*a1+2*a0*a2)*omega+4*a0*a1*d1+a0*a0*d2)/PI_LO
    return M

def parse_number(z):return F(int(z["numerator"]),int(z["denominator"]))

def decimal_out(x,up,places=45):
    scale=10**places
    k=x.numerator*scale//x.denominator
    if up and F(k,scale)<x:k+=1
    sign="-" if k<0 else "";k=abs(k)
    return f"{sign}{k//scale}.{k%scale:0{places}d}"

def audit(e7):
    data=e7.read_bytes()
    assert hashlib.sha256(data).hexdigest()==E7_SHA,"E7 provenance hash differs"
    rows=[]
    with zipfile.ZipFile(e7) as arc:
        for kind,name in SOURCES.items():
            bands=json.loads(arc.read(name))["mesh"]
            prev=F(0)
            for x in bands:
                l,h=F(x["from"]),F(x["to"])
                assert l==prev and l<h and x["cells"]>0
                old=parse_number(x["f2_absolute_bound"] if kind=="G11" else x["f2_bound"])
                independent=m2_independent(l,h,kind)
                assert independent<=old, f"independent M exceeds B3: {kind} {l}-{h}"
                rows.append({"component":kind,"band":[str(l),str(h)],
                  "independent_M_upper":decimal_out(independent,True),
                  "historical_M_upper":decimal_out(old,True),
                  "verified_margin_lower":decimal_out(old-independent,False),
                  "encloses_independent_calculus_bound":True})
                prev=h
            assert prev==160
    assert len(rows)==11
    assert sum(x["component"]=="G11" for x in rows)==6
    assert sum(x["component"]=="G22" for x in rows)==5
    return {"schema":"NUMERIA.V23.INDEPENDENT_ANALYTIC_DERIVATIVE_BOUNDS.1",
      "status":"ELEVEN_SECOND_DERIVATIVE_MAJORANTS_INDEPENDENTLY_REDERIVED_CONDITIONALLY",
      "provenance_e7_sha256":E7_SHA,"bands_count":11,"rows":rows,
      "bounds_sources":["Omega=psi(1/4+i r/2).real-log(pi) positive-term series",
                        "Integral comparison for Omega series tail",
                        "Derivative identities and rational envelopes for Omega' and Omega''",
                        "Compact cosine/sine Fourier integral bounds at removable poles",
                        "Rational denominator differentiation away from poles",
                        "Product rule for (A_k)^2*Omega/pi",
                        "Euler-Maclaurin alternating gamma remainder, Machin pi and atanh log"],
      "analytic_inputs_not_formalized":["Real Gamma/digamma reflection and duplication identities",
                       "Euler-Maclaurin alternating remainder for harmonic numbers",
                       "Classical Fourier representation of compactly supported cosine/sine modes",
                       "Alternating series and integral-comparison inequalities"],
      "not_certified":["Lean proof of all calculus identities and assumptions",
                       "Recomputed full nodal B5 transcendental intervals or LDL matrices",
                       "Universal Weil positivity or Riemann Hypothesis"]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--e7",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    try:j=audit(a.e7)
    except Exception as e:
        print("FAIL",type(e).__name__,e)
        return 2
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(j,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("PASS: independent rational second-derivative upper bounds for 11/11 B3 bands")
    return 0

if __name__=="__main__":raise SystemExit(main())
