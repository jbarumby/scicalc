#!/usr/bin/env python3
"""Independent Arb GL60, 200-mode, 50-panel quadrature; NO RH claim.
Values are interval enclosures of the QUADRATURE SUM only, not the integral.
All 20,100 upper-triangular entries are emitted as exact dyadic balls.
"""
import argparse, gzip, hashlib, json, math, os, time
from fractions import Fraction
from pathlib import Path
from flint import arb, acb, arb_mat, ctx

SCHEMA = "NUMERIA.ES8P226.ARB_200MODE_50PANELS.1"

def bernoulli_even(m=24):
    b=[Fraction(1)]
    for n in range(1,2*m+1):
        b.append(-sum((Fraction(math.comb(n+1,k))*b[k] for k in range(n)),Fraction(0))/(n+1))
    return [arb(b[2*k].numerator)/arb(b[2*k].denominator) for k in range(1,m)]

def sphere_anchor(n,x):
    pref=x**n
    for i in range(1,n+1): pref/=arb(2*i+1)
    a=arb(1);s=arb(0);xx=x*x
    for k in range(700):
        r=Fraction(161**2,2*(k+1)*(2*n+2*k+3))
        if r<=Fraction(1,2) and a.upper() < arb(10)**-205:
            return pref*(s+arb(0,(2*a.upper()).upper())),k
        s+=a if k%2==0 else -a
        a*=xx/(2*(k+1)*(2*n+2*k+3))
    raise ArithmeticError("Bessel series remainder not bounded")

def sphere_all(x):
    if not (x.lower()>0 and x.upper()<161):
        raise ValueError("outside Bessel argument domain")
    j=[None]*402
    j[400],k0=sphere_anchor(400,x)
    j[401],k1=sphere_anchor(401,x)
    for n in range(400,0,-1):
        j[n-1]=(2*n+1)*j[n]/x-j[n+1]
    return j[:400],(k0,k1)

class Symbol:
    def __init__(self,T):
        self.pi=arb.pi(); self.l2=arb(2).log();self.l3=arb(3).log()
        self.l4=2*self.l2
        self.a2=2*self.l2/arb(2).sqrt()
        self.a3=2*self.l3/arb(3).sqrt()
        self.a4=self.l2
        self.beta=(arb(T)/(2*self.pi)).log()-arb(1)/T-self.a2-self.a3-self.a4
        self.B=bernoulli_even(24)
        E=(arb(4*math.factorial(47))/arb(480**48)).upper()
        self.remainder=arb(0,E)
    def digamma(self,t):
        z=acb(arb(1)/4,t/2);w=z+80;inv=1/w;p=inv*inv
        v=w.log()-inv/2
        for k,b in enumerate(self.B,1):
            v-=b*p/(2*k);p*=inv*inv
        for k in range(80): v-=1/(z+k)
        return v.real+self.remainder
    def kernel(self,t):
        return (self.digamma(t)-self.pi.log()
                -self.a2*(t*self.l2).cos()
                -self.a3*(t*self.l3).cos()
                -self.a4*(t*self.l4).cos()
                -self.beta)/self.pi

def ball(x):
    m,e=x.mid().man_exp();r,k=x.rad().upper().man_exp()
    return [[str(m),str(e)],[str(r),str(k)]]

def integrate(sector,start,stop,bits):
    if sector not in ("even","odd") or stop-start!=50 or start%50:
        raise ValueError("invalid sector or 50-panel geometry")
    T=200 if sector=="even" else 150
    if start<0 or stop>T*2: raise ValueError("range beyond T")
    if bits<768: raise ValueError("precision below policy")
    ctx.prec=bits
    st=time.monotonic()
    roots=[arb.legendre_p_root(60,k,weight=True) for k in range(60)]
    assert all(roots[k][0].lower()>roots[k+1][0].upper() for k in range(59))
    assert sum((w for x,w in roots),arb(0)).contains(2)
    symbol=Symbol(T)
    L=arb(4)/5
    modes=[2*i+(sector=="odd") for i in range(200)]
    amps=[(-1)**(n//2)*(2*L*(2*n+1)).sqrt() for n in modes]
    accum=arb_mat(200,200); sample=None; lowcheck=0
    kmin=1000;kmax=0
    quarter=arb(1)/4
    for panel in range(start,stop):
        center=arb(2*panel+1)/4
        V=[];WV=[]
        for root,w in roots:
            t=center+quarter*root
            j,ks=sphere_all(L*t)
            kmin=min(kmin,*ks);kmax=max(kmax,*ks)
            if lowcheck<12:
                assert j[0].overlaps((L*t).sin()/(L*t))
                assert j[1].overlaps((L*t).sin()/(L*t)**2-(L*t).cos()/(L*t))
                lowcheck+=1
            u=[amps[i]*j[n] for i,n in enumerate(modes)]
            wt=quarter*w*symbol.kernel(t)
            V.extend(u);WV.extend([wt*a for a in u])
        J=arb_mat(60,200,V);WJ=arb_mat(60,200,WV)
        P=J.transpose()*WJ
        if panel==0: sample=[[i,j,ball(P[i,j])] for i,j in [(0,0),(0,1),(1,1),(10,120),(199,199)]]
        accum+=P
        print(f"verified panel {panel+1}/{stop}  elapsed {time.monotonic()-st:.1f}s",flush=True)
    entries=[[i,j,ball(accum[i,j])] for i in range(200) for j in range(i,200)]
    return {"schema":SCHEMA,"backend":"python-flint/Arb",
      "flint_version":__import__("flint").__version__,
      "precision_bits":bits,"sector":sector,"Tsharp":T,
      "L":"4/5","panel_range":[start,stop],
      "mode_count":200,"node_count":(stop-start)*60,
      "source_method":"independent arb.legendre_p_root + Binet-Stirling + 400/401 Bessel anchor",
      "gauss_degree":60,"quadrature_analytic_remainder_added":False,
      "beta_and_poles_added":False,"Riemann_hypothesis_proved":False,
      "status":"QUADRATURE_SUM_ENCLOSURES_ONLY",
      "Bessel_anchor_terms_minmax":[kmin,kmax],
      "first_panel_sample":sample,"entries_count":len(entries),
      "C_quadrature_upper_triangle_dyads":entries,
      "elapsed_seconds_diagnostic":round(time.monotonic()-st,1),
      "workflow_commit":os.environ.get("GITHUB_SHA","unknown")}

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--sector",choices=["even","odd"],default="even")
    ap.add_argument("--start",type=int,default=0)
    ap.add_argument("--stop",type=int,default=50)
    ap.add_argument("--bits",type=int,default=800)
    ap.add_argument("--out",default="P226_EVEN_P000_050.json.gz")
    x=ap.parse_args()
    d=integrate(x.sector,x.start,x.stop,x.bits)
    assert d["entries_count"]==20100
    raw=json.dumps(d,separators=(",",":")).encode()
    Path(x.out).write_bytes(gzip.compress(raw,compresslevel=9,mtime=0))
    print("SUCCESS",x.sector,x.start,x.stop,
          d["entries_count"],"entries",
          "SHA256",hashlib.sha256(Path(x.out).read_bytes()).hexdigest(),flush=True)
