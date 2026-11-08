#!/usr/bin/env python3
"""ES8-P89 strict Arb intervals: poles + prime powers for 8 smooth bump translates.
No claim about the archimedean component, full matrix positivity, or RH.
"""
import argparse,hashlib,json,platform,sys
from pathlib import Path
from flint import acb,arb,ctx
import flint
NAMES=(1,2,3,5,7,11,13,17)
def primepowers(k):
    d={}
    for p in range(2,k+1):
        if any(p%j==0 for j in range(2,int(p**0.5)+1)):continue
        n=p
        while n<=k:
            assert n not in d
            d[n]=p
            n*=p
    return d

def run(shard,shards,bits=256,tol="1e-12"):
    if not (0<=shard<shards and bits>=160): raise ValueError("bad args")
    ctx.prec=bits
    assert ctx.prec==bits
    h=arb(2)/5
    centers=[arb(n).log() for n in NAMES]
    def fmt(b):
        return [str(y) for x in (b.mid(),b.rad()) for y in x.man_exp()]
    def bump(z,analytic):
        real=z.real
        if bool(real.lower()>-h) and bool(real.upper()<h):
            return (1-acb(1)/(1-(z/h)**2)).exp()
        if bool(real.upper() < -h) or bool(real.lower()>h):
            return acb(0)
        if analytic:return acb(arb("nan"))
        if not bool(z.imag==0):return acb(arb("nan"))
        return acb(arb("0.5 +/- 0.5"))
    def integrate(f,a,b,tag):
        val=acb.integral(f,a,b,abs_tol=arb(tol),rel_tol=arb("1e-25"),
                         deg_limit=150,eval_limit=200000,depth_limit=80)
        if not val.is_finite():raise ArithmeticError("nonfinite: "+tag)
        if not bool(val.imag.contains(0)):raise ArithmeticError("nonreal: "+tag)
        if not bool(val.real.rad().upper()<arb("5e-9")):
            raise ArithmeticError("inaccurate: "+tag+" "+str(val.real))
        return val.real
    U=integrate(lambda z,a:bump(z,a)*(z/2).exp(),-h,h,"moment +")
    V=integrate(lambda z,a:bump(z,a)*(-z/2).exp(),-h,h,"moment -")
    if not bool((U-V).contains(0)):raise ArithmeticError("unequal moments")
    R0=None
    count={"integrated_correlations":0,"small_shift":0,"disjoint_zero":0}
    def corr(s,tag):
        nonlocal R0
        if bool(s.upper()<0):s=-s
        if bool(s.lower()>=2*h):
            count["disjoint_zero"]+=1
            return arb(0)
        if bool(s.lower()<=0):
            if not bool(s.abs_upper()<arb("1e-25")):raise ArithmeticError("sign ambiguous: "+tag)
            if R0 is None:
                R0=integrate(lambda z,a:bump(z,a)**2,-h,h,"R0")
            count["small_shift"]+=1
            # |R(s)-R(0)| <= 2 |s| since ||b'||_1=2 and ||b||_inf=1.
            return R0+arb(0,2*s.abs_upper())
        if not bool(s.upper()<2*h):
            width=(2*h-s).abs_upper()
            if not bool(width<arb("1e-12")):raise ArithmeticError("support ambiguous: "+tag)
            return arb(0,width)
        length=2*h-s
        def f(z,analytic):
            y=-h+length*z
            return length*bump(y,analytic)*bump(y+s,analytic)
        count["integrated_correlations"]+=1
        return integrate(f,arb(0),arb(1),tag)
    if not bool(arb(40).log()>centers[-1]+2*h):
        raise ArithmeticError("primes beyond 39 may contribute")
    powers=primepowers(39)
    records=[]
    idx=0
    for i in range(8):
        for j in range(i,8):
            if idx%shards!=shard:
                idx+=1;continue
            idx+=1
            delta=centers[j]-centers[i]
            P=U*V*((delta/2).exp()+(-delta/2).exp())
            S=arb(0)
            for n,p in powers.items():
                x=arb(n).log()
                S+=-arb(p).log()/arb(n).sqrt()*(corr(x+delta,f"{i},{j} {n}+")+corr(x-delta,f"{i},{j} {n}-"))
            records.append({"i":i,"j":j,"poles":fmt(P),"primes":fmt(S)})
            print("COMPLETED",i,j,flush=True)
    return {"status":"CERTIFIED_POLE_PRIME_COMPONENTS_ONLY","shard":shard,"shards":shards,
            "precision_bits":ctx.prec,"python_flint":flint.__version__,
            "center_integers":list(NAMES),"h":"2/5","records":records,
            "prime_powers_considered":sorted(powers),"evaluations":count,
            "scientific_limit":"Archimedean integrals not enclosed: no full positivity/RH claim",
            "environment":{"platform":platform.platform(),"python":sys.version}}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--shard",type=int,required=True)
    p.add_argument("--shards",type=int,default=8)
    p.add_argument("--bits",type=int,default=256)
    p.add_argument("--out",required=True)
    p.add_argument("--tol",default="1e-12")
    a=p.parse_args()
    out=Path(a.out)
    if out.exists():raise FileExistsError(out)
    result=run(a.shard,a.shards,a.bits,a.tol)
    out.write_text(json.dumps(result,indent=2)+"\n")
