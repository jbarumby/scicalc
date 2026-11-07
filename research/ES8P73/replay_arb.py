#!/usr/bin/env python3
"""Fresh ES8-P73 Arb sample replay. No historical coefficient is used to calculate a new one."""
import argparse
import hashlib
import json
import platform
import sys
import time
from fractions import Fraction as Q
from pathlib import Path

def endpoints(enc):
    m, e, r, f = map(int, enc)
    if r < 0: raise ValueError("negative radius")
    def power(n): return Q(2**n) if n >= 0 else Q(1, 2**-n)
    c, rad = m*power(e), r*power(f)
    return c-rad, c+rad

def relation(new, old):
    a, b = endpoints(new)
    c, d = endpoints(old)
    if b < c or d < a: return "DISJOINT"
    if c <= a and b <= d: return "REPLAY_SUBSET_HISTORICAL"
    if a <= c and d <= b: return "HISTORICAL_SUBSET_REPLAY"
    return "OVERLAP"

def selected_pairs(data, shard):
    pairs = [tuple(p) for p in data["pairs"]]
    if len(pairs) != 19 or len(set(pairs)) != 19: raise ValueError("unexpected sample")
    if set(data["entries"]) != {f"{i},{j}" for i,j in pairs}: raise ValueError("fixture incomplete")
    return [p for n,p in enumerate(pairs) if n % 4 == shard]

def run(args):
    try:
        import flint
        from flint import arb, acb, ctx
    except ImportError as exc:
        raise RuntimeError("python-flint unavailable: not an Arb replay") from exc
    ctx.prec = args.bits
    data = json.loads(Path(args.historical).read_text(encoding="utf-8"))
    if data["degrees"] != [14,22,14]: raise ValueError("wrong degrees")
    expected_idx = [(block,p) for block,n in enumerate([14,22,14]) for p in range(n)]
    pairs = selected_pairs(data, args.shard)
    out = Path(args.out)
    if out.exists(): raise FileExistsError(f"Refuse overwrite: {out}")
    t = arb(2).log()
    a, R, pi = arb(9)/20, arb(30), arb.pi()
    # Mode supports: left, center, right
    pieces = [(-t/2,a-t/2), (arb(0),t-a), (t/2,a-t/2)]
    def oddfactorial(k):
        ans=arb(1)
        for j in range(1,2*k+2,2): ans*=j
        return ans
    def spherical_j(k,z):
        return z**k / oddfactorial(k) * (-z*z/4).hypgeom_0f1(arb(k)+arb(3)/2)
    def spherical_i(k,z):
        return z**k / oddfactorial(k) * (z*z/4).hypgeom_0f1(arb(k)+arb(3)/2)
    def psi_even(z):
        base=acb(arb(1)/4)
        iz=acb(0,1)*z/2
        return ((base+iz).digamma()+(base-iz).digamma())/2
    # Pole stress test for complex analytic continuation
    for sgn in (-1,1):
        if psi_even(acb(0,sgn*arb(1)/2)).is_finite():
            raise ArithmeticError("psi pole must be non-finite")
    psi_R=acb(arb(1)/4,R/2).digamma().real
    m_R=psi_R-pi.log()
    uv={}
    for i,(block,p) in enumerate(expected_idx):
        c,h=pieces[block]
        coef=((arb(2*p+1)*h/2).sqrt())*2*spherical_i(p,h/2)
        uv[i]=((c/2).exp()*coef, (-c/2).exp()*coef*((-1)**p))
    def enc(v):
        m,e=v.mid().man_exp()
        r,f=v.rad().man_exp()
        return [int(m),int(e),int(r),int(f)]
    counts={"evaluations":0,"nonfinite":0}
    records=[]
    for i,j in pairs:
        clock=time.monotonic()
        I,p=expected_idx[i]; J,q=expected_idx[j]
        ci,hi=pieces[I];cj,hj=pieces[J]
        D=ci-cj
        ui,vi=uv[i];uj,vj=uv[j]
        poles=ui*vj+vi*uj
        primes=t/arb(2).sqrt() if {I,J}=={0,2} and p==q and I!=J else arb(0)
        # A is exactly zero for opposite parity on the same subinterval.
        if I==J and (p+q)%2:
            A=arb(0)
        else:
            even=(p+q)%2==0
            phase=(-1)**(((p-q)//2) if even else ((p-q+1)//2))
            def integrand(z, analytic):
                tr=(z*D).cos() if even else (z*D).sin()
                val=spherical_j(p,hi*z)*spherical_j(q,hj*z)*(psi_R-psi_even(z))*tr
                counts["evaluations"]+=1
                if not val.is_finite(): counts["nonfinite"]+=1
                return val
            total=acb(0)
            for left,right in [(0,5),(5,15),(15,30)]:
                term=acb.integral(integrand,acb(left),acb(right),
                                  abs_tol=arb(10)**(-args.abs_exp),
                                  rel_tol=arb(10)**(-args.rel_exp))
                if not term.is_finite():
                    raise ArithmeticError(f"non-finite integral for {i},{j} on {left},{right}")
                total+=term
            A=arb(2)/pi*((hi*hj*arb(2*p+1)*arb(2*q+1)).sqrt())*phase*total.real
        value=(m_R if i==j else arb(0))+poles-primes-A
        new=enc(value)
        old=data["entries"][f"{i},{j}"]
        check=relation(new,old)
        records.append({"pair":[i,j],"new_ball":new,"old_ball":old,
                        "comparison":check,"seconds":round(time.monotonic()-clock,2)})
        print(f"{i},{j}: {check}",flush=True)
        if check=="DISJOINT": raise AssertionError(f"interval disjoint at {i},{j}")
    result={"status":"ARB_EXECUTED","scientific_scope":"sample_interval_consistency_not_RH",
            "shard":args.shard,"bits":args.bits,"python_flint":flint.__version__,
            "python":sys.version,"platform":platform.platform(),
            "historical_sample_sha256":hashlib.sha256(Path(args.historical).read_bytes()).hexdigest(),
            "evaluations":counts,"records":records}
    out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--historical",default="historical_sample.json")
    ap.add_argument("--out",default="REPLAY_ARB_ES8P73.json")
    ap.add_argument("--bits",type=int,default=256)
    ap.add_argument("--abs-exp",type=int,default=34)
    ap.add_argument("--rel-exp",type=int,default=36)
    ap.add_argument("--shard",type=int,choices=range(4),default=0)
    run(ap.parse_args())
