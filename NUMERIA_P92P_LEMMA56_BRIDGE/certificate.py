"""P92-P: second exact implementation of the two-index Schur certificate.

Analytic proof: README.md.  Tests only certify specific rational cases, not Lemmas
5 or 6 in full generality.
"""
from fractions import Fraction as F
from random import Random
from pathlib import Path
import argparse, json, hashlib

SCHEMA = "NUMERIA.P92P.LEMMA56.EXACT_POTENTIAL_SCHUR.1"

def arctan_bounds(z, n=8):
    z=F(z)
    v=sum(((-1 if j%2 else 1)*z**(2*j+1)/F(2*j+1) for j in range(n)),F(0))
    nxt=z**(2*n+1)/(2*n+1)
    return ((v-nxt,v) if n%2 else (v,v+nxt))

def pi_bounds():
    a,b=arctan_bounds(F(1,5))
    c,d=arctan_bounds(F(1,239))
    return 16*a-4*d,16*b-4*c

def gaps(xs):
    return [min(abs(xs[i]-xs[j]) for j in range(len(xs)) if j!=i)
            for i in range(len(xs))]

def energy(xs,ds,i):
    return sum((ds[j]/(xs[i]-xs[j])**2 for j in range(len(xs)) if j!=i),F(0))

def potential(xs,ds,i):
    return sum((ds[j]/(xs[i]-xs[j]) for j in range(len(xs)) if j!=i),F(0))

def ldlt(a):
    """Independent dense fraction-only LDL^T of a symmetric matrix."""
    a=[[F(v) for v in row] for row in a]
    n=len(a)
    assert all(len(row)==n for row in a)
    assert all(a[i][j]==a[j][i] for i in range(n) for j in range(n))
    p=[]
    for k in range(n):
        q=a[k][k]
        p.append(q)
        if q<=0:break
        for i in range(k+1,n):
            for j in range(i,n):
                v=a[i][j]-a[i][k]*a[j][k]/q
                a[i][j]=a[j][i]=v
    return p

def audit(xs,s,t,kappa,do_ldlt):
    xs=[F(x) for x in xs]
    if len(xs)<2 or any(xs[i]>=xs[i+1] for i in range(len(xs)-1)):
        raise ValueError("need strictly increasing points")
    if not 0<=s<t<len(xs):
        raise ValueError("need s<t")
    ds=gaps(xs)
    h=xs[t]-xs[s]
    d=ds[s]+ds[t]
    es,et=energy(xs,ds,s),energy(xs,ds,t)
    ts,tt=potential(xs,ds,s),potential(xs,ds,t)
    others=[i for i in range(len(xs)) if i not in (s,t)]
    y=[1/(xs[i]-xs[s])-1/(xs[i]-xs[t]) for i in others]
    direct_h2=sum((h*h*ds[i]/((xs[i]-xs[s])**2*(xs[i]-xs[t])**2)
                   for i in others),F(0))
    norm_y=sum((ds[i]*yy*yy for i,yy in zip(others,y)),F(0))
    ov=sum((ds[i]/((xs[i]-xs[s])*(xs[i]-xs[t])) for i in others),F(0))
    assert direct_h2==norm_y
    assert ov==(ts-tt)/h+d/h**2
    assert norm_y+3*d/h**2==es+et-2*(ts-tt)/h
    budget=kappa*(1/ds[s]+1/ds[t])-3*d/h**2
    defect_s=kappa/ds[s]-es
    defect_t=kappa/ds[t]-et
    slack=budget-norm_y
    assert slack==defect_s+defect_t+2*(ts-tt)/h
    assert slack>0
    pivot=None
    if do_ldlt:
        n=len(others)+1
        a=[[F(0) for _ in range(n)] for _ in range(n)]
        for j,i in enumerate(others):
            a[j][j]=1/ds[i]
            a[j][-1]=a[-1][j]=y[j]
        a[-1][-1]=budget
        pivots=ldlt(a)
        assert len(pivots)==n and min(pivots)>0 and pivots[-1]==slack
        pivot=str(min(pivots))
    return {"n":len(xs),"s":s,"t":t,"index_separation":t-s,
            "potential_difference":str(ts-tt),
            "h":str(h),"energy_slack_s":str(defect_s),
            "energy_slack_t":str(defect_t),
            "h_squared_kernel":str(norm_y),"rhs_budget":str(budget),
            "slack":str(slack),"ldlt_min_pivot":pivot,
            "monotone_potential_suffices":ts>=tt}

def certify():
    pi_lo,pi_hi=pi_bounds()
    assert F(3)<pi_lo<pi_hi<F(22,7)
    c=pi_lo**2/3
    main=audit([0,8000,8300,32300,32500,32510,33510,33515],0,6,c,True)
    assert main["monotone_potential_suffices"] and F(main["potential_difference"])>F(47,100)
    simple=audit([0,7,8,16,17,32,33],0,5,c,True)
    assert F(simple["potential_difference"])==F(31333,261800)
    uniform=audit(list(range(13)),2,9,c,True)
    assert not uniform["monotone_potential_suffices"]
    cases=[main,simple,uniform]
    rand=Random(20261122)
    for k in range(60):
        n=3+k%10
        x=[0]
        for _ in range(n-1):
            x.append(x[-1]+rand.choice([1,2,3,5,7,13,21,80,300]))
        s=rand.randrange(n-1)
        t=rand.randrange(s+1,n)
        cases.append(audit(x,s,t,c,k<24))
    assert len(cases)==63
    assert sum(c["ldlt_min_pivot"] is not None for c in cases)==27
    return {"schema":SCHEMA,"pi_lower":str(pi_lo),"pi_upper":str(pi_hi),
            "universal_algebraic_bridge_proved":True,
            "rational_schur_equivalence_proved":True,
            "sufficient_potential_monotonicity_proved_conditional_on_energy_lemma":True,
            "all_configurations_lemma5_proved":False,
            "all_configurations_lemma6_independently_proved":False,
            "Riemann_Hypothesis_proved":False,
            "exact_fraction_case_count":63,
            "independent_dense_ldlt_count":27,
            "cases":cases}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    o=certify()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("P92P PASS: 63 exact rational cases and 27 dense LDL")
    print("SHA256",hashlib.sha256(a.out.read_bytes()).hexdigest())
    print("General L5/L6 and RH remain OPEN.")

if __name__=="__main__":
    main()
